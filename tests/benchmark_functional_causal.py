"""
Functional-causal benchmark: L1 vs SuSiE vs FINEMAP.

Selects causal variants from eQTL-annotated variants (realistic).
In real biology, causal variants are enriched for functional annotations.
This is where L1's annotation integration should show its advantage.

Usage:
    python tests/benchmark_functional_causal.py 2>/dev/null
"""

import csv, time, json, os, subprocess, tempfile, shutil
import numpy as np
from scipy import stats as sp_stats

import graphgwas.config as cfg
cfg._N_SAMPLES_DETECTED = False

from graphgwas.db import GraphGWASConnection
from graphgwas.config import detect_n_samples
from graphgwas.finemapping_v2 import (
    _load_locus_variants, _compute_association_stats, _compute_ld_matrix,
    FinemapCandidate,
)
from graphgwas.genotype import get_all_indices, get_phenotype_values, build_dosage

OUTDIR = "/mnt/data/GraphGWAS/results/benchmark_v2/tool_comparison"
os.makedirs(OUTDIR, exist_ok=True)

FINEMAP_BIN = "/home/jfmao/bin/finemap"
EQTL_CACHE = "/mnt/data/GraphGWAS/data/annotations/gtex_chr22_enhanced_cache.json"

WINDOW = 50000
N_REPS = 30

print("Loading eQTL annotation cache...", flush=True)
with open(EQTL_CACHE) as f:
    EQTL_DATA = json.load(f)
print(f"  {len(EQTL_DATA)} annotations loaded")


def _simulate_functional_causal(conn, all_idx, eqtl_variants_in_locus,
                                  beta=0.3, h2_target=0.05, seed=42):
    """Simulate phenotype with causal variant chosen from eQTL variants."""
    rng = np.random.default_rng(seed)

    if not eqtl_variants_in_locus:
        return None, None

    # Pick a random eQTL variant as causal
    causal = rng.choice(eqtl_variants_in_locus)
    dosage = causal["dosage"].copy()
    dosage = np.where(np.isnan(dosage), np.nanmean(dosage), dosage)

    af = causal.get("af_total", 0)
    if af < 0.02 or af > 0.98:
        return None, None

    genetic = beta * dosage
    var_g = np.var(genetic)
    if var_g < 1e-10:
        return None, None
    var_e = var_g * (1 / h2_target - 1)
    noise = rng.normal(0, np.sqrt(var_e), len(all_idx))
    phenotype = genetic + noise

    return phenotype, causal


def _store_phenotype(conn, phenotype, all_idx, label):
    """Store simulated phenotype on Sample nodes."""
    conn.execute_write("""
        UNWIND $data AS row
        MATCH (s:Sample {packed_index: row.idx})
        SET s.simulated_phenotype = row.val
    """, {"data": [{"idx": int(all_idx[i]), "val": float(phenotype[i])}
                   for i in range(len(all_idx))]})


def _fast_l1_finemap(variants, pheno, r2_smooth=0.3, alpha=0.5,
                      use_annotations=True, credible_set_coverage=0.95):
    """Fast L1 fine-mapping with file-based annotations."""
    n_var = len(variants)
    if n_var < 3:
        return []

    z_stats = _compute_association_stats(variants, pheno)
    R = _compute_ld_matrix(variants)

    unique_stats = np.zeros(n_var)
    n_neighbors = np.zeros(n_var, dtype=int)
    for i in range(n_var):
        neighbors = np.where(R[i] > r2_smooth)[0]
        neighbors = neighbors[neighbors != i]
        n_neighbors[i] = len(neighbors)
        if len(neighbors) > 0:
            ld_contribution = np.sum(R[i, neighbors] * z_stats[neighbors])
            unique_stats[i] = z_stats[i] - ld_contribution / len(neighbors)
        else:
            unique_stats[i] = z_stats[i]

    if unique_stats.max() > 0:
        unique_stats = unique_stats / unique_stats.max() * z_stats.max()
    unique_stats = np.maximum(unique_stats, 0)

    z_func = np.zeros(n_var)
    if use_annotations:
        for i, v in enumerate(variants):
            vid = v["variantId"]
            eqtl = EQTL_DATA.get(vid)
            if eqtl:
                score = eqtl.get("composite", 0)
                n_tissues = eqtl.get("n_tissues", 0)
                if 1 <= n_tissues <= 3:
                    score *= 1.5
                elif n_tissues >= 8:
                    score *= 0.3
                z_func[i] = np.log1p(score)

        if z_func.max() > 0:
            z_func = z_func / z_func.max() * z_stats.max()

    if use_annotations and z_func.max() > 0:
        combined = alpha * unique_stats + (1 - alpha) * z_func
    else:
        combined = unique_stats

    exp_scores = np.exp(combined - combined.max())
    pip = exp_scores / exp_scores.sum()

    order = np.argsort(-pip)
    cumsum = np.cumsum(pip[order])
    in_cs = np.zeros(n_var, dtype=bool)
    for k, idx in enumerate(order):
        in_cs[idx] = True
        if cumsum[k] >= credible_set_coverage:
            break

    candidates = []
    for i in range(n_var):
        candidates.append(FinemapCandidate(
            variant_id=variants[i]["variantId"],
            chr="chr22", pos=variants[i]["pos"],
            ref=variants[i].get("ref", ""),
            alt=variants[i].get("alt", ""),
            af=variants[i].get("af_total", 0),
            z_stat=float(z_stats[i]),
            z_functional=float(z_func[i]),
            unique_stat=float(unique_stats[i]),
            combined_score=float(combined[i]),
            pip=float(pip[i]),
            in_credible_set=bool(in_cs[i]),
            annotations=[], n_ld_neighbors=int(n_neighbors[i]),
        ))
    candidates.sort(key=lambda c: -c.combined_score)
    return candidates


def _build_finemap_inputs(variants, pheno, tmpdir, n_samples):
    n_var = len(variants)
    D = np.array([v["dosage"] for v in variants], dtype=np.float64)
    row_means = np.nanmean(D, axis=1, keepdims=True)
    D = np.where(np.isnan(D), row_means, D)
    y = pheno.copy()
    valid = ~np.isnan(y)
    y_clean = np.where(np.isnan(y), 0, y)
    D_v, y_v = D[:, valid], y_clean[valid]
    n = int(valid.sum())

    y_c = y_v - y_v.mean()
    D_c = D_v - D_v.mean(axis=1, keepdims=True)
    var_g = np.var(D_v, axis=1)
    safe_var = np.where(var_g > 1e-10, var_g, 1.0)
    beta = (D_c @ y_c) / n / safe_var
    predicted = D_c * beta[:, np.newaxis]
    mse = np.sum((y_c - predicted)**2, axis=1) / (n - 2)
    se = np.sqrt(mse / (n * safe_var))
    safe_se = np.where(se > 1e-10, se, 1.0)
    keep = var_g > 1e-10
    keep_idx = np.where(keep)[0]

    with open(f"{tmpdir}/data.z", "w") as f:
        f.write("rsid chromosome position allele1 allele2 maf beta se\n")
        for i in keep_idx:
            v = variants[i]
            rsid = v["variantId"].replace(":", "_").replace("/", "_")
            maf = min(v.get("af_total", 0.1), 1 - v.get("af_total", 0.1))
            ref = v.get("ref", "A")[:1] or "A"
            alt = v.get("alt", "G")[:1] or "G"
            f.write(f"{rsid} 22 {v['pos']} {ref} {alt} {max(maf,0.001):.6f} "
                    f"{beta[i]:.6f} {safe_se[i]:.6f}\n")

    stds = D_v.std(axis=1, keepdims=True)
    stds = np.where(stds < 1e-10, 1.0, stds)
    D_std = (D_v - D_v.mean(axis=1, keepdims=True)) / stds
    R = (D_std @ D_std.T) / D_v.shape[1]
    np.fill_diagonal(R, 1.0)
    R_keep = R[np.ix_(keep, keep)]

    with open(f"{tmpdir}/data.ld", "w") as f:
        for i in range(R_keep.shape[0]):
            f.write(" ".join(f"{R_keep[i,j]:.6f}" for j in range(R_keep.shape[1])) + "\n")

    with open(f"{tmpdir}/data.master", "w") as f:
        f.write("z;ld;snp;config;cred;log;n_samples\n")
        f.write(f"{tmpdir}/data.z;{tmpdir}/data.ld;{tmpdir}/data.snp;{tmpdir}/data.config;"
                f"{tmpdir}/data.cred;{tmpdir}/data.log;{n_samples}\n")

    return f"{tmpdir}/data.master"


def _run_finemap(master_path, tmpdir, causal_vid):
    causal_rsid = causal_vid.replace(":", "_").replace("/", "_")
    try:
        proc = subprocess.run(
            [FINEMAP_BIN, "--sss", "--in-files", master_path, "--n-causal-snps", "5"],
            capture_output=True, text=True, timeout=120, cwd=tmpdir)

        snp_path = f"{tmpdir}/data.snp"
        if not os.path.exists(snp_path) or os.path.getsize(snp_path) == 0:
            return -1, 0.0, 0, False

        pips = {}
        with open(snp_path) as f:
            f.readline()
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 11:
                    pips[parts[1]] = float(parts[10])

        if not pips:
            return -1, 0.0, 0, False

        causal_pip = pips.get(causal_rsid, 0.0)
        rank = sum(1 for p in pips.values() if p >= causal_pip) if causal_pip > 0 else -1

        sorted_pips = sorted(pips.items(), key=lambda x: -x[1])
        cs = set()
        cumsum = 0
        for rsid, pip in sorted_pips:
            cs.add(rsid)
            cumsum += pip
            if cumsum >= 0.95:
                break
        return rank, causal_pip, len(cs), causal_rsid in cs

    except Exception:
        return -1, 0.0, 0, False


def _run_susie(variants, pheno, all_idx, causal_idx):
    n_var = len(variants)
    X = np.zeros((len(all_idx), n_var))
    for i, v in enumerate(variants):
        d = v["dosage"]
        X[:, i] = np.where(np.isnan(d), np.nanmean(d), d)
    y = np.where(np.isnan(pheno), np.nanmean(pheno), pheno)

    tmp = tempfile.mkdtemp()
    np.savetxt(f"{tmp}/X.csv", X, delimiter=",")
    np.savetxt(f"{tmp}/y.csv", y, delimiter=",")

    r_script = f'''
    library(susieR)
    X <- as.matrix(read.csv("{tmp}/X.csv", header=FALSE))
    y <- scan("{tmp}/y.csv")
    fit <- tryCatch(susie(X, y, L=5, verbose=FALSE), error=function(e) NULL)
    if(!is.null(fit)) {{
        pips <- fit$pip
        ci <- {causal_idx + 1}
        cat(sprintf("PIP=%.6f\\n", pips[ci]))
        cat(sprintf("RANK=%d\\n", sum(pips >= pips[ci])))
        cs <- susie_get_cs(fit)
        in_cs <- FALSE; cs_size <- 0
        for(i in seq_along(cs$cs)) if(ci %in% cs$cs[[i]]) {{
            in_cs <- TRUE; cs_size <- length(cs$cs[[i]]); break
        }}
        cat(sprintf("IN_CS=%s\\n", in_cs))
        cat(sprintf("CS_SIZE=%d\\n", cs_size))
    }} else cat("FAILED\\n")
    '''
    try:
        proc = subprocess.run(
            ["R", "--no-save", "--no-restore", "-e", r_script],
            capture_output=True, text=True, timeout=120)
        pip = rank = cs_size = 0
        in_cs = False
        for line in proc.stdout.split("\n"):
            if line.startswith("PIP="):
                pip = float(line.split("=")[1])
            elif line.startswith("RANK="):
                rank = int(line.split("=")[1])
            elif line.startswith("IN_CS="):
                in_cs = line.split("=")[1].strip() == "TRUE"
            elif line.startswith("CS_SIZE="):
                cs_size = int(line.split("=")[1])
        return rank, pip, cs_size, in_cs
    except Exception:
        return -1, 0.0, 0, False
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def run_benchmark():
    results = []
    t0_total = time.time()

    # Find loci with eQTL-annotated variants
    # Scan chr22 in 50kb windows and find regions with >= 20 eQTL variants
    eqtl_positions = sorted([int(k.split(":")[1]) for k in EQTL_DATA.keys()
                              if k.startswith("chr22:")])
    print(f"eQTL positions on chr22: {len(eqtl_positions)}")

    # Group into 50kb windows
    windows = {}
    for pos in eqtl_positions:
        wkey = (pos // WINDOW) * WINDOW
        windows.setdefault(wkey, []).append(pos)

    # Select windows with many eQTL variants
    good_windows = [(k, len(v)) for k, v in windows.items() if len(v) >= 10]
    good_windows.sort(key=lambda x: -x[1])
    print(f"Windows with >=10 eQTLs: {len(good_windows)}")
    centers = [w[0] + WINDOW // 2 for w in good_windows[:N_REPS]]

    print(f"FUNCTIONAL-CAUSAL BENCHMARK: causal selected from eQTL variants")
    print(f"beta=0.3, h2=0.05, {len(centers)} windows, chr22")
    print("=" * 90)

    with GraphGWASConnection("bolt://localhost:7688", "neo4j", "graphgwas") as conn:
        detect_n_samples(conn)
        n_samples = cfg.N_SAMPLES
        all_idx = get_all_indices(conn)
        print(f"N_SAMPLES={n_samples}")

        for rep, center in enumerate(centers):
            seed = 200 + rep

            # Load variants in locus
            variants = _load_locus_variants(conn, "chr22", center, WINDOW, all_idx)
            n_var = len(variants)

            # Find eQTL-annotated variants in this locus
            eqtl_vars = [v for v in variants if v["variantId"] in EQTL_DATA
                          and 0.02 < v.get("af_total", 0) < 0.98]

            if len(eqtl_vars) < 5 or n_var < 20:
                print(f"  Rep {rep}: skip (eqtl_vars={len(eqtl_vars)}, nv={n_var})")
                continue

            # Simulate with eQTL variant as causal
            pheno, causal_v = _simulate_functional_causal(
                conn, all_idx, eqtl_vars, beta=0.3, h2_target=0.05, seed=seed)
            if pheno is None:
                print(f"  Rep {rep}: sim failed")
                continue

            causal_vid = causal_v["variantId"]
            causal_af = causal_v.get("af_total", 0)
            causal_idx = next((i for i, v in enumerate(variants) if v["variantId"] == causal_vid), -1)

            if causal_idx < 0:
                print(f"  Rep {rep}: causal not found in variants")
                continue

            # Store phenotype
            _store_phenotype(conn, pheno, all_idx, f"fc_rep{rep}")

            row = {"rep": rep, "n_var": n_var, "n_eqtl_in_locus": len(eqtl_vars),
                   "causal_af": causal_af, "causal_pos": causal_v["pos"],
                   "causal_eqtl_score": EQTL_DATA[causal_vid].get("composite", 0)}

            # L1 stat
            t0 = time.time()
            l1s = _fast_l1_finemap(variants, pheno, use_annotations=False)
            row["l1s_time"] = time.time() - t0
            row["l1s_rank"] = next((i+1 for i, c in enumerate(l1s) if c.variant_id == causal_vid), -1)
            row["l1s_pip"] = next((c.pip for c in l1s if c.variant_id == causal_vid), 0)
            row["l1s_cs"] = sum(1 for c in l1s if c.in_credible_set)
            row["l1s_in_cs"] = any(c.variant_id == causal_vid and c.in_credible_set for c in l1s)

            # L1 annotated
            t0 = time.time()
            l1a = _fast_l1_finemap(variants, pheno, use_annotations=True)
            row["l1a_time"] = time.time() - t0
            row["l1a_rank"] = next((i+1 for i, c in enumerate(l1a) if c.variant_id == causal_vid), -1)
            row["l1a_pip"] = next((c.pip for c in l1a if c.variant_id == causal_vid), 0)
            row["l1a_cs"] = sum(1 for c in l1a if c.in_credible_set)
            row["l1a_in_cs"] = any(c.variant_id == causal_vid and c.in_credible_set for c in l1a)

            # FINEMAP
            tmpdir = tempfile.mkdtemp()
            try:
                t0 = time.time()
                master = _build_finemap_inputs(variants, pheno, tmpdir, n_samples)
                fm_rank, fm_pip, fm_cs, fm_in_cs = _run_finemap(master, tmpdir, causal_vid)
                row["fm_time"] = time.time() - t0
                row["fm_rank"] = fm_rank
                row["fm_pip"] = fm_pip
                row["fm_cs"] = fm_cs
                row["fm_in_cs"] = fm_in_cs
            finally:
                shutil.rmtree(tmpdir, ignore_errors=True)

            # SuSiE
            t0 = time.time()
            su_rank, su_pip, su_cs, su_in_cs = _run_susie(variants, pheno, all_idx, causal_idx)
            row["su_time"] = time.time() - t0
            row["su_rank"] = su_rank
            row["su_pip"] = su_pip
            row["su_cs"] = su_cs
            row["su_in_cs"] = su_in_cs

            results.append(row)
            print(f"  Rep {rep}: nv={n_var} eqtl={len(eqtl_vars)} AF={causal_af:.2f} | "
                  f"L1s=#{row['l1s_rank']} L1a=#{row['l1a_rank']} "
                  f"FM=#{row['fm_rank']} Su=#{row['su_rank']}", flush=True)

    total_time = time.time() - t0_total

    print(f"\n{'='*90}")
    print(f"FUNCTIONAL-CAUSAL BENCHMARK — {len(results)} reps, {total_time:.0f}s")
    print("=" * 90)

    for name, pre in [("L1 (stat)", "l1s"), ("L1 (annot)", "l1a"),
                       ("FINEMAP", "fm"), ("SuSiE", "su")]:
        ranks = [r[f"{pre}_rank"] for r in results if r.get(f"{pre}_rank", -1) > 0]
        pips = [r[f"{pre}_pip"] for r in results if r.get(f"{pre}_rank", -1) > 0]
        in_cs = [r.get(f"{pre}_in_cs", False) for r in results]
        css = [r[f"{pre}_cs"] for r in results if r.get(f"{pre}_cs", 0) > 0]
        times = [r[f"{pre}_time"] for r in results if f"{pre}_time" in r]
        rank1 = sum(1 for r in ranks if r == 1)
        n = len(results)

        print(f"\n{name}:")
        print(f"  Valid: {len(ranks)}/{n}")
        if ranks:
            print(f"  Mean rank: {np.mean(ranks):.1f} | Median: {np.median(ranks):.0f}")
            print(f"  Rank #1: {rank1}/{n} ({rank1/n*100:.0f}%)")
            print(f"  Mean PIP: {np.mean(pips):.4f}")
        print(f"  Coverage: {sum(in_cs)}/{n} ({np.mean(in_cs)*100:.0f}%)")
        if css:
            print(f"  Mean CS: {np.mean(css):.1f}")
        if times:
            print(f"  Mean time: {np.mean(times):.2f}s")

    # Head-to-head
    print(f"\n{'='*90}")
    print("HEAD-TO-HEAD")
    print("=" * 90)
    for a_name, a_pre, b_name, b_pre in [
        ("L1(annot)", "l1a", "SuSiE", "su"),
        ("L1(annot)", "l1a", "FINEMAP", "fm"),
        ("L1(stat)", "l1s", "SuSiE", "su"),
        ("L1(stat)", "l1s", "FINEMAP", "fm"),
        ("L1(annot)", "l1a", "L1(stat)", "l1s"),
    ]:
        a_wins = b_wins = ties = 0
        for r in results:
            a_r = r.get(f"{a_pre}_rank", -1)
            b_r = r.get(f"{b_pre}_rank", -1)
            if a_r > 0 and b_r > 0:
                if a_r < b_r: a_wins += 1
                elif b_r < a_r: b_wins += 1
                else: ties += 1
        t = a_wins + b_wins + ties
        if t > 0:
            print(f"  {a_name} vs {b_name}: {a_wins}:{b_wins}:{ties} "
                  f"[{a_wins/t*100:.0f}% win]")

    with open(f"{OUTDIR}/functional_causal_comparison.json", "w") as f:
        json.dump(results, f, indent=2, default=str)

    print(f"\nSaved to {OUTDIR}/")
    return results


if __name__ == "__main__":
    run_benchmark()
