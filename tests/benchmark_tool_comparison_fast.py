"""
Fast L1 vs SuSiE vs FINEMAP benchmark.

Avoids slow Neo4j annotation queries by:
1. Using pure statistical L1 (no annotation overhead)
2. Running annotation-enhanced L1 with file-based cache
3. SuSiE and FINEMAP get the same data

Usage:
    python tests/benchmark_tool_comparison_fast.py 2>/dev/null
"""

import csv, time, json, os, subprocess, tempfile, shutil
import numpy as np
from scipy import stats as sp_stats

import graphgwas.config as cfg
cfg._N_SAMPLES_DETECTED = False

from graphgwas.db import GraphGWASConnection
from graphgwas.config import detect_n_samples
from graphgwas.simulate import simulate_F1_single_causal_in_ld, load_simulation_as_phenotype
from graphgwas.finemapping_v2 import (
    _load_locus_variants, _compute_association_stats, _compute_ld_matrix,
    FinemapCandidate,
)
from graphgwas.genotype import get_all_indices, get_phenotype_values

OUTDIR = "/mnt/data/GraphGWAS/results/benchmark_v2/tool_comparison"
os.makedirs(OUTDIR, exist_ok=True)

FINEMAP_BIN = "/home/jfmao/bin/finemap"
EQTL_CACHE = "/mnt/data/GraphGWAS/data/annotations/gtex_chr22_enhanced_cache.json"

WINDOW = 50000
N_REPS = 20
CENTERS = [17000000, 19000000, 21000000, 23000000, 25000000,
           27000000, 29000000, 31000000, 33000000, 35000000,
           37000000, 39000000, 41000000, 43000000, 45000000,
           17500000, 20000000, 24000000, 28000000, 36000000]

# Load eQTL cache once
print("Loading eQTL annotation cache...", flush=True)
with open(EQTL_CACHE) as f:
    EQTL_DATA = json.load(f)
print(f"  {len(EQTL_DATA)} variant annotations loaded")


def _fast_l1_finemap(variants, pheno, r2_smooth=0.3, alpha=0.5,
                      use_annotations=True, credible_set_coverage=0.95):
    """Fast L1 fine-mapping without Neo4j annotation queries."""
    n_var = len(variants)
    if n_var < 3:
        return []

    # Step 1: Association statistics
    z_stats = _compute_association_stats(variants, pheno)

    # Step 2: LD matrix
    R = _compute_ld_matrix(variants)

    # Step 3: LD deconvolution
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

    # Step 4: Annotation scores from file cache
    z_func = np.zeros(n_var)
    if use_annotations:
        for i, v in enumerate(variants):
            vid = v["variantId"]
            # Try to match variant ID to eQTL cache
            # variantId format: chr22:pos:ref:alt or chromosome22:pos:ref/alt
            eqtl = EQTL_DATA.get(vid)
            if eqtl is None:
                # Try alternate format
                pos = v["pos"]
                ref = v.get("ref", "")
                alt = v.get("alt", "")
                for fmt in [f"chr22:{pos}:{ref}:{alt}", f"chr22:{pos}:{ref}/{alt}"]:
                    eqtl = EQTL_DATA.get(fmt)
                    if eqtl:
                        break

            if eqtl:
                score = eqtl.get("composite", 0)
                n_tissues = eqtl.get("n_tissues", 0)
                # Tissue specificity: 1-2 tissues = regulatory, 8+ = LD proxy
                if 1 <= n_tissues <= 3:
                    score *= 1.5  # boost tissue-specific
                elif n_tissues >= 8:
                    score *= 0.3  # penalize ubiquitous (likely LD proxy)
                z_func[i] = np.log1p(score)

        if z_func.max() > 0:
            z_func = z_func / z_func.max() * z_stats.max()

    # Step 5: Combined score (Bayesian: log_posterior = unique_stats + prior)
    if use_annotations and z_func.max() > 0:
        combined = alpha * unique_stats + (1 - alpha) * z_func
    else:
        combined = unique_stats

    # Step 6: PIP
    exp_scores = np.exp(combined - combined.max())
    pip = exp_scores / exp_scores.sum()

    # Step 7: Credible set
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
            chr="chr22",
            pos=variants[i]["pos"],
            ref=variants[i].get("ref", ""),
            alt=variants[i].get("alt", ""),
            af=variants[i].get("af_total", 0),
            z_stat=float(z_stats[i]),
            z_functional=float(z_func[i]),
            unique_stat=float(unique_stats[i]),
            combined_score=float(combined[i]),
            pip=float(pip[i]),
            in_credible_set=bool(in_cs[i]),
            annotations=[],
            n_ld_neighbors=int(n_neighbors[i]),
        ))
    candidates.sort(key=lambda c: -c.combined_score)
    return candidates


def _build_finemap_inputs(variants, pheno, tmpdir, n_samples):
    """Build FINEMAP input files."""
    n_var = len(variants)
    D = np.array([v["dosage"] for v in variants], dtype=np.float64)
    row_means = np.nanmean(D, axis=1, keepdims=True)
    D = np.where(np.isnan(D), row_means, D)

    y = pheno.copy()
    valid = ~np.isnan(y)
    y_clean = np.where(np.isnan(y), 0, y)
    D_v = D[:, valid]
    y_v = y_clean[valid]
    n = int(valid.sum())

    y_centered = y_v - y_v.mean()
    D_centered = D_v - D_v.mean(axis=1, keepdims=True)
    var_g = np.var(D_v, axis=1)
    safe_var = np.where(var_g > 1e-10, var_g, 1.0)
    cov_gy = (D_centered @ y_centered) / n
    beta = cov_gy / safe_var
    predicted = D_centered * beta[:, np.newaxis]
    residuals = y_centered[np.newaxis, :] - predicted
    mse = np.sum(residuals**2, axis=1) / (n - 2)
    se = np.sqrt(mse / (n * safe_var))
    safe_se = np.where(se > 1e-10, se, 1.0)

    keep = var_g > 1e-10
    keep_idx = np.where(keep)[0]
    n_written = int(keep.sum())

    z_path = os.path.join(tmpdir, "data.z")
    with open(z_path, "w") as f:
        f.write("rsid chromosome position allele1 allele2 maf beta se\n")
        for i in keep_idx:
            v = variants[i]
            rsid = v["variantId"].replace(":", "_").replace("/", "_")
            pos = v["pos"]
            ref = v.get("ref", "A")[:1] or "A"
            alt = v.get("alt", "G")[:1] or "G"
            maf = min(v.get("af_total", 0.1), 1 - v.get("af_total", 0.1))
            maf = max(maf, 0.001)
            f.write(f"{rsid} 22 {pos} {ref} {alt} {maf:.6f} {beta[i]:.6f} {safe_se[i]:.6f}\n")

    # LD correlation matrix (signed, not r²)
    stds = D_v.std(axis=1, keepdims=True)
    stds = np.where(stds < 1e-10, 1.0, stds)
    D_std = (D_v - D_v.mean(axis=1, keepdims=True)) / stds
    R = (D_std @ D_std.T) / D_v.shape[1]
    np.fill_diagonal(R, 1.0)
    R_keep = R[np.ix_(keep, keep)]

    ld_path = os.path.join(tmpdir, "data.ld")
    with open(ld_path, "w") as f:
        for i in range(R_keep.shape[0]):
            f.write(" ".join(f"{R_keep[i, j]:.6f}" for j in range(R_keep.shape[1])) + "\n")

    master_path = os.path.join(tmpdir, "data.master")
    with open(master_path, "w") as f:
        f.write("z;ld;snp;config;cred;log;n_samples\n")
        f.write(f"{z_path};{ld_path};{tmpdir}/data.snp;{tmpdir}/data.config;"
                f"{tmpdir}/data.cred;{tmpdir}/data.log;{n_samples}\n")

    # Return mapping from keep_idx to rsid for result parsing
    rsid_map = {}
    for i in keep_idx:
        rsid = variants[i]["variantId"].replace(":", "_").replace("/", "_")
        rsid_map[rsid] = variants[i]["variantId"]

    return master_path, n_written, rsid_map


def _run_finemap(master_path, tmpdir, causal_vid):
    """Run FINEMAP and parse results."""
    causal_rsid = causal_vid.replace(":", "_").replace("/", "_")
    try:
        proc = subprocess.run(
            [FINEMAP_BIN, "--sss", "--in-files", master_path,
             "--n-causal-snps", "5"],
            capture_output=True, text=True, timeout=120, cwd=tmpdir)

        snp_path = os.path.join(tmpdir, "data.snp")
        if not os.path.exists(snp_path) or os.path.getsize(snp_path) == 0:
            return -1, 0.0, 0, False, proc.stderr[:200]

        pips = {}
        with open(snp_path) as f:
            header = f.readline()
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 11:
                    pips[parts[1]] = float(parts[10])

        if not pips:
            return -1, 0.0, 0, False, "empty"

        causal_pip = pips.get(causal_rsid, 0.0)
        rank = sum(1 for p in pips.values() if p >= causal_pip) if causal_pip > 0 else -1

        cred_path = os.path.join(tmpdir, "data.cred")
        in_cs = False
        cs_size = 0
        if os.path.exists(cred_path) and os.path.getsize(cred_path) > 0:
            with open(cred_path) as f:
                header = f.readline()
                for line in f:
                    parts = line.strip().split()
                    if causal_rsid in parts:
                        in_cs = True
                    cs_size = max(cs_size, len(parts) - 1)

        return rank, causal_pip, cs_size, in_cs, ""

    except subprocess.TimeoutExpired:
        return -1, 0.0, 0, False, "timeout"
    except Exception as e:
        return -1, 0.0, 0, False, str(e)[:100]


def _run_susie(variants, pheno, all_idx, causal_idx):
    """Run SuSiE via R."""
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

    has_finemap = os.path.exists(FINEMAP_BIN)
    has_susie = shutil.which("R") is not None

    print(f"Tools: L1(stat+annot) FINEMAP={has_finemap} SuSiE(R)={has_susie}")
    print(f"{N_REPS} replicates, window={WINDOW}bp, chr22")
    print("=" * 90)

    with GraphGWASConnection("bolt://localhost:7688", "neo4j", "graphgwas") as conn:
        detect_n_samples(conn)
        n_samples = cfg.N_SAMPLES
        all_idx = get_all_indices(conn)
        print(f"N_SAMPLES={n_samples}, n_idx={len(all_idx)}")

        for rep in range(N_REPS):
            center = CENTERS[rep % len(CENTERS)]
            seed = 42 + rep

            try:
                sim = simulate_F1_single_causal_in_ld(
                    conn, "chr22", locus_center=center,
                    beta=0.5, h2_target=0.1, seed=seed)
            except Exception as e:
                print(f"  Rep {rep}: sim failed ({e})")
                continue

            causal = sim.causal_variants[0]
            load_simulation_as_phenotype(conn, sim, f"tc_rep{rep}")
            pheno = get_phenotype_values(conn, all_idx)

            # Load variants once
            variants = _load_locus_variants(conn, "chr22", causal.pos, WINDOW, all_idx)
            n_var = len(variants)
            causal_idx = next((i for i, v in enumerate(variants)
                               if v["variantId"] == causal.variant_id), -1)

            if causal_idx < 0 or n_var < 10:
                print(f"  Rep {rep}: skip (causal_idx={causal_idx}, nv={n_var})")
                continue

            row = {"rep": rep, "n_var": n_var, "causal_af": causal.af,
                   "causal_pos": causal.pos}

            # === L1 (pure statistical, no annotations) ===
            t0 = time.time()
            l1_stat = _fast_l1_finemap(variants, pheno, use_annotations=False)
            row["l1_stat_time"] = time.time() - t0
            row["l1_stat_rank"] = next((i+1 for i, c in enumerate(l1_stat)
                                         if c.variant_id == causal.variant_id), -1)
            row["l1_stat_pip"] = next((c.pip for c in l1_stat
                                        if c.variant_id == causal.variant_id), 0)
            row["l1_stat_cs"] = sum(1 for c in l1_stat if c.in_credible_set)
            row["l1_stat_in_cs"] = any(c.variant_id == causal.variant_id
                                        and c.in_credible_set for c in l1_stat)

            # === L1 (with eQTL annotations from file cache) ===
            t0 = time.time()
            l1_annot = _fast_l1_finemap(variants, pheno, use_annotations=True)
            row["l1_annot_time"] = time.time() - t0
            row["l1_annot_rank"] = next((i+1 for i, c in enumerate(l1_annot)
                                          if c.variant_id == causal.variant_id), -1)
            row["l1_annot_pip"] = next((c.pip for c in l1_annot
                                         if c.variant_id == causal.variant_id), 0)
            row["l1_annot_cs"] = sum(1 for c in l1_annot if c.in_credible_set)
            row["l1_annot_in_cs"] = any(c.variant_id == causal.variant_id
                                         and c.in_credible_set for c in l1_annot)

            # === FINEMAP ===
            if has_finemap:
                tmpdir = tempfile.mkdtemp()
                try:
                    t0 = time.time()
                    master_path, n_written, rsid_map = _build_finemap_inputs(
                        variants, pheno, tmpdir, n_samples)
                    fm_rank, fm_pip, fm_cs, fm_in_cs, fm_err = _run_finemap(
                        master_path, tmpdir, causal.variant_id)
                    row["finemap_time"] = time.time() - t0
                    row["finemap_rank"] = fm_rank
                    row["finemap_pip"] = fm_pip
                    row["finemap_cs"] = fm_cs
                    row["finemap_in_cs"] = fm_in_cs
                    if fm_err:
                        row["finemap_err"] = fm_err
                finally:
                    shutil.rmtree(tmpdir, ignore_errors=True)

            # === SuSiE ===
            if has_susie:
                t0 = time.time()
                su_rank, su_pip, su_cs, su_in_cs = _run_susie(
                    variants, pheno, all_idx, causal_idx)
                row["susie_time"] = time.time() - t0
                row["susie_rank"] = su_rank
                row["susie_pip"] = su_pip
                row["susie_cs"] = su_cs
                row["susie_in_cs"] = su_in_cs

            results.append(row)

            # Progress
            parts = [f"Rep {rep}: nv={n_var} AF={causal.af:.2f}"]
            parts.append(f"L1s=#{row['l1_stat_rank']}")
            parts.append(f"L1a=#{row['l1_annot_rank']}")
            if "finemap_rank" in row:
                parts.append(f"FM=#{row['finemap_rank']}" + (f"[{row.get('finemap_err','')}]" if row.get('finemap_err') else ""))
            if "susie_rank" in row:
                parts.append(f"Su=#{row['susie_rank']}")
            print("  " + " | ".join(parts), flush=True)

    total_time = time.time() - t0_total

    # Summary
    print(f"\n{'='*90}")
    print(f"TOOL COMPARISON — {len(results)} replicates, {total_time:.0f}s total")
    print("=" * 90)

    methods = [
        ("L1 (statistical)", "l1_stat"),
        ("L1 (annotated)", "l1_annot"),
        ("FINEMAP v1.4.2", "finemap"),
        ("SuSiE", "susie"),
    ]
    for name, pre in methods:
        if not any(f"{pre}_rank" in r for r in results):
            continue
        ranks = [r[f"{pre}_rank"] for r in results if r.get(f"{pre}_rank", -1) > 0]
        pips = [r[f"{pre}_pip"] for r in results if r.get(f"{pre}_rank", -1) > 0]
        in_cs = [r.get(f"{pre}_in_cs", False) for r in results]
        css = [r[f"{pre}_cs"] for r in results if r.get(f"{pre}_cs", 0) > 0]
        times = [r[f"{pre}_time"] for r in results if f"{pre}_time" in r]
        rank1 = sum(1 for r in ranks if r == 1)
        n_total = len(results)

        print(f"\n{name}:")
        print(f"  Valid results: {len(ranks)}/{n_total}")
        if ranks:
            print(f"  Mean rank: {np.mean(ranks):.1f} | Median rank: {np.median(ranks):.0f}")
            print(f"  Rank #1: {rank1}/{n_total} ({rank1/n_total*100:.0f}%)")
            print(f"  Mean PIP: {np.mean(pips):.4f}")
        print(f"  Coverage: {sum(in_cs)}/{n_total} ({np.mean(in_cs)*100:.0f}%)")
        if css:
            print(f"  Mean CS size: {np.mean(css):.1f}")
        if times:
            print(f"  Mean time: {np.mean(times):.2f}s")

    # Head-to-head
    print(f"\n{'='*90}")
    print("HEAD-TO-HEAD (per replicate)")
    print("=" * 90)

    for a_name, a_pre, b_name, b_pre in [
        ("L1(annot)", "l1_annot", "SuSiE", "susie"),
        ("L1(annot)", "l1_annot", "FINEMAP", "finemap"),
        ("L1(stat)", "l1_stat", "SuSiE", "susie"),
    ]:
        a_wins = b_wins = ties = 0
        for r in results:
            a_rank = r.get(f"{a_pre}_rank", -1)
            b_rank = r.get(f"{b_pre}_rank", -1)
            if a_rank > 0 and b_rank > 0:
                if a_rank < b_rank:
                    a_wins += 1
                elif b_rank < a_rank:
                    b_wins += 1
                else:
                    ties += 1
        total = a_wins + b_wins + ties
        if total > 0:
            print(f"  {a_name} vs {b_name}: {a_wins}:{b_wins}:{ties} (win:loss:tie) "
                  f"[{a_wins/total*100:.0f}% win rate]")

    # Save
    with open(f"{OUTDIR}/l1_vs_susie_vs_finemap_v2.json", "w") as f:
        json.dump(results, f, indent=2, default=str)
    with open(f"{OUTDIR}/l1_vs_susie_vs_finemap_v2.tsv", "w", newline="") as f:
        if results:
            keys = sorted(set().union(*(r.keys() for r in results)))
            w = csv.DictWriter(f, fieldnames=keys, delimiter="\t")
            w.writeheader()
            w.writerows(results)

    print(f"\nResults saved to {OUTDIR}/")
    return results


if __name__ == "__main__":
    run_benchmark()
