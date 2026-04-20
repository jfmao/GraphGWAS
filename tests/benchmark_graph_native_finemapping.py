"""
Benchmark: Graph-Native Fine-Mapping Methods (GRSD, HBP, CLGF)
vs L1, SuSiE, FINEMAP.

Three scenarios:
  1. Single-locus random causal (20 reps) — all methods on equal footing
  2. Single-locus functional causal (20 reps) — causal is eQTL variant
  3. Multi-locus pathway (5 reps) — CLGF vs per-locus methods

Usage:
    python tests/benchmark_graph_native_finemapping.py 2>/dev/null
"""

import csv, time, json, os, subprocess, tempfile, shutil
import numpy as np

import graphgwas.config as cfg
cfg._N_SAMPLES_DETECTED = False

from graphgwas.db import GraphGWASConnection
from graphgwas.config import detect_n_samples
from graphgwas.simulate import (
    simulate_F1_single_causal_in_ld, load_simulation_as_phenotype,
    simulate_F6_multi_locus_pathway,
)
from graphgwas.finemapping_v2 import (
    _load_locus_variants, _compute_association_stats, _compute_ld_matrix,
    _ld_deconvolve, _softmax, _build_credible_set, _load_eqtl_cache,
    _compute_ld_correlation, _compute_z_scores,
    graph_regularized_sparse_finemap,
    hierarchical_bp_finemap,
    cross_locus_graph_finemap,
    FinemapCandidate,
)
from graphgwas.genotype import get_all_indices, get_phenotype_values

OUTDIR = "/mnt/data/GraphGWAS/results/benchmark_v2/graph_native"
os.makedirs(OUTDIR, exist_ok=True)
FINEMAP_BIN = "/home/jfmao/bin/finemap"
EQTL_PATH = "/mnt/data/GraphGWAS/data/annotations/gtex_chr22_enhanced_cache.json"

WINDOW = 50000
CENTERS = [17000000, 19000000, 21000000, 23000000, 25000000,
           27000000, 29000000, 31000000, 33000000, 35000000,
           37000000, 39000000, 41000000, 43000000, 45000000,
           17500000, 20000000, 24000000, 28000000, 36000000]


# ---- Fast L1 (no Neo4j annotation queries) ----
def _fast_l1(variants, pheno, eqtl_cache, alpha=0.9):
    z_stats = _compute_association_stats(variants, pheno)
    R_sq = _compute_ld_matrix(variants)
    unique, n_nb = _ld_deconvolve(z_stats, R_sq, 0.3)
    z_func = np.zeros(len(variants))
    if eqtl_cache:
        for i, v in enumerate(variants):
            eq = eqtl_cache.get(v["variantId"])
            if eq:
                s = eq.get("composite", 0)
                nt = eq.get("n_tissues", 0)
                if 1 <= nt <= 3: s *= 1.5
                elif nt >= 8: s *= 0.3
                z_func[i] = np.log1p(s)
        if z_func.max() > 0:
            z_func = z_func / z_func.max() * z_stats.max()
    combined = alpha * unique + (1 - alpha) * z_func
    pip = _softmax(combined)
    return pip


# ---- FINEMAP wrapper ----
def _run_finemap(variants, pheno, tmpdir, n_samples, causal_vid):
    n_var = len(variants)
    D = np.array([v["dosage"] for v in variants], dtype=np.float64)
    D = np.where(np.isnan(D), np.nanmean(D, axis=1, keepdims=True), D)
    y = np.where(np.isnan(pheno), 0, pheno)
    valid = ~np.isnan(pheno)
    Dv, yv = D[:, valid], y[valid]
    n = int(valid.sum())
    yc = yv - yv.mean()
    Dc = Dv - Dv.mean(axis=1, keepdims=True)
    vg = np.var(Dv, axis=1)
    sv = np.where(vg > 1e-10, vg, 1.0)
    beta = (Dc @ yc) / n / sv
    pred = Dc * beta[:, np.newaxis]
    mse = np.sum((yc - pred)**2, axis=1) / (n - 2)
    se = np.sqrt(mse / (n * sv))
    sse = np.where(se > 1e-10, se, 1.0)
    keep = vg > 1e-10

    with open(f"{tmpdir}/d.z", "w") as f:
        f.write("rsid chromosome position allele1 allele2 maf beta se\n")
        for i in np.where(keep)[0]:
            v = variants[i]
            r = v["variantId"].replace(":", "_").replace("/", "_")
            m = min(v.get("af_total", .1), 1 - v.get("af_total", .1))
            f.write(f"{r} 22 {v['pos']} A G {max(m,.001):.6f} {beta[i]:.6f} {sse[i]:.6f}\n")

    stds = Dv.std(axis=1, keepdims=True)
    stds = np.where(stds < 1e-10, 1.0, stds)
    Ds = (Dv - Dv.mean(axis=1, keepdims=True)) / stds
    R = (Ds @ Ds.T) / Dv.shape[1]
    np.fill_diagonal(R, 1.0)
    Rk = R[np.ix_(keep, keep)]
    with open(f"{tmpdir}/d.ld", "w") as f:
        for i in range(Rk.shape[0]):
            f.write(" ".join(f"{Rk[i,j]:.6f}" for j in range(Rk.shape[1])) + "\n")
    with open(f"{tmpdir}/d.master", "w") as f:
        f.write("z;ld;snp;config;cred;log;n_samples\n")
        f.write(f"{tmpdir}/d.z;{tmpdir}/d.ld;{tmpdir}/d.snp;{tmpdir}/d.config;"
                f"{tmpdir}/d.cred;{tmpdir}/d.log;{n_samples}\n")

    rsid = causal_vid.replace(":", "_").replace("/", "_")
    try:
        subprocess.run([FINEMAP_BIN, "--sss", "--in-files", f"{tmpdir}/d.master",
                        "--n-causal-snps", "5"],
                       capture_output=True, text=True, timeout=120, cwd=tmpdir)
        pips = {}
        with open(f"{tmpdir}/d.snp") as f:
            f.readline()
            for line in f:
                p = line.strip().split()
                if len(p) >= 11: pips[p[1]] = float(p[10])
        if not pips: return -1, 0.0
        cp = pips.get(rsid, 0.0)
        rank = sum(1 for p in pips.values() if p >= cp) if cp > 0 else -1
        return rank, cp
    except Exception:
        return -1, 0.0


# ---- SuSiE wrapper ----
def _run_susie(variants, pheno, all_idx, causal_idx):
    nv = len(variants)
    X = np.zeros((len(all_idx), nv))
    for i, v in enumerate(variants):
        d = v["dosage"]
        X[:, i] = np.where(np.isnan(d), np.nanmean(d), d)
    y = np.where(np.isnan(pheno), np.nanmean(pheno), pheno)
    tmp = tempfile.mkdtemp()
    np.savetxt(f"{tmp}/X.csv", X, delimiter=",")
    np.savetxt(f"{tmp}/y.csv", y, delimiter=",")
    rs = f'''library(susieR)
X<-as.matrix(read.csv("{tmp}/X.csv",header=FALSE));y<-scan("{tmp}/y.csv")
fit<-tryCatch(susie(X,y,L=5,verbose=FALSE),error=function(e) NULL)
if(!is.null(fit)){{pips<-fit$pip;ci<-{causal_idx+1}
cat(sprintf("PIP=%.6f\\n",pips[ci]));cat(sprintf("RANK=%d\\n",sum(pips>=pips[ci])))
}}else cat("FAILED\\n")'''
    try:
        proc = subprocess.run(["R", "--no-save", "--no-restore", "-e", rs],
                               capture_output=True, text=True, timeout=120)
        pip = rank = 0
        for l in proc.stdout.split("\n"):
            if l.startswith("PIP="): pip = float(l.split("=")[1])
            elif l.startswith("RANK="): rank = int(l.split("=")[1])
        return rank, pip
    except Exception:
        return -1, 0.0
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def _eval_method(pip_array, variants, causal_vid):
    """Get rank and PIP for causal variant."""
    cidx = next((i for i, v in enumerate(variants) if v["variantId"] == causal_vid), -1)
    if cidx < 0:
        return -1, 0.0
    causal_pip = pip_array[cidx]
    rank = int(np.sum(pip_array >= causal_pip))
    return rank, float(causal_pip)


def run_scenario_1(conn, all_idx, eqtl_cache, n_reps=20):
    """Single-locus random causal: all methods compared."""
    print(f"\n{'='*90}")
    print("SCENARIO 1: Single-locus random causal (β=0.3, h²=0.05)")
    print("="*90)

    results = []
    for rep in range(n_reps):
        center = CENTERS[rep % len(CENTERS)]
        try:
            sim = simulate_F1_single_causal_in_ld(
                conn, "chr22", center, beta=0.3, h2_target=0.05, seed=500 + rep)
        except Exception as e:
            print(f"  Rep {rep}: sim failed ({e})")
            continue

        causal = sim.causal_variants[0]
        load_simulation_as_phenotype(conn, sim, f"gn_s1_{rep}")
        pheno = get_phenotype_values(conn, all_idx)

        variants = _load_locus_variants(conn, "chr22", causal.pos, WINDOW, all_idx)
        cidx = next((i for i, v in enumerate(variants) if v["variantId"] == causal.variant_id), -1)
        if cidx < 0 or len(variants) < 10:
            continue

        row = {"rep": rep, "nv": len(variants), "af": causal.af}

        # L1
        t0 = time.time()
        l1_pip = _fast_l1(variants, pheno, eqtl_cache)
        row["l1_time"] = time.time() - t0
        row["l1_rank"], row["l1_pip"] = _eval_method(l1_pip, variants, causal.variant_id)

        # GRSD
        t0 = time.time()
        grsd = graph_regularized_sparse_finemap(
            conn, "chr22", causal.pos, window=WINDOW, verbose=False)
        row["grsd_time"] = time.time() - t0
        row["grsd_rank"] = next((i+1 for i, c in enumerate(grsd) if c.variant_id == causal.variant_id), -1)
        row["grsd_pip"] = next((c.pip for c in grsd if c.variant_id == causal.variant_id), 0)

        # HBP
        t0 = time.time()
        hbp = hierarchical_bp_finemap(
            conn, "chr22", causal.pos, window=WINDOW,
            eqtl_cache=eqtl_cache, verbose=False)
        row["hbp_time"] = time.time() - t0
        row["hbp_rank"] = next((i+1 for i, c in enumerate(hbp) if c.variant_id == causal.variant_id), -1)
        row["hbp_pip"] = next((c.pip for c in hbp if c.variant_id == causal.variant_id), 0)

        # FINEMAP
        td = tempfile.mkdtemp()
        try:
            t0 = time.time()
            row["fm_rank"], row["fm_pip"] = _run_finemap(
                variants, pheno, td, cfg.N_SAMPLES, causal.variant_id)
            row["fm_time"] = time.time() - t0
        finally:
            shutil.rmtree(td, ignore_errors=True)

        # SuSiE
        t0 = time.time()
        row["su_rank"], row["su_pip"] = _run_susie(variants, pheno, all_idx, cidx)
        row["su_time"] = time.time() - t0

        results.append(row)
        print(f"  Rep {rep}: nv={row['nv']} | L1=#{row['l1_rank']} GRSD=#{row['grsd_rank']} "
              f"HBP=#{row['hbp_rank']} FM=#{row['fm_rank']} Su=#{row['su_rank']}", flush=True)

    return results


def run_scenario_2(conn, all_idx, eqtl_cache, n_reps=20):
    """Single-locus functional causal: causal is eQTL variant."""
    print(f"\n{'='*90}")
    print("SCENARIO 2: Single-locus functional causal (eQTL variant)")
    print("="*90)

    # Find loci rich in eQTL variants
    eqtl_positions = sorted([int(k.split(":")[1]) for k in eqtl_cache if k.startswith("chr22:")])
    windows = {}
    for pos in eqtl_positions:
        wk = (pos // WINDOW) * WINDOW
        windows.setdefault(wk, []).append(pos)
    good = [(k, len(v)) for k, v in windows.items() if len(v) >= 10]
    good.sort(key=lambda x: -x[1])
    centers = [w[0] + WINDOW // 2 for w in good[:n_reps]]

    results = []
    for rep, center in enumerate(centers):
        variants = _load_locus_variants(conn, "chr22", center, WINDOW, all_idx)
        eqtl_vars = [v for v in variants if v["variantId"] in eqtl_cache
                      and 0.02 < v.get("af_total", 0) < 0.98]
        if len(eqtl_vars) < 5 or len(variants) < 20:
            continue

        rng = np.random.default_rng(700 + rep)
        causal_v = rng.choice(eqtl_vars)
        dosage = causal_v["dosage"].copy()
        dosage = np.where(np.isnan(dosage), np.nanmean(dosage), dosage)
        af = causal_v.get("af_total", 0)
        if af < 0.02 or af > 0.98:
            continue

        beta = 0.3
        genetic = beta * dosage
        var_g = np.var(genetic)
        if var_g < 1e-10:
            continue
        var_e = var_g * (1 / 0.05 - 1)
        pheno = genetic + rng.normal(0, np.sqrt(var_e), len(all_idx))

        # Store phenotype
        conn.execute_write("""
            UNWIND $data AS row
            MATCH (s:Sample {packed_index: row.idx})
            SET s.gwas_value = row.val, s.is_case = null, s.is_control = null
        """, {"data": [{"idx": int(all_idx[i]), "val": float(pheno[i])} for i in range(len(all_idx))]})

        causal_vid = causal_v["variantId"]
        cidx = next((i for i, v in enumerate(variants) if v["variantId"] == causal_vid), -1)
        if cidx < 0:
            continue

        row = {"rep": rep, "nv": len(variants), "af": af,
               "n_eqtl": len(eqtl_vars)}

        # L1
        t0 = time.time()
        l1_pip = _fast_l1(variants, pheno, eqtl_cache)
        row["l1_time"] = time.time() - t0
        row["l1_rank"], row["l1_pip"] = _eval_method(l1_pip, variants, causal_vid)

        # GRSD
        t0 = time.time()
        grsd = graph_regularized_sparse_finemap(
            conn, "chr22", causal_v["pos"], window=WINDOW, verbose=False)
        row["grsd_time"] = time.time() - t0
        row["grsd_rank"] = next((i+1 for i, c in enumerate(grsd) if c.variant_id == causal_vid), -1)
        row["grsd_pip"] = next((c.pip for c in grsd if c.variant_id == causal_vid), 0)

        # HBP
        t0 = time.time()
        hbp = hierarchical_bp_finemap(
            conn, "chr22", causal_v["pos"], window=WINDOW,
            eqtl_cache=eqtl_cache, verbose=False)
        row["hbp_time"] = time.time() - t0
        row["hbp_rank"] = next((i+1 for i, c in enumerate(hbp) if c.variant_id == causal_vid), -1)
        row["hbp_pip"] = next((c.pip for c in hbp if c.variant_id == causal_vid), 0)

        # FINEMAP
        td = tempfile.mkdtemp()
        try:
            t0 = time.time()
            row["fm_rank"], row["fm_pip"] = _run_finemap(
                variants, pheno, td, cfg.N_SAMPLES, causal_vid)
            row["fm_time"] = time.time() - t0
        finally:
            shutil.rmtree(td, ignore_errors=True)

        # SuSiE
        t0 = time.time()
        row["su_rank"], row["su_pip"] = _run_susie(variants, pheno, all_idx, cidx)
        row["su_time"] = time.time() - t0

        results.append(row)
        print(f"  Rep {rep}: nv={row['nv']} eqtl={row['n_eqtl']} | "
              f"L1=#{row['l1_rank']} GRSD=#{row['grsd_rank']} "
              f"HBP=#{row['hbp_rank']} FM=#{row['fm_rank']} Su=#{row['su_rank']}", flush=True)

    return results


def run_scenario_3(conn, all_idx, eqtl_cache, n_reps=5):
    """Multi-locus pathway: CLGF vs per-locus L1."""
    print(f"\n{'='*90}")
    print("SCENARIO 3: Multi-locus pathway (CLGF vs per-locus)")
    print("="*90)

    results = []
    for rep in range(n_reps):
        try:
            sim = simulate_F6_multi_locus_pathway(
                conn, "chr22", n_pathway_loci=3, n_null_loci=3,
                beta=0.3, h2_target=0.10, seed=800 + rep)
        except Exception as e:
            print(f"  Rep {rep}: sim failed ({e})")
            continue

        load_simulation_as_phenotype(conn, sim, f"gn_s3_{rep}")
        pheno = get_phenotype_values(conn, all_idx)

        # Build loci list from causal variant positions
        loci = [{"chr": cv.chr, "lead_pos": cv.pos, "window": WINDOW}
                for cv in sim.causal_variants]

        # CLGF (cross-locus)
        t0 = time.time()
        clgf_results = cross_locus_graph_finemap(
            conn, loci, pheno, all_idx,
            eqtl_cache=eqtl_cache, verbose=False)
        clgf_time = time.time() - t0

        # Per-locus L1 for comparison
        for cv in sim.causal_variants:
            key = f"{cv.chr}:{cv.pos}"
            variants = _load_locus_variants(conn, cv.chr, cv.pos, WINDOW, all_idx)
            if len(variants) < 3:
                continue

            row = {"rep": rep, "causal_vid": cv.variant_id,
                   "in_pathway": cv.pathway is not None,
                   "pathway": cv.pathway or "none",
                   "nv": len(variants), "af": cv.af}

            # L1
            l1_pip = _fast_l1(variants, pheno, eqtl_cache)
            row["l1_rank"], row["l1_pip"] = _eval_method(l1_pip, variants, cv.variant_id)

            # CLGF
            if key in clgf_results:
                clgf_cands = clgf_results[key]
                row["clgf_rank"] = next((i+1 for i, c in enumerate(clgf_cands)
                                          if c.variant_id == cv.variant_id), -1)
                row["clgf_pip"] = next((c.pip for c in clgf_cands
                                         if c.variant_id == cv.variant_id), 0)
            else:
                row["clgf_rank"] = -1
                row["clgf_pip"] = 0

            row["clgf_time"] = clgf_time / len(loci)  # amortized

            results.append(row)
            pw_tag = "PW" if cv.pathway else "null"
            print(f"  Rep {rep} [{pw_tag}]: L1=#{row['l1_rank']} CLGF=#{row['clgf_rank']}", flush=True)

    return results


def print_summary(name, results, methods):
    """Print summary table for a scenario."""
    n = len(results)
    print(f"\n--- {name} ({n} data points) ---")
    for mname, pre in methods:
        ranks = [r[f"{pre}_rank"] for r in results if r.get(f"{pre}_rank", -1) > 0]
        pips = [r[f"{pre}_pip"] for r in results if r.get(f"{pre}_rank", -1) > 0]
        times = [r[f"{pre}_time"] for r in results if f"{pre}_time" in r]
        r1 = sum(1 for r in ranks if r == 1)
        if ranks:
            print(f"  {mname:15s}: rank={np.mean(ranks):5.1f} (#{1} {r1}/{n}={r1/n*100:.0f}%) "
                  f"PIP={np.mean(pips):.3f}"
                  + (f" time={np.mean(times):.2f}s" if times else ""))

    # Head-to-head
    method_list = [m for m in methods]
    for i, (an, ap) in enumerate(method_list):
        for bn, bp in method_list[i+1:]:
            w = l = t = 0
            for r in results:
                ar, br = r.get(f"{ap}_rank", -1), r.get(f"{bp}_rank", -1)
                if ar > 0 and br > 0:
                    if ar < br: w += 1
                    elif br < ar: l += 1
                    else: t += 1
            tot = w + l + t
            if tot > 0:
                print(f"    {an} vs {bn}: {w}:{l}:{t}")


def run_benchmark():
    print("Loading eQTL cache...", flush=True)
    eqtl_cache = _load_eqtl_cache(EQTL_PATH)
    print(f"  {len(eqtl_cache)} annotations loaded")

    with GraphGWASConnection("bolt://localhost:7688", "neo4j", "graphgwas") as conn:
        detect_n_samples(conn)
        all_idx = get_all_indices(conn)
        print(f"N_SAMPLES={cfg.N_SAMPLES}")

        s1 = run_scenario_1(conn, all_idx, eqtl_cache, n_reps=20)
        s2 = run_scenario_2(conn, all_idx, eqtl_cache, n_reps=20)
        s3 = run_scenario_3(conn, all_idx, eqtl_cache, n_reps=5)

    print(f"\n{'='*90}")
    print("GRAPH-NATIVE FINE-MAPPING BENCHMARK — RESULTS")
    print("="*90)

    methods_5 = [("L1", "l1"), ("GRSD", "grsd"), ("HBP", "hbp"),
                  ("FINEMAP", "fm"), ("SuSiE", "su")]
    methods_2 = [("L1", "l1"), ("CLGF", "clgf")]

    print_summary("Scenario 1: Random causal", s1, methods_5)
    print_summary("Scenario 2: Functional causal", s2, methods_5)
    print_summary("Scenario 3: Multi-locus (all)", s3, methods_2)

    # Scenario 3 stratified
    s3_pw = [r for r in s3 if r.get("in_pathway")]
    s3_null = [r for r in s3 if not r.get("in_pathway")]
    if s3_pw:
        print_summary("Scenario 3: Pathway loci only", s3_pw, methods_2)
    if s3_null:
        print_summary("Scenario 3: Null loci only", s3_null, methods_2)

    # Save
    all_data = {"scenario_1": s1, "scenario_2": s2, "scenario_3": s3}
    with open(f"{OUTDIR}/graph_native_benchmark.json", "w") as f:
        json.dump(all_data, f, indent=2, default=str)
    print(f"\nSaved to {OUTDIR}/graph_native_benchmark.json")


if __name__ == "__main__":
    run_benchmark()
