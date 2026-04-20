"""
Head-to-head: Fast HBP vs SuSiE vs FINEMAP (50 replicates).

Pre-caches graph structure once, then runs fast HBP (~0.1s/locus).
Two scenarios: random causal + functional causal.

Usage:
    python tests/benchmark_hbp_vs_susie_finemap.py 2>/dev/null
"""

import time, json, os, subprocess, tempfile, shutil
import numpy as np

import graphgwas.config as cfg
cfg._N_SAMPLES_DETECTED = False

from graphgwas.db import GraphGWASConnection
from graphgwas.config import detect_n_samples
from graphgwas.simulate import simulate_F1_single_causal_in_ld, load_simulation_as_phenotype
from graphgwas.finemapping_v2 import (
    _load_locus_variants, _compute_association_stats, _compute_ld_matrix,
    _ld_deconvolve, _softmax, _load_eqtl_cache,
    _cache_chr_graph_structure, fast_hbp_finemap,
)
from graphgwas.genotype import get_all_indices, get_phenotype_values

OUTDIR = "/mnt/data/GraphGWAS/results/benchmark_v2/hbp_h2h"
os.makedirs(OUTDIR, exist_ok=True)
FINEMAP_BIN = "/home/jfmao/bin/finemap"
EQTL_PATH = "/mnt/data/GraphGWAS/data/annotations/gtex_chr22_enhanced_cache.json"
WINDOW = 50000

CENTERS = [17000000, 19000000, 21000000, 23000000, 25000000,
           27000000, 29000000, 31000000, 33000000, 35000000,
           37000000, 39000000, 41000000, 43000000, 45000000,
           17500000, 20000000, 24000000, 28000000, 36000000,
           18000000, 22000000, 26000000, 30000000, 34000000,
           38000000, 40000000, 42000000, 44000000, 46000000,
           16500000, 18500000, 20500000, 22500000, 24500000,
           26500000, 28500000, 30500000, 32500000, 34500000,
           36500000, 38500000, 40500000, 42500000, 44500000,
           17200000, 19200000, 21200000, 23200000, 25200000]


def _fast_l1(variants, pheno, eqtl_cache, alpha=0.9):
    z = _compute_association_stats(variants, pheno)
    R = _compute_ld_matrix(variants)
    u, n = _ld_deconvolve(z, R, 0.3)
    zf = np.zeros(len(variants))
    if eqtl_cache:
        for i, v in enumerate(variants):
            eq = eqtl_cache.get(v["variantId"])
            if eq:
                s = eq.get("composite", 0)
                nt = eq.get("n_tissues", 0)
                if 1 <= nt <= 3: s *= 1.5
                elif nt >= 8: s *= 0.3
                zf[i] = np.log1p(s)
        if zf.max() > 0:
            zf = zf / zf.max() * z.max()
    return _softmax(alpha * u + (1 - alpha) * zf)


def _run_finemap(variants, pheno, tmpdir, n_samples, causal_vid):
    nv = len(variants)
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
        return (sum(1 for p in pips.values() if p >= cp) if cp > 0 else -1), cp
    except Exception:
        return -1, 0.0


def _run_susie(variants, pheno, all_idx, cidx):
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
if(!is.null(fit)){{pips<-fit$pip;ci<-{cidx+1}
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


def _eval(pip_array, variants, causal_vid):
    cidx = next((i for i, v in enumerate(variants) if v["variantId"] == causal_vid), -1)
    if cidx < 0: return -1, 0.0
    return int(np.sum(pip_array >= pip_array[cidx])), float(pip_array[cidx])


def run_scenario(name, conn, all_idx, eqtl_cache, graph_cache, n_reps,
                  beta, h2, seed_base, functional=False):
    print(f"\n{'='*90}")
    print(f"{name} ({n_reps} reps, β={beta}, h²={h2})")
    print("="*90)

    results = []

    if functional:
        # Find eQTL-rich loci
        eqtl_pos = sorted([int(k.split(":")[1]) for k in eqtl_cache if k.startswith("chr22:")])
        windows = {}
        for pos in eqtl_pos:
            wk = (pos // WINDOW) * WINDOW
            windows.setdefault(wk, []).append(pos)
        good = [(k, len(v)) for k, v in windows.items() if len(v) >= 10]
        good.sort(key=lambda x: -x[1])
        locus_centers = [w[0] + WINDOW // 2 for w in good[:n_reps]]
    else:
        locus_centers = CENTERS[:n_reps]

    for rep in range(min(n_reps, len(locus_centers))):
        center = locus_centers[rep]
        seed = seed_base + rep

        if functional:
            # Load variants and pick eQTL causal
            variants = _load_locus_variants(conn, "chr22", center, WINDOW, all_idx)
            eqtl_vars = [v for v in variants if v["variantId"] in eqtl_cache
                          and 0.02 < v.get("af_total", 0) < 0.98]
            if len(eqtl_vars) < 5 or len(variants) < 20:
                continue
            rng = np.random.default_rng(seed)
            causal_v = rng.choice(eqtl_vars)
            dosage = causal_v["dosage"].copy()
            dosage = np.where(np.isnan(dosage), np.nanmean(dosage), dosage)
            if causal_v.get("af_total", 0) < 0.02:
                continue
            genetic = beta * dosage
            var_g = np.var(genetic)
            if var_g < 1e-10:
                continue
            var_e = var_g * (1 / h2 - 1)
            pheno = genetic + rng.normal(0, np.sqrt(var_e), len(all_idx))
            conn.execute_write("""
                UNWIND $data AS row
                MATCH (s:Sample {packed_index: row.idx})
                SET s.gwas_value = row.val, s.is_case = null, s.is_control = null
            """, {"data": [{"idx": int(all_idx[i]), "val": float(pheno[i])} for i in range(len(all_idx))]})
            causal_vid = causal_v["variantId"]
            causal_pos = causal_v["pos"]
            causal_af = causal_v.get("af_total", 0)
        else:
            try:
                sim = simulate_F1_single_causal_in_ld(
                    conn, "chr22", center, beta=beta, h2_target=h2, seed=seed)
            except Exception as e:
                print(f"  Rep {rep}: sim failed ({e})")
                continue
            causal_vid = sim.causal_variants[0].variant_id
            causal_pos = sim.causal_variants[0].pos
            causal_af = sim.causal_variants[0].af
            load_simulation_as_phenotype(conn, sim, f"h2h_{seed}")
            pheno = get_phenotype_values(conn, all_idx)
            variants = _load_locus_variants(conn, "chr22", causal_pos, WINDOW, all_idx)

        cidx = next((i for i, v in enumerate(variants) if v["variantId"] == causal_vid), -1)
        if cidx < 0 or len(variants) < 10:
            continue

        nv = len(variants)
        row = {"rep": rep, "nv": nv, "af": causal_af}

        # L1
        t0 = time.time()
        l1_pip = _fast_l1(variants, pheno, eqtl_cache)
        row["l1_time"] = time.time() - t0
        row["l1_rank"], row["l1_pip"] = _eval(l1_pip, variants, causal_vid)

        # Fast HBP
        t0 = time.time()
        hbp = fast_hbp_finemap(variants, pheno, graph_cache,
                                eqtl_cache=eqtl_cache)
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
        print(f"  Rep {rep}: nv={nv} AF={causal_af:.2f} | "
              f"L1=#{row['l1_rank']} HBP=#{row['hbp_rank']}({row['hbp_time']:.2f}s) "
              f"FM=#{row['fm_rank']} Su=#{row['su_rank']}", flush=True)

    return results


def print_results(name, results):
    n = len(results)
    print(f"\n--- {name} ({n} reps) ---")
    for mname, pre in [("L1", "l1"), ("HBP", "hbp"), ("FINEMAP", "fm"), ("SuSiE", "su")]:
        ranks = [r[f"{pre}_rank"] for r in results if r.get(f"{pre}_rank", -1) > 0]
        pips = [r[f"{pre}_pip"] for r in results if r.get(f"{pre}_rank", -1) > 0]
        times = [r[f"{pre}_time"] for r in results if f"{pre}_time" in r]
        r1 = sum(1 for r in ranks if r == 1)
        if ranks:
            print(f"  {mname:10s}: rank={np.mean(ranks):5.1f} (#{1} {r1}/{n}={r1/n*100:.0f}%) "
                  f"PIP={np.mean(pips):.3f} time={np.mean(times):.2f}s")

    # Head-to-head
    pairs = [("HBP", "hbp", "SuSiE", "su"), ("HBP", "hbp", "FINEMAP", "fm"),
             ("HBP", "hbp", "L1", "l1"), ("L1", "l1", "SuSiE", "su")]
    for an, ap, bn, bp in pairs:
        w = l = t = 0
        for r in results:
            ar, br = r.get(f"{ap}_rank", -1), r.get(f"{bp}_rank", -1)
            if ar > 0 and br > 0:
                if ar < br: w += 1
                elif br < ar: l += 1
                else: t += 1
        tot = w + l + t
        if tot > 0:
            print(f"    {an} vs {bn}: {w}:{l}:{t} (win:loss:tie)")


def main():
    print("Loading caches...", flush=True)
    eqtl_cache = _load_eqtl_cache(EQTL_PATH)
    print(f"  eQTL: {len(eqtl_cache)} annotations")

    with GraphGWASConnection("bolt://localhost:7688", "neo4j", "graphgwas") as conn:
        detect_n_samples(conn)
        all_idx = get_all_indices(conn)
        print(f"  N_SAMPLES={cfg.N_SAMPLES}")

        print("Caching chr22 graph structure (variant→gene→pathway)...", flush=True)
        t0 = time.time()
        graph_cache = _cache_chr_graph_structure(conn, "chr22")
        print(f"  {len(graph_cache)} variants with annotations ({time.time()-t0:.1f}s)")

        # Scenario A: Random causal, strong signal
        s_strong = run_scenario(
            "RANDOM CAUSAL, STRONG (β=0.5, h²=0.10)",
            conn, all_idx, eqtl_cache, graph_cache,
            n_reps=30, beta=0.5, h2=0.10, seed_base=1000)

        # Scenario B: Random causal, weak signal
        s_weak = run_scenario(
            "RANDOM CAUSAL, WEAK (β=0.2, h²=0.02)",
            conn, all_idx, eqtl_cache, graph_cache,
            n_reps=30, beta=0.2, h2=0.02, seed_base=2000)

        # Scenario C: Functional causal (eQTL variant), weak signal
        s_func = run_scenario(
            "FUNCTIONAL CAUSAL (eQTL), WEAK (β=0.2, h²=0.02)",
            conn, all_idx, eqtl_cache, graph_cache,
            n_reps=30, beta=0.2, h2=0.02, seed_base=3000,
            functional=True)

    print(f"\n{'='*90}")
    print("HBP vs SuSiE vs FINEMAP — HEAD-TO-HEAD RESULTS")
    print("="*90)

    print_results("Strong signal (random causal)", s_strong)
    print_results("Weak signal (random causal)", s_weak)
    print_results("Weak signal (functional causal)", s_func)

    all_data = {"strong": s_strong, "weak": s_weak, "functional": s_func}
    with open(f"{OUTDIR}/hbp_h2h_50rep.json", "w") as f:
        json.dump(all_data, f, indent=2, default=str)
    print(f"\nSaved to {OUTDIR}/")


if __name__ == "__main__":
    main()
