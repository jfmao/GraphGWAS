"""
Null simulation FPR benchmark.

Generates pure-noise phenotypes (no causal variants) and runs all fine-mapping
methods on random loci. Measures:
  - Max PIP across the locus (should stay low under null)
  - Credible set size (should be LARGE under null — no concentration)
  - FPR at PIP > 0.5 and PIP > 0.9 thresholds

A well-calibrated method has:
  - Mean max PIP ≈ 1/n_variants (uniform)
  - Credible set ≈ 95% of variants (no concentration)
  - FPR ≤ 5% at PIP > 0.5

Usage:
    python tests/benchmark_null_simulations.py 2>/dev/null
"""

import csv, time, json, os, subprocess, tempfile, shutil
import numpy as np

import graphgwas.config as cfg
cfg._N_SAMPLES_DETECTED = False

from graphgwas.db import GraphGWASConnection
from graphgwas.config import detect_n_samples
from graphgwas.finemapping_v2 import (
    _load_locus_variants, _compute_association_stats, _compute_ld_matrix,
    _ld_deconvolve, _softmax, _build_credible_set, _load_eqtl_cache,
    _cache_chr_graph_structure, fast_hbp_finemap,
)
from graphgwas.genotype import get_all_indices, get_phenotype_values

OUTDIR = "/mnt/data/GraphGWAS/results/benchmark_v2/null_calibration"
os.makedirs(OUTDIR, exist_ok=True)
FINEMAP_BIN = "/home/jfmao/bin/finemap"

WINDOW = 50000
N_REPS = 100  # 100 null simulations

# Random loci across yeast chr1-16
np.random.seed(42)
LOCI = []
chr_lens = {  # approximate yeast chromosome lengths in bp
    1: 230218, 2: 813184, 3: 316620, 4: 1531933, 5: 576874,
    6: 270161, 7: 1090940, 8: 562643, 9: 439888, 10: 745751,
    11: 666816, 12: 1078177, 13: 924431, 14: 784333,
    15: 1091291, 16: 948066,
}
rng = np.random.default_rng(42)
chrs = list(chr_lens.keys())
for _ in range(N_REPS):
    c = rng.choice(chrs)
    pos = int(rng.integers(WINDOW + 1, chr_lens[c] - WINDOW))
    LOCI.append((f"chromosome{c}", pos))


def _fast_l1(variants, pheno, alpha=0.9):
    z = _compute_association_stats(variants, pheno)
    R = _compute_ld_matrix(variants)
    u, n = _ld_deconvolve(z, R, 0.3)
    return _softmax(alpha * u)


def _run_finemap(variants, pheno, tmpdir, n_samples):
    """Run FINEMAP and return max PIP, credible set size."""
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
    if keep.sum() < 5:
        return -1, -1, 0

    with open(f"{tmpdir}/d.z", "w") as f:
        f.write("rsid chromosome position allele1 allele2 maf beta se\n")
        for i in np.where(keep)[0]:
            v = variants[i]
            r = v["variantId"].replace(":", "_").replace("/", "_")
            m = min(v.get("af_total", .1), 1 - v.get("af_total", .1))
            f.write(f"{r} 1 {v['pos']} A G {max(m,.001):.6f} {beta[i]:.6f} {sse[i]:.6f}\n")

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

    try:
        subprocess.run([FINEMAP_BIN, "--sss", "--in-files", f"{tmpdir}/d.master",
                        "--n-causal-snps", "5"],
                       capture_output=True, text=True, timeout=120, cwd=tmpdir)
        pips = []
        with open(f"{tmpdir}/d.snp") as f:
            f.readline()
            for line in f:
                p = line.strip().split()
                if len(p) >= 11:
                    pips.append(float(p[10]))
        if not pips:
            return -1, -1, 0
        pips = np.array(pips)
        # Credible set from cumulative PIP
        sp = np.sort(pips)[::-1]
        cs = (np.cumsum(sp) < 0.95).sum() + 1
        return float(pips.max()), float(np.median(pips)), int(cs)
    except Exception:
        return -1, -1, 0


def _run_susie(variants, pheno, all_idx):
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
if(!is.null(fit)){{pips<-fit$pip
cat(sprintf("MAXPIP=%.6f\\n",max(pips)))
cat(sprintf("MEDPIP=%.6f\\n",median(pips)))
cs<-susie_get_cs(fit)
total_cs<-0
for(i in seq_along(cs$cs)) total_cs<-total_cs+length(cs$cs[[i]])
cat(sprintf("CS_SIZE=%d\\n",total_cs))
}}else cat("FAILED\\n")'''
    try:
        proc = subprocess.run(["R", "--no-save", "--no-restore", "-e", rs],
                               capture_output=True, text=True, timeout=120)
        max_pip = med_pip = -1
        cs_size = 0
        for l in proc.stdout.split("\n"):
            if l.startswith("MAXPIP="): max_pip = float(l.split("=")[1])
            elif l.startswith("MEDPIP="): med_pip = float(l.split("=")[1])
            elif l.startswith("CS_SIZE="): cs_size = int(l.split("=")[1])
        return max_pip, med_pip, cs_size
    except Exception:
        return -1, -1, 0
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def run_benchmark():
    print("=" * 90)
    print("NULL SIMULATION BENCHMARK — Type I Error Rate Calibration")
    print("=" * 90)
    print(f"{N_REPS} null simulations on yeast chr1-16, window={WINDOW}bp")
    print("Phenotype = pure noise (no causal variant)")
    print()

    with GraphGWASConnection("bolt://localhost:7688", "neo4j", "graphgwas") as conn:
        detect_n_samples(conn)
        all_idx = get_all_indices(conn)
        n_samples = cfg.N_SAMPLES
        print(f"N_SAMPLES={n_samples}")

        # Cache graph structure for all chromosomes
        print("\nCaching graph structure...", flush=True)
        graph_caches = {}
        for i in range(1, 17):
            chr_name = f"chromosome{i}"
            graph_caches[chr_name] = _cache_chr_graph_structure(conn, chr_name)
        print(f"  {sum(len(v) for v in graph_caches.values())} annotated variants")

        results = []
        rng = np.random.default_rng(42)

        for rep, (chr_name, lead_pos) in enumerate(LOCI):
            # Generate pure-noise phenotype
            pheno = rng.normal(0, 1, len(all_idx))

            # Store as gwas_value
            conn.execute_write("""
                UNWIND $data AS row
                MATCH (s:Sample {packed_index: row.idx})
                SET s.gwas_value = row.val, s.is_case = null, s.is_control = null
            """, {"data": [{"idx": int(all_idx[i]), "val": float(pheno[i])}
                            for i in range(len(all_idx))]})

            # Load locus variants
            variants = _load_locus_variants(conn, chr_name, lead_pos, WINDOW, all_idx)
            n_var = len(variants)
            if n_var < 10:
                continue

            row = {"rep": rep, "chr": chr_name, "lead_pos": lead_pos, "n_var": n_var}

            # L1
            t0 = time.time()
            l1_pip = _fast_l1(variants, pheno)
            row["l1_time"] = time.time() - t0
            row["l1_max_pip"] = float(l1_pip.max())
            row["l1_med_pip"] = float(np.median(l1_pip))
            l1_cs = _build_credible_set(l1_pip, 0.95)
            row["l1_cs_size"] = int(l1_cs.sum())

            # HBP
            t0 = time.time()
            hbp = fast_hbp_finemap(variants, pheno, graph_caches.get(chr_name, {}),
                                    chr_name=chr_name)
            row["hbp_time"] = time.time() - t0
            hbp_pips = np.array([c.pip for c in hbp])
            row["hbp_max_pip"] = float(hbp_pips.max())
            row["hbp_med_pip"] = float(np.median(hbp_pips))
            row["hbp_cs_size"] = sum(1 for c in hbp if c.in_credible_set)

            # FINEMAP
            td = tempfile.mkdtemp()
            try:
                t0 = time.time()
                fm_max, fm_med, fm_cs = _run_finemap(variants, pheno, td, n_samples)
                row["fm_time"] = time.time() - t0
                row["fm_max_pip"] = fm_max
                row["fm_med_pip"] = fm_med
                row["fm_cs_size"] = fm_cs
            finally:
                shutil.rmtree(td, ignore_errors=True)

            # SuSiE
            t0 = time.time()
            su_max, su_med, su_cs = _run_susie(variants, pheno, all_idx)
            row["su_time"] = time.time() - t0
            row["su_max_pip"] = su_max
            row["su_med_pip"] = su_med
            row["su_cs_size"] = su_cs

            results.append(row)

            if (rep + 1) % 10 == 0:
                print(f"  Rep {rep+1}/{N_REPS}: nv={n_var} | "
                      f"L1={row['l1_max_pip']:.3f} HBP={row['hbp_max_pip']:.3f} "
                      f"FM={row['fm_max_pip']:.3f} Su={row['su_max_pip']:.3f}",
                      flush=True)

    # Analysis
    print(f"\n{'='*90}")
    print("CALIBRATION RESULTS")
    print("="*90)

    methods = [("L1", "l1"), ("HBP", "hbp"), ("FINEMAP", "fm"), ("SuSiE", "su")]
    summary = {}

    print(f"\n{'Method':<10s} {'Mean MaxPIP':>12s} {'P95 MaxPIP':>12s} "
          f"{'FPR>0.5':>10s} {'FPR>0.9':>10s} {'Mean CS':>10s} {'CS Frac':>10s}")
    print("-" * 80)

    for name, pre in methods:
        max_pips = [r[f"{pre}_max_pip"] for r in results
                     if r.get(f"{pre}_max_pip", -1) > 0]
        cs_sizes = [r[f"{pre}_cs_size"] for r in results
                     if r.get(f"{pre}_cs_size", -1) > 0]
        n_var_arr = [r["n_var"] for r in results
                      if r.get(f"{pre}_max_pip", -1) > 0]
        if not max_pips:
            continue
        mean_pip = np.mean(max_pips)
        p95_pip = np.percentile(max_pips, 95)
        fpr_05 = np.mean([1 if p > 0.5 else 0 for p in max_pips])
        fpr_09 = np.mean([1 if p > 0.9 else 0 for p in max_pips])
        mean_cs = np.mean(cs_sizes) if cs_sizes else 0
        cs_frac = np.mean([cs_sizes[i] / n_var_arr[i] for i in range(len(cs_sizes))]) if cs_sizes else 0

        print(f"{name:<10s} {mean_pip:>12.4f} {p95_pip:>12.4f} "
              f"{fpr_05:>10.1%} {fpr_09:>10.1%} {mean_cs:>10.1f} {cs_frac:>10.1%}")

        summary[name] = {
            "n_reps": len(max_pips),
            "mean_max_pip": float(mean_pip),
            "p95_max_pip": float(p95_pip),
            "fpr_at_0.5": float(fpr_05),
            "fpr_at_0.9": float(fpr_09),
            "mean_cs_size": float(mean_cs),
            "mean_cs_fraction": float(cs_frac),
        }

    print(f"\nInterpretation:")
    print(f"  - Mean MaxPIP should be ~1/n_var (uniform under null)")
    print(f"  - FPR>0.5 should be ≤ 5% (well-calibrated null)")
    print(f"  - FPR>0.9 should be ~0% (high-confidence calls reliable)")
    print(f"  - Mean CS Frac should be near 95% (full uncertainty)")

    # Save
    with open(f"{OUTDIR}/null_calibration.json", "w") as f:
        json.dump({"summary": summary, "results": results}, f, indent=2, default=str)
    print(f"\nSaved to {OUTDIR}/")


if __name__ == "__main__":
    run_benchmark()
