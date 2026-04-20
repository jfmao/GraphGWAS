"""
PIP Calibration Benchmark.

Generates 200 simulations with known causal variants at varying h² (0.02-0.20).
For each method, bins variants by PIP and computes TDR (True Discovery Rate)
per bin. A well-calibrated method has TDR ≈ midpoint(PIP bin).

Standard Figure 2 in every fine-mapping paper.

Usage:
    python tests/benchmark_pip_calibration.py 2>/dev/null
"""

import csv, time, json, os, subprocess, tempfile, shutil
import numpy as np

import graphgwas.config as cfg
cfg._N_SAMPLES_DETECTED = False

from graphgwas.db import GraphGWASConnection
from graphgwas.config import detect_n_samples
from graphgwas.simulate import simulate_F1_single_causal_in_ld, load_simulation_as_phenotype
from graphgwas.finemapping_v2 import (
    _load_locus_variants, _compute_association_stats, _compute_ld_matrix,
    _ld_deconvolve, _softmax, _build_credible_set,
    _cache_chr_graph_structure, fast_hbp_finemap,
)
from graphgwas.genotype import get_all_indices, get_phenotype_values

OUTDIR = "/mnt/data/GraphGWAS/results/benchmark_v2/pip_calibration"
os.makedirs(OUTDIR, exist_ok=True)
FINEMAP_BIN = "/home/jfmao/bin/finemap"

WINDOW = 50000
N_REPS_PER_H2 = 50
H2_LEVELS = [0.02, 0.05, 0.10, 0.20]  # 4 levels × 50 reps = 200 sims

# Yeast loci spread across chr1-16
np.random.seed(42)
chr_lens = {1:230218, 2:813184, 3:316620, 4:1531933, 5:576874,
            6:270161, 7:1090940, 8:562643, 9:439888, 10:745751,
            11:666816, 12:1078177, 13:924431, 14:784333,
            15:1091291, 16:948066}

rng = np.random.default_rng(42)
LOCI = []
for _ in range(N_REPS_PER_H2):
    c = rng.choice(list(chr_lens.keys()))
    pos = int(rng.integers(WINDOW + 1, chr_lens[c] - WINDOW))
    LOCI.append((f"chromosome{c}", pos))


def _fast_l1(variants, pheno, alpha=0.9):
    z = _compute_association_stats(variants, pheno)
    R = _compute_ld_matrix(variants)
    u, n = _ld_deconvolve(z, R, 0.3)
    return _softmax(alpha * u)


def _run_finemap_pips(variants, pheno, tmpdir, n_samples):
    """Run FINEMAP and return PIPs aligned to variants list."""
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
        return None

    # Map variant ID -> index in original list
    rsid_to_idx = {}
    with open(f"{tmpdir}/d.z", "w") as f:
        f.write("rsid chromosome position allele1 allele2 maf beta se\n")
        for i in np.where(keep)[0]:
            v = variants[i]
            r = v["variantId"].replace(":", "_").replace("/", "_")
            rsid_to_idx[r] = i
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
        pips = np.zeros(nv)
        with open(f"{tmpdir}/d.snp") as f:
            f.readline()
            for line in f:
                p = line.strip().split()
                if len(p) >= 11:
                    rsid = p[1]
                    pip = float(p[10])
                    if rsid in rsid_to_idx:
                        pips[rsid_to_idx[rsid]] = pip
        return pips
    except Exception:
        return None


def _run_susie_pips(variants, pheno, all_idx):
    nv = len(variants)
    X = np.zeros((len(all_idx), nv))
    for i, v in enumerate(variants):
        d = v["dosage"]
        X[:, i] = np.where(np.isnan(d), np.nanmean(d), d)
    y = np.where(np.isnan(pheno), np.nanmean(pheno), pheno)
    tmp = tempfile.mkdtemp()
    np.savetxt(f"{tmp}/X.csv", X, delimiter=",")
    np.savetxt(f"{tmp}/y.csv", y, delimiter=",")
    np.savetxt(f"{tmp}/pip_out.csv", np.zeros(nv), delimiter=",")  # init
    rs = f'''library(susieR)
X<-as.matrix(read.csv("{tmp}/X.csv",header=FALSE));y<-scan("{tmp}/y.csv")
fit<-tryCatch(susie(X,y,L=5,verbose=FALSE),error=function(e) NULL)
if(!is.null(fit)) write.table(fit$pip,file="{tmp}/pip_out.csv",row.names=FALSE,col.names=FALSE)
'''
    try:
        subprocess.run(["R", "--no-save", "--no-restore", "-e", rs],
                       capture_output=True, text=True, timeout=120)
        pips = np.loadtxt(f"{tmp}/pip_out.csv")
        return pips if pips.shape == (nv,) else None
    except Exception:
        return None
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def run_calibration():
    print("=" * 90)
    print("PIP CALIBRATION BENCHMARK — TDR vs PIP")
    print("=" * 90)
    print(f"{N_REPS_PER_H2} reps × {len(H2_LEVELS)} h² levels = {N_REPS_PER_H2*len(H2_LEVELS)} simulations")
    print(f"h² levels: {H2_LEVELS}")
    print()

    with GraphGWASConnection("bolt://localhost:7688", "neo4j", "graphgwas") as conn:
        detect_n_samples(conn)
        all_idx = get_all_indices(conn)
        n_samples = cfg.N_SAMPLES
        print(f"N_SAMPLES={n_samples}")

        print("Caching graph structure...", flush=True)
        graph_caches = {}
        for i in range(1, 17):
            graph_caches[f"chromosome{i}"] = _cache_chr_graph_structure(conn, f"chromosome{i}")
        print(f"  {sum(len(v) for v in graph_caches.values())} annotated variants\n")

        # Collect (pip, is_causal) pairs across all simulations
        # Format: {method: [(pip, is_causal), ...]}
        collected = {"l1": [], "hbp": [], "fm": [], "su": []}

        for h2 in H2_LEVELS:
            print(f"\n{'='*90}")
            print(f"h² = {h2}, β=0.3")
            print("="*90)

            for rep, (chr_name, lead_pos) in enumerate(LOCI):
                seed = int(h2 * 1000) + rep + 9000
                try:
                    sim = simulate_F1_single_causal_in_ld(
                        conn, chr_name, lead_pos, beta=0.3, h2_target=h2, seed=seed)
                except Exception:
                    continue

                causal_vid = sim.causal_variants[0].variant_id
                load_simulation_as_phenotype(conn, sim, f"calib_{seed}")
                pheno = get_phenotype_values(conn, all_idx)

                variants = _load_locus_variants(conn, chr_name, sim.causal_variants[0].pos,
                                                  WINDOW, all_idx)
                if len(variants) < 10:
                    continue

                cidx = next((i for i, v in enumerate(variants)
                              if v["variantId"] == causal_vid), -1)
                if cidx < 0:
                    continue

                # is_causal vector
                is_causal = np.zeros(len(variants), dtype=bool)
                is_causal[cidx] = True

                # L1
                l1_pip = _fast_l1(variants, pheno)
                for i in range(len(variants)):
                    collected["l1"].append((float(l1_pip[i]), bool(is_causal[i])))

                # HBP
                hbp = fast_hbp_finemap(variants, pheno, graph_caches.get(chr_name, {}),
                                        chr_name=chr_name)
                for c in hbp:
                    is_c = c.variant_id == causal_vid
                    collected["hbp"].append((float(c.pip), bool(is_c)))

                # FINEMAP
                td = tempfile.mkdtemp()
                try:
                    fm_pips = _run_finemap_pips(variants, pheno, td, n_samples)
                    if fm_pips is not None:
                        for i in range(len(variants)):
                            collected["fm"].append((float(fm_pips[i]), bool(is_causal[i])))
                finally:
                    shutil.rmtree(td, ignore_errors=True)

                # SuSiE
                su_pips = _run_susie_pips(variants, pheno, all_idx)
                if su_pips is not None:
                    for i in range(len(variants)):
                        collected["su"].append((float(su_pips[i]), bool(is_causal[i])))

                if (rep + 1) % 10 == 0:
                    print(f"  Rep {rep+1}/{N_REPS_PER_H2} done", flush=True)

    # Compute calibration in PIP bins
    print(f"\n{'='*90}")
    print("PIP CALIBRATION TABLES")
    print("="*90)

    bins = [(0.0, 0.05), (0.05, 0.1), (0.1, 0.2), (0.2, 0.3), (0.3, 0.5),
            (0.5, 0.7), (0.7, 0.9), (0.9, 1.01)]

    summary = {}
    for name, pre in [("L1", "l1"), ("HBP", "hbp"), ("FINEMAP", "fm"), ("SuSiE", "su")]:
        data = collected[pre]
        if not data:
            continue
        print(f"\n--- {name} ---")
        print(f"  {'PIP bin':<12s} {'N variants':>10s} {'N causal':>10s} {'TDR':>8s}")
        print(f"  {'-'*46}")
        bin_results = []
        for lo, hi in bins:
            in_bin = [(p, c) for p, c in data if lo <= p < hi]
            n = len(in_bin)
            n_caus = sum(1 for _, c in in_bin if c)
            tdr = n_caus / n if n > 0 else 0
            print(f"  {lo:.2f}-{hi:.2f}    {n:>10d} {n_caus:>10d} {tdr:>8.4f}")
            bin_results.append({
                "bin_lo": lo, "bin_hi": hi, "n": n, "n_causal": n_caus,
                "tdr": tdr, "expected": (lo + hi) / 2,
            })
        summary[name] = bin_results

    # Save
    with open(f"{OUTDIR}/pip_calibration.json", "w") as f:
        json.dump({"summary": summary, "n_collected": {k: len(v) for k, v in collected.items()}},
                  f, indent=2, default=str)
    print(f"\nSaved to {OUTDIR}/pip_calibration.json")


if __name__ == "__main__":
    run_calibration()
