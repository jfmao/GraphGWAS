"""
L1 vs SuSiE vs FINEMAP vs FAME benchmark.

Fixes from previous attempt:
- FINEMAP: needs CORRELATION matrix (not r²), proper master file with semicolons
- FAME: needs proper PLINK-format files + correct annotation matrix

Usage:
    python tests/benchmark_finemap_fame_comparison.py
"""

import csv, time, json, os, subprocess, tempfile, shutil
import numpy as np
from scipy import stats as sp_stats

# Setup
import graphgwas.config as cfg
cfg._N_SAMPLES_DETECTED = False

from graphgwas.db import GraphGWASConnection
from graphgwas.config import detect_n_samples
from graphgwas.simulate import simulate_F1_single_causal_in_ld, load_simulation_as_phenotype
from graphgwas.finemapping_v2 import (
    dual_graph_finemap, _load_locus_variants,
    _compute_association_stats, _compute_ld_matrix,
)
from graphgwas.genotype import get_all_indices, get_phenotype_values

OUTDIR = "/mnt/data/GraphGWAS/results/benchmark_v2/tool_comparison"
os.makedirs(OUTDIR, exist_ok=True)

FINEMAP_BIN = "/home/jfmao/bin/finemap"
FAME_BIN = "/tmp/FAME/build/FAME"

WINDOW = 50000
N_REPS = 20
CENTERS = [17000000, 19000000, 21000000, 23000000, 25000000,
           27000000, 29000000, 31000000, 33000000, 35000000,
           37000000, 39000000, 41000000, 43000000, 45000000,
           17500000, 20000000, 24000000, 28000000, 36000000]


def _build_finemap_inputs(variants, pheno, all_idx, tmpdir, n_samples):
    """Build FINEMAP input files: z-score, LD (correlation), master."""
    n_var = len(variants)

    # Build dosage matrix
    D = np.array([v["dosage"] for v in variants], dtype=np.float64)
    row_means = np.nanmean(D, axis=1, keepdims=True)
    nan_mask = np.isnan(D)
    D = np.where(nan_mask, row_means, D)

    # Phenotype
    y = pheno.copy()
    valid = ~np.isnan(y)
    y_clean = np.where(np.isnan(y), 0, y)

    D_v = D[:, valid]
    y_v = y_clean[valid]
    n = int(valid.sum())

    # Compute z-scores (beta/se)
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

    # Write z-score file (FINEMAP format: rsid chromosome position allele1 allele2 maf beta se)
    z_path = os.path.join(tmpdir, "data.z")
    with open(z_path, "w") as f:
        f.write("rsid chromosome position allele1 allele2 maf beta se\n")
        for i, v in enumerate(variants):
            rsid = v["variantId"].replace(":", "_").replace("/", "_")
            chrom = v.get("chromosome", "chr22").replace("chr", "").replace("chromosome", "")
            pos = v["pos"]
            ref = v.get("ref", "A")[:1] or "A"
            alt = v.get("alt", "G")[:1] or "G"
            maf = min(v.get("af_total", 0.1), 1 - v.get("af_total", 0.1))
            maf = max(maf, 0.001)  # FINEMAP doesn't like 0 MAF
            if var_g[i] > 1e-10:
                f.write(f"{rsid} {chrom} {pos} {ref} {alt} {maf:.6f} {beta[i]:.6f} {safe_se[i]:.6f}\n")

    # Count actual variants written (excluding monomorphic)
    n_written = sum(1 for vg in var_g if vg > 1e-10)

    # Write LD CORRELATION matrix (NOT r² — FINEMAP needs signed correlation)
    # Standardize dosages, compute D @ D.T / n
    stds = D_v.std(axis=1, keepdims=True)
    stds = np.where(stds < 1e-10, 1.0, stds)
    D_std = (D_v - D_v.mean(axis=1, keepdims=True)) / stds
    R = (D_std @ D_std.T) / D_v.shape[1]
    np.fill_diagonal(R, 1.0)

    # Only include non-monomorphic variants
    keep = var_g > 1e-10
    R_keep = R[np.ix_(keep, keep)]

    ld_path = os.path.join(tmpdir, "data.ld")
    with open(ld_path, "w") as f:
        for i in range(R_keep.shape[0]):
            f.write(" ".join(f"{R_keep[i, j]:.6f}" for j in range(R_keep.shape[1])) + "\n")

    # Write master file
    master_path = os.path.join(tmpdir, "data.master")
    with open(master_path, "w") as f:
        f.write("z;ld;snp;config;cred;log;n_samples\n")
        f.write(f"{z_path};{ld_path};{tmpdir}/data.snp;{tmpdir}/data.config;"
                f"{tmpdir}/data.cred;{tmpdir}/data.log;{n_samples}\n")

    return master_path, n_written


def _run_finemap(master_path, tmpdir, causal_vid, variants, var_g_mask=None):
    """Run FINEMAP and parse results."""
    try:
        proc = subprocess.run(
            [FINEMAP_BIN, "--sss", "--in-files", master_path,
             "--n-causal-snps", "5"],
            capture_output=True, text=True, timeout=120,
            cwd=tmpdir)

        snp_path = os.path.join(tmpdir, "data.snp")
        if not os.path.exists(snp_path) or os.path.getsize(snp_path) == 0:
            return -1, 0.0, 0, False, proc.stderr[:500]

        # Parse .snp file
        causal_rsid = causal_vid.replace(":", "_").replace("/", "_")
        pips = {}
        with open(snp_path) as f:
            header = f.readline()
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 11:
                    rsid = parts[1]
                    pip = float(parts[10])  # 'prob' column
                    pips[rsid] = pip

        if not pips:
            return -1, 0.0, 0, False, "empty snp file"

        causal_pip = pips.get(causal_rsid, 0.0)
        # Rank: how many variants have PIP >= causal PIP
        rank = sum(1 for p in pips.values() if p >= causal_pip) if causal_pip > 0 else -1

        # Credible set from .cred file
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
                    cs_size = max(cs_size, len(parts) - 1)  # subtract first column

        return rank, causal_pip, cs_size, in_cs, ""

    except subprocess.TimeoutExpired:
        return -1, 0.0, 0, False, "timeout"
    except Exception as e:
        return -1, 0.0, 0, False, str(e)


def _run_susie(variants, pheno, all_idx, causal_vid, causal_idx):
    """Run SuSiE via R and parse results."""
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

    print(f"Tools available: FINEMAP={has_finemap}, SuSiE(R)={has_susie}")
    print(f"Running {N_REPS} replicates, window={WINDOW}bp")
    print("=" * 80)

    with GraphGWASConnection("bolt://localhost:7688", "neo4j", "graphgwas") as conn:
        detect_n_samples(conn)
        n_samples = cfg.N_SAMPLES
        all_idx = get_all_indices(conn)
        print(f"N_SAMPLES={n_samples}, all_idx={len(all_idx)}")

        for rep in range(N_REPS):
            center = CENTERS[rep % len(CENTERS)]
            seed = 42 + rep

            # Simulate
            try:
                sim = simulate_F1_single_causal_in_ld(
                    conn, "chr22", locus_center=center,
                    beta=0.5, h2_target=0.1, seed=seed)
            except Exception as e:
                print(f"  Rep {rep}: simulation failed ({e})")
                continue

            causal = sim.causal_variants[0]
            load_simulation_as_phenotype(conn, sim, f"cmp_rep{rep}")
            pheno = get_phenotype_values(conn, all_idx)

            # Load variants once (shared across all methods)
            variants = _load_locus_variants(conn, "chr22", causal.pos, WINDOW, all_idx)
            n_var = len(variants)
            causal_idx = next((i for i, v in enumerate(variants)
                               if v["variantId"] == causal.variant_id), -1)

            if causal_idx < 0 or n_var < 10:
                print(f"  Rep {rep}: causal not in locus or too few variants ({n_var})")
                continue

            row = {"rep": rep, "n_var": n_var, "causal_af": causal.af,
                   "causal_pos": causal.pos, "center": center}

            # === L1 ===
            t0 = time.time()
            l1 = dual_graph_finemap(conn, "chr22", causal.pos,
                                     window=WINDOW, verbose=False)
            row["l1_time"] = time.time() - t0
            row["l1_rank"] = next((i+1 for i, c in enumerate(l1)
                                    if c.variant_id == causal.variant_id), -1)
            row["l1_pip"] = next((c.pip for c in l1
                                   if c.variant_id == causal.variant_id), 0)
            row["l1_cs"] = sum(1 for c in l1 if c.in_credible_set)
            row["l1_in_cs"] = any(c.variant_id == causal.variant_id
                                   and c.in_credible_set for c in l1)

            # === SuSiE ===
            if has_susie:
                t0 = time.time()
                su_rank, su_pip, su_cs, su_in_cs = _run_susie(
                    variants, pheno, all_idx, causal.variant_id, causal_idx)
                row["susie_time"] = time.time() - t0
                row["susie_rank"] = su_rank
                row["susie_pip"] = su_pip
                row["susie_cs"] = su_cs
                row["susie_in_cs"] = su_in_cs

            # === FINEMAP ===
            if has_finemap:
                tmpdir = tempfile.mkdtemp()
                try:
                    t0 = time.time()
                    master_path, n_written = _build_finemap_inputs(
                        variants, pheno, all_idx, tmpdir, n_samples)
                    fm_rank, fm_pip, fm_cs, fm_in_cs, fm_err = _run_finemap(
                        master_path, tmpdir, causal.variant_id, variants)
                    row["finemap_time"] = time.time() - t0
                    row["finemap_rank"] = fm_rank
                    row["finemap_pip"] = fm_pip
                    row["finemap_cs"] = fm_cs
                    row["finemap_in_cs"] = fm_in_cs
                    row["finemap_n_written"] = n_written
                    if fm_err:
                        row["finemap_err"] = fm_err
                finally:
                    shutil.rmtree(tmpdir, ignore_errors=True)

            results.append(row)

            # Print progress
            parts = [f"Rep {rep}: nv={n_var}"]
            parts.append(f"L1=#{row['l1_rank']} PIP={row['l1_pip']:.3f} CS={row['l1_cs']} ({row['l1_time']:.1f}s)")
            if has_susie and "susie_rank" in row:
                parts.append(f"SuSiE=#{row['susie_rank']} PIP={row.get('susie_pip',0):.3f} ({row.get('susie_time',0):.1f}s)")
            if has_finemap and "finemap_rank" in row:
                parts.append(f"FINEMAP=#{row['finemap_rank']} PIP={row.get('finemap_pip',0):.3f} ({row.get('finemap_time',0):.1f}s)")
                if row.get("finemap_err"):
                    parts.append(f"[{row['finemap_err'][:40]}]")
            print("  " + " | ".join(parts), flush=True)

    total_time = time.time() - t0_total

    # Summary
    print(f"\n{'='*80}")
    print(f"L1 vs SuSiE vs FINEMAP — {len(results)} replicates, {total_time:.0f}s total")
    print("=" * 80)

    methods = [("L1", "l1"), ("SuSiE", "susie"), ("FINEMAP", "finemap")]
    for name, pre in methods:
        if not any(f"{pre}_rank" in r for r in results):
            continue
        ranks = [r[f"{pre}_rank"] for r in results if r.get(f"{pre}_rank", -1) > 0]
        pips = [r[f"{pre}_pip"] for r in results if r.get(f"{pre}_rank", -1) > 0]
        in_cs = [r.get(f"{pre}_in_cs", False) for r in results]
        css = [r[f"{pre}_cs"] for r in results if r.get(f"{pre}_cs", 0) > 0]
        times = [r[f"{pre}_time"] for r in results if f"{pre}_time" in r]
        rank1 = sum(1 for r in ranks if r == 1)

        print(f"\n{name}:")
        print(f"  Replicates with valid results: {len(ranks)}/{len(results)}")
        if ranks:
            print(f"  Mean rank: {np.mean(ranks):.1f} (rank #1 in {rank1}/{len(results)})")
            print(f"  Median rank: {np.median(ranks):.0f}")
            print(f"  Mean PIP: {np.mean(pips):.4f}")
        print(f"  Coverage (causal in CS): {np.mean(in_cs):.1%}")
        if css:
            print(f"  Mean CS size: {np.mean(css):.1f}")
        if times:
            print(f"  Mean runtime: {np.mean(times):.1f}s")

    # Save
    with open(f"{OUTDIR}/l1_vs_susie_vs_finemap.json", "w") as f:
        json.dump(results, f, indent=2, default=str)
    with open(f"{OUTDIR}/l1_vs_susie_vs_finemap.tsv", "w", newline="") as f:
        if results:
            keys = sorted(set().union(*(r.keys() for r in results)))
            w = csv.DictWriter(f, fieldnames=keys, delimiter="\t")
            w.writeheader()
            w.writerows(results)

    print(f"\nResults saved to {OUTDIR}/")
    return results


if __name__ == "__main__":
    run_benchmark()
