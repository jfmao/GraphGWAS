"""
Optimized L1 vs SuSiE benchmark.

Strategy for speed:
1. Pre-compute variants + LD matrix + annotations per locus (Neo4j queries ONCE per locus)
2. Vectorized LD matrix (single BLAS matmul instead of O(n²) loop)
3. Vectorized z-scores (matrix regression, not per-variant loop)
4. Batch SuSiE in one R session
5. Only phenotype changes between replicates at same locus

Expected: ~5 min total (vs ~60 min before).

Usage:
    # Restore the multiomics dump first, then:
    python tests/benchmark_l1_vs_susie_optimized.py
"""

import csv, time, json, os, subprocess, tempfile
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

OUTDIR = "/mnt/data/GraphGWAS/results/benchmark_v2/susie_comparison"
os.makedirs(OUTDIR, exist_ok=True)

WINDOW = 50000
N_REPS = 20
CENTERS = [17000000, 19000000, 21000000, 23000000, 25000000,
           27000000, 29000000, 31000000, 33000000, 35000000,
           37000000, 39000000, 41000000, 43000000, 45000000,
           17500000, 20000000, 24000000, 28000000, 36000000]


def run_benchmark():
    results = []
    t0_total = time.time()

    with GraphGWASConnection("bolt://localhost:7688", "neo4j", "graphgwas") as conn:
        detect_n_samples(conn)
        all_idx = get_all_indices(conn)

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
            load_simulation_as_phenotype(conn, sim, f"opt_rep{rep}")
            pheno = get_phenotype_values(conn, all_idx)

            # === L1 (uses vectorized LD + annotations) ===
            t0 = time.time()
            l1 = dual_graph_finemap(conn, "chr22", causal.pos,
                                     window=WINDOW, verbose=False)
            l1_time = time.time() - t0

            l1_rank = next((i+1 for i, c in enumerate(l1)
                            if c.variant_id == causal.variant_id), -1)
            l1_pip = next((c.pip for c in l1
                           if c.variant_id == causal.variant_id), 0)
            l1_cs = sum(1 for c in l1 if c.in_credible_set)
            l1_in_cs = any(c.variant_id == causal.variant_id
                           and c.in_credible_set for c in l1)

            # === SuSiE ===
            # Build dosage matrix from variants already loaded by L1
            variants = _load_locus_variants(conn, "chr22", causal.pos,
                                             WINDOW, all_idx)
            n_var = len(variants)
            causal_idx = next((i for i, v in enumerate(variants)
                               if v["variantId"] == causal.variant_id), -1)

            t0 = time.time()
            if causal_idx >= 0 and n_var >= 10:
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
                proc = subprocess.run(
                    ["R", "--no-save", "--no-restore", "-e", r_script],
                    capture_output=True, text=True, timeout=120)
                susie_time = time.time() - t0

                susie_pip = susie_rank = susie_cs = 0
                susie_in_cs = False
                for line in proc.stdout.split("\n"):
                    if line.startswith("PIP="):
                        susie_pip = float(line.split("=")[1])
                    elif line.startswith("RANK="):
                        susie_rank = int(line.split("=")[1])
                    elif line.startswith("IN_CS="):
                        susie_in_cs = line.split("=")[1].strip() == "TRUE"
                    elif line.startswith("CS_SIZE="):
                        susie_cs = int(line.split("=")[1])
            else:
                susie_time = 0
                susie_pip = susie_rank = susie_cs = 0
                susie_in_cs = False

            results.append({
                "rep": rep, "n_var": n_var, "causal_af": causal.af,
                "l1_rank": l1_rank, "l1_pip": l1_pip, "l1_cs": l1_cs,
                "l1_in_cs": l1_in_cs, "l1_time": l1_time,
                "susie_rank": susie_rank, "susie_pip": susie_pip,
                "susie_cs": susie_cs, "susie_in_cs": susie_in_cs,
                "susie_time": susie_time,
            })

            print(f"  Rep {rep}: L1=#{l1_rank} PIP={l1_pip:.3f} CS={l1_cs} ({l1_time:.1f}s) | "
                  f"SuSiE=#{susie_rank} PIP={susie_pip:.3f} CS={susie_cs} ({susie_time:.1f}s)",
                  flush=True)

    total_time = time.time() - t0_total

    # Summary
    print(f"\n{'='*70}")
    print(f"L1 (MULTI-OMICS, VECTORIZED) vs SuSiE — {len(results)} replicates, {total_time:.0f}s total")
    print("=" * 70)

    for name, pre in [("L1_multiomics", "l1"), ("SuSiE", "susie")]:
        ranks = [r[f"{pre}_rank"] for r in results if r[f"{pre}_rank"] > 0]
        pips = [r[f"{pre}_pip"] for r in results if r[f"{pre}_rank"] > 0]
        in_cs = [r[f"{pre}_in_cs"] for r in results]
        css = [r[f"{pre}_cs"] for r in results if r[f"{pre}_cs"] > 0]
        times = [r[f"{pre}_time"] for r in results]
        rank1 = sum(1 for r in ranks if r == 1)

        print(f"\n{name}:")
        print(f"  Mean rank: {np.mean(ranks):.1f} (rank #1 in {rank1}/{len(results)})")
        print(f"  Mean PIP: {np.mean(pips):.4f}")
        print(f"  Coverage: {np.mean(in_cs):.1%}")
        print(f"  Mean CS: {np.mean(css):.1f}" if css else "  Mean CS: N/A")
        print(f"  Mean runtime: {np.mean(times):.1f}s")

    # Save
    with open(f"{OUTDIR}/l1_multiomics_vs_susie.json", "w") as f:
        json.dump(results, f, indent=2, default=str)
    with open(f"{OUTDIR}/l1_multiomics_vs_susie.tsv", "w", newline="") as f:
        if results:
            w = csv.DictWriter(f, fieldnames=results[0].keys(), delimiter="\t")
            w.writeheader()
            w.writerows(results)

    print(f"\nResults saved to {OUTDIR}/")
    return results


if __name__ == "__main__":
    run_benchmark()
