"""Validate population structure correction on yeast data.

Compares lambda GC and hit counts across four methods:
  1. None (naive regression)
  2. PCA covariates (Level 1)
  3. GRAMMAR residualization (Level 2)
  4. Graph-native kinship (Level 3)

Expected: lambda GC drops from ~3-8 (none) to ~1.0-1.1 (corrected).

Usage: python tests/yeast_popstruct_validation.py
"""

import glob
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src", "python"))

DATA_DIR = os.path.join(os.path.dirname(__file__), "data", "yeast")


def main():
    t0 = time.time()
    print("=" * 70)
    print("  Population Structure Correction Validation (Yeast)")
    print("=" * 70)

    from graphgwas.popstruct import (compute_grm, compute_pca,
                                      load_distance_matrix,
                                      corrected_gwas, lambda_gc)
    from graphgwas.sv import fast_gwas

    # --- Load genotypes ---
    print("\n--- Loading data ---")
    from pandas_plink import read_plink
    bim, fam, geno_dask = read_plink(os.path.join(DATA_DIR, "1011GWAS_matrix"))
    geno = geno_dask.compute().astype(np.float64)
    sample_ids = list(fam["iid"].values)

    af = np.nanmean(geno, axis=1) / 2
    maf = np.minimum(af, 1 - af)
    geno = geno[maf >= 0.01]
    print(f"  Genotypes: {geno.shape[0]:,} × {geno.shape[1]}")

    # --- Load phenotypes ---
    growth_dir = os.path.join(DATA_DIR, "Phenotypes_8391Traits", "GrowthTraits")
    phenotypes = {}
    for f in sorted(glob.glob(os.path.join(growth_dir, "*.phen"))):
        name = os.path.basename(f).replace(".phen", "")
        data = {}
        with open(f) as fh:
            fh.readline()
            for line in fh:
                parts = line.strip().split("\t")
                if len(parts) >= 2:
                    try:
                        data[parts[0]] = float(parts[1])
                    except ValueError:
                        pass
        if len(data) > 100:
            phenotypes[name] = data

    # Select test traits
    test_traits = ["YPETHANOL_48h", "YPDCAFEIN40_48h", "YPD14_48h",
                   "YPDKCL2M_48h", "YPDHU_48h", "YPDBENOMYL500_40h",
                   "YPDCUSO410MM_48h", "SCCuSO405mM_38h",
                   "YPDANISO10_48h", "YPDCAFEIN50_48h"]
    test_traits = [t for t in test_traits if t in phenotypes]
    print(f"  Test traits: {len(test_traits)}")

    # --- Align data for first trait (to compute GRM once) ---
    trait0 = test_traits[0]
    shared = [s for s in sample_ids if s in phenotypes[trait0]]
    id_to_idx = {s: i for i, s in enumerate(sample_ids)}
    cols = np.array([id_to_idx[s] for s in shared])
    geno_aligned = geno[:, cols]
    n_sam = len(shared)
    print(f"  Shared samples: {n_sam}")

    # --- Compute GRM once (reused across all methods and traits) ---
    print(f"\n{'='*70}")
    print(f"  Step 1: Compute GRM")
    print(f"{'='*70}")
    grm = compute_grm(geno_aligned, verbose=True)

    # --- Compute graph-native kinship once ---
    print(f"\n{'='*70}")
    print(f"  Step 2: Compute graph-native kinship")
    print(f"{'='*70}")
    from graphgwas.popstruct import graph_kinship
    gk = graph_kinship(geno_aligned, verbose=True)

    # --- Also try the precomputed distance matrix ---
    dist_path = os.path.join(DATA_DIR, "1011DistanceMatrixBasedOnSNPs.tab.gz")
    if os.path.exists(dist_path):
        print(f"\n  Loading precomputed distance matrix...")
        kinship_precomp, kin_ids = load_distance_matrix(dist_path, sample_ids=shared)
    else:
        kinship_precomp = None

    # =================================================================
    # Run all four methods on each trait
    # =================================================================
    print(f"\n{'='*70}")
    print(f"  Step 3: GWAS with 4 correction methods × {len(test_traits)} traits")
    print(f"{'='*70}")

    methods = ["none", "pca", "grammar", "graph"]
    all_results = {m: [] for m in methods}

    for trait in test_traits:
        shared_t = [s for s in sample_ids if s in phenotypes[trait]]
        cols_t = np.array([id_to_idx[s] for s in shared_t])
        geno_t = geno[:, cols_t]
        pheno_t = np.array([phenotypes[trait][s] for s in shared_t])

        # Need kinship aligned to this trait's samples
        # Quick: recompute from the aligned genotypes (small yeast data)
        grm_t = compute_grm(geno_t, verbose=False)

        print(f"\n  --- {trait} ({len(shared_t)} samples) ---")

        for method in methods:
            if method == "graph":
                gk_t = graph_kinship(geno_t, verbose=False)
                kin = gk_t["kinship_combined"]
            else:
                kin = grm_t

            result = corrected_gwas(geno_t, pheno_t, kinship=kin,
                                    method=method, n_pcs=10, verbose=False)
            all_results[method].append({
                "trait": trait,
                "lambda_gc": result["lambda_gc"],
                "n_bonf": result["n_bonferroni"],
                "n_sug": result["n_suggestive"],
                "min_p": result["min_p"],
            })
            print(f"    {method:<10s} λ_GC={result['lambda_gc']:.3f}  "
                  f"bonf={result['n_bonferroni']:>6,}  "
                  f"sug={result['n_suggestive']:>6,}  "
                  f"min_p={result['min_p']:.2e}")

    # =================================================================
    # Summary
    # =================================================================
    print(f"\n{'='*70}")
    print(f"  SUMMARY")
    print(f"{'='*70}")

    print(f"\n  Mean Lambda GC across {len(test_traits)} traits:")
    print(f"  {'Method':<12} {'Mean λ_GC':>10} {'Median λ_GC':>12} "
          f"{'Mean Bonf':>10} {'Mean Sug':>10}")
    print(f"  {'-'*58}")
    for method in methods:
        results = all_results[method]
        lgcs = [r["lambda_gc"] for r in results]
        bonfs = [r["n_bonf"] for r in results]
        sugs = [r["n_sug"] for r in results]
        print(f"  {method:<12} {np.mean(lgcs):>10.3f} {np.median(lgcs):>12.3f} "
              f"{np.mean(bonfs):>10,.0f} {np.mean(sugs):>10,.0f}")

    # Reduction ratios
    none_lgcs = [r["lambda_gc"] for r in all_results["none"]]
    for method in ["pca", "grammar", "graph"]:
        method_lgcs = [r["lambda_gc"] for r in all_results[method]]
        reduction = (1 - np.mean(method_lgcs) / np.mean(none_lgcs)) * 100
        print(f"\n  {method}: {reduction:.1f}% reduction in λ_GC vs naive")

    elapsed = time.time() - t0
    print(f"\n  Total time: {elapsed:.0f}s ({elapsed/60:.1f} min)")


if __name__ == "__main__":
    main()
