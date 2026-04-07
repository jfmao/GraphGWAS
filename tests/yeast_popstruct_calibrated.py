"""Test GRAMMAR+ calibration and auto method on yeast.

Compares all 7 methods: none, pca, grammar, grammar+, graph, graph+, auto.

Usage: python tests/yeast_popstruct_calibrated.py
"""

import glob, os, sys, time
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src", "python"))
DATA_DIR = os.path.join(os.path.dirname(__file__), "data", "yeast")


def main():
    t0 = time.time()
    print("=" * 70)
    print("  GRAMMAR+ Calibration & Auto Method Validation")
    print("=" * 70)

    from graphgwas.popstruct import compute_grm, corrected_gwas

    # Load data
    from pandas_plink import read_plink
    bim, fam, geno_dask = read_plink(os.path.join(DATA_DIR, "1011GWAS_matrix"))
    geno = geno_dask.compute().astype(np.float64)
    sample_ids = list(fam["iid"].values)
    af = np.nanmean(geno, axis=1) / 2
    maf = np.minimum(af, 1 - af)
    geno = geno[maf >= 0.01]
    print(f"  Genotypes: {geno.shape[0]:,} × {geno.shape[1]}")

    # Load phenotypes
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

    test_traits = ["YPETHANOL_48h", "YPDCAFEIN40_48h", "YPD14_48h",
                   "YPDCUSO410MM_48h", "YPDBENOMYL500_40h"]
    test_traits = [t for t in test_traits if t in phenotypes]

    methods = ["none", "pca", "grammar", "grammar+", "graph", "graph+", "auto"]

    all_results = {m: [] for m in methods}

    for trait in test_traits:
        shared = [s for s in sample_ids if s in phenotypes[trait]]
        id_to_idx = {s: i for i, s in enumerate(sample_ids)}
        cols = np.array([id_to_idx[s] for s in shared])
        geno_t = geno[:, cols]
        pheno_t = np.array([phenotypes[trait][s] for s in shared])

        grm_t = compute_grm(geno_t, verbose=False)

        print(f"\n  === {trait} ({len(shared)} samples) ===")
        for method in methods:
            result = corrected_gwas(geno_t, pheno_t, kinship=grm_t,
                                    method=method, n_pcs=10, verbose=False)
            all_results[method].append(result)
            print(f"    {method:<12s} λ_GC={result['lambda_gc']:.3f}  "
                  f"bonf={result['n_bonferroni']:>6,}  "
                  f"sug={result['n_suggestive']:>6,}  "
                  f"min_p={result['min_p']:.2e}")

    # Summary
    print(f"\n{'='*70}")
    print(f"  SUMMARY: Mean λ_GC across {len(test_traits)} traits")
    print(f"{'='*70}")
    print(f"  {'Method':<12s} {'Mean λ_GC':>10} {'|λ-1|':>8} {'Mean Bonf':>10} {'Mean Sug':>10}")
    print(f"  {'-'*54}")
    for method in methods:
        lgcs = [r["lambda_gc"] for r in all_results[method]]
        bonfs = [r["n_bonferroni"] for r in all_results[method]]
        sugs = [r["n_suggestive"] for r in all_results[method]]
        mean_lgc = np.mean(lgcs)
        deviation = abs(mean_lgc - 1.0)
        print(f"  {method:<12s} {mean_lgc:>10.3f} {deviation:>8.3f} "
              f"{np.mean(bonfs):>10,.0f} {np.mean(sugs):>10,.0f}")

    elapsed = time.time() - t0
    print(f"\n  Time: {elapsed:.0f}s ({elapsed/60:.1f} min)")


if __name__ == "__main__":
    main()
