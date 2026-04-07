"""Yeast SV validation: SNPs + PAV + CNV + Frameshift combined GWAS.

Compares: SNP-only vs SNP+SV heritability and GWAS power.
Validates the 14.3% heritability improvement reported by Loegler et al. 2025.

Usage:
    python tests/yeast_sv_validation.py
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
    print("  Yeast SV Validation: SNPs + PAV + CNV + Frameshift")
    print("=" * 70)

    from graphgwas.sv import (load_pav_matrix, load_cnv_matrix,
                               load_frameshift_matrix,
                               combine_genotype_matrices,
                               sv_gwas, sv_heritability, fast_gwas)

    # --- Load SNP genotypes ---
    print("\n--- Loading SNP genotypes ---")
    from pandas_plink import read_plink
    bim, fam, geno_dask = read_plink(os.path.join(DATA_DIR, "1011GWAS_matrix"))
    snp_geno = geno_dask.compute().astype(np.float64)
    snp_sample_ids = list(fam["iid"].values)

    # MAF filter
    af = np.nanmean(snp_geno, axis=1) / 2
    maf = np.minimum(af, 1 - af)
    snp_geno = snp_geno[maf >= 0.01]
    print(f"  SNPs: {snp_geno.shape[0]:,} × {snp_geno.shape[1]}")

    # --- Load SV matrices ---
    print("\n--- Loading SV matrices ---")
    pav = load_pav_matrix(os.path.join(DATA_DIR, "genesMatrix_PresenceAbsence.tab.gz"),
                          sample_ids=snp_sample_ids)
    cnv = load_cnv_matrix(os.path.join(DATA_DIR, "genesMatrix_CopyNumber.tab.gz"),
                          sample_ids=snp_sample_ids)
    lof = load_frameshift_matrix(os.path.join(DATA_DIR, "genesMatrix_Frameshift.tab.gz"),
                                  sample_ids=snp_sample_ids)

    # --- Combine ---
    print("\n--- Combining SNP + SV matrices ---")
    combined = combine_genotype_matrices(snp_geno, [pav, cnv, lof], snp_sample_ids)

    # --- Load growth phenotypes ---
    print("\n--- Loading phenotypes ---")
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
    print(f"  {len(phenotypes)} growth traits")

    # =================================================================
    # Test 1: SV-GWAS on top traits
    # =================================================================
    print(f"\n{'='*70}")
    print(f"  Test 1: SV-GWAS (top 10 traits)")
    print(f"{'='*70}")

    # Pick the most heritable traits from our previous run
    top_traits = ["YPDCUSO410MM_48h", "YPDCAFEIN50_48h", "YPDCAFEIN40_48h",
                  "YPDBENOMYL500_40h", "SCCuSO405mM_38h", "YPD14_48h",
                  "YPDKCL2M_48h", "YPETHANOL_48h", "YPDANISO10_48h",
                  "YPDHU_48h"]
    top_traits = [t for t in top_traits if t in phenotypes]

    for trait in top_traits[:5]:
        print(f"\n  --- {trait} ---")
        result = sv_gwas(combined, phenotypes[trait])

    # =================================================================
    # Test 2: SNP-only vs SNP+SV heritability
    # =================================================================
    print(f"\n{'='*70}")
    print(f"  Test 2: Heritability — SNPs only vs SNPs+SVs")
    print(f"{'='*70}")

    h2_comparison = {}
    for trait in top_traits:
        print(f"\n  --- {trait} ---")
        h2 = sv_heritability(combined, phenotypes[trait])
        if "error" not in h2:
            h2_comparison[trait] = h2

    # Summary
    if h2_comparison:
        print(f"\n{'='*70}")
        print(f"  Heritability Summary")
        print(f"{'='*70}")
        print(f"  {'Trait':<30s} {'h²_SNP':>8} {'h²_SV':>8} {'h²_All':>8} {'Δh²':>8} {'%Increase':>10}")
        print(f"  {'-'*80}")
        sv_increases = []
        for trait, h2 in h2_comparison.items():
            delta = h2["h2_sv_contribution"]
            pct = h2["h2_sv_pct_increase"]
            sv_increases.append(pct)
            print(f"  {trait:<30s} {h2['h2_snp_only']:>8.4f} {h2['h2_sv_only']:>8.4f} "
                  f"{h2['h2_combined']:>8.4f} {delta:>+8.4f} {pct:>+10.1f}%")

        mean_increase = np.mean(sv_increases)
        print(f"\n  Mean SV heritability increase: {mean_increase:+.1f}%")
        print(f"  (Published: +14.3% from Loegler et al. 2025)")

    # =================================================================
    # Test 3: SV pleiotropy vs SNP pleiotropy
    # =================================================================
    print(f"\n{'='*70}")
    print(f"  Test 3: SV vs SNP Pleiotropy (20 traits)")
    print(f"{'='*70}")

    # Use first 20 traits for speed
    trait_subset = {t: phenotypes[t] for t in sorted(phenotypes.keys())[:20]}

    from graphgwas.sv import sv_pleiotropy
    pleio = sv_pleiotropy(combined, trait_subset)

    # =================================================================
    # Summary
    # =================================================================
    elapsed = time.time() - t0
    print(f"\n{'='*70}")
    print(f"  FINAL SUMMARY")
    print(f"{'='*70}")
    print(f"  Total time: {elapsed:.0f}s ({elapsed/60:.1f} min)")
    print(f"  Combined matrix: {combined['genotype_matrix'].shape[0]:,} variants "
          f"({combined['n_snps']:,} SNPs + {combined['n_svs']:,} SVs)")
    if h2_comparison:
        mean_snp = np.mean([h["h2_snp_only"] for h in h2_comparison.values()])
        mean_all = np.mean([h["h2_combined"] for h in h2_comparison.values()])
        print(f"  Mean h² (SNPs only): {mean_snp:.4f}")
        print(f"  Mean h² (SNPs+SVs):  {mean_all:.4f}")
        print(f"  Mean SV increase:    {np.mean(sv_increases):+.1f}%")


if __name__ == "__main__":
    main()
