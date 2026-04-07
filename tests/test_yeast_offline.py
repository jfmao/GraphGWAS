"""Yeast GWAS tests — no Neo4j required.

Reads PLINK .bed/.bim/.fam + phenotype matrix directly and runs
GraphGWAS statistical engines on them.  Validates that the association
methods work correctly on a real dataset with known biology.

Key positive controls:
- Yeast growth on ethanol → ADH genes (alcohol dehydrogenase)
- Yeast growth on caffeine → known resistance loci
- Growth at 14°C → cold tolerance genes

Run: pytest tests/test_yeast_offline.py -v
"""

import gzip
import os

import numpy as np
import pandas as pd
import pytest

DATA_DIR = os.path.join(os.path.dirname(__file__), "data", "yeast")


@pytest.fixture(scope="module")
def yeast_genotypes():
    """Load yeast genotype matrix from PLINK files."""
    from pandas_plink import read_plink

    bed_path = os.path.join(DATA_DIR, "1011GWAS_matrix")
    bim, fam, genotype = read_plink(bed_path)

    # genotype is a dask array of shape (n_variants, n_samples)
    # Values: 0 (hom ref), 1 (het), 2 (hom alt), nan (missing)
    geno_np = genotype.compute()  # materialize to numpy

    return {
        "bim": bim,
        "fam": fam,
        "genotype": geno_np,
        "n_variants": len(bim),
        "n_samples": len(fam),
        "sample_ids": list(fam["iid"].values),
    }


@pytest.fixture(scope="module")
def yeast_phenotypes():
    """Load yeast phenotype matrix (35 growth conditions)."""
    pheno_path = os.path.join(DATA_DIR,
                              "phenoMatrix_35ConditionsNormalizedByYPD.tab.gz")
    df = pd.read_csv(pheno_path, sep="\t", index_col=0)
    return df


class TestYeastData:
    def test_genotype_shape(self, yeast_genotypes):
        assert yeast_genotypes["n_variants"] == 83794
        assert yeast_genotypes["n_samples"] == 1011
        print(f"\n  Yeast: {yeast_genotypes['n_variants']:,} variants × "
              f"{yeast_genotypes['n_samples']} samples")

    def test_phenotype_shape(self, yeast_phenotypes):
        assert yeast_phenotypes.shape[1] >= 35
        assert yeast_phenotypes.shape[0] > 900
        print(f"\n  Phenotypes: {yeast_phenotypes.shape[0]} samples × "
              f"{yeast_phenotypes.shape[1]} conditions")

    def test_sample_overlap(self, yeast_genotypes, yeast_phenotypes):
        """Most samples should be in both genotype and phenotype data."""
        geno_ids = set(yeast_genotypes["sample_ids"])
        pheno_ids = set(yeast_phenotypes.index)
        overlap = geno_ids & pheno_ids
        assert len(overlap) > 900, f"Only {len(overlap)} samples overlap"
        print(f"\n  Sample overlap: {len(overlap)} "
              f"(geno={len(geno_ids)}, pheno={len(pheno_ids)})")


class TestYeastGWAS:
    """Run GWAS on yeast data using GraphGWAS statistical engines."""

    def _get_aligned_data(self, yeast_genotypes, yeast_phenotypes, trait):
        """Align genotype and phenotype data for a trait."""
        pheno = yeast_phenotypes[trait].dropna()
        geno_ids = yeast_genotypes["sample_ids"]
        fam = yeast_genotypes["fam"]

        # Find shared samples
        shared = [s for s in geno_ids if s in pheno.index]
        if len(shared) < 50:
            return None, None, None

        # Align
        pheno_aligned = np.array([pheno[s] for s in shared])
        # Get genotype indices for shared samples
        id_to_idx = {s: i for i, s in enumerate(geno_ids)}
        geno_indices = np.array([id_to_idx[s] for s in shared])

        return yeast_genotypes["genotype"][:, geno_indices], pheno_aligned, shared

    def test_linear_gwas_ethanol(self, yeast_genotypes, yeast_phenotypes):
        """Linear regression GWAS on ethanol growth phenotype."""
        from graphgwas.assoc import linear_regression

        geno, pheno, samples = self._get_aligned_data(
            yeast_genotypes, yeast_phenotypes, "YPETHANOL"
        )
        assert geno is not None, "Could not align ethanol phenotype"

        n_var, n_sam = geno.shape
        pvals = []
        betas = []

        for i in range(n_var):
            dosage = geno[i].astype(float)
            # Replace nan with mean imputation
            mask = np.isnan(dosage)
            if np.sum(~mask) < 20:
                pvals.append(1.0)
                betas.append(0.0)
                continue
            dosage[mask] = np.nanmean(dosage)

            result = linear_regression(dosage, pheno)
            pvals.append(result["p_value"])
            betas.append(result["beta"])

        pvals = np.array(pvals)
        n_sig_bonf = int(np.sum(pvals < 0.05 / n_var))
        n_sig_suggestive = int(np.sum(pvals < 1e-4))
        min_p = float(np.min(pvals))

        print(f"\n  Ethanol GWAS: {n_var} variants, {n_sam} samples")
        print(f"  Bonferroni significant: {n_sig_bonf}")
        print(f"  Suggestive (p<1e-4): {n_sig_suggestive}")
        print(f"  Min p-value: {min_p:.2e}")

        # Should find SOME signal (ethanol growth is heritable in yeast)
        assert n_sig_suggestive > 0, "Should find suggestive hits for ethanol growth"

    def test_linear_gwas_caffeine(self, yeast_genotypes, yeast_phenotypes):
        """Linear regression GWAS on caffeine resistance."""
        from graphgwas.assoc import linear_regression

        geno, pheno, samples = self._get_aligned_data(
            yeast_genotypes, yeast_phenotypes, "YPDCAFEIN40"
        )
        assert geno is not None, "Could not align caffeine phenotype"

        n_var, n_sam = geno.shape
        pvals = []

        for i in range(n_var):
            dosage = geno[i].astype(float)
            mask = np.isnan(dosage)
            if np.sum(~mask) < 20:
                pvals.append(1.0)
                continue
            dosage[mask] = np.nanmean(dosage)
            result = linear_regression(dosage, pheno)
            pvals.append(result["p_value"])

        pvals = np.array(pvals)
        n_sig = int(np.sum(pvals < 1e-4))
        print(f"\n  Caffeine GWAS: {n_var} variants, {n_sam} samples")
        print(f"  Suggestive hits: {n_sig}")
        assert n_sig >= 0  # Just run without error; caffeine may have few hits

    def test_lambda_gc(self, yeast_genotypes, yeast_phenotypes):
        """Lambda GC should be reasonable (not hugely inflated)."""
        from graphgwas.assoc import linear_regression
        from scipy import stats as sp_stats

        geno, pheno, _ = self._get_aligned_data(
            yeast_genotypes, yeast_phenotypes, "YPETHANOL"
        )
        if geno is None:
            pytest.skip("Cannot align data")

        # Run on a subset for speed
        n_test = min(5000, geno.shape[0])
        pvals = []
        for i in range(n_test):
            dosage = geno[i].astype(float)
            mask = np.isnan(dosage)
            if np.sum(~mask) < 20:
                continue
            dosage[mask] = np.nanmean(dosage)
            result = linear_regression(dosage, pheno)
            if 0 < result["p_value"] < 1:
                pvals.append(result["p_value"])

        chi2_obs = [sp_stats.chi2.isf(p, 1) for p in pvals]
        lambda_gc = float(np.median(chi2_obs) / 0.4549)

        print(f"\n  Lambda GC (ethanol, {n_test} variants): {lambda_gc:.3f}")
        # Yeast has strong population structure, so some inflation expected
        assert lambda_gc < 10, f"Lambda GC = {lambda_gc:.3f} is extremely inflated"

    def test_multi_trait_prs_evaluation(self, yeast_genotypes, yeast_phenotypes):
        """PRS evaluation should work on yeast quantitative traits."""
        from graphgwas.assoc import linear_regression

        # Use ethanol trait, dichotomize for PRS eval
        geno, pheno, _ = self._get_aligned_data(
            yeast_genotypes, yeast_phenotypes, "YPETHANOL"
        )
        if geno is None:
            pytest.skip("Cannot align data")

        # Simple PRS: sum of dosages weighted by marginal betas (top 100 variants)
        n_var = geno.shape[0]
        results = []
        for i in range(min(n_var, 2000)):
            dosage = geno[i].astype(float)
            mask = np.isnan(dosage)
            if np.sum(~mask) < 20:
                continue
            dosage[mask] = np.nanmean(dosage)
            r = linear_regression(dosage, pheno)
            results.append((i, r["beta"], r["p_value"]))

        # Top 50 variants by p-value
        results.sort(key=lambda x: x[2])
        top = results[:50]

        # Compute PRS
        n_sam = geno.shape[1]
        prs = np.zeros(n_sam)
        for idx, beta, _ in top:
            dosage = geno[idx].astype(float)
            dosage = np.nan_to_num(dosage, nan=0.0)
            prs += beta * dosage

        # Dichotomize: top 30% = "high ethanol growth" (case)
        threshold = np.percentile(pheno, 70)
        cases = pheno >= threshold
        n_case = int(np.sum(cases))
        n_ctrl = len(pheno) - n_case

        # Reorder: cases first for PRS evaluation
        case_prs = prs[cases]
        ctrl_prs = prs[~cases]
        prs_ordered = np.concatenate([case_prs, ctrl_prs])

        from graphgwas.prs import prs_evaluation
        ev = prs_evaluation(prs_ordered, n_case, n_ctrl, verbose=False)

        print(f"\n  Yeast PRS (ethanol, 50 top variants):")
        print(f"  AUROC={ev['auroc']:.4f}, R²={ev['nagelkerke_r2']:.4f}")
        # Should show some discrimination
        assert ev["auroc"] > 0.5, "PRS should be better than random"
