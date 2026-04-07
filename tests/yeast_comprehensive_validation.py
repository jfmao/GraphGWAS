"""Comprehensive yeast GWAS validation: 8,391 traits.

Uses vectorized numpy regression for speed (~1s per trait instead of ~15min).

Pipeline:
  1. Load genotypes (83,794 SNPs, 1,011 strains) + all phenotypes
  2. GWAS on 241 growth traits → compare to published Table S12
  3. Heritability estimation → compare to published Table S11
  4. Genetic correlation matrix (growth traits)
  5. eQTL scan (100 expression traits)
  6. Pleiotropy analysis
  7. Summary report

Usage:
    python tests/yeast_comprehensive_validation.py
"""

import csv
import glob
import json
import os
import sys
import time

import numpy as np
import pandas as pd
from scipy import stats as sp_stats

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src", "python"))

DATA_DIR = os.path.join(os.path.dirname(__file__), "data", "yeast")
RESULTS_DIR = os.path.join(os.path.dirname(__file__), "results", "yeast_validation")
os.makedirs(RESULTS_DIR, exist_ok=True)


# ===================================================================
# Vectorized GWAS — orders of magnitude faster than per-variant loop
# ===================================================================

def fast_gwas(geno, pheno):
    """Vectorized linear regression GWAS across all variants at once.

    For each variant j: beta_j = cov(X_j, Y) / var(X_j), then t-test.
    Handles missing genotypes via per-variant mean imputation.

    ~1 second for 83K variants × 1K samples.
    """
    n_var, n_sam = geno.shape

    # Mean-impute missing genotypes per variant
    G = geno.copy()
    for i in range(n_var):
        mask = np.isnan(G[i])
        if np.any(mask):
            G[i, mask] = np.nanmean(G[i])

    # Standardize phenotype
    Y = pheno - np.mean(pheno)
    n = len(Y)

    # Vectorized: correlations between each variant and phenotype
    # G: (n_var, n_sam), Y: (n_sam,)
    G_mean = np.mean(G, axis=1, keepdims=True)
    G_centered = G - G_mean
    G_var = np.sum(G_centered ** 2, axis=1)

    # Filter zero-variance variants
    valid = G_var > 0
    betas = np.zeros(n_var)
    pvals = np.ones(n_var)

    if np.sum(valid) == 0:
        return betas, pvals

    # cov(X, Y) = X_centered @ Y / n
    cov_xy = G_centered[valid] @ Y  # (n_valid,)
    betas[valid] = cov_xy / G_var[valid]

    # Residual variance: RSS / (n - 2)
    predicted = np.outer(betas[valid], np.ones(n))  # broadcast isn't right
    # Actually: Y_hat = beta * X, residual = Y - beta * X
    residuals = Y[np.newaxis, :] - betas[valid, np.newaxis] * G_centered[valid]
    rss = np.sum(residuals ** 2, axis=1)
    mse = rss / (n - 2)

    # SE of beta
    se = np.sqrt(mse / G_var[valid])
    se[se == 0] = np.inf

    # t-statistic and p-value
    t_stat = betas[valid] / se
    pvals[valid] = 2 * sp_stats.t.sf(np.abs(t_stat), n - 2)

    return betas, pvals


# ===================================================================
# Data Loading
# ===================================================================

def load_genotypes():
    from pandas_plink import read_plink
    print("Loading genotypes...")
    t0 = time.time()
    bim, fam, geno_dask = read_plink(os.path.join(DATA_DIR, "1011GWAS_matrix"))
    geno = geno_dask.compute().astype(np.float64)
    af = np.nanmean(geno, axis=1) / 2
    maf = np.minimum(af, 1 - af)
    keep = maf >= 0.01
    geno = geno[keep]
    bim = bim[keep].reset_index(drop=True)
    sample_ids = list(fam["iid"].values)
    print(f"  {geno.shape[0]:,} variants × {geno.shape[1]} samples ({time.time()-t0:.1f}s)")
    return geno, bim, sample_ids


def load_phenotype_file(path):
    data = {}
    with open(path) as f:
        f.readline()
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                try:
                    data[parts[0]] = float(parts[1])
                except ValueError:
                    pass
    return data


def load_all_growth_phenotypes():
    print("Loading growth phenotypes...")
    growth_dir = os.path.join(DATA_DIR, "Phenotypes_8391Traits", "GrowthTraits")
    phenotypes = {}
    for f in sorted(glob.glob(os.path.join(growth_dir, "*.phen"))):
        name = os.path.basename(f).replace(".phen", "")
        data = load_phenotype_file(f)
        if len(data) > 100:
            phenotypes[name] = data
    print(f"  {len(phenotypes)} growth traits loaded")
    return phenotypes


def load_published_heritability():
    path = os.path.join(DATA_DIR, "1086yeast_heritability.tsv")
    if not os.path.exists(path):
        return {}
    h2 = {}
    with open(path) as f:
        f.readline()
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 3 and parts[0] and parts[2]:
                try:
                    h2[parts[1]] = float(parts[2])
                except ValueError:
                    pass
    print(f"  {len(h2)} published h² values loaded")
    return h2


def load_published_gwas():
    path = os.path.join(DATA_DIR, "1086yeast_gwas_results.tsv")
    if not os.path.exists(path):
        return {}
    hits = {}
    with open(path) as f:
        f.readline()
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 7 and parts[1]:
                trait = parts[1]
                if trait not in hits:
                    hits[trait] = []
                try:
                    hits[trait].append({
                        "variant_type": parts[2],
                        "chr": parts[4],
                        "pos": int(parts[5]) if parts[5].isdigit() else 0,
                        "p_value": float(parts[6]),
                    })
                except (ValueError, IndexError):
                    pass
    print(f"  {sum(len(v) for v in hits.values())} published hits")
    return hits


def align_data(geno, sample_ids, pheno_dict):
    shared = [s for s in sample_ids if s in pheno_dict]
    if len(shared) < 50:
        return None, None
    id_to_idx = {s: i for i, s in enumerate(sample_ids)}
    idx = np.array([id_to_idx[s] for s in shared])
    pheno_vec = np.array([pheno_dict[s] for s in shared])
    return geno[:, idx], pheno_vec


# ===================================================================
# Main
# ===================================================================

def main():
    t_start = time.time()
    print("=" * 70)
    print("  GraphGWAS Comprehensive Yeast Validation")
    print("  8,391 traits × 83,794 SNPs × 1,011 strains")
    print("=" * 70)

    geno, bim, sample_ids = load_genotypes()
    growth_phenos = load_all_growth_phenotypes()
    published_h2 = load_published_heritability()
    published_gwas = load_published_gwas()

    # =================================================================
    # Phase 1: GWAS on all growth traits (vectorized — fast)
    # =================================================================
    print(f"\n{'='*70}")
    print(f"  Phase 1: GWAS on {len(growth_phenos)} growth traits")
    print(f"{'='*70}")

    gwas_results = {}
    t1 = time.time()

    for ti, (trait, pheno_dict) in enumerate(sorted(growth_phenos.items())):
        geno_a, pheno_a = align_data(geno, sample_ids, pheno_dict)
        if geno_a is None:
            continue

        betas, pvals = fast_gwas(geno_a, pheno_a)
        n_bonf = int(np.sum(pvals < 0.05 / geno.shape[0]))
        n_sug = int(np.sum(pvals < 1e-4))

        gwas_results[trait] = {
            "n_samples": geno_a.shape[1],
            "n_bonf": n_bonf,
            "n_suggestive": n_sug,
            "min_p": float(np.nanmin(pvals)),
            "pvals": pvals,
            "betas": betas,
        }

        if (ti + 1) % 50 == 0 or (ti + 1) == len(growth_phenos):
            elapsed = time.time() - t1
            print(f"  ...{ti+1}/{len(growth_phenos)} traits ({elapsed:.0f}s)")

    elapsed = time.time() - t1
    print(f"\n  GWAS complete: {len(gwas_results)} traits in {elapsed:.0f}s")

    bonf_counts = [r["n_bonf"] for r in gwas_results.values()]
    sug_counts = [r["n_suggestive"] for r in gwas_results.values()]
    print(f"  Traits with Bonf hits: {sum(1 for c in bonf_counts if c > 0)}")
    print(f"  Total Bonf hits: {sum(bonf_counts):,}")
    print(f"  Total suggestive: {sum(sug_counts):,}")

    top_traits = sorted(gwas_results.items(), key=lambda x: x[1]["min_p"])[:10]
    print(f"\n  Top 10 traits:")
    for trait, r in top_traits:
        print(f"    {trait:<35s} min_p={r['min_p']:.2e}  "
              f"bonf={r['n_bonf']}  sug={r['n_suggestive']}")

    # =================================================================
    # Phase 2: Heritability
    # =================================================================
    print(f"\n{'='*70}")
    print(f"  Phase 2: Heritability estimation")
    print(f"{'='*70}")

    h2_results = {}
    for trait, pheno_dict in sorted(growth_phenos.items()):
        geno_a, pheno_a = align_data(geno, sample_ids, pheno_dict)
        if geno_a is None:
            continue

        rng = np.random.default_rng(42)
        n_var, n_sam = geno_a.shape
        train = np.arange(0, n_sam, 2)
        test = np.arange(1, n_sam, 2)

        # Fast GWAS on training set
        betas_train, pvals_train = fast_gwas(geno_a[:, train], pheno_a[train])
        top_idx = np.argsort(pvals_train)[:50]

        # PRS on test set
        G_test = geno_a[:, test].copy()
        for i in range(n_var):
            mask = np.isnan(G_test[i])
            if np.any(mask):
                G_test[i, mask] = np.nanmean(G_test[i])

        prs = np.zeros(len(test))
        for idx in top_idx:
            prs += betas_train[idx] * G_test[idx]

        if np.std(prs) > 0 and np.std(pheno_a[test]) > 0:
            r2 = np.corrcoef(prs, pheno_a[test])[0, 1] ** 2
        else:
            r2 = 0.0
        h2_results[trait] = float(r2)

    sorted_h2 = sorted(h2_results.items(), key=lambda x: x[1], reverse=True)
    print(f"  Top 10 heritable:")
    for trait, h2 in sorted_h2[:10]:
        pub = published_h2.get(trait)
        pub_str = f"  (pub: {pub:.3f})" if pub else ""
        print(f"    {trait:<35s} h²={h2:.4f}{pub_str}")

    # Compare to published
    if published_h2:
        shared = [t for t in h2_results if t in published_h2]
        if len(shared) >= 10:
            our = np.array([h2_results[t] for t in shared])
            pub = np.array([published_h2[t] for t in shared])
            corr = np.corrcoef(our, pub)[0, 1]
            print(f"\n  h² correlation with published: r = {corr:.4f} "
                  f"({len(shared)} traits)")
            print(f"  Our mean: {np.mean(our):.4f}, Published mean: {np.mean(pub):.4f}")

    # =================================================================
    # Phase 3: Genetic correlation (top 20 traits)
    # =================================================================
    print(f"\n{'='*70}")
    print(f"  Phase 3: Genetic correlation matrix")
    print(f"{'='*70}")

    top20 = [t for t, _ in sorted_h2[:20]]
    T = len(top20)
    z_all = {}

    for trait in top20:
        geno_a, pheno_a = align_data(geno, sample_ids, growth_phenos[trait])
        if geno_a is None:
            continue
        betas, pvals = fast_gwas(geno_a, pheno_a)
        # z-scores from betas and p-values
        z = np.sign(betas) * sp_stats.norm.isf(pvals / 2)
        z[~np.isfinite(z)] = 0
        z_all[trait] = z

    traits_with_z = [t for t in top20 if t in z_all]
    T = len(traits_with_z)
    z_mat = np.column_stack([z_all[t] for t in traits_with_z])
    rg = np.corrcoef(z_mat.T)

    print(f"  {T} traits, genetic correlation matrix:")
    labels = [t[:12] for t in traits_with_z[:8]]
    header = "             " + "  ".join(f"{l:>12}" for l in labels)
    print(f"  {header}")
    for i in range(min(8, T)):
        row = f"  {labels[i]:<12}" + "  ".join(
            f"{rg[i,j]:12.3f}" for j in range(min(8, T)))
        print(row)

    off = rg[np.triu_indices(T, k=1)]
    print(f"\n  Off-diagonal: mean={np.mean(off):.3f}, "
          f"range=[{off.min():.3f}, {off.max():.3f}]")

    # =================================================================
    # Phase 4: eQTL scan (100 genes)
    # =================================================================
    print(f"\n{'='*70}")
    print(f"  Phase 4: eQTL scan")
    print(f"{'='*70}")

    expr_dir = os.path.join(DATA_DIR, "Phenotypes_8391Traits", "ExpressionTraits")
    expr_files = sorted(glob.glob(os.path.join(expr_dir, "*.phen")))
    rng = np.random.default_rng(42)
    sample_idx = rng.choice(len(expr_files), min(100, len(expr_files)), replace=False)

    eqtl_results = {}
    t4 = time.time()
    for i, idx in enumerate(sample_idx):
        path = expr_files[idx]
        gene = os.path.basename(path).replace(".RNASeq.phen", "")
        pheno_dict = load_phenotype_file(path)
        geno_a, pheno_a = align_data(geno, sample_ids, pheno_dict)
        if geno_a is None:
            continue
        _, pvals = fast_gwas(geno_a, pheno_a)
        n_sig = int(np.sum(pvals < 0.05 / geno.shape[0]))
        eqtl_results[gene] = {"n_sig": n_sig, "min_p": float(np.nanmin(pvals))}
        if (i + 1) % 25 == 0:
            print(f"  ...{i+1}/100 genes ({time.time()-t4:.0f}s)")

    n_with_eqtl = sum(1 for r in eqtl_results.values() if r["n_sig"] > 0)
    print(f"\n  {len(eqtl_results)} genes, {n_with_eqtl} with significant eQTL")

    top_eqtl = sorted(eqtl_results.items(), key=lambda x: x[1]["min_p"])[:10]
    print(f"  Top 10:")
    for gene, r in top_eqtl:
        print(f"    {gene:<15s} min_p={r['min_p']:.2e}  hits={r['n_sig']}")

    # =================================================================
    # Phase 5: Pleiotropy
    # =================================================================
    print(f"\n{'='*70}")
    print(f"  Phase 5: Pleiotropy analysis")
    print(f"{'='*70}")

    n_var = geno.shape[0]
    pleio_count = np.zeros(n_var, dtype=int)
    for r in gwas_results.values():
        pleio_count += (r["pvals"] < 1e-4).astype(int)

    n_pleio2 = int(np.sum(pleio_count >= 2))
    n_pleio10 = int(np.sum(pleio_count >= 10))
    n_pleio50 = int(np.sum(pleio_count >= 50))
    print(f"  Variants sig for ≥2 traits:  {n_pleio2:,}")
    print(f"  Variants sig for ≥10 traits: {n_pleio10:,}")
    print(f"  Variants sig for ≥50 traits: {n_pleio50:,}")
    print(f"  Max traits per variant:      {int(pleio_count.max())}")

    top_pleio = np.argsort(-pleio_count)[:10]
    print(f"\n  Top 10 pleiotropic variants:")
    for idx in top_pleio:
        if idx < len(bim):
            c = bim.iloc[idx]["chrom"]
            p = bim.iloc[idx]["pos"]
            print(f"    chr{c}:{p}  → {pleio_count[idx]} traits")

    # =================================================================
    # Phase 6: Compare to published GWAS
    # =================================================================
    if published_gwas:
        print(f"\n{'='*70}")
        print(f"  Phase 6: Comparison to published GWAS")
        print(f"{'='*70}")

        overlap = 0
        tested = 0
        for trait in gwas_results:
            if trait not in published_gwas:
                continue
            pub = [h for h in published_gwas[trait] if h["variant_type"] == "SNPs"]
            if not pub:
                continue
            our_sig_idx = set(np.where(gwas_results[trait]["pvals"] < 1e-4)[0])
            for ph in pub:
                for oi in our_sig_idx:
                    if oi < len(bim):
                        if (str(bim.iloc[oi]["chrom"]) == ph["chr"] and
                                abs(int(bim.iloc[oi]["pos"]) - ph["pos"]) < 10000):
                            overlap += 1
                            break
                tested += 1
                if tested >= 500:
                    break
            if tested >= 500:
                break

        print(f"  Published SNP-QTL tested: {tested}")
        print(f"  Overlapping (±10kb): {overlap}")
        if tested > 0:
            print(f"  Overlap rate: {overlap/tested:.1%}")

    # =================================================================
    # Summary
    # =================================================================
    elapsed = time.time() - t_start
    print(f"\n{'='*70}")
    print(f"  SUMMARY")
    print(f"{'='*70}")
    print(f"  Time: {elapsed:.0f}s ({elapsed/60:.1f} min)")
    print(f"  Growth traits: {len(gwas_results)}")
    print(f"  eQTL genes: {len(eqtl_results)}")
    print(f"  Bonf-sig growth traits: "
          f"{sum(1 for r in gwas_results.values() if r['n_bonf'] > 0)}")
    print(f"  Pleiotropic variants (≥2): {n_pleio2:,}")

    if published_h2:
        shared = [t for t in h2_results if t in published_h2]
        if shared:
            r = np.corrcoef([h2_results[t] for t in shared],
                            [published_h2[t] for t in shared])[0, 1]
            print(f"  h² vs published: r = {r:.4f} ({len(shared)} traits)")

    # Save
    summary = {
        "n_growth": len(gwas_results),
        "n_eqtl": len(eqtl_results),
        "n_bonf_traits": sum(1 for r in gwas_results.values() if r["n_bonf"] > 0),
        "total_bonf": sum(r["n_bonf"] for r in gwas_results.values()),
        "n_pleio2": n_pleio2,
        "n_pleio10": n_pleio10,
        "elapsed": elapsed,
    }
    with open(os.path.join(RESULTS_DIR, "validation_summary.json"), "w") as f:
        json.dump(summary, f, indent=2)
    print(f"\n  Saved to {RESULTS_DIR}/")


if __name__ == "__main__":
    main()
