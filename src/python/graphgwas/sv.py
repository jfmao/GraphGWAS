"""Structural Variant (SV) support for GraphGWAS.

Extends GraphGWAS to handle structural variants: gene presence/absence (PAV),
copy number variation (CNV), frameshifts (LoF), and arbitrary SV genotype
matrices.  SVs are 3× more likely to be causal than SNPs and 2× more
pleiotropic (Loegler et al. Nature 2025).

The key insight: most SVs are biallelic (carrier/non-carrier) and use the
same dosage-based statistics as SNPs.  CNVs with >2 copy states use a
generalised dosage.  The graph naturally unifies SNPs and SVs as Variant
nodes with a `variant_type` property.

Components:
    S1. load_pav_matrix       — gene presence/absence → dosage matrix
    S2. load_cnv_matrix       — gene copy number → dosage matrix
    S3. load_frameshift_matrix — gene frameshifts → binary matrix
    S4. sv_gwas               — unified GWAS across SNPs + SVs
    S5. sv_heritability       — heritability from combined SNP+SV matrix
    S6. sv_pleiotropy         — SV-specific pleiotropy analysis
"""

from __future__ import annotations

import gzip
import os
import time

import numpy as np
from scipy import stats as sp_stats


# ===================================================================
# S1. Load PAV Matrix
# ===================================================================

def load_pav_matrix(path: str, sample_ids: list[str] | None = None,
                    verbose: bool = True) -> dict:
    """Load gene presence/absence variation matrix.

    Format: TSV with strain IDs as rows, gene IDs as columns.
    Values: 0 (absent) or 1 (present).

    Args:
        path: path to .tab.gz or .tsv file.
        sample_ids: if provided, align to these sample IDs.
        verbose: print progress.

    Returns:
        dict with dosage_matrix (n_genes × n_samples), gene_ids, sample_ids.
    """
    return _load_gene_matrix(path, "PAV", sample_ids, verbose)


# ===================================================================
# S2. Load CNV Matrix
# ===================================================================

def load_cnv_matrix(path: str, sample_ids: list[str] | None = None,
                    verbose: bool = True) -> dict:
    """Load gene copy number variation matrix.

    Values: integer copy number (0, 1, 2, 3, ...).

    Args:
        path: path to .tab.gz or .tsv file.
        sample_ids: align to these sample IDs.
        verbose: print progress.

    Returns:
        dict with dosage_matrix (n_genes × n_samples), gene_ids, sample_ids.
    """
    return _load_gene_matrix(path, "CNV", sample_ids, verbose)


# ===================================================================
# S3. Load Frameshift Matrix
# ===================================================================

def load_frameshift_matrix(path: str, sample_ids: list[str] | None = None,
                           verbose: bool = True) -> dict:
    """Load gene frameshift (loss-of-function) matrix.

    Format: TSV with gene metadata columns + strain columns.
    Values: 0 (intact) or 1 (frameshift).

    The frameshift matrix has a different format: rows are genes with
    metadata (CODE, Gene, SGDorder, CH, START, STOP) then strain columns.

    Args:
        path: path to .tab.gz or .tsv file.
        sample_ids: align to these sample IDs.
        verbose: print progress.

    Returns:
        dict with dosage_matrix, gene_ids, gene_info, sample_ids.
    """
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as f:
        header = f.readline().strip().split("\t")

        # Find where strain columns start (after metadata columns)
        meta_cols = []
        strain_start = 0
        for i, col in enumerate(header):
            if col in ("CODE", "Gene", "SGDorder", "CH", "START", "STOP"):
                meta_cols.append(col)
                strain_start = i + 1
            else:
                break
        # If no recognized meta columns, strain start is at column 6
        if not meta_cols:
            strain_start = 6
            meta_cols = header[:6]

        strain_ids = header[strain_start:]

        genes = []
        gene_info = []
        dosages = []

        for line in f:
            parts = line.strip().split("\t")
            if len(parts) <= strain_start:
                continue

            # Gene metadata
            info = {meta_cols[i]: parts[i] for i in range(min(len(meta_cols), len(parts)))}
            gene_id = info.get("Gene", parts[1] if len(parts) > 1 else f"gene_{len(genes)}")
            genes.append(gene_id)
            gene_info.append(info)

            # Dosage values
            vals = []
            for v in parts[strain_start:]:
                try:
                    vals.append(float(v))
                except ValueError:
                    vals.append(np.nan)
            dosages.append(vals)

    dosage_matrix = np.array(dosages, dtype=np.float64)

    # Align to sample_ids if provided
    if sample_ids is not None:
        strain_idx_map = {s: i for i, s in enumerate(strain_ids)}
        shared = [s for s in sample_ids if s in strain_idx_map]
        if shared:
            col_idx = np.array([strain_idx_map[s] for s in shared])
            dosage_matrix = dosage_matrix[:, col_idx]
            strain_ids = shared

    # Filter zero-variance genes
    var = np.nanvar(dosage_matrix, axis=1)
    keep = var > 0
    dosage_matrix = dosage_matrix[keep]
    genes = [g for g, k in zip(genes, keep) if k]
    gene_info = [g for g, k in zip(gene_info, keep) if k]

    if verbose:
        print(f"  Frameshift matrix: {dosage_matrix.shape[0]} genes × "
              f"{dosage_matrix.shape[1]} samples (LoF)")

    return {
        "dosage_matrix": dosage_matrix,
        "variant_ids": [f"LoF_{g}" for g in genes],
        "gene_ids": genes,
        "gene_info": gene_info,
        "sample_ids": strain_ids,
        "variant_type": "frameshift",
    }


def _load_gene_matrix(path: str, variant_type: str,
                      sample_ids: list[str] | None = None,
                      verbose: bool = True) -> dict:
    """Generic loader for strain × gene matrices (PAV, CNV)."""
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as f:
        header = f.readline().strip().split("\t")
        gene_ids = header[1:]  # first column is strain ID

        strain_ids = []
        rows = []
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 2:
                continue
            strain_ids.append(parts[0])
            vals = []
            for v in parts[1:]:
                try:
                    vals.append(float(v))
                except ValueError:
                    vals.append(np.nan)
            rows.append(vals)

    # Matrix: strains × genes → transpose to genes × strains (variant × sample)
    raw = np.array(rows, dtype=np.float64)  # (n_strains, n_genes)
    dosage_matrix = raw.T  # (n_genes, n_strains)

    # Align to sample_ids if provided
    if sample_ids is not None:
        strain_idx_map = {s: i for i, s in enumerate(strain_ids)}
        shared = [s for s in sample_ids if s in strain_idx_map]
        if shared:
            col_idx = np.array([strain_idx_map[s] for s in shared])
            dosage_matrix = dosage_matrix[:, col_idx]
            strain_ids = shared

    # Filter: remove zero-variance genes
    var = np.nanvar(dosage_matrix, axis=1)
    keep = var > 0
    dosage_matrix = dosage_matrix[keep]
    filtered_genes = [g for g, k in zip(gene_ids, keep) if k]

    if verbose:
        print(f"  {variant_type} matrix: {dosage_matrix.shape[0]} genes × "
              f"{dosage_matrix.shape[1]} samples")
        if variant_type == "CNV":
            max_cn = int(np.nanmax(dosage_matrix))
            print(f"    Copy number range: [0, {max_cn}]")

    return {
        "dosage_matrix": dosage_matrix,
        "variant_ids": [f"{variant_type}_{g}" for g in filtered_genes],
        "gene_ids": filtered_genes,
        "sample_ids": strain_ids if isinstance(strain_ids, list) else list(strain_ids),
        "variant_type": variant_type,
    }


# ===================================================================
# S4. Unified GWAS: SNPs + SVs
# ===================================================================

def combine_genotype_matrices(snp_geno: np.ndarray,
                              sv_matrices: list[dict],
                              snp_sample_ids: list[str],
                              verbose: bool = True) -> dict:
    """Combine SNP genotype matrix with SV matrices.

    Aligns all matrices to a common sample set and stacks vertically.

    Args:
        snp_geno: SNP genotype matrix (n_snps × n_samples).
        sv_matrices: list of SV dicts from load_*_matrix().
        snp_sample_ids: sample IDs for the SNP matrix.
        verbose: print progress.

    Returns:
        dict with combined matrix, variant IDs, variant types, sample IDs.
    """
    # Find shared samples across all matrices
    shared = set(snp_sample_ids)
    for sv in sv_matrices:
        shared &= set(sv["sample_ids"])
    shared = sorted(shared)

    if not shared:
        raise ValueError("No shared samples across SNP and SV matrices")

    # Align SNP matrix
    snp_idx = {s: i for i, s in enumerate(snp_sample_ids)}
    snp_cols = np.array([snp_idx[s] for s in shared])
    snp_aligned = snp_geno[:, snp_cols]

    # Stack: SNPs first, then each SV type
    matrices = [snp_aligned]
    variant_ids = [f"SNP_{i}" for i in range(snp_geno.shape[0])]
    variant_types = ["SNP"] * snp_geno.shape[0]

    for sv in sv_matrices:
        sv_idx = {s: i for i, s in enumerate(sv["sample_ids"])}
        sv_cols = np.array([sv_idx[s] for s in shared])
        sv_aligned = sv["dosage_matrix"][:, sv_cols]
        matrices.append(sv_aligned)
        variant_ids.extend(sv["variant_ids"])
        variant_types.extend([sv["variant_type"]] * sv["dosage_matrix"].shape[0])

    combined = np.vstack(matrices)

    if verbose:
        n_snp = snp_aligned.shape[0]
        n_sv = combined.shape[0] - n_snp
        print(f"\n  Combined matrix: {combined.shape[0]:,} variants × "
              f"{combined.shape[1]} samples")
        print(f"    SNPs: {n_snp:,}")
        for sv in sv_matrices:
            n = sv["dosage_matrix"].shape[0]
            print(f"    {sv['variant_type']}: {n:,}")
        print(f"    Total SVs: {n_sv:,}")

    return {
        "genotype_matrix": combined,
        "variant_ids": variant_ids,
        "variant_types": np.array(variant_types),
        "sample_ids": shared,
        "n_snps": snp_aligned.shape[0],
        "n_svs": combined.shape[0] - snp_aligned.shape[0],
    }


def fast_gwas(geno: np.ndarray, pheno: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Vectorized linear regression GWAS. ~1s for 90K variants × 1K samples."""
    n_var, n_sam = geno.shape

    G = geno.copy()
    # Vectorized imputation: replace NaN with row mean
    row_means = np.nanmean(G, axis=1)
    nan_mask = np.isnan(G)
    # Use broadcasting for imputation
    for i in range(n_var):
        if np.any(nan_mask[i]):
            G[i, nan_mask[i]] = row_means[i] if not np.isnan(row_means[i]) else 0.0

    Y = pheno - np.mean(pheno)
    n = len(Y)

    G_mean = np.mean(G, axis=1, keepdims=True)
    G_centered = G - G_mean
    G_var = np.sum(G_centered ** 2, axis=1)

    valid = G_var > 0
    betas = np.zeros(n_var)
    pvals = np.ones(n_var)

    if np.sum(valid) == 0:
        return betas, pvals

    cov_xy = G_centered[valid] @ Y
    betas[valid] = cov_xy / G_var[valid]

    residuals = Y[np.newaxis, :] - betas[valid, np.newaxis] * G_centered[valid]
    rss = np.sum(residuals ** 2, axis=1)
    mse = rss / (n - 2)
    se = np.sqrt(mse / G_var[valid])
    se[se == 0] = np.inf

    t_stat = betas[valid] / se
    pvals[valid] = 2 * sp_stats.t.sf(np.abs(t_stat), n - 2)

    return betas, pvals


def sv_gwas(combined: dict, pheno_dict: dict,
            verbose: bool = True) -> dict:
    """Run GWAS on combined SNP+SV matrix.

    Args:
        combined: output of combine_genotype_matrices().
        pheno_dict: {sample_id: phenotype_value}.
        verbose: print progress.

    Returns:
        dict with betas, pvals, per-variant-type breakdown.
    """
    # Align phenotype to combined samples
    shared = [s for s in combined["sample_ids"] if s in pheno_dict]
    if len(shared) < 50:
        return {"error": "Too few samples"}

    sam_idx = {s: i for i, s in enumerate(combined["sample_ids"])}
    cols = np.array([sam_idx[s] for s in shared])
    geno = combined["genotype_matrix"][:, cols]
    pheno = np.array([pheno_dict[s] for s in shared])

    betas, pvals = fast_gwas(geno, pheno)

    n_var = len(pvals)
    types = combined["variant_types"]
    bonf_threshold = 0.05 / n_var

    # Breakdown by variant type
    breakdown = {}
    for vtype in np.unique(types):
        mask = types == vtype
        n_total = int(np.sum(mask))
        n_bonf = int(np.sum(pvals[mask] < bonf_threshold))
        n_sug = int(np.sum(pvals[mask] < 1e-4))
        min_p = float(np.nanmin(pvals[mask])) if n_total > 0 else 1.0

        # Fraction significant (enrichment metric)
        frac_sig = n_sug / n_total if n_total > 0 else 0

        breakdown[vtype] = {
            "n_tested": n_total,
            "n_bonferroni": n_bonf,
            "n_suggestive": n_sug,
            "min_p": min_p,
            "fraction_significant": frac_sig,
        }

    if verbose:
        print(f"\n  SV-GWAS Results ({len(shared)} samples, {n_var:,} variants):")
        print(f"  {'Type':<12} {'Tested':>8} {'Bonf':>8} {'Sug':>8} "
              f"{'%Sig':>8} {'Min P':>12}")
        print(f"  {'-'*60}")
        for vtype, info in sorted(breakdown.items()):
            print(f"  {vtype:<12} {info['n_tested']:>8,} {info['n_bonferroni']:>8,} "
                  f"{info['n_suggestive']:>8,} {info['fraction_significant']:>8.1%} "
                  f"{info['min_p']:>12.2e}")

    return {
        "betas": betas,
        "pvals": pvals,
        "variant_types": types,
        "breakdown": breakdown,
        "n_samples": len(shared),
        "bonf_threshold": bonf_threshold,
    }


# ===================================================================
# S5. SV-Aware Heritability
# ===================================================================

def sv_heritability(combined: dict, pheno_dict: dict,
                    n_top: int = 50, verbose: bool = True) -> dict:
    """Estimate heritability from combined SNP+SV matrix.

    Compares: h² from SNPs only vs h² from SNPs+SVs.
    The difference = SV contribution to heritability.

    Args:
        combined: output of combine_genotype_matrices().
        pheno_dict: {sample_id: value}.
        n_top: top variants for PRS.
        verbose: print progress.

    Returns:
        dict with h2_snp_only, h2_snp_sv, h2_sv_contribution.
    """
    shared = [s for s in combined["sample_ids"] if s in pheno_dict]
    if len(shared) < 50:
        return {"error": "Too few samples"}

    sam_idx = {s: i for i, s in enumerate(combined["sample_ids"])}
    cols = np.array([sam_idx[s] for s in shared])
    geno = combined["genotype_matrix"][:, cols]
    pheno = np.array([pheno_dict[s] for s in shared])
    types = combined["variant_types"]
    n_sam = len(shared)

    rng = np.random.default_rng(42)
    train = np.arange(0, n_sam, 2)
    test = np.arange(1, n_sam, 2)

    def _prs_r2(geno_sub, pheno_vec, train_idx, test_idx, top_k):
        betas_train, pvals_train = fast_gwas(geno_sub[:, train_idx], pheno_vec[train_idx])
        top = np.argsort(pvals_train)[:top_k]

        G_test = geno_sub[:, test_idx].copy()
        row_means = np.nanmean(G_test, axis=1)
        for i in range(G_test.shape[0]):
            mask = np.isnan(G_test[i])
            if np.any(mask):
                G_test[i, mask] = row_means[i] if not np.isnan(row_means[i]) else 0.0

        prs = np.zeros(len(test_idx))
        for idx in top:
            prs += betas_train[idx] * G_test[idx]

        if np.std(prs) == 0 or np.std(pheno_vec[test_idx]) == 0:
            return 0.0
        return float(np.corrcoef(prs, pheno_vec[test_idx])[0, 1] ** 2)

    # h² from SNPs only
    snp_mask = types == "SNP"
    h2_snp = _prs_r2(geno[snp_mask], pheno, train, test, n_top)

    # h² from SNPs + all SVs
    h2_all = _prs_r2(geno, pheno, train, test, n_top)

    # h² from SVs only
    sv_mask = types != "SNP"
    if np.sum(sv_mask) > 0:
        h2_sv = _prs_r2(geno[sv_mask], pheno, train, test, min(n_top, int(np.sum(sv_mask))))
    else:
        h2_sv = 0.0

    sv_contribution = h2_all - h2_snp

    if verbose:
        print(f"  h² (SNPs only):   {h2_snp:.4f}")
        print(f"  h² (SVs only):    {h2_sv:.4f}")
        print(f"  h² (SNPs + SVs):  {h2_all:.4f}")
        print(f"  SV contribution:  {sv_contribution:+.4f} "
              f"({sv_contribution/max(h2_snp,0.001)*100:+.1f}%)")

    return {
        "h2_snp_only": h2_snp,
        "h2_sv_only": h2_sv,
        "h2_combined": h2_all,
        "h2_sv_contribution": sv_contribution,
        "h2_sv_pct_increase": float(sv_contribution / max(h2_snp, 0.001) * 100),
    }


# ===================================================================
# S6. SV Pleiotropy Analysis
# ===================================================================

def sv_pleiotropy(combined: dict, pheno_dicts: dict[str, dict],
                  p_threshold: float = 1e-4,
                  verbose: bool = True) -> dict:
    """Compare pleiotropy between SNPs and SVs.

    For each variant type, count how many traits each variant is
    significant for.  SVs should be more pleiotropic than SNPs
    (Loegler et al.: 2.82 vs 1.45 traits per QTL).

    Args:
        combined: output of combine_genotype_matrices().
        pheno_dicts: {trait_name: {sample_id: value}}.
        p_threshold: significance threshold.
        verbose: print progress.

    Returns:
        dict with per-type pleiotropy statistics.
    """
    types = combined["variant_types"]
    n_var = len(types)
    pleio_count = np.zeros(n_var, dtype=int)

    t0 = time.time()
    n_traits = 0
    for trait, pheno_dict in pheno_dicts.items():
        result = sv_gwas(combined, pheno_dict, verbose=False)
        if "error" in result:
            continue
        pleio_count += (result["pvals"] < p_threshold).astype(int)
        n_traits += 1

    # Per-type statistics
    type_stats = {}
    for vtype in np.unique(types):
        mask = types == vtype
        counts = pleio_count[mask]
        sig = counts > 0
        multi = counts >= 2

        type_stats[vtype] = {
            "n_variants": int(np.sum(mask)),
            "n_significant_any": int(np.sum(sig)),
            "n_pleiotropic": int(np.sum(multi)),
            "mean_traits_per_qtl": float(np.mean(counts[sig])) if np.sum(sig) > 0 else 0,
            "max_traits": int(np.max(counts)) if len(counts) > 0 else 0,
            "fraction_pleiotropic": float(np.sum(multi) / max(np.sum(sig), 1)),
        }

    if verbose:
        elapsed = time.time() - t0
        print(f"\n  SV Pleiotropy ({n_traits} traits, {elapsed:.0f}s):")
        print(f"  {'Type':<12} {'Variants':>10} {'Any sig':>10} "
              f"{'Pleiotropic':>12} {'Mean traits':>12} {'%Pleio':>8}")
        print(f"  {'-'*68}")
        for vtype, s in sorted(type_stats.items()):
            print(f"  {vtype:<12} {s['n_variants']:>10,} {s['n_significant_any']:>10,} "
                  f"{s['n_pleiotropic']:>12,} {s['mean_traits_per_qtl']:>12.2f} "
                  f"{s['fraction_pleiotropic']:>8.1%}")

    return {
        "type_stats": type_stats,
        "n_traits": n_traits,
        "pleiotropy_counts": pleio_count,
    }
