"""MPAT — Message-Passing Association Test (gene-level).

Aggregates per-variant z-scores within genes using graph topology (HAS_CONSEQUENCE),
corrects for LD correlation, and produces gene-level p-values with an analytical null.

Two modes:
- Directed MPAT: T = w^T z / sqrt(w^T Sigma w) ~ N(0,1). Preserves effect direction.
- Undirected (MAGMA-like): Q = z^T Sigma^{-1} z ~ weighted chi-squared.
"""

from __future__ import annotations

import numpy as np
from scipy import stats as sp_stats

from . import config as _cfg
from .db import GraphGWASConnection
from .genotype import (
    build_dosage,
    get_phenotype_indices,
    build_count_table,
    compute_z_score,
)


def build_gene_ld_matrix(conn: GraphGWASConnection, gene_symbol: str,
                         indices: np.ndarray) -> tuple[np.ndarray, list[dict]]:
    """Build within-gene LD correlation matrix from gt_packed.

    Args:
        conn: database connection.
        gene_symbol: gene symbol to query.
        indices: sample packed_indices to use for LD computation.

    Returns:
        (Sigma: k x k correlation matrix, variants: list of variant dicts with z_score).
    """
    # Get variants in this gene via HAS_CONSEQUENCE
    result = conn.execute_read(
        """
        MATCH (v:Variant)-[:HAS_CONSEQUENCE]->(g:Gene {symbol: $gene})
        RETURN v.variantId AS variantId, v.pos AS pos, v.gt_packed AS gt_packed,
               v.af_total AS af_total
        ORDER BY v.pos
        """,
        {"gene": gene_symbol},
    )
    variants = [dict(r) for r in result]

    if len(variants) < 2:
        return np.array([[1.0]]) if variants else np.array([[]]), variants

    k = len(variants)

    # Build dosage matrix (k variants x N samples)
    dosage_matrix = np.zeros((k, len(indices)))
    for i, v in enumerate(variants):
        if v["gt_packed"] is not None:
            d = build_dosage(v["gt_packed"], indices, _cfg.N_SAMPLES)
            # Replace NaN with mean for correlation computation
            mean_d = np.nanmean(d)
            d = np.where(np.isnan(d), mean_d, d)
            dosage_matrix[i] = d

    # Correlation matrix (k x k)
    # Use signed Pearson r (not r²) for directed MPAT
    Sigma = np.corrcoef(dosage_matrix)
    # Fix NaN correlations (monomorphic variants)
    Sigma = np.nan_to_num(Sigma, nan=0.0)
    np.fill_diagonal(Sigma, 1.0)

    return Sigma, variants


def mpat_gene_test(z_scores: np.ndarray, ld_matrix: np.ndarray,
                   method: str = "directed",
                   weights: np.ndarray | None = None) -> dict:
    """MPAT gene-level test.

    Args:
        z_scores: per-variant z-scores (length k).
        ld_matrix: k x k LD correlation matrix (signed Pearson r).
        method: 'directed' (signed sum, normal null) or 'undirected' (sum of squares, chi-sq null).
        weights: optional per-variant weights (default: 1/sqrt(k)).

    Returns:
        dict: {z_gene, p_value, n_variants, method}.
    """
    k = len(z_scores)
    if k == 0:
        return {"z_gene": 0.0, "p_value": 1.0, "n_variants": 0, "method": method}

    if weights is None:
        weights = np.ones(k) / np.sqrt(k)

    if method == "directed":
        # T = w^T z, Var(T) = w^T Sigma w
        T = weights @ z_scores
        var_T = weights @ ld_matrix @ weights
        if var_T <= 0:
            return {"z_gene": 0.0, "p_value": 1.0, "n_variants": k, "method": method}
        Z_gene = T / np.sqrt(var_T)
        p_value = float(2 * sp_stats.norm.sf(abs(Z_gene)))
        return {
            "z_gene": float(Z_gene), "p_value": p_value,
            "n_variants": k, "method": "mpat_directed",
        }

    elif method == "undirected":
        # Q = z^T z (without LD correction: overly anti-conservative)
        # With LD correction via Liu et al. (2009) moment-matching:
        # Q ~ a * chi2(d) + b
        eigenvalues = np.linalg.eigvalsh(ld_matrix)
        eigenvalues = eigenvalues[eigenvalues > 1e-6]

        if len(eigenvalues) == 0:
            return {"z_gene": 0.0, "p_value": 1.0, "n_variants": k, "method": method}

        Q = float(np.sum(z_scores ** 2))

        # Liu et al. (2009) moment-matching approximation
        c1 = np.sum(eigenvalues)
        c2 = np.sum(eigenvalues ** 2)
        c3 = np.sum(eigenvalues ** 3)
        c4 = np.sum(eigenvalues ** 4)

        s1 = c3 / (c2 ** 1.5)
        s2 = c4 / (c2 ** 2)

        if abs(s1) < 1e-10:
            # Approximate as scaled chi-squared
            a = c2 / c1
            d = c1 ** 2 / c2
            p_value = float(sp_stats.chi2.sf(Q / a, d))
        else:
            # Liu et al. (2009) approximation
            if s1 ** 2 > s2:
                a_val = 1.0 / (s1 - np.sqrt(s1 ** 2 - s2))
                delta = s1 * a_val ** 3 - a_val ** 2
                l_val = a_val ** 2 - 2 * delta
            else:
                a_val = 1.0 / s1
                delta = 0
                l_val = c2 ** 3 / c3 ** 2

            mu_Q = c1
            sigma_Q = np.sqrt(2 * c2)

            # Standardize Q and compute p-value
            Q_star = (Q - mu_Q) / sigma_Q
            Q_chi2 = Q_star * np.sqrt(2 * l_val) + l_val
            p_value = float(sp_stats.chi2.sf(max(Q_chi2, 0), l_val))

        return {
            "z_gene": float(np.sqrt(Q)), "p_value": p_value,
            "n_variants": k, "n_effective": float(len(eigenvalues)),
            "method": "mpat_undirected",
        }

    else:
        raise ValueError(f"Unknown MPAT method: {method}")


def mpat_scan(conn: GraphGWASConnection, chr: str | None = None,
              genes: list[str] | None = None,
              method: str = "directed",
              min_variants: int = 2,
              verbose: bool = True) -> list[dict]:
    """Run MPAT across genes in a region or gene list.

    Args:
        conn: database connection.
        chr: chromosome to scan (all genes on this chr).
        genes: specific gene symbols to test.
        method: 'directed' or 'undirected'.
        min_variants: minimum variants per gene to test.

    Returns:
        list of gene-level result dicts sorted by p-value.
    """
    case_idx, ctrl_idx = get_phenotype_indices(conn)
    all_idx = np.concatenate([case_idx, ctrl_idx])

    # Get gene list
    if genes:
        gene_list = genes
    elif chr:
        result = conn.execute_read(
            """
            MATCH (v:Variant)-[:HAS_CONSEQUENCE]->(g:Gene)
            WHERE v.chr = $chr
            RETURN DISTINCT g.symbol AS symbol
            ORDER BY symbol
            """,
            {"chr": chr},
        )
        gene_list = [r["symbol"] for r in result if r["symbol"]]
    else:
        raise ValueError("Provide either chr or genes parameter")

    if verbose:
        print(f"MPAT scan: {len(gene_list)} genes, method={method}")

    results = []
    for i, gene in enumerate(gene_list):
        # Build LD matrix and get variants
        Sigma, variants = build_gene_ld_matrix(conn, gene, all_idx)

        if len(variants) < min_variants:
            continue

        # Compute z-scores for these variants
        z_scores = []
        for v in variants:
            if v["gt_packed"] is None:
                z_scores.append(0.0)
                continue
            table = build_count_table(v["gt_packed"], case_idx, ctrl_idx, _cfg.N_SAMPLES)
            z = compute_z_score(table, len(case_idx), len(ctrl_idx))
            z_scores.append(z)

        z_arr = np.array(z_scores)

        # Run MPAT
        result = mpat_gene_test(z_arr, Sigma, method=method)
        result["gene"] = gene
        result["chr"] = chr or ""
        results.append(result)

        if verbose and (i + 1) % 100 == 0:
            print(f"  ...{i + 1}/{len(gene_list)} genes tested", flush=True)

    # Sort by p-value
    results.sort(key=lambda r: r["p_value"])

    if verbose:
        sig = sum(1 for r in results if r["p_value"] < 0.05 / max(len(results), 1))
        print(f"Done: {len(results)} genes tested, {sig} Bonferroni-significant")

    return results
