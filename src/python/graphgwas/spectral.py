"""Graph spectral phenotype decomposition.

Generalizes PCA-based population structure correction by decomposing
the phenotype signal in the frequency domain of a rare-variant-based
sample similarity graph. Low-frequency components capture genetic signal;
high-frequency components are noise/environment.

PCA uses eigenvectors of the common-variant GRM (a special case).
Spectral filtering generalizes to rare-variant similarity, pathway-based
similarity, and any graph topology.
"""

from __future__ import annotations

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import eigsh

from . import config as _cfg
from .db import GraphGWASConnection
from .genotype import (
    build_carrier_set,
    get_all_indices,
    get_phenotype_indices,
    variant_iterator,
)


def build_sample_similarity_graph(conn: GraphGWASConnection,
                                  af_threshold: float = 0.01,
                                  chromosomes: list[str] | None = None,
                                  verbose: bool = True) -> np.ndarray:
    """Build sample-sample similarity matrix from shared rare variant carriage.

    For each rare variant (AF < threshold), the carrier set defines which
    samples share this variant. The outer product of carrier sets, accumulated
    over all rare variants, gives the similarity matrix.

    This captures rare-variant-based population structure that common-variant
    PCA misses.

    Returns: (N, N) float64 similarity matrix where N = number of non-excluded samples.
    """
    all_idx = get_all_indices(conn)
    n = len(all_idx)

    if verbose:
        print(f"Building rare-variant similarity matrix: {n} samples, AF < {af_threshold}")

    similarity = np.zeros((n, n), dtype=np.float64)

    if chromosomes is None:
        result = conn.execute_read(
            "MATCH (v:Variant) RETURN DISTINCT v.chr AS c ORDER BY c"
        )
        chromosomes = [r["c"] for r in result]

    n_variants_used = 0
    for chrom in chromosomes:
        if verbose:
            print(f"  Processing {chrom}...", end="", flush=True)

        chrom_count = 0
        for v in variant_iterator(conn, chrom):
            af = v.get("af_total", 1.0)
            if af >= af_threshold or af <= 0:
                continue
            gt_packed = v["gt_packed"]
            if gt_packed is None:
                continue

            # Carrier set as float for outer product
            carriers = build_carrier_set(gt_packed, all_idx, _cfg.N_SAMPLES).astype(np.float64)
            n_carriers = np.sum(carriers)
            if n_carriers < 2:
                continue

            # Weight by 1/AF (upweight ultra-rare variants)
            weight = 1.0 / af

            # Rank-1 update: similarity += weight * carriers * carriers^T
            # Vectorized: update entire K×K submatrix at once
            carrier_idx = np.where(carriers > 0)[0]
            similarity[np.ix_(carrier_idx, carrier_idx)] += weight

            chrom_count += 1

        n_variants_used += chrom_count
        if verbose:
            print(f" {chrom_count} rare variants")

    # Normalize: divide by number of variants used
    if n_variants_used > 0:
        similarity /= n_variants_used

    if verbose:
        print(f"Similarity matrix: {n}×{n}, {n_variants_used} rare variants used")

    return similarity


def spectral_decompose(similarity_matrix: np.ndarray,
                       n_components: int = 20,
                       normalized: bool = True,
                       verbose: bool = True) -> tuple[np.ndarray, np.ndarray]:
    """Compute graph Laplacian eigendecomposition.

    Args:
        similarity_matrix: (N, N) similarity matrix.
        n_components: number of eigenvectors to compute.
        normalized: if True, use normalized Laplacian L_sym = I - D^{-1/2} W D^{-1/2}.

    Returns:
        (eigenvalues[k], eigenvectors[N, k]) — smallest eigenvalues first.
    """
    N = similarity_matrix.shape[0]
    n_components = min(n_components, N - 1)

    # Degree matrix
    D = np.diag(np.sum(similarity_matrix, axis=1))

    if normalized:
        # Normalized Laplacian: L_sym = I - D^{-1/2} W D^{-1/2}
        D_inv_sqrt = np.diag(1.0 / np.sqrt(np.maximum(np.diag(D), 1e-10)))
        L = np.eye(N) - D_inv_sqrt @ similarity_matrix @ D_inv_sqrt
    else:
        # Unnormalized Laplacian: L = D - W
        L = D - similarity_matrix

    # Eigendecomposition (smallest eigenvalues = smoothest graph signals)
    L_sparse = csr_matrix(L)
    eigenvalues, eigenvectors = eigsh(L_sparse, k=n_components, which="SM")

    # Sort by eigenvalue (ascending = smoothest first)
    order = np.argsort(eigenvalues)
    eigenvalues = eigenvalues[order]
    eigenvectors = eigenvectors[:, order]

    if verbose:
        print(f"Spectral decomposition: {n_components} components")
        print(f"  Eigenvalue range: [{eigenvalues[0]:.4f}, {eigenvalues[-1]:.4f}]")
        # Spectral gap: difference between 1st and 2nd eigenvalue
        if len(eigenvalues) > 1:
            print(f"  Spectral gap (λ₂ - λ₁): {eigenvalues[1] - eigenvalues[0]:.4f}")

    return eigenvalues, eigenvectors


def filter_phenotype_signal(phenotype: np.ndarray,
                            eigenvectors: np.ndarray,
                            eigenvalues: np.ndarray,
                            n_keep: int = 10,
                            verbose: bool = True) -> np.ndarray:
    """Decompose phenotype into graph frequencies, keep low-frequency signal.

    Low-frequency components (small eigenvalues) = smooth on the graph =
    genetic signal. High-frequency components = noise/environment.

    Args:
        phenotype: (N,) phenotype vector (1=case, 0=control or continuous).
        eigenvectors: (N, k) from spectral_decompose.
        eigenvalues: (k,) from spectral_decompose.
        n_keep: number of low-frequency components to retain.

    Returns:
        (N,) filtered phenotype vector (denoised genetic signal).
    """
    n_keep = min(n_keep, eigenvectors.shape[1])

    # Project phenotype onto eigenvectors
    V = eigenvectors[:, :n_keep]  # (N, n_keep)
    coefficients = V.T @ phenotype  # (n_keep,)

    # Reconstruct from low-frequency components only
    filtered = V @ coefficients  # (N,)

    if verbose:
        # Variance explained
        total_var = np.var(phenotype)
        filtered_var = np.var(filtered)
        residual_var = np.var(phenotype - filtered)
        pct_explained = (filtered_var / total_var * 100) if total_var > 0 else 0
        print(f"Spectral filtering: kept {n_keep} components")
        print(f"  Variance explained: {pct_explained:.1f}%")
        print(f"  Original var: {total_var:.4f}, Filtered var: {filtered_var:.4f}")

    return filtered


def spectral_gwas_correction(conn: GraphGWASConnection,
                             n_components: int = 10,
                             af_threshold: float = 0.01,
                             chromosomes: list[str] | None = None,
                             verbose: bool = True) -> dict:
    """Full spectral correction pipeline.

    1. Build rare-variant similarity matrix
    2. Eigendecompose the Laplacian
    3. Filter the phenotype signal
    4. Return filtered phenotype + spectral covariates (eigenvectors)

    The spectral covariates can be used as alternatives to PCA in
    association tests, capturing rare-variant-based structure.
    """
    # Build similarity
    sim = build_sample_similarity_graph(
        conn, af_threshold=af_threshold,
        chromosomes=chromosomes, verbose=verbose
    )

    # Eigendecompose
    eigenvalues, eigenvectors = spectral_decompose(
        sim, n_components=n_components, verbose=verbose
    )

    # Get phenotype (supports both binary and quantitative)
    from .genotype import get_phenotype_values
    all_idx = get_all_indices(conn)
    phenotype = get_phenotype_values(conn, all_idx)
    phenotype = np.where(np.isnan(phenotype), np.nanmean(phenotype), phenotype)

    # Filter
    filtered_phenotype = filter_phenotype_signal(
        phenotype, eigenvectors, eigenvalues,
        n_keep=n_components, verbose=verbose
    )

    return {
        "similarity_matrix": sim,
        "eigenvalues": eigenvalues,
        "eigenvectors": eigenvectors,
        "phenotype_original": phenotype,
        "phenotype_filtered": filtered_phenotype,
        "spectral_covariates": eigenvectors,  # use columns as covariates
        "n_components": n_components,
    }
