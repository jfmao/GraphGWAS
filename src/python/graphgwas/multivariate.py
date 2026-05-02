"""Multi-trait genetic architecture in a graph-native world.

Estimates genetic covariance, genetic correlation, bivariate heritability,
and the G-matrix by projecting multiple phenotypes onto the graph Fourier
basis of a sample similarity graph.  Decomposes cross-trait architecture
by biological scale (variant → gene → pathway), identifies pleiotropic
genes via direct graph queries, and provides spectral coherence between
trait pairs.

Classical tools (bivariate GREML, cross-trait LDSC) produce a single scalar
r_G.  Graph-native methods produce multi-resolution decomposition: which
pathways mediate the genetic correlation, which genes are pleiotropic, and
whether the correlation is driven by LD or biology.

Components:
    M1. spectral_genetic_covariance  — Cov_G, r_G, h²_biv from graph Fourier
    M2. g_matrix                     — T×T genetic (co)variance, one eigendecomp
    M3. multiresolution_correlation  — r_G at variant / gene / pathway scales
    M4. spectral_coherence           — frequency-resolved genetic vs env covariance
    M5. flow_covariance              — per-pathway genetic covariance via max-flow
    M6. pleiotropic_genes            — direct graph query for shared genes/pathways
    M7. multivariate_report          — integrated multi-method report
"""

from __future__ import annotations

import numpy as np

from .db import GraphGWASConnection


# ===================================================================
# Phenotype loading helpers
# ===================================================================

def load_trait_vectors(conn: GraphGWASConnection,
                      trait_names: list[str],
                      verbose: bool = True) -> tuple[np.ndarray, list[str], np.ndarray]:
    """Load multiple phenotype trait vectors aligned to a common sample set.

    Queries pheno_<trait> properties on Sample nodes.  Returns only samples
    that have non-null values for ALL requested traits.

    Args:
        conn: database connection.
        trait_names: list of trait names (without 'pheno_' prefix).
        verbose: print progress.

    Returns:
        (trait_matrix, trait_names, sample_indices) where:
        - trait_matrix: (N, T) float64 array, samples × traits
        - trait_names: ordered list of T trait names
        - sample_indices: (N,) int32 array of packed_index values
    """
    # Build RETURN clause for each trait
    safe_names = [t.replace(" ", "_").replace("-", "_") for t in trait_names]
    return_parts = [f"s.pheno_{sn} AS trait_{i}" for i, sn in enumerate(safe_names)]

    result = conn.execute_read(
        f"""
        MATCH (s:Sample)
        WHERE s.packed_index IS NOT NULL
        RETURN s.packed_index AS idx, {', '.join(return_parts)}
        ORDER BY s.packed_index
        """,
    )

    records = list(result)
    T = len(trait_names)

    # Filter to samples with all traits present
    valid_rows = []
    for rec in records:
        vals = []
        skip = False
        for i in range(T):
            v = rec[f"trait_{i}"]
            if v is None:
                skip = True
                break
            try:
                vals.append(float(v))
            except (ValueError, TypeError):
                skip = True
                break
        if not skip:
            valid_rows.append((rec["idx"], vals))

    if not valid_rows:
        if verbose:
            print(f"No samples found with all traits: {trait_names}")
        return np.empty((0, T)), trait_names, np.empty(0, dtype=np.int32)

    indices = np.array([r[0] for r in valid_rows], dtype=np.int32)
    matrix = np.array([r[1] for r in valid_rows], dtype=np.float64)

    if verbose:
        print(f"Loaded {len(indices)} samples × {T} traits: {trait_names}")

    return matrix, trait_names, indices


def _standardise(x: np.ndarray) -> np.ndarray:
    """Standardise to zero mean, unit variance.  Returns zeros if no variance."""
    s = np.std(x)
    if s == 0:
        return np.zeros_like(x)
    return (x - np.mean(x)) / s


# ===================================================================
# M1. Spectral Genetic Covariance / Correlation
# ===================================================================

def spectral_genetic_covariance(conn: GraphGWASConnection,
                                trait1: str, trait2: str,
                                af_threshold: float = 0.01,
                                n_components: int = 20,
                                chromosomes: list[str] | None = None,
                                verbose: bool = True) -> dict:
    """Genetic covariance and correlation between two traits.

    Projects both phenotypes onto the low-frequency eigenvectors of the
    sample similarity graph.  Their inner product in that subspace is the
    genetic covariance.

    Generalises bivariate GREML: GREML uses GRM eigenvectors; this uses
    eigenvectors of any graph (rare-variant, pathway-mediated, etc.).

    Args:
        conn: database connection.
        trait1, trait2: trait names (pheno_<trait> properties on Sample nodes).
        af_threshold: AF cutoff for rare-variant similarity graph.
        n_components: number of Laplacian eigenvectors.
        chromosomes: restrict to these chromosomes.
        verbose: print progress.

    Returns:
        dict with cov_g, r_g, h2_biv, per-component cross-spectrum,
        h2_trait1, h2_trait2.
    """
    from .spectral import build_sample_similarity_graph, spectral_decompose

    # Load both traits
    matrix, names, indices = load_trait_vectors(conn, [trait1, trait2],
                                                verbose=verbose)
    if matrix.shape[0] < 20:
        return _null_bivariate(trait1, trait2, "Too few samples with both traits")

    y1 = _standardise(matrix[:, 0])
    y2 = _standardise(matrix[:, 1])
    N = len(y1)

    # Build similarity graph
    if verbose:
        print(f"\nBuilding sample similarity graph (AF < {af_threshold})...")

    # We need the similarity graph for the samples that have BOTH traits.
    # Reuse the full graph builder, then slice to our sample indices.
    W_full = build_sample_similarity_graph(conn, af_threshold=af_threshold,
                                           chromosomes=chromosomes,
                                           verbose=verbose)

    # The full graph is built for all phenotyped samples (case+control).
    # Map our indices into the full graph's ordering.
    from .genotype import get_phenotype_indices
    case_idx, ctrl_idx = get_phenotype_indices(conn)
    all_pheno_idx = np.concatenate([case_idx, ctrl_idx])

    # Build index mapping: packed_index → position in full similarity matrix
    full_pos = {int(idx): pos for pos, idx in enumerate(all_pheno_idx)}
    our_positions = np.array([full_pos[int(i)] for i in indices
                              if int(i) in full_pos])

    if len(our_positions) < 20:
        return _null_bivariate(trait1, trait2, "Too few overlap with phenotyped samples")

    # Slice similarity matrix to our samples
    W = W_full[np.ix_(our_positions, our_positions)]
    # Also slice phenotypes to matching samples
    mask = np.array([int(i) in full_pos for i in indices])
    y1 = _standardise(y1[mask])
    y2 = _standardise(y2[mask])
    N = len(y1)

    # Eigendecompose
    n_comp = min(n_components, N - 1)
    eigenvalues, eigenvectors = spectral_decompose(W, n_components=n_comp,
                                                    verbose=verbose)

    # Graph Fourier coefficients for both traits
    coeff1 = eigenvectors.T @ y1  # (K,)
    coeff2 = eigenvectors.T @ y2  # (K,)

    # Cross-spectrum and auto-spectra
    cross_spectrum = coeff1 * coeff2       # per-component
    auto1 = coeff1 ** 2
    auto2 = coeff2 ** 2

    # Find genetic subspace cutoff via eigenvalue gap
    if len(eigenvalues) > 2:
        gaps = np.diff(eigenvalues)
        cutoff = int(np.argmax(gaps)) + 1
    else:
        cutoff = max(1, n_comp // 2)
    cutoff = max(1, min(cutoff, len(eigenvalues)))

    # Genetic quantities (low-frequency subspace)
    cov_g = float(np.sum(cross_spectrum[:cutoff]))
    var_g1 = float(np.sum(auto1[:cutoff]))
    var_g2 = float(np.sum(auto2[:cutoff]))

    # Genetic correlation
    denom = np.sqrt(var_g1 * var_g2)
    r_g = float(cov_g / denom) if denom > 0 else 0.0
    r_g = max(-1.0, min(1.0, r_g))  # clamp to [-1, 1]

    # Heritabilities
    total_var1 = float(np.sum(auto1))
    total_var2 = float(np.sum(auto2))
    h2_1 = float(var_g1 / total_var1) if total_var1 > 0 else 0.0
    h2_2 = float(var_g2 / total_var2) if total_var2 > 0 else 0.0

    # Bivariate heritability
    cov_p = float(np.cov(y1, y2)[0, 1])
    h2_biv = float(cov_g / cov_p) if abs(cov_p) > 1e-10 else 0.0

    if verbose:
        print("\n=== Spectral Genetic Covariance ===")
        print(f"  Traits: {trait1} × {trait2}")
        print(f"  N = {N} samples, K = {cutoff} genetic components")
        print(f"  Cov_G  = {cov_g:.6f}")
        print(f"  r_G    = {r_g:.4f}")
        print(f"  h²_biv = {h2_biv:.4f}")
        print(f"  h²({trait1}) = {h2_1:.4f}")
        print(f"  h²({trait2}) = {h2_2:.4f}")

    return {
        "method": "spectral",
        "trait1": trait1,
        "trait2": trait2,
        "cov_g": cov_g,
        "r_g": r_g,
        "h2_bivariate": h2_biv,
        "h2_trait1": h2_1,
        "h2_trait2": h2_2,
        "cov_p": cov_p,
        "cross_spectrum": cross_spectrum.tolist(),
        "eigenvalues": eigenvalues.tolist(),
        "genetic_cutoff": cutoff,
        "n_samples": N,
    }


# ===================================================================
# M2. G-Matrix
# ===================================================================

def g_matrix(conn: GraphGWASConnection,
             trait_names: list[str],
             af_threshold: float = 0.01,
             n_components: int = 20,
             chromosomes: list[str] | None = None,
             verbose: bool = True) -> dict:
    """Compute the T×T genetic variance-covariance matrix (G-matrix).

    From a single eigendecomposition of the sample similarity graph,
    project all T traits onto the low-frequency subspace.  The G-matrix
    is G = Ŷ_genetic^T × Ŷ_genetic, where Ŷ_genetic is K×T.

    O(N² + T×N×K) vs GCTA's O(T² × N² × M).

    Args:
        conn: database connection.
        trait_names: list of T trait names.
        af_threshold: AF cutoff for similarity graph.
        n_components: number of eigenvectors.
        chromosomes: restrict to these chromosomes.
        verbose: print progress.

    Returns:
        dict with g_matrix (T×T), genetic_correlations (T×T),
        heritabilities (T,), trait_names, eigenvalues.
    """
    from .spectral import build_sample_similarity_graph, spectral_decompose

    T = len(trait_names)
    if T < 2:
        raise ValueError("G-matrix requires at least 2 traits")

    # Load all traits
    matrix, names, indices = load_trait_vectors(conn, trait_names,
                                                verbose=verbose)
    if matrix.shape[0] < 20:
        return {"error": "Too few samples", "trait_names": trait_names}

    N = matrix.shape[0]

    # Standardise each trait
    Y = np.column_stack([_standardise(matrix[:, t]) for t in range(T)])

    # Build similarity graph
    W_full = build_sample_similarity_graph(conn, af_threshold=af_threshold,
                                           chromosomes=chromosomes,
                                           verbose=verbose)

    # Map to phenotyped sample ordering
    from .genotype import get_phenotype_indices
    case_idx, ctrl_idx = get_phenotype_indices(conn)
    all_pheno_idx = np.concatenate([case_idx, ctrl_idx])
    full_pos = {int(idx): pos for pos, idx in enumerate(all_pheno_idx)}
    our_positions = np.array([full_pos[int(i)] for i in indices
                              if int(i) in full_pos])

    if len(our_positions) < 20:
        return {"error": "Too few overlap", "trait_names": trait_names}

    W = W_full[np.ix_(our_positions, our_positions)]
    mask = np.array([int(i) in full_pos for i in indices])
    Y = np.column_stack([_standardise(Y[mask, t]) for t in range(T)])
    N = Y.shape[0]

    # Eigendecompose (single decomposition for all traits)
    n_comp = min(n_components, N - 1)
    eigenvalues, eigenvectors = spectral_decompose(W, n_components=n_comp,
                                                    verbose=verbose)

    # Project all traits onto eigenvectors: Ŷ = U^T Y  (K × T)
    Y_hat = eigenvectors.T @ Y  # (K, T)

    # Genetic subspace cutoff
    if len(eigenvalues) > 2:
        gaps = np.diff(eigenvalues)
        cutoff = int(np.argmax(gaps)) + 1
    else:
        cutoff = max(1, n_comp // 2)
    cutoff = max(1, min(cutoff, len(eigenvalues)))

    # G-matrix from genetic subspace
    Y_genetic = Y_hat[:cutoff, :]   # (K_genetic, T)
    G = Y_genetic.T @ Y_genetic     # (T, T)

    # Heritabilities (diagonal)
    Y_total = Y_hat                  # (K, T)
    total_var = np.sum(Y_total ** 2, axis=0)  # (T,)
    genetic_var = np.diag(G)                   # (T,)
    heritabilities = np.where(total_var > 0,
                              genetic_var / total_var,
                              0.0)

    # Genetic correlation matrix
    D = np.sqrt(np.maximum(genetic_var, 1e-10))
    R = G / np.outer(D, D)
    np.fill_diagonal(R, 1.0)
    R = np.clip(R, -1.0, 1.0)

    # G-matrix eigendecomposition (directions of max genetic variation)
    g_eigenvalues, g_eigenvectors = np.linalg.eigh(G)
    order = np.argsort(-g_eigenvalues)  # descending
    g_eigenvalues = g_eigenvalues[order]
    g_eigenvectors = g_eigenvectors[:, order]

    if verbose:
        print(f"\n=== G-Matrix ({T} traits, {N} samples) ===")
        print(f"  Genetic subspace: {cutoff} / {n_comp} components")
        print("\n  Heritabilities:")
        for t, h2 in zip(trait_names, heritabilities):
            print(f"    h²({t}) = {h2:.4f}")
        print("\n  Genetic Correlation Matrix:")
        header = "         " + "  ".join(f"{t[:8]:>8}" for t in trait_names)
        print(header)
        for i, t in enumerate(trait_names):
            row = f"  {t[:8]:<8}" + "  ".join(f"{R[i,j]:8.4f}"
                                               for j in range(T))
            print(row)
        print(f"\n  G-matrix eigenvalues: "
              f"{', '.join(f'{v:.4f}' for v in g_eigenvalues[:5])}")

    return {
        "method": "g_matrix",
        "trait_names": trait_names,
        "g_matrix": G.tolist(),
        "genetic_correlations": R.tolist(),
        "heritabilities": heritabilities.tolist(),
        "n_samples": N,
        "genetic_cutoff": cutoff,
        "g_eigenvalues": g_eigenvalues.tolist(),
        "g_eigenvectors": g_eigenvectors.tolist(),
        "eigenvalues": eigenvalues.tolist(),
    }


# ===================================================================
# M3. Multi-Resolution Genetic Correlation
# ===================================================================

def multiresolution_correlation(conn: GraphGWASConnection,
                                trait1: str, trait2: str,
                                af_threshold: float = 0.01,
                                n_components: int = 20,
                                verbose: bool = True) -> dict:
    """Genetic correlation at variant, gene, and pathway scales.

    Builds three similarity graphs (variant-level, gene-level, pathway-level)
    and computes r_G on each.  The differences decompose:

        r_G_variant − r_G_gene     = LD-driven correlation (not true pleiotropy)
        r_G_gene − r_G_pathway     = gene-level pleiotropy outside known pathways
        r_G_pathway                = pathway-mediated correlation

    Args:
        conn: database connection.
        trait1, trait2: trait names.
        af_threshold: AF cutoff.
        n_components: eigenvectors per graph.
        verbose: print progress.

    Returns:
        dict with r_g at each scale and decomposition.
    """
    from .spectral import build_sample_similarity_graph, spectral_decompose
    from .heritability import _build_gene_similarity, _build_pathway_similarity

    # Load traits
    matrix, names, indices = load_trait_vectors(conn, [trait1, trait2],
                                                verbose=verbose)
    if matrix.shape[0] < 20:
        return _null_bivariate(trait1, trait2, "Too few samples")

    y1 = _standardise(matrix[:, 0])
    y2 = _standardise(matrix[:, 1])

    def _compute_rg(W, y1, y2, n_comp, label):
        """Compute r_G from a similarity matrix."""
        if np.max(W) == 0:
            return 0.0
        k = min(n_comp, W.shape[0] - 1)
        if k < 2:
            return 0.0
        evals, evecs = spectral_decompose(W, n_components=k, verbose=False)
        c1 = evecs.T @ y1
        c2 = evecs.T @ y2
        # Cutoff
        if len(evals) > 2:
            cut = int(np.argmax(np.diff(evals))) + 1
        else:
            cut = max(1, k // 2)
        cut = max(1, min(cut, len(evals)))
        cov_g = np.sum(c1[:cut] * c2[:cut])
        var1 = np.sum(c1[:cut] ** 2)
        var2 = np.sum(c2[:cut] ** 2)
        denom = np.sqrt(var1 * var2)
        rg = float(cov_g / denom) if denom > 0 else 0.0
        if verbose:
            print(f"  r_G_{label} = {rg:.4f} (cutoff={cut})")
        return max(-1.0, min(1.0, rg))

    if verbose:
        print("\n=== Multi-Resolution Genetic Correlation ===")
        print(f"  Traits: {trait1} × {trait2}")

    # Variant-level
    if verbose:
        print("\n--- Variant-level ---")
    W_var = build_sample_similarity_graph(conn, af_threshold=af_threshold,
                                          verbose=verbose)
    # Slice to our samples
    from .genotype import get_phenotype_indices
    case_idx, ctrl_idx = get_phenotype_indices(conn)
    all_pheno = np.concatenate([case_idx, ctrl_idx])
    full_pos = {int(idx): pos for pos, idx in enumerate(all_pheno)}
    pos = np.array([full_pos[int(i)] for i in indices if int(i) in full_pos])
    mask = np.array([int(i) in full_pos for i in indices])
    y1m, y2m = _standardise(y1[mask]), _standardise(y2[mask])
    W_v = W_var[np.ix_(pos, pos)]
    rg_variant = _compute_rg(W_v, y1m, y2m, n_components, "variant")

    # Gene-level
    if verbose:
        print("\n--- Gene-level ---")
    W_gene = _build_gene_similarity(conn, af_threshold=af_threshold,
                                    verbose=verbose)
    W_g = W_gene[np.ix_(pos, pos)]
    rg_gene = _compute_rg(W_g, y1m, y2m, n_components, "gene")

    # Pathway-level
    if verbose:
        print("\n--- Pathway-level ---")
    W_pw = _build_pathway_similarity(conn, af_threshold=af_threshold,
                                     verbose=verbose)
    W_p = W_pw[np.ix_(pos, pos)]
    rg_pathway = _compute_rg(W_p, y1m, y2m, n_components, "pathway")

    # Decomposition
    rg_ld_driven = rg_variant - rg_gene
    rg_gene_specific = rg_gene - rg_pathway

    if verbose:
        print("\n--- Decomposition ---")
        print(f"  r_G_variant  = {rg_variant:.4f}")
        print(f"  r_G_gene     = {rg_gene:.4f}")
        print(f"  r_G_pathway  = {rg_pathway:.4f}")
        print("  ---")
        print(f"  LD-driven (variant − gene)     = {rg_ld_driven:.4f}")
        print(f"  Gene-specific (gene − pathway) = {rg_gene_specific:.4f}")
        print(f"  Pathway-mediated               = {rg_pathway:.4f}")

    return {
        "method": "multiresolution",
        "trait1": trait1,
        "trait2": trait2,
        "r_g_variant": rg_variant,
        "r_g_gene": rg_gene,
        "r_g_pathway": rg_pathway,
        "r_g_ld_driven": rg_ld_driven,
        "r_g_gene_specific": rg_gene_specific,
        "decomposition": {
            "ld_driven": rg_ld_driven,
            "gene_specific": rg_gene_specific,
            "pathway_mediated": rg_pathway,
        },
    }


# ===================================================================
# M4. Spectral Coherence
# ===================================================================

def spectral_coherence(conn: GraphGWASConnection,
                       trait1: str, trait2: str,
                       af_threshold: float = 0.01,
                       n_components: int = 20,
                       chromosomes: list[str] | None = None,
                       verbose: bool = True) -> dict:
    """Frequency-resolved genetic vs environmental covariance.

    Spectral coherence C(λ_k) = S₁₂²(λ_k) / (S₁₁(λ_k) × S₂₂(λ_k))
    at each graph frequency.

    High coherence at low frequencies = genetic correlation.
    High coherence at high frequencies = shared environmental factors.
    Low coherence = independent traits.

    Args:
        conn: database connection.
        trait1, trait2: trait names.
        af_threshold: AF cutoff.
        n_components: eigenvectors.
        chromosomes: restrict to these chromosomes.
        verbose: print progress.

    Returns:
        dict with per-frequency coherence, genetic vs environmental
        coherence summary, eigenvalues.
    """
    from .spectral import build_sample_similarity_graph, spectral_decompose

    matrix, names, indices = load_trait_vectors(conn, [trait1, trait2],
                                                verbose=verbose)
    if matrix.shape[0] < 20:
        return _null_bivariate(trait1, trait2, "Too few samples")

    y1 = _standardise(matrix[:, 0])
    y2 = _standardise(matrix[:, 1])

    # Build and slice similarity graph
    W_full = build_sample_similarity_graph(conn, af_threshold=af_threshold,
                                           chromosomes=chromosomes,
                                           verbose=verbose)
    from .genotype import get_phenotype_indices
    case_idx, ctrl_idx = get_phenotype_indices(conn)
    all_pheno = np.concatenate([case_idx, ctrl_idx])
    full_pos = {int(idx): pos for pos, idx in enumerate(all_pheno)}
    pos = np.array([full_pos[int(i)] for i in indices if int(i) in full_pos])
    mask = np.array([int(i) in full_pos for i in indices])
    y1 = _standardise(y1[mask])
    y2 = _standardise(y2[mask])
    W = W_full[np.ix_(pos, pos)]
    N = len(y1)

    n_comp = min(n_components, N - 1)
    eigenvalues, eigenvectors = spectral_decompose(W, n_components=n_comp,
                                                    verbose=verbose)

    c1 = eigenvectors.T @ y1
    c2 = eigenvectors.T @ y2

    # Per-frequency quantities
    cross = c1 * c2
    auto1 = c1 ** 2
    auto2 = c2 ** 2

    # Coherence (avoid division by zero)
    denom = auto1 * auto2
    coherence = np.where(denom > 1e-20, cross ** 2 / denom, 0.0)

    # Signed coherence (preserves direction)
    signed_coherence = np.where(denom > 1e-20, cross / np.sqrt(denom), 0.0)

    # Genetic vs environmental summary
    if len(eigenvalues) > 2:
        gaps = np.diff(eigenvalues)
        cutoff = int(np.argmax(gaps)) + 1
    else:
        cutoff = max(1, n_comp // 2)
    cutoff = max(1, min(cutoff, len(eigenvalues)))

    mean_genetic_coherence = float(np.mean(coherence[:cutoff]))
    mean_env_coherence = float(np.mean(coherence[cutoff:]))if cutoff < len(coherence) else 0.0
    mean_genetic_signed = float(np.mean(signed_coherence[:cutoff]))
    mean_env_signed = float(np.mean(signed_coherence[cutoff:])) if cutoff < len(coherence) else 0.0

    if verbose:
        print(f"\n=== Spectral Coherence: {trait1} × {trait2} ===")
        print(f"  N = {N}, K = {n_comp} components, cutoff = {cutoff}")
        print(f"  Mean genetic coherence  (λ < τ): {mean_genetic_coherence:.4f} "
              f"(signed: {mean_genetic_signed:+.4f})")
        print(f"  Mean environmental coh. (λ ≥ τ): {mean_env_coherence:.4f} "
              f"(signed: {mean_env_signed:+.4f})")
        print("\n  Per-component coherence (first 10):")
        for k in range(min(10, len(coherence))):
            tag = "G" if k < cutoff else "E"
            print(f"    λ_{k}={eigenvalues[k]:.4f}  C={coherence[k]:.4f}  "
                  f"signed={signed_coherence[k]:+.4f}  [{tag}]")

    return {
        "method": "spectral_coherence",
        "trait1": trait1,
        "trait2": trait2,
        "coherence": coherence.tolist(),
        "signed_coherence": signed_coherence.tolist(),
        "eigenvalues": eigenvalues.tolist(),
        "genetic_cutoff": cutoff,
        "mean_genetic_coherence": mean_genetic_coherence,
        "mean_env_coherence": mean_env_coherence,
        "mean_genetic_signed": mean_genetic_signed,
        "mean_env_signed": mean_env_signed,
        "n_samples": N,
    }


# ===================================================================
# M5. Flow Covariance Decomposition
# ===================================================================

def flow_covariance(conn: GraphGWASConnection,
                    trait1: str, trait2: str,
                    af_threshold: float = 0.05,
                    verbose: bool = True) -> dict:
    """Per-pathway genetic covariance via max-flow overlap.

    Builds flow networks for each trait separately, computes pathway flow
    vectors, then measures their cosine similarity as pathway-level
    genetic correlation.  Per-pathway contribution = f1_p × f2_p.

    Args:
        conn: database connection.
        trait1, trait2: trait names (must have been activated in turn).
        af_threshold: AF threshold for flow network.
        verbose: print progress.

    Note:
        Requires running phenotype activation for each trait in sequence.
        Flow vectors are computed from the current active phenotype.
        For a true two-trait analysis, run flow for each trait and
        pass the results to flow_covariance_from_vectors().

    Returns:
        dict with r_g_flow, per-pathway covariance, cosine similarity.
    """
    if verbose:
        print(f"\n=== Flow Covariance: {trait1} × {trait2} ===")
        print("  NOTE: requires sequential phenotype activation.")
        print("  Use flow_covariance_from_vectors() with pre-computed flows.")

    return {
        "method": "flow_covariance",
        "trait1": trait1,
        "trait2": trait2,
        "note": "Use flow_covariance_from_vectors() with pre-computed flow results.",
    }


def flow_covariance_from_vectors(flows_trait1: list[dict],
                                 flows_trait2: list[dict],
                                 trait1: str = "trait1",
                                 trait2: str = "trait2",
                                 verbose: bool = True) -> dict:
    """Compute pathway-level genetic covariance from pre-computed flow results.

    Args:
        flows_trait1: output of compute_pathway_flows() for trait 1.
        flows_trait2: output of compute_pathway_flows() for trait 2.
        trait1, trait2: trait names for labelling.
        verbose: print progress.

    Returns:
        dict with r_g_flow (cosine similarity), per-pathway contributions.
    """
    # Build flow vectors indexed by pathway
    f1 = {r["pathway"]: r["flow_value"] for r in flows_trait1}
    f2 = {r["pathway"]: r["flow_value"] for r in flows_trait2}

    # Union of pathways
    all_pathways = sorted(set(f1.keys()) | set(f2.keys()))
    if not all_pathways:
        return _null_bivariate(trait1, trait2, "No pathways in flow results")

    v1 = np.array([f1.get(p, 0.0) for p in all_pathways])
    v2 = np.array([f2.get(p, 0.0) for p in all_pathways])

    # Cosine similarity = pathway-level genetic correlation
    norm1 = np.linalg.norm(v1)
    norm2 = np.linalg.norm(v2)
    if norm1 > 0 and norm2 > 0:
        r_g_flow = float(np.dot(v1, v2) / (norm1 * norm2))
    else:
        r_g_flow = 0.0

    # Per-pathway contribution
    total_product = float(np.dot(v1, v2))
    contributions = []
    for p, a, b in zip(all_pathways, v1, v2):
        product = float(a * b)
        contributions.append({
            "pathway": p,
            "flow_trait1": float(a),
            "flow_trait2": float(b),
            "covariance_contribution": product,
            "fraction": float(product / total_product) if abs(total_product) > 1e-10 else 0.0,
        })

    contributions.sort(key=lambda c: abs(c["covariance_contribution"]), reverse=True)

    # Shared vs unique pathways
    shared = [p for p in all_pathways if f1.get(p, 0) > 0 and f2.get(p, 0) > 0]
    unique_t1 = [p for p in all_pathways if f1.get(p, 0) > 0 and f2.get(p, 0) == 0]
    unique_t2 = [p for p in all_pathways if f1.get(p, 0) == 0 and f2.get(p, 0) > 0]

    if verbose:
        print(f"\n=== Flow Covariance: {trait1} × {trait2} ===")
        print(f"  r_G_flow (cosine) = {r_g_flow:.4f}")
        print(f"  Pathways: {len(shared)} shared, "
              f"{len(unique_t1)} {trait1}-only, {len(unique_t2)} {trait2}-only")
        print("  Top pathway contributions:")
        for c in contributions[:5]:
            print(f"    {c['pathway']}: cov={c['covariance_contribution']:.4f} "
                  f"({c['fraction']*100:.1f}%)")

    return {
        "method": "flow_covariance",
        "trait1": trait1,
        "trait2": trait2,
        "r_g_flow": r_g_flow,
        "pathway_contributions": contributions,
        "n_shared_pathways": len(shared),
        "n_unique_trait1": len(unique_t1),
        "n_unique_trait2": len(unique_t2),
    }


# ===================================================================
# M6. Pleiotropic Gene / Pathway Discovery
# ===================================================================

def pleiotropic_genes(conn: GraphGWASConnection,
                      trait1_run_id: str, trait2_run_id: str,
                      p_threshold: float = 5e-8,
                      verbose: bool = True) -> dict:
    """Find genes with significant associations for both traits.

    Directly queries the graph for genes where AssociationResult nodes
    from both traits pass through via FOR_VARIANT → HAS_CONSEQUENCE.

    Distinguishes:
    - Biological pleiotropy: same variant → same gene → both traits
    - Mediated pleiotropy: different variants → same gene → both traits
    - Pathway pleiotropy: different genes → same pathway → both traits

    Args:
        conn: database connection.
        trait1_run_id: GWAS run ID for trait 1.
        trait2_run_id: GWAS run ID for trait 2.
        p_threshold: significance threshold.
        verbose: print progress.

    Returns:
        dict with pleiotropic genes, pathways, and pleiotropy type.
    """
    # Find shared genes (variant → gene)
    result = conn.execute_read(
        """
        MATCH (ar1:AssociationResult)-[:FOR_VARIANT]->(v1:Variant)
              -[:HAS_CONSEQUENCE]->(g:Gene)
        WHERE ar1.run_id = $run1 AND ar1.p_value < $pval
        WITH g, collect(DISTINCT v1.variantId) AS variants_t1
        MATCH (ar2:AssociationResult)-[:FOR_VARIANT]->(v2:Variant)
              -[:HAS_CONSEQUENCE]->(g)
        WHERE ar2.run_id = $run2 AND ar2.p_value < $pval
        WITH g, variants_t1, collect(DISTINCT v2.variantId) AS variants_t2
        OPTIONAL MATCH (g)-[:IN_PATHWAY]->(p:Pathway)
        RETURN g.symbol AS gene,
               variants_t1, variants_t2,
               collect(DISTINCT p.name) AS pathways
        ORDER BY size(variants_t1) + size(variants_t2) DESC
        """,
        {"run1": trait1_run_id, "run2": trait2_run_id, "pval": p_threshold},
    )

    genes = []
    shared_variants_total = set()
    for rec in result:
        v1_set = set(rec["variants_t1"])
        v2_set = set(rec["variants_t2"])
        shared = v1_set & v2_set
        shared_variants_total.update(shared)

        if shared:
            pleiotropy_type = "biological"  # same variant hits both
        else:
            pleiotropy_type = "mediated"    # different variants, same gene

        genes.append({
            "gene": rec["gene"],
            "pleiotropy_type": pleiotropy_type,
            "n_variants_trait1": len(v1_set),
            "n_variants_trait2": len(v2_set),
            "n_shared_variants": len(shared),
            "shared_variants": sorted(shared),
            "pathways": [p for p in rec["pathways"] if p],
        })

    # Find shared pathways (gene → pathway)
    result_pw = conn.execute_read(
        """
        MATCH (ar1:AssociationResult)-[:FOR_VARIANT]->(:Variant)
              -[:HAS_CONSEQUENCE]->(:Gene)-[:IN_PATHWAY]->(p:Pathway)
        WHERE ar1.run_id = $run1 AND ar1.p_value < $pval
        WITH p, collect(DISTINCT ar1.variant_id) AS v1s
        MATCH (ar2:AssociationResult)-[:FOR_VARIANT]->(:Variant)
              -[:HAS_CONSEQUENCE]->(:Gene)-[:IN_PATHWAY]->(p)
        WHERE ar2.run_id = $run2 AND ar2.p_value < $pval
        RETURN p.name AS pathway, size(v1s) AS n_v1,
               count(DISTINCT ar2.variant_id) AS n_v2
        ORDER BY n_v1 + n_v2 DESC
        """,
        {"run1": trait1_run_id, "run2": trait2_run_id, "pval": p_threshold},
    )

    pathways = [{"pathway": r["pathway"],
                 "n_variants_trait1": r["n_v1"],
                 "n_variants_trait2": r["n_v2"]}
                for r in result_pw if r["pathway"]]

    n_biological = sum(1 for g in genes if g["pleiotropy_type"] == "biological")
    n_mediated = sum(1 for g in genes if g["pleiotropy_type"] == "mediated")

    if verbose:
        print("\n=== Pleiotropic Genes ===")
        print(f"  Runs: {trait1_run_id} × {trait2_run_id}")
        print(f"  p < {p_threshold:.0e}")
        print(f"  Pleiotropic genes: {len(genes)} "
              f"({n_biological} biological, {n_mediated} mediated)")
        print(f"  Shared pathways: {len(pathways)}")
        print(f"  Shared variants: {len(shared_variants_total)}")
        for g in genes[:10]:
            print(f"    {g['gene']} ({g['pleiotropy_type']}): "
                  f"{g['n_variants_trait1']}+{g['n_variants_trait2']} variants, "
                  f"{g['n_shared_variants']} shared")

    return {
        "method": "pleiotropic_genes",
        "pleiotropic_genes": genes,
        "pleiotropic_pathways": pathways,
        "n_biological_pleiotropy": n_biological,
        "n_mediated_pleiotropy": n_mediated,
        "n_shared_variants": len(shared_variants_total),
        "n_shared_pathways": len(pathways),
    }


# ===================================================================
# M6b. Gene Pleiotropy Strength Index (PSI)
# ===================================================================

def gene_pleiotropy_score(conn: GraphGWASConnection,
                          run_ids: list[str],
                          trait_labels: list[str] | None = None,
                          p_threshold: float = 5e-8,
                          verbose: bool = True) -> dict:
    """Compute a continuous Pleiotropy Strength Index (PSI) for every gene.

    PSI(gene) = breadth × magnitude × mechanism_bonus

    Where:
        breadth          = n_traits_significant / T
        magnitude        = mean(-log10(p)) across significant traits
        mechanism_bonus  = 1 + shared_variant_ratio
                           (biological pleiotropy through shared variants
                            gets up to 2× boost over mediated pleiotropy)

    High PSI = gene strongly affects many traits through shared biology.

    Args:
        conn: database connection.
        run_ids: list of GWAS run IDs (one per trait).
        trait_labels: optional human-readable labels (defaults to run_ids).
        p_threshold: significance threshold per trait.
        verbose: print progress.

    Returns:
        dict with ranked gene list, each containing PSI score and
        per-trait breakdown.
    """
    T = len(run_ids)
    if T < 2:
        raise ValueError("PSI requires at least 2 GWAS runs")

    if trait_labels is None:
        trait_labels = run_ids

    if verbose:
        print("=== Gene Pleiotropy Strength Index ===")
        print(f"  {T} traits, p < {p_threshold:.0e}")

    # Query: for each gene, gather per-run stats
    # Build a UNION query across all runs
    run_list = ", ".join(f"'{r}'" for r in run_ids)
    result = conn.execute_read(
        f"""
        MATCH (ar:AssociationResult)-[:FOR_VARIANT]->(v:Variant)
              -[:HAS_CONSEQUENCE]->(g:Gene)
        WHERE ar.run_id IN [{run_list}] AND ar.p_value < $pval
        WITH g.symbol AS gene, ar.run_id AS run_id,
             collect(DISTINCT v.variantId) AS variants,
             min(ar.p_value) AS best_p,
             avg(ar.p_value_log10) AS mean_log10p,
             count(DISTINCT ar) AS n_hits
        RETURN gene, run_id, variants, best_p, mean_log10p, n_hits
        ORDER BY gene, run_id
        """,
        {"pval": p_threshold},
    )

    # Aggregate per gene across traits
    gene_data: dict[str, dict] = {}
    for rec in result:
        gene = rec["gene"]
        run_id = rec["run_id"]
        if gene not in gene_data:
            gene_data[gene] = {
                "per_trait": {},
                "all_variants": set(),
            }
        variants = set(rec["variants"])
        gene_data[gene]["per_trait"][run_id] = {
            "variants": variants,
            "best_p": rec["best_p"],
            "mean_log10p": rec["mean_log10p"],
            "n_hits": rec["n_hits"],
        }
        gene_data[gene]["all_variants"].update(variants)

    # Compute PSI for each gene
    scored_genes = []
    for gene, data in gene_data.items():
        per_trait = data["per_trait"]
        n_traits_sig = len(per_trait)

        # Breadth: fraction of traits where gene is significant
        breadth = n_traits_sig / T

        # Magnitude: mean -log10(best_p) across significant traits
        log10ps = [info["mean_log10p"] for info in per_trait.values()
                   if info["mean_log10p"] is not None]
        magnitude = float(np.mean(log10ps)) if log10ps else 0.0

        # Mechanism bonus: proportion of variants shared across traits
        # For each pair of traits, count shared variants
        trait_variants = [info["variants"] for info in per_trait.values()]
        if len(trait_variants) >= 2:
            # Union and intersection across all trait pairs
            all_union = set()
            all_intersect_count = 0
            n_pairs = 0
            for i in range(len(trait_variants)):
                for j in range(i + 1, len(trait_variants)):
                    shared = trait_variants[i] & trait_variants[j]
                    union = trait_variants[i] | trait_variants[j]
                    all_intersect_count += len(shared)
                    all_union.update(union)
                    n_pairs += 1
            total_variants = len(all_union)
            avg_shared = all_intersect_count / n_pairs if n_pairs > 0 else 0
            shared_ratio = avg_shared / total_variants if total_variants > 0 else 0
        else:
            shared_ratio = 0.0
            total_variants = len(data["all_variants"])

        mechanism_bonus = 1.0 + shared_ratio  # range [1.0, 2.0]

        # PSI = breadth × magnitude × mechanism_bonus
        psi = breadth * magnitude * mechanism_bonus

        # Determine dominant pleiotropy type
        if shared_ratio > 0.5:
            pleiotropy_type = "biological"
        elif shared_ratio > 0:
            pleiotropy_type = "mixed"
        else:
            pleiotropy_type = "mediated"

        # Per-trait detail
        trait_detail = []
        for rid, label in zip(run_ids, trait_labels):
            info = per_trait.get(rid)
            if info:
                trait_detail.append({
                    "trait": label,
                    "run_id": rid,
                    "best_p": info["best_p"],
                    "mean_log10p": info["mean_log10p"],
                    "n_variants": len(info["variants"]),
                    "significant": True,
                })
            else:
                trait_detail.append({
                    "trait": label,
                    "run_id": rid,
                    "significant": False,
                })

        scored_genes.append({
            "gene": gene,
            "psi": float(psi),
            "breadth": float(breadth),
            "magnitude": float(magnitude),
            "mechanism_bonus": float(mechanism_bonus),
            "shared_variant_ratio": float(shared_ratio),
            "pleiotropy_type": pleiotropy_type,
            "n_traits_significant": n_traits_sig,
            "n_total_variants": total_variants,
            "per_trait": trait_detail,
        })

    # Sort by PSI descending
    scored_genes.sort(key=lambda g: g["psi"], reverse=True)

    # Assign ranks
    for i, g in enumerate(scored_genes):
        g["rank"] = i + 1

    if verbose:
        print(f"  Genes with ≥1 significant trait: {len(scored_genes)}")
        multi = sum(1 for g in scored_genes if g["n_traits_significant"] >= 2)
        print(f"  Genes significant for ≥2 traits: {multi}")
        print("\n  Top genes by PSI:")
        print(f"  {'Rank':<6}{'Gene':<15}{'PSI':<10}{'Breadth':<10}"
              f"{'Magnitude':<12}{'Mech':<8}{'Type':<12}{'Traits'}")
        print(f"  {'-'*85}")
        for g in scored_genes[:20]:
            traits_str = ",".join(t["trait"][:6] for t in g["per_trait"]
                                  if t["significant"])
            print(f"  {g['rank']:<6}{g['gene']:<15}{g['psi']:<10.3f}"
                  f"{g['breadth']:<10.2f}{g['magnitude']:<12.2f}"
                  f"{g['mechanism_bonus']:<8.2f}{g['pleiotropy_type']:<12}"
                  f"{traits_str}")

    return {
        "method": "gene_pleiotropy_score",
        "genes": scored_genes,
        "n_genes_scored": len(scored_genes),
        "n_multi_trait": sum(1 for g in scored_genes
                             if g["n_traits_significant"] >= 2),
        "n_traits": T,
        "trait_labels": trait_labels,
        "p_threshold": p_threshold,
    }


# ===================================================================
# M7. Integrated Report
# ===================================================================

def multivariate_report(conn: GraphGWASConnection,
                        trait1: str, trait2: str,
                        af_threshold: float = 0.01,
                        n_components: int = 20,
                        verbose: bool = True) -> dict:
    """Run spectral covariance + coherence for two traits.

    Args:
        conn: database connection.
        trait1, trait2: trait names.
        af_threshold: AF cutoff.
        n_components: eigenvectors.
        verbose: print progress.

    Returns:
        Unified dict with all estimates.
    """
    if verbose:
        print("=" * 60)
        print(f"  Multivariate Report: {trait1} × {trait2}")
        print("=" * 60)

    estimates = {}

    # M1: Spectral genetic covariance
    if verbose:
        print("\n" + "-" * 60)
    estimates["spectral"] = spectral_genetic_covariance(
        conn, trait1, trait2, af_threshold=af_threshold,
        n_components=n_components, verbose=verbose,
    )

    # M4: Spectral coherence
    if verbose:
        print("\n" + "-" * 60)
    estimates["coherence"] = spectral_coherence(
        conn, trait1, trait2, af_threshold=af_threshold,
        n_components=n_components, verbose=verbose,
    )

    # Summary
    if verbose:
        print("\n" + "=" * 60)
        print("  Summary")
        print("=" * 60)
        sc = estimates["spectral"]
        coh = estimates["coherence"]
        print(f"  r_G                = {sc.get('r_g', 0):.4f}")
        print(f"  Cov_G              = {sc.get('cov_g', 0):.6f}")
        print(f"  h²_biv             = {sc.get('h2_bivariate', 0):.4f}")
        print(f"  h²({trait1})       = {sc.get('h2_trait1', 0):.4f}")
        print(f"  h²({trait2})       = {sc.get('h2_trait2', 0):.4f}")
        print(f"  Genetic coherence  = {coh.get('mean_genetic_coherence', 0):.4f}")
        print(f"  Env. coherence     = {coh.get('mean_env_coherence', 0):.4f}")

    return {
        "estimates": estimates,
        "summary": {
            "r_g": estimates["spectral"].get("r_g", 0),
            "cov_g": estimates["spectral"].get("cov_g", 0),
            "h2_bivariate": estimates["spectral"].get("h2_bivariate", 0),
            "h2_trait1": estimates["spectral"].get("h2_trait1", 0),
            "h2_trait2": estimates["spectral"].get("h2_trait2", 0),
            "genetic_coherence": estimates["coherence"].get("mean_genetic_coherence", 0),
            "env_coherence": estimates["coherence"].get("mean_env_coherence", 0),
        },
    }


# ===================================================================
# Helpers
# ===================================================================

def _null_bivariate(trait1: str, trait2: str, reason: str = "") -> dict:
    return {
        "trait1": trait1,
        "trait2": trait2,
        "r_g": 0.0,
        "cov_g": 0.0,
        "error": reason or "Insufficient data",
    }
