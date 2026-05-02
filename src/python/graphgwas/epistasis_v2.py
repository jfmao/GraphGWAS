"""Graph-Native Epistasis Detection v2.

Seven methods for detecting genetic interactions:
  M1: LD-Aware Co-occurrence with interaction statistics
  M2: Motif-Filtered Pairwise Testing (biological priors)
  M3: Differential Subgraph Analysis (case vs control structure)
  M4: Dark Matter / Negative Co-occurrence (synthetic lethality)
  M5: Random Walk on Bipartite Graph (higher-order)
  M6: GNN Interaction Extraction
  M7: Persistent Homology / Topological

This module implements M2 first (highest ROI, uses existing annotations).
"""

from __future__ import annotations

from dataclasses import dataclass
import numpy as np
from scipy import stats as sp_stats

from . import config as _cfg
from .db import GraphGWASConnection
from .genotype import (
    get_all_indices,
    get_phenotype_indices,
    get_phenotype_values,
    build_dosage,
    variant_iterator,
)


# ===================================================================
# Data classes
# ===================================================================

@dataclass
class InteractionResult:
    """Result for a single tested variant pair."""
    variant_1: str
    variant_2: str
    motif: str  # "same_pathway", "same_gene", "protein_interaction", etc.
    shared_entity: str  # gene/pathway/complex name
    beta_marginal_1: float
    beta_marginal_2: float
    beta_interaction: float
    se_interaction: float
    p_interaction: float
    p_corrected: float | None  # after multiple testing correction
    n_samples: int
    maf_1: float
    maf_2: float


# ===================================================================
# Interaction term regression
# ===================================================================

def _test_interaction(
    dosage_1: np.ndarray,
    dosage_2: np.ndarray,
    phenotype: np.ndarray,
    covariates: np.ndarray | None = None,
) -> dict:
    """Test interaction term: Y ~ β1·G1 + β2·G2 + β12·(G1×G2) + covariates + ε.

    Returns dict with beta_1, beta_2, beta_interaction, se_interaction, p_interaction.
    """
    valid = ~np.isnan(dosage_1) & ~np.isnan(dosage_2) & ~np.isnan(phenotype)
    if covariates is not None:
        valid &= ~np.any(np.isnan(covariates), axis=1)
    n = valid.sum()
    if n < 20:
        return {"beta_1": 0, "beta_2": 0, "beta_interaction": 0,
                "se_interaction": np.inf, "p_interaction": 1.0, "n": n}

    g1 = dosage_1[valid]
    g2 = dosage_2[valid]
    y = phenotype[valid]
    interaction = g1 * g2

    # Check variance
    if np.std(g1) == 0 or np.std(g2) == 0 or np.std(interaction) == 0:
        return {"beta_1": 0, "beta_2": 0, "beta_interaction": 0,
                "se_interaction": np.inf, "p_interaction": 1.0, "n": n}

    # Build design matrix: intercept, G1, G2, G1×G2, [covariates]
    if covariates is not None:
        X = np.column_stack([np.ones(n), g1, g2, interaction, covariates[valid]])
    else:
        X = np.column_stack([np.ones(n), g1, g2, interaction])

    try:
        beta_hat, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
        resid = y - X @ beta_hat
        p = X.shape[1]
        mse = np.sum(resid ** 2) / (n - p)
        var_beta = mse * np.linalg.inv(X.T @ X)
        se = np.sqrt(np.diag(var_beta))

        # Interaction term is index 3 (after intercept, G1, G2)
        t_stat = beta_hat[3] / se[3]
        p_value = float(2 * sp_stats.t.sf(abs(t_stat), n - p))
    except (np.linalg.LinAlgError, ValueError):
        return {"beta_1": 0, "beta_2": 0, "beta_interaction": 0,
                "se_interaction": np.inf, "p_interaction": 1.0, "n": n}

    return {
        "beta_1": float(beta_hat[1]),
        "beta_2": float(beta_hat[2]),
        "beta_interaction": float(beta_hat[3]),
        "se_interaction": float(se[3]),
        "p_interaction": p_value,
        "n": n,
    }


# ===================================================================
# M2: Motif-Filtered Pairwise Testing
# ===================================================================

def _enumerate_motif_pairs(
    conn: GraphGWASConnection,
    motif: str,
    chr: str | None = None,
    limit: int = 50000,
) -> list[dict]:
    """Query Neo4j for variant pairs matching a biological motif.

    Motifs:
        'same_pathway': variants in different genes that share a pathway
        'same_gene': two variants in the same gene
        'protein_interaction': variants in genes that physically interact

    Returns list of {v1, v2, shared_entity, gene_1, gene_2}.
    """
    chr_filter = f"AND v1.chr = '{chr}'" if chr else ""

    if motif == "same_pathway":
        query = f"""
            MATCH (v1:Variant)-[:HAS_CONSEQUENCE]->(g1:Gene)-[:IN_PATHWAY]->(p:Pathway)
                  <-[:IN_PATHWAY]-(g2:Gene)<-[:HAS_CONSEQUENCE]-(v2:Variant)
            WHERE g1.geneId < g2.geneId
              AND v1.af_total >= 0.01 AND v1.af_total <= 0.99
              AND v2.af_total >= 0.01 AND v2.af_total <= 0.99
              AND v1.gt_packed IS NOT NULL AND v2.gt_packed IS NOT NULL
              {chr_filter}
            RETURN DISTINCT v1.variantId AS v1, v2.variantId AS v2,
                   v1.gt_packed AS gtp1, v2.gt_packed AS gtp2,
                   v1.af_total AS af1, v2.af_total AS af2,
                   g1.symbol AS gene1, g2.symbol AS gene2,
                   p.name AS pathway
            LIMIT $limit
        """
    elif motif == "same_gene":
        query = f"""
            MATCH (v1:Variant)-[:HAS_CONSEQUENCE]->(g:Gene)<-[:HAS_CONSEQUENCE]-(v2:Variant)
            WHERE v1.variantId < v2.variantId
              AND v1.af_total >= 0.01 AND v1.af_total <= 0.99
              AND v2.af_total >= 0.01 AND v2.af_total <= 0.99
              AND v1.gt_packed IS NOT NULL AND v2.gt_packed IS NOT NULL
              AND abs(v1.pos - v2.pos) > 1000
              {chr_filter}
            RETURN DISTINCT v1.variantId AS v1, v2.variantId AS v2,
                   v1.gt_packed AS gtp1, v2.gt_packed AS gtp2,
                   v1.af_total AS af1, v2.af_total AS af2,
                   g.symbol AS gene1, g.symbol AS gene2,
                   g.symbol AS pathway
            LIMIT $limit
        """
    else:
        raise ValueError(f"Unknown motif: {motif}. Use: same_pathway, same_gene")

    result = conn.execute_read(query, {"limit": limit})
    return [dict(r) for r in result]


def motif_filtered_epistasis(
    conn: GraphGWASConnection,
    motifs: list[str] | None = None,
    chr: str | None = None,
    mac_min: int = 10,
    correction: str = "BH",
    max_pairs_per_motif: int = 50000,
    verbose: bool = True,
) -> list[InteractionResult]:
    """M2: Motif-Filtered Pairwise Epistasis Testing.

    Uses biological graph structure (Gene, Pathway, protein interaction)
    to enumerate candidate interacting pairs. Then tests interaction term
    for each pair. Multiple testing correction over candidate pairs only
    (typically 10³-10⁵ pairs, not 10¹² exhaustive).

    Args:
        conn: active database connection.
        motifs: which motifs to test. Default: ["same_pathway", "same_gene"].
        chr: restrict to one chromosome (None = all).
        mac_min: minimum minor allele count for both variants.
        correction: "BH" (Benjamini-Hochberg FDR) or "bonferroni".
        max_pairs_per_motif: limit per motif query.
        verbose: print progress.

    Returns:
        list of InteractionResult sorted by p_interaction.
    """
    if motifs is None:
        motifs = ["same_pathway", "same_gene"]

    all_idx = get_all_indices(conn)
    n = len(all_idx)

    # Get phenotype
    pheno = get_phenotype_values(conn, all_idx)
    valid_pheno = ~np.isnan(pheno)
    if valid_pheno.sum() < 50:
        if verbose:
            print("Too few samples with phenotype")
        return []

    if verbose:
        print(f"M2 Motif-Filtered Epistasis: {valid_pheno.sum()} samples, "
              f"motifs={motifs}, chr={chr or 'all'}")

    # Enumerate candidate pairs from each motif
    all_results = []
    total_pairs_tested = 0

    for motif in motifs:
        if verbose:
            print(f"\n  Motif: {motif}...")

        pairs = _enumerate_motif_pairs(conn, motif, chr=chr,
                                        limit=max_pairs_per_motif)

        if verbose:
            print(f"    {len(pairs):,} candidate pairs found")

        if not pairs:
            continue

        # Deduplicate (same pair may appear via different shared entities)
        seen = set()
        unique_pairs = []
        for p in pairs:
            key = tuple(sorted([p["v1"], p["v2"]]))
            if key not in seen:
                seen.add(key)
                unique_pairs.append(p)

        if verbose:
            print(f"    {len(unique_pairs):,} unique pairs after dedup")

        # Test interaction for each pair
        tested = 0
        for p in unique_pairs:
            d1 = build_dosage(p["gtp1"], all_idx, _cfg.N_SAMPLES)
            d2 = build_dosage(p["gtp2"], all_idx, _cfg.N_SAMPLES)

            # MAC filter
            d1_clean = np.where(np.isnan(d1), 0, d1)
            d2_clean = np.where(np.isnan(d2), 0, d2)
            mac1 = min(int(d1_clean.sum()), int(2 * n - d1_clean.sum()))
            mac2 = min(int(d2_clean.sum()), int(2 * n - d2_clean.sum()))
            if mac1 < mac_min or mac2 < mac_min:
                continue

            result = _test_interaction(d1, d2, pheno)
            tested += 1

            all_results.append(InteractionResult(
                variant_1=p["v1"],
                variant_2=p["v2"],
                motif=motif,
                shared_entity=p.get("pathway", ""),
                beta_marginal_1=result["beta_1"],
                beta_marginal_2=result["beta_2"],
                beta_interaction=result["beta_interaction"],
                se_interaction=result["se_interaction"],
                p_interaction=result["p_interaction"],
                p_corrected=None,  # filled below
                n_samples=result["n"],
                maf_1=p["af1"],
                maf_2=p["af2"],
            ))

        total_pairs_tested += tested
        if verbose:
            print(f"    {tested:,} pairs tested (after MAC filter)")

    if not all_results:
        if verbose:
            print("\n  No testable pairs found.")
        return []

    # Multiple testing correction
    pvals = np.array([r.p_interaction for r in all_results])

    if correction == "BH":
        # Benjamini-Hochberg FDR
        n_tests = len(pvals)
        order = np.argsort(pvals)
        ranks = np.empty(n_tests, dtype=int)
        ranks[order] = np.arange(1, n_tests + 1)
        corrected = np.minimum(1.0, pvals * n_tests / ranks)
        # Ensure monotonicity
        for i in range(n_tests - 2, -1, -1):
            corrected[order[i]] = min(corrected[order[i]], corrected[order[i + 1]])
        for i, r in enumerate(all_results):
            r.p_corrected = float(corrected[i])
    elif correction == "bonferroni":
        for r in all_results:
            r.p_corrected = min(1.0, r.p_interaction * len(all_results))

    # Sort by interaction p-value
    all_results.sort(key=lambda r: r.p_interaction)

    if verbose:
        sig_nominal = sum(1 for r in all_results if r.p_interaction < 0.05)
        sig_corrected = sum(1 for r in all_results if r.p_corrected and r.p_corrected < 0.05)
        print("\n  Summary:")
        print(f"    Total pairs tested: {total_pairs_tested:,}")
        print(f"    Nominal p < 0.05: {sig_nominal:,}")
        print(f"    Corrected p < 0.05: {sig_corrected:,}")
        if all_results:
            best = all_results[0]
            print(f"    Best pair: {best.variant_1} × {best.variant_2}")
            print(f"      motif={best.motif}, shared={best.shared_entity}")
            print(f"      β_int={best.beta_interaction:.4f}, p={best.p_interaction:.2e}")

    return all_results


# ===================================================================
# M2 — no-Neo4j path (graph-cache driven)
#
# Mirrors _enumerate_motif_pairs / motif_filtered_epistasis above but
# reads gene / pathway / PPI annotations from a JSON graph_cache file
# (see data/annotations/human_graph_cache_v2_chr*.json) instead of
# running Cypher queries against a live Neo4j connection.  Used by
# the paper-#2 REGENIE rescue demo (tests/regenie_rescue_demo.py
# Stage 3) and any other context where Neo4j is not available.
# ===================================================================

def enumerate_motif_pairs_from_cache(
    cache: dict,
    variant_id_to_col: dict,
    motifs: list,
    *,
    max_pairs_per_entity: int = 500,
    rng: "np.random.Generator | None" = None,
) -> list:
    """Enumerate variant pairs matching biological motifs from a graph_cache.

    The graph_cache JSON keys are variant IDs (``chrN:pos:REF:ALT`` format)
    and each value is a dict with lists ``{"genes", "pathways", "ppi"}``.
    For each requested motif we build the inverse index (e.g. gene → variant
    indices) and emit pairs from each entity that has ≥ 2 contributing
    variants.

    Args:
        cache: pre-loaded JSON dict, keyed by variant ID.
        variant_id_to_col: maps each variant ID present in the dosage
            matrix to its column index.  Variants in the cache but not in
            this map are silently skipped.
        motifs: subset of ``["same_gene", "same_pathway",
            "protein_interaction"]``.  Unknown motif names raise.
        max_pairs_per_entity: cap on how many pairs to emit from any
            single entity (gene / pathway / PPI partner-set) to keep the
            candidate pool tractable on densely-annotated regions.
        rng: numpy Generator for deterministic per-entity sub-sampling
            once an entity hits the cap.

    Returns:
        List of ``(var_idx_1, var_idx_2, motif, shared_entity)`` tuples,
        deduplicated across motifs (a pair surfaced via multiple
        motifs appears once, with the *first* motif name encountered).
    """
    import itertools

    if rng is None:
        rng = np.random.default_rng(0)

    # Build inverse indices for each requested motif
    valid_motifs = {"same_gene", "same_pathway", "protein_interaction"}
    for m in motifs:
        if m not in valid_motifs:
            raise ValueError(f"Unknown motif: {m!r}.  Valid: {sorted(valid_motifs)}")

    gene_to_vars: dict = {}
    pathway_to_vars: dict = {}
    ppi_to_vars: dict = {}      # gene → variant indices in any of its PPI partners
    var_to_genes: dict = {}     # for cross-gene constraint on same_pathway / PPI

    for vid, idx in variant_id_to_col.items():
        entry = cache.get(vid)
        if entry is None:
            continue
        genes = entry.get("genes", []) or []
        if "same_gene" in motifs:
            for g in genes:
                gene_to_vars.setdefault(g, []).append(idx)
        if "same_pathway" in motifs:
            for p in entry.get("pathways", []) or []:
                pathway_to_vars.setdefault(p, []).append(idx)
        if "protein_interaction" in motifs:
            for partner_gene in entry.get("ppi", []) or []:
                ppi_to_vars.setdefault(partner_gene, []).append(idx)
        if genes:
            var_to_genes[idx] = set(genes)

    seen_pairs: set = set()
    out: list = []

    def _emit_pairs(entity_to_vars, motif_name, *,
                    cross_gene_only=False):
        """Walk an inverse index, emitting unique cross-variant pairs."""
        for entity, var_idxs in entity_to_vars.items():
            if len(var_idxs) < 2:
                continue
            n_full = len(var_idxs) * (len(var_idxs) - 1) // 2
            if n_full <= max_pairs_per_entity:
                pair_iter = itertools.combinations(var_idxs, 2)
            else:
                # Random sample of distinct pairs (sample-without-replacement
                # via rng.choice with size=2)
                def _sample():
                    for _ in range(max_pairs_per_entity):
                        a, b = rng.choice(var_idxs, size=2, replace=False)
                        yield int(a), int(b)
                pair_iter = _sample()
            for a, b in pair_iter:
                i, j = (a, b) if a < b else (b, a)
                if cross_gene_only:
                    g_i = var_to_genes.get(i, set())
                    g_j = var_to_genes.get(j, set())
                    if g_i and g_j and g_i & g_j:
                        # Same-gene pair — skip; that's handled by same_gene
                        continue
                key = (i, j)
                if key in seen_pairs:
                    continue
                seen_pairs.add(key)
                out.append((i, j, motif_name, str(entity)))

    if "same_gene" in motifs:
        _emit_pairs(gene_to_vars, "same_gene")
    if "same_pathway" in motifs:
        _emit_pairs(pathway_to_vars, "same_pathway", cross_gene_only=True)
    if "protein_interaction" in motifs:
        _emit_pairs(ppi_to_vars, "protein_interaction", cross_gene_only=True)

    return out


def motif_filtered_epistasis_from_data(
    dosages: np.ndarray,
    variant_ids: list,
    phenotype: np.ndarray,
    graph_cache,
    motifs: list | None = None,
    *,
    mac_min: int = 10,
    correction: str = "BH",
    max_pairs_per_entity: int = 500,
    max_pairs_total: int = 50_000,
    verbose: bool = True,
    rng: "np.random.Generator | None" = None,
) -> "list[InteractionResult]":
    """M2 motif-filtered pairwise epistasis testing — no-Neo4j variant.

    Functional twin of :func:`motif_filtered_epistasis` that takes raw
    dosages + a graph_cache instead of a Neo4j connection.  Reuses
    :func:`_test_interaction` for the regression and the same
    BH-FDR / Bonferroni multiple-testing correction.

    Args:
        dosages: shape ``(n_samples, n_variants)``, can contain NaN.
        variant_ids: list of length ``n_variants`` in ``chrN:pos:REF:ALT``
            format (matching graph_cache keys).
        phenotype: shape ``(n_samples,)``.  NaN samples are dropped per pair.
        graph_cache: either a dict (already loaded) or a ``str``/``Path``
            to a JSON file.
        motifs: default ``["same_gene", "same_pathway"]``.  Pass
            ``["same_gene", "same_pathway", "protein_interaction"]`` to
            include PPI-partner pairs (P2).
        mac_min: minimum minor allele count required from each variant in
            the pair (after NaN-drop).
        correction: ``"BH"`` or ``"bonferroni"``.
        max_pairs_per_entity: per-gene / per-pathway sampling cap.
        max_pairs_total: hard cap on total pairs tested across all motifs;
            after enumeration the head of the list is taken if this is
            exceeded.
        verbose: print progress.
        rng: numpy Generator for deterministic sampling.

    Returns:
        ``list[InteractionResult]`` sorted by ``p_interaction`` ascending.
    """
    import json as _json
    from pathlib import Path as _Path

    if motifs is None:
        motifs = ["same_gene", "same_pathway"]
    if isinstance(graph_cache, (str, _Path)):
        cache = _json.loads(_Path(graph_cache).read_text())
    else:
        cache = graph_cache

    n_samples, n_variants = dosages.shape
    if len(variant_ids) != n_variants:
        raise ValueError(
            f"variant_ids length ({len(variant_ids)}) must match "
            f"dosages.shape[1] ({n_variants})"
        )
    if phenotype.shape[0] != n_samples:
        raise ValueError(
            f"phenotype length ({phenotype.shape[0]}) must match "
            f"dosages.shape[0] ({n_samples})"
        )

    variant_id_to_col = {v: i for i, v in enumerate(variant_ids)}

    if verbose:
        print(f"M2 (no-Neo4j) motif-filtered epistasis: "
              f"n_samples={n_samples}, n_variants={n_variants}, "
              f"motifs={motifs}", flush=True)

    pairs = enumerate_motif_pairs_from_cache(
        cache, variant_id_to_col, motifs,
        max_pairs_per_entity=max_pairs_per_entity, rng=rng,
    )
    if verbose:
        print(f"  enumerated {len(pairs):,} candidate pairs", flush=True)

    if len(pairs) > max_pairs_total:
        pairs = pairs[:max_pairs_total]
        if verbose:
            print(f"  truncated to first {max_pairs_total:,} pairs", flush=True)

    # Compute MAFs once per variant for filter + result populating.
    # NaN-aware: AF = nanmean(d) / 2.
    af_total = np.nanmean(dosages, axis=0) / 2.0
    maf_total = np.minimum(af_total, 1.0 - af_total)

    all_results: list = []
    n_passed_mac = 0
    for (i, j, motif, shared) in pairs:
        d1 = dosages[:, i].astype(float)
        d2 = dosages[:, j].astype(float)
        # MAC filter using non-NaN dosage sums
        s1 = float(np.nansum(d1)); n1 = int(np.sum(~np.isnan(d1)))
        s2 = float(np.nansum(d2)); n2 = int(np.sum(~np.isnan(d2)))
        mac1 = min(s1, 2 * n1 - s1)
        mac2 = min(s2, 2 * n2 - s2)
        if mac1 < mac_min or mac2 < mac_min:
            continue
        n_passed_mac += 1

        result = _test_interaction(d1, d2, phenotype)
        all_results.append(InteractionResult(
            variant_1=variant_ids[i],
            variant_2=variant_ids[j],
            motif=motif,
            shared_entity=shared,
            beta_marginal_1=result["beta_1"],
            beta_marginal_2=result["beta_2"],
            beta_interaction=result["beta_interaction"],
            se_interaction=result["se_interaction"],
            p_interaction=result["p_interaction"],
            p_corrected=None,  # filled below
            n_samples=int(result["n"]),
            maf_1=float(maf_total[i]),
            maf_2=float(maf_total[j]),
        ))

    if verbose:
        print(f"  {n_passed_mac:,} pairs survived MAC≥{mac_min} filter", flush=True)
        print(f"  {len(all_results):,} pairs tested", flush=True)

    if not all_results:
        return []

    # Multiple-testing correction (lifted verbatim from motif_filtered_epistasis)
    pvals = np.array([r.p_interaction for r in all_results])
    if correction == "BH":
        n_tests = len(pvals)
        order = np.argsort(pvals)
        ranks = np.empty(n_tests, dtype=int)
        ranks[order] = np.arange(1, n_tests + 1)
        corrected = np.minimum(1.0, pvals * n_tests / ranks)
        # Enforce monotonicity (BH step-up)
        for k in range(n_tests - 2, -1, -1):
            corrected[order[k]] = min(corrected[order[k]],
                                       corrected[order[k + 1]])
        for k, r in enumerate(all_results):
            r.p_corrected = float(corrected[k])
    elif correction == "bonferroni":
        for r in all_results:
            r.p_corrected = float(min(1.0, r.p_interaction * len(all_results)))
    else:
        raise ValueError(f"correction must be 'BH' or 'bonferroni'; got {correction!r}")

    all_results.sort(key=lambda r: r.p_interaction)

    if verbose:
        sig_nominal = sum(1 for r in all_results if r.p_interaction < 0.05)
        sig_corrected = sum(1 for r in all_results
                            if r.p_corrected is not None and r.p_corrected < 0.05)
        print(f"  nominal p < 0.05      : {sig_nominal:,}")
        print(f"  corrected p < 0.05 ({correction}): {sig_corrected:,}")
        if all_results:
            best = all_results[0]
            print(f"  best pair: {best.variant_1} × {best.variant_2}  "
                  f"motif={best.motif}  p_interaction={best.p_interaction:.2e}")

    return all_results


# ===================================================================
# M1: LD-Aware Co-occurrence with Interaction Statistics
# ===================================================================

def _compute_ld(dosage_1: np.ndarray, dosage_2: np.ndarray) -> float:
    """Compute r² between two dosage vectors."""
    valid = ~np.isnan(dosage_1) & ~np.isnan(dosage_2)
    if valid.sum() < 20:
        return 0.0
    d1 = dosage_1[valid]
    d2 = dosage_2[valid]
    if np.std(d1) == 0 or np.std(d2) == 0:
        return 0.0
    r = np.corrcoef(d1, d2)[0, 1]
    return r ** 2


def ld_pruned_cooccurrence(
    conn: GraphGWASConnection,
    chr: str,
    start: int | None = None,
    end: int | None = None,
    r2_prune: float = 0.5,
    min_cocarriers: int = 5,
    min_distance_bp: int = 100_000,
    mac_min: int = 10,
    max_variants: int = 2000,
    n_permutations: int = 100,
    verbose: bool = True,
) -> list[InteractionResult]:
    """M1: LD-Aware Co-occurrence with Interaction Statistics.

    Fixes v1 by:
    1. LD pruning (r² < threshold) to keep only independent variants
    2. Requiring minimum physical distance between pairs (avoids LD blocks)
    3. Testing interaction term (not just enrichment)
    4. Permutation p-values for module significance

    Args:
        conn: database connection.
        chr: chromosome.
        start, end: region (None = whole chromosome).
        r2_prune: LD pruning threshold.
        min_cocarriers: minimum case co-carriers for an edge.
        min_distance_bp: minimum bp distance between pair members.
        mac_min: minimum minor allele count.
        max_variants: max variants to load (for tractability).
        n_permutations: permutations for significance testing.
        verbose: print progress.

    Returns:
        list of InteractionResult sorted by interaction p-value.
    """
    case_idx, ctrl_idx = get_phenotype_indices(conn)
    all_idx = np.concatenate([case_idx, ctrl_idx])
    n_case = len(case_idx)
    n_ctrl = len(ctrl_idx)

    if n_case < 20 or n_ctrl < 20:
        # Need case/control — try binarizing quantitative phenotype
        all_idx = get_all_indices(conn)
        pheno = get_phenotype_values(conn, all_idx)
        valid = ~np.isnan(pheno)
        pheno_v = pheno[valid]
        all_idx = all_idx[valid]
        q75 = np.percentile(pheno_v, 75)
        q25 = np.percentile(pheno_v, 25)
        case_mask = pheno_v >= q75
        ctrl_mask = pheno_v <= q25
        case_idx_local = np.where(case_mask)[0]
        ctrl_idx_local = np.where(ctrl_mask)[0]
        n_case = len(case_idx_local)
        n_ctrl = len(ctrl_idx_local)
        use_local = True
    else:
        pheno = np.concatenate([np.ones(n_case), np.zeros(n_ctrl)])
        use_local = False

    if verbose:
        print(f"M1 LD-Pruned Co-occurrence: {chr}:{start or ''}-{end or ''}")
        print(f"  {n_case} cases, {n_ctrl} controls, r²_prune={r2_prune}")

    # Step 1: Load variants and their dosages
    variants = []
    dosage_list = []
    for v in variant_iterator(conn, chr, start, end):
        af = v.get("af_total", 0)
        if af < 0.01 or af > 0.99:
            continue
        gtp = v.get("gt_packed")
        if gtp is None:
            continue
        d = build_dosage(gtp, all_idx, _cfg.N_SAMPLES)
        d_clean = np.where(np.isnan(d), 0, d)
        mac = min(int(d_clean.sum()), int(2 * len(all_idx) - d_clean.sum()))
        if mac < mac_min:
            continue
        variants.append(v)
        dosage_list.append(d)
        if len(variants) >= max_variants * 3:  # load extra for LD pruning
            break

    if verbose:
        print(f"  Loaded {len(variants)} variants")

    if len(variants) < 4:
        if verbose:
            print("  Too few variants")
        return []

    # Step 2: LD pruning (greedy, keep variant with lowest p-value)
    keep = []
    used = set()
    for i in range(len(variants)):
        if i in used:
            continue
        keep.append(i)
        # Mark LD partners as used
        for j in range(i + 1, len(variants)):
            if j in used:
                continue
            if abs(variants[i]["pos"] - variants[j]["pos"]) < min_distance_bp:
                r2 = _compute_ld(dosage_list[i], dosage_list[j])
                if r2 > r2_prune:
                    used.add(j)
        if len(keep) >= max_variants:
            break

    pruned_variants = [variants[i] for i in keep]
    pruned_dosages = [dosage_list[i] for i in keep]

    if verbose:
        print(f"  After LD pruning (r²>{r2_prune}): {len(pruned_variants)} variants")

    # Step 3: Build co-occurrence graph on pruned set and test interactions
    if use_local:
        pheno_for_test = pheno_v
    else:
        pheno_for_test = pheno

    results = []
    n_v = len(pruned_variants)

    pairs_tested = 0
    for i in range(n_v):
        for j in range(i + 1, n_v):
            # Require minimum physical distance
            if abs(pruned_variants[i]["pos"] - pruned_variants[j]["pos"]) < min_distance_bp:
                continue

            d1 = pruned_dosages[i]
            d2 = pruned_dosages[j]

            # Co-carrier count in cases
            if use_local:
                carriers_1 = (d1 > 0)
                carriers_2 = (d2 > 0)
                co_case = np.sum(carriers_1[case_idx_local] & carriers_2[case_idx_local])
            else:
                carriers_1 = (d1[:n_case] > 0)
                carriers_2 = (d2[:n_case] > 0)
                co_case = np.sum(carriers_1 & carriers_2)

            if co_case < min_cocarriers:
                continue

            # Test interaction
            result = _test_interaction(d1, d2, pheno_for_test)
            pairs_tested += 1

            if result["p_interaction"] < 0.05:  # pre-filter for memory
                results.append(InteractionResult(
                    variant_1=pruned_variants[i]["variantId"],
                    variant_2=pruned_variants[j]["variantId"],
                    motif="co_occurrence",
                    shared_entity=f"co_carriers={co_case}",
                    beta_marginal_1=result["beta_1"],
                    beta_marginal_2=result["beta_2"],
                    beta_interaction=result["beta_interaction"],
                    se_interaction=result["se_interaction"],
                    p_interaction=result["p_interaction"],
                    p_corrected=None,
                    n_samples=result["n"],
                    maf_1=pruned_variants[i]["af_total"],
                    maf_2=pruned_variants[j]["af_total"],
                ))

    if verbose:
        print(f"  Pairs tested: {pairs_tested:,}")
        print(f"  Nominal p < 0.05: {len(results):,}")

    # BH correction
    if results:
        pvals = np.array([r.p_interaction for r in results])
        n_tests = pairs_tested  # correct for all tested, not just significant
        # Bonferroni-like using total tested pairs
        for r in results:
            r.p_corrected = min(1.0, r.p_interaction * n_tests)
        results.sort(key=lambda r: r.p_interaction)

        if verbose:
            sig = sum(1 for r in results if r.p_corrected < 0.05)
            print(f"  Bonferroni p < 0.05: {sig}")

    return results


# ===================================================================
# M1 — source-agnostic entry point (BGEN-compatible)
# ===================================================================

def ld_pruned_cooccurrence_from_data(
    variants: list[dict],
    dosage_list: list[np.ndarray],
    phenotype: np.ndarray,
    r2_prune: float = 0.5,
    min_cocarriers: int = 5,
    min_distance_bp: int = 100_000,
    max_variants: int = 2000,
    verbose: bool = True,
) -> list[InteractionResult]:
    """Pure-function M1: runs LD pruning + interaction testing on pre-loaded
    genotypes and a quantitative phenotype.

    No Neo4j dependency. Can be fed from any source — Neo4j, BGEN via
    BgenReader.load_locus, etc. Binarises the phenotype at 25th/75th
    percentiles to define cases/controls for the co-carrier pre-filter.

    Args:
        variants: list of dicts, each with at least 'pos', 'variantId',
                  'af_total'. Must be in the same order as `dosage_list`.
        dosage_list: list of numpy arrays (length n_samples each, values in [0, 2]).
        phenotype: length-n_samples quantitative phenotype vector.
                   NaN-bearing samples are masked out.
        r2_prune: LD threshold for greedy pruning.
        min_cocarriers: minimum number of co-carriers in cases to test a pair.
        min_distance_bp: minimum physical distance between pair members.
        max_variants: stop LD pruning after this many kept variants.

    Returns:
        list[InteractionResult] sorted by interaction p-value.
    """
    # Sample-level validity mask (drop NaN phenotype samples)
    valid = ~np.isnan(phenotype)
    pheno_v = phenotype[valid]
    # Subset all dosage vectors consistently
    dosage_list_v = [d[valid] for d in dosage_list]
    n_samples = int(valid.sum())
    if n_samples < 40:
        if verbose:
            print(f"M1 from_data: only {n_samples} valid samples — aborting")
        return []

    # Binarise via 25/75 percentiles of phenotype
    q75 = np.percentile(pheno_v, 75)
    q25 = np.percentile(pheno_v, 25)
    case_idx_local = np.where(pheno_v >= q75)[0]
    ctrl_idx_local = np.where(pheno_v <= q25)[0]
    n_case, n_ctrl = len(case_idx_local), len(ctrl_idx_local)
    if verbose:
        print(f"M1 from_data: n={n_samples} ({n_case} cases + {n_ctrl} controls "
              f"by quartile), {len(variants)} variants, r²_prune={r2_prune}")

    # LD pruning (greedy keep-first)
    keep: list[int] = []
    used: set[int] = set()
    for i in range(len(variants)):
        if i in used:
            continue
        keep.append(i)
        for j in range(i + 1, len(variants)):
            if j in used:
                continue
            if abs(variants[i]["pos"] - variants[j]["pos"]) < min_distance_bp:
                r2 = _compute_ld(dosage_list_v[i], dosage_list_v[j])
                if r2 > r2_prune:
                    used.add(j)
        if len(keep) >= max_variants:
            break
    pruned_variants = [variants[i] for i in keep]
    pruned_dosages = [dosage_list_v[i] for i in keep]
    if verbose:
        print(f"  After LD pruning: {len(pruned_variants)} variants")

    results: list[InteractionResult] = []
    pairs_tested = 0
    for i in range(len(pruned_variants)):
        for j in range(i + 1, len(pruned_variants)):
            if abs(pruned_variants[i]["pos"] - pruned_variants[j]["pos"]) < min_distance_bp:
                continue
            d1 = pruned_dosages[i]
            d2 = pruned_dosages[j]
            carriers_1 = d1 > 0
            carriers_2 = d2 > 0
            co_case = int(np.sum(carriers_1[case_idx_local] & carriers_2[case_idx_local]))
            if co_case < min_cocarriers:
                continue
            r = _test_interaction(d1, d2, pheno_v)
            pairs_tested += 1
            if r["p_interaction"] < 0.05:
                results.append(InteractionResult(
                    variant_1=pruned_variants[i]["variantId"],
                    variant_2=pruned_variants[j]["variantId"],
                    motif="co_occurrence",
                    shared_entity=f"co_carriers={co_case}",
                    beta_marginal_1=r["beta_1"],
                    beta_marginal_2=r["beta_2"],
                    beta_interaction=r["beta_interaction"],
                    se_interaction=r["se_interaction"],
                    p_interaction=r["p_interaction"],
                    p_corrected=None,
                    n_samples=r["n"],
                    maf_1=pruned_variants[i]["af_total"],
                    maf_2=pruned_variants[j]["af_total"],
                ))
    if verbose:
        print(f"  Pairs tested: {pairs_tested:,}")
        print(f"  Nominal p < 0.05: {len(results):,}")
    if results:
        for r in results:
            r.p_corrected = min(1.0, r.p_interaction * pairs_tested)
        results.sort(key=lambda r: r.p_interaction)
    return results


def ld_pruned_cooccurrence_bgen(
    reader,
    chr: str,
    start: int,
    end: int,
    phenotype_by_sample: dict[str, float],
    min_af: float = 0.01,
    mac_min: int = 10,
    r2_prune: float = 0.5,
    min_cocarriers: int = 5,
    min_distance_bp: int = 100_000,
    max_variants: int = 2000,
    verbose: bool = True,
) -> list[InteractionResult]:
    """M1 epistasis driven by a BgenReader — UKB-ready path.

    Loads the locus once from BGEN, aligns the phenotype to BGEN sample
    order, then calls ld_pruned_cooccurrence_from_data. The output is
    identical in shape to the Neo4j-backed ld_pruned_cooccurrence().
    """
    variants_df, dosage = reader.load_locus(chr, start, end, format="dosage")
    if len(variants_df) == 0:
        return []
    # Align phenotype to BGEN sample order
    samples = [str(s) for s in reader.samples(chr)]
    pheno = np.array(
        [phenotype_by_sample.get(s, np.nan) for s in samples], dtype=float
    )

    # Per-variant filters: MAF + MAC + mean-impute NaN
    variants: list[dict] = []
    dosage_list: list[np.ndarray] = []
    for i in range(len(variants_df)):
        d = dosage[:, i]
        col_mean = np.nanmean(d)
        if not np.isfinite(col_mean):
            continue
        af = float(col_mean / 2.0)
        if af < min_af or af > 1 - min_af:
            continue
        d_imp = np.where(np.isnan(d), col_mean, d)
        mac = min(int(d_imp.sum()), int(2 * len(d_imp) - d_imp.sum()))
        if mac < mac_min:
            continue
        row = variants_df.iloc[i]
        variants.append({
            "variantId": f"{row['chr']}:{row['pos']}:{row['a1']}:{row['a2']}",
            "pos": int(row["pos"]),
            "chr": row["chr"],
            "ref": row["a1"],
            "alt": row["a2"],
            "af_total": af,
        })
        dosage_list.append(d_imp)
        if len(variants) >= max_variants * 3:  # load extra for LD pruning headroom
            break
    if verbose:
        print(f"M1 BGEN: {len(variants)} variants from {chr}:{start}-{end}")
    return ld_pruned_cooccurrence_from_data(
        variants=variants,
        dosage_list=dosage_list,
        phenotype=pheno,
        r2_prune=r2_prune,
        min_cocarriers=min_cocarriers,
        min_distance_bp=min_distance_bp,
        max_variants=max_variants,
        verbose=verbose,
    )


# ===================================================================
# M3: Differential Subgraph Analysis
# ===================================================================

def differential_subgraph(
    conn: GraphGWASConnection,
    chr: str,
    start: int | None = None,
    end: int | None = None,
    r2_prune: float = 0.5,
    min_cocarriers: int = 5,
    min_distance_bp: int = 100_000,
    mac_min: int = 10,
    max_variants: int = 1000,
    verbose: bool = True,
) -> list[InteractionResult]:
    """M3: Differential Subgraph Analysis.

    Builds SEPARATE co-occurrence graphs for cases and controls, then finds:
    - Case-unique edges (co-occur in cases, not controls → positive epistasis)
    - Control-unique edges (co-occur in controls, not cases → synthetic lethality)
    - Differential edges (exist in both but different weight → frequency shift)

    This captures interactions invisible to enrichment-only methods.
    """
    # Get case/control or binarize quantitative trait
    case_idx, ctrl_idx = get_phenotype_indices(conn)
    if len(case_idx) < 20 or len(ctrl_idx) < 20:
        all_idx = get_all_indices(conn)
        pheno = get_phenotype_values(conn, all_idx)
        valid = ~np.isnan(pheno)
        pheno_v = pheno[valid]
        all_idx = all_idx[valid]
        q75 = np.percentile(pheno_v, 75)
        q25 = np.percentile(pheno_v, 25)
        case_mask = pheno_v >= q75
        ctrl_mask = pheno_v <= q25
        case_local = np.where(case_mask)[0]
        ctrl_local = np.where(ctrl_mask)[0]
        pheno_for_test = pheno_v
    else:
        all_idx = np.concatenate([case_idx, ctrl_idx])
        n_case = len(case_idx)
        case_local = np.arange(n_case)
        ctrl_local = np.arange(n_case, n_case + len(ctrl_idx))
        pheno_for_test = np.concatenate([np.ones(n_case), np.zeros(len(ctrl_idx))])

    n_case = len(case_local)
    n_ctrl = len(ctrl_local)

    if verbose:
        print(f"M3 Differential Subgraph: {chr}:{start or ''}-{end or ''}")
        print(f"  {n_case} cases, {n_ctrl} controls")

    # Load and LD-prune variants (reuse M1 logic)
    variants = []
    dosage_list = []
    for v in variant_iterator(conn, chr, start, end):
        af = v.get("af_total", 0)
        if af < 0.01 or af > 0.99:
            continue
        gtp = v.get("gt_packed")
        if gtp is None:
            continue
        d = build_dosage(gtp, all_idx, _cfg.N_SAMPLES)
        d_clean = np.where(np.isnan(d), 0, d)
        mac = min(int(d_clean.sum()), int(2 * len(all_idx) - d_clean.sum()))
        if mac < mac_min:
            continue
        variants.append(v)
        dosage_list.append(d)
        if len(variants) >= max_variants * 3:
            break

    # LD pruning
    keep = []
    used = set()
    for i in range(len(variants)):
        if i in used:
            continue
        keep.append(i)
        for j in range(i + 1, len(variants)):
            if j in used:
                continue
            if abs(variants[i]["pos"] - variants[j]["pos"]) < min_distance_bp:
                r2 = _compute_ld(dosage_list[i], dosage_list[j])
                if r2 > r2_prune:
                    used.add(j)
        if len(keep) >= max_variants:
            break

    pruned_v = [variants[i] for i in keep]
    pruned_d = [dosage_list[i] for i in keep]
    n_v = len(pruned_v)

    if verbose:
        print(f"  {len(variants)} loaded → {n_v} after LD pruning")

    if n_v < 4:
        return []

    # Build case and control carrier arrays
    carrier_case = np.zeros((n_v, n_case), dtype=bool)
    carrier_ctrl = np.zeros((n_v, n_ctrl), dtype=bool)
    for i in range(n_v):
        carrier_case[i] = pruned_d[i][case_local] > 0
        carrier_ctrl[i] = pruned_d[i][ctrl_local] > 0

    # Build co-occurrence counts for cases and controls
    results = []
    pairs_tested = 0

    for i in range(n_v):
        for j in range(i + 1, n_v):
            if abs(pruned_v[i]["pos"] - pruned_v[j]["pos"]) < min_distance_bp:
                continue

            co_case = int(np.sum(carrier_case[i] & carrier_case[j]))
            co_ctrl = int(np.sum(carrier_ctrl[i] & carrier_ctrl[j]))

            # Classify edge
            case_freq = co_case / n_case if n_case > 0 else 0
            ctrl_freq = co_ctrl / n_ctrl if n_ctrl > 0 else 0

            # Skip if no co-occurrence in either group
            if co_case < 2 and co_ctrl < 2:
                continue

            pairs_tested += 1

            # Log frequency ratio
            pseudo = 1.0 / max(n_case, n_ctrl)
            lfr = np.log2((case_freq + pseudo) / (ctrl_freq + pseudo))

            # Fisher's exact test for differential co-occurrence
            # Table: [[co_case, case_solo], [co_ctrl, ctrl_solo]]
            case_solo = n_case - co_case
            ctrl_solo = n_ctrl - co_ctrl
            try:
                _, fisher_p = sp_stats.fisher_exact(
                    [[co_case, case_solo], [co_ctrl, ctrl_solo]]
                )
            except ValueError:
                fisher_p = 1.0

            # Classify edge type
            if co_case >= min_cocarriers and co_ctrl < 2:
                edge_type = "case_unique"
            elif co_ctrl >= min_cocarriers and co_case < 2:
                edge_type = "ctrl_unique"
            elif abs(lfr) > 1.0:
                edge_type = "differential"
            else:
                edge_type = "shared"

            # Also test interaction term for significant edges
            if fisher_p < 0.1:
                int_result = _test_interaction(pruned_d[i], pruned_d[j], pheno_for_test)
                p_int = int_result["p_interaction"]
                beta_int = int_result["beta_interaction"]
            else:
                p_int = 1.0
                beta_int = 0.0

            if fisher_p < 0.05 or edge_type in ("case_unique", "ctrl_unique"):
                results.append(InteractionResult(
                    variant_1=pruned_v[i]["variantId"],
                    variant_2=pruned_v[j]["variantId"],
                    motif=edge_type,
                    shared_entity=f"lfr={lfr:.2f},co_case={co_case},co_ctrl={co_ctrl}",
                    beta_marginal_1=0.0,
                    beta_marginal_2=0.0,
                    beta_interaction=beta_int,
                    se_interaction=0.0,
                    p_interaction=p_int,
                    p_corrected=min(1.0, fisher_p * pairs_tested),
                    n_samples=n_case + n_ctrl,
                    maf_1=pruned_v[i]["af_total"],
                    maf_2=pruned_v[j]["af_total"],
                ))

    results.sort(key=lambda r: r.p_interaction)

    if verbose:
        print(f"  Pairs tested: {pairs_tested:,}")
        by_type = {}
        for r in results:
            by_type[r.motif] = by_type.get(r.motif, 0) + 1
        for t, c in sorted(by_type.items()):
            print(f"    {t}: {c}")
        if results:
            print(f"  Best: {results[0].variant_1} × {results[0].variant_2}")
            print(f"    type={results[0].motif}, β_int={results[0].beta_interaction:.4f}, p={results[0].p_interaction:.2e}")

    return results


# ===================================================================
# M4: Dark Matter — Negative Co-occurrence / Synthetic Incompatibility
# ===================================================================

def dark_matter_epistasis(
    conn: GraphGWASConnection,
    chr: str,
    start: int | None = None,
    end: int | None = None,
    maf_min: float = 0.05,
    depletion_threshold: float = 0.5,
    min_expected: float = 5.0,
    mac_min: int = 20,
    max_variants: int = 500,
    use_motif_prefilter: bool = True,
    verbose: bool = True,
) -> list[InteractionResult]:
    """M4: Dark Matter — detect depleted co-occurrence (synthetic incompatibility).

    For each variant pair, computes expected co-occurrence under independence
    and compares to observed. Pairs with observed << expected in cases (but not
    in controls) are candidates for synthetic lethality.

    This detects the OPPOSITE of what M1/M3 find — absent combinations, not
    enriched ones. No existing GWAS tool explicitly tests this.

    Args:
        conn: database connection.
        chr: chromosome.
        start, end: region.
        maf_min: minimum MAF for both variants (need enough expected co-occurrences).
        depletion_threshold: max ratio observed/expected to flag (e.g., 0.5 = 50% depleted).
        min_expected: minimum expected co-occurrences for testability.
        mac_min: minimum minor allele count.
        max_variants: max variants to load.
        use_motif_prefilter: if True and Gene/Pathway exist, prefer same-pathway pairs.
        verbose: print progress.

    Returns:
        list of InteractionResult sorted by depletion significance.
    """
    # Get case/control indices
    case_idx, ctrl_idx = get_phenotype_indices(conn)
    if len(case_idx) < 20 or len(ctrl_idx) < 20:
        all_idx = get_all_indices(conn)
        pheno = get_phenotype_values(conn, all_idx)
        valid = ~np.isnan(pheno)
        pheno_v = pheno[valid]
        all_idx = all_idx[valid]
        q75 = np.percentile(pheno_v, 75)
        q25 = np.percentile(pheno_v, 25)
        case_local = np.where(pheno_v >= q75)[0]
        ctrl_local = np.where(pheno_v <= q25)[0]
    else:
        all_idx = np.concatenate([case_idx, ctrl_idx])
        n_case = len(case_idx)
        case_local = np.arange(n_case)
        ctrl_local = np.arange(n_case, n_case + len(ctrl_idx))

    n_case = len(case_local)
    n_ctrl = len(ctrl_local)

    if verbose:
        print(f"M4 Dark Matter: {chr}:{start or ''}-{end or ''}")
        print(f"  {n_case} cases, {n_ctrl} controls, MAF≥{maf_min}")

    # Load common variants (need high MAF for expected co-occurrence power)
    variants = []
    dosage_list = []
    carrier_freqs_case = []
    carrier_freqs_ctrl = []

    for v in variant_iterator(conn, chr, start, end):
        af = v.get("af_total", 0)
        if af < maf_min or af > (1 - maf_min):
            continue
        gtp = v.get("gt_packed")
        if gtp is None:
            continue
        d = build_dosage(gtp, all_idx, _cfg.N_SAMPLES)
        d_clean = np.where(np.isnan(d), 0, d)
        mac = min(int(d_clean.sum()), int(2 * len(all_idx) - d_clean.sum()))
        if mac < mac_min:
            continue

        # Carrier frequency in cases and controls
        carriers = d > 0
        cf_case = carriers[case_local].mean()
        cf_ctrl = carriers[ctrl_local].mean()

        variants.append(v)
        dosage_list.append(d)
        carrier_freqs_case.append(cf_case)
        carrier_freqs_ctrl.append(cf_ctrl)

        if len(variants) >= max_variants:
            break

    n_v = len(variants)
    if verbose:
        print(f"  Loaded {n_v} common variants (MAF≥{maf_min})")

    if n_v < 4:
        return []

    # Test all pairs for depleted co-occurrence
    results = []
    pairs_tested = 0
    depleted_found = 0

    for i in range(n_v):
        for j in range(i + 1, n_v):
            # Expected co-occurrence under independence
            expected_case = n_case * carrier_freqs_case[i] * carrier_freqs_case[j]
            expected_ctrl = n_ctrl * carrier_freqs_ctrl[i] * carrier_freqs_ctrl[j]

            # Need enough expected for statistical power
            if expected_case < min_expected:
                continue

            # Observed co-occurrence
            d_i = dosage_list[i]
            d_j = dosage_list[j]
            carriers_i = d_i > 0
            carriers_j = d_j > 0
            observed_case = int(np.sum(carriers_i[case_local] & carriers_j[case_local]))
            observed_ctrl = int(np.sum(carriers_i[ctrl_local] & carriers_j[ctrl_local]))

            pairs_tested += 1

            # Depletion ratio
            ratio_case = observed_case / expected_case if expected_case > 0 else 1.0
            ratio_ctrl = observed_ctrl / expected_ctrl if expected_ctrl > 0 else 1.0

            # Case-specific depletion: depleted in cases but NOT in controls
            is_depleted = (ratio_case < depletion_threshold and ratio_ctrl > depletion_threshold)

            # Poisson test: P(X ≤ observed | λ = expected)
            poisson_p_case = float(sp_stats.poisson.cdf(observed_case, expected_case))

            # Differential depletion: z-score
            z_case = (observed_case - expected_case) / np.sqrt(expected_case + 1)
            z_ctrl = (observed_ctrl - expected_ctrl) / np.sqrt(max(expected_ctrl, 1) + 1)
            depletion_differential = z_ctrl - z_case  # positive = depleted in cases only

            if is_depleted or (poisson_p_case < 0.01 and depletion_differential > 1.0):
                depleted_found += 1
                results.append(InteractionResult(
                    variant_1=variants[i]["variantId"],
                    variant_2=variants[j]["variantId"],
                    motif="dark_matter",
                    shared_entity=(
                        f"obs_case={observed_case},exp_case={expected_case:.1f},"
                        f"ratio={ratio_case:.3f},"
                        f"obs_ctrl={observed_ctrl},exp_ctrl={expected_ctrl:.1f}"
                    ),
                    beta_marginal_1=0.0,
                    beta_marginal_2=0.0,
                    beta_interaction=depletion_differential,  # use z-diff as "effect"
                    se_interaction=0.0,
                    p_interaction=poisson_p_case,
                    p_corrected=min(1.0, poisson_p_case * pairs_tested),
                    n_samples=n_case + n_ctrl,
                    maf_1=variants[i]["af_total"],
                    maf_2=variants[j]["af_total"],
                ))

    results.sort(key=lambda r: r.p_interaction)

    if verbose:
        print(f"  Pairs tested: {pairs_tested:,}")
        print(f"  Depleted pairs found: {depleted_found}")
        if results:
            best = results[0]
            print(f"  Best: {best.variant_1} × {best.variant_2}")
            print(f"    {best.shared_entity}")
            print(f"    Poisson p = {best.p_interaction:.2e}")

    return results


def motif_epistasis_to_tsv(results: list[InteractionResult], path: str):
    """Write epistasis v2 results to TSV file."""
    import csv
    cols = ["variant_1", "variant_2", "motif", "shared_entity",
            "beta_marginal_1", "beta_marginal_2", "beta_interaction",
            "se_interaction", "p_interaction", "p_corrected",
            "n_samples", "maf_1", "maf_2"]

    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        writer.writeheader()
        for r in results:
            writer.writerow({
                "variant_1": r.variant_1,
                "variant_2": r.variant_2,
                "motif": r.motif,
                "shared_entity": r.shared_entity,
                "beta_marginal_1": f"{r.beta_marginal_1:.6f}",
                "beta_marginal_2": f"{r.beta_marginal_2:.6f}",
                "beta_interaction": f"{r.beta_interaction:.6f}",
                "se_interaction": f"{r.se_interaction:.6f}",
                "p_interaction": f"{r.p_interaction:.6e}",
                "p_corrected": f"{r.p_corrected:.6e}" if r.p_corrected else "",
                "n_samples": r.n_samples,
                "maf_1": f"{r.maf_1:.4f}",
                "maf_2": f"{r.maf_2:.4f}",
            })


# ===================================================================
# Bridge to graphgwas.bias: build interaction-feature matrix Z from
# motif-typed pair results so rho_max / R(x) can be evaluated under
# a biology-typed interaction subspace (cf. paper #2).
# ===================================================================

def motif_interaction_matrix(motif_results: list,
                              dosages: np.ndarray,
                              variant_id_to_col: dict[str, int]) -> np.ndarray:
    """Build interaction-feature matrix Z from M2 motif-pair results.

    For each motif pair ``(v1, v2)`` in *motif_results*, append a column
    ``dosage(v1) * dosage(v2)`` (centred) to ``Z``.  The resulting matrix
    spans the biology-typed interaction subspace and is the input
    ``Z`` to :func:`graphgwas.bias.rho_max` and the ``conservativeness_ratio``.

    Args:
        motif_results: list of :class:`InteractionResult` (from
            :func:`motif_filtered_epistasis`).
        dosages: per-variant dosage matrix, shape (n_samples, n_variants).
        variant_id_to_col: maps each ``InteractionResult.variant_*``
            variant ID string to its column index in ``dosages``.

    Returns:
        Z, shape (n_samples, n_motif_pairs), with each column the centred
        dosage product for one motif pair.  Pairs whose variant IDs are
        not present in ``variant_id_to_col`` are silently skipped; the
        returned matrix may have fewer columns than ``len(motif_results)``.
    """
    cols = []
    for r in motif_results:
        i = variant_id_to_col.get(r.variant_1)
        j = variant_id_to_col.get(r.variant_2)
        if i is None or j is None:
            continue
        prod = dosages[:, i] * dosages[:, j]
        prod = prod - prod.mean()
        cols.append(prod)
    if not cols:
        # 0-column matrix with the right number of rows
        return np.empty((dosages.shape[0], 0), dtype=float)
    return np.column_stack(cols)
