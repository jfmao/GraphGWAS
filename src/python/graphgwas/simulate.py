"""Ground-truth simulation framework for GraphGWAS v2 benchmarks.

Simulates phenotypes with known genetic architectures:
- Epistasis scenarios: S1 (pure interaction), S2 (additive+interaction),
  S3 (synthetic lethality), S4 (pathway epistasis), S5 (multi-way)
- Fine-mapping scenarios: F1 (single causal in LD), F2 (two independent causal),
  F3 (causal+proxy), F4 (annotated causal), F5 (multi-signal)

Each simulator returns phenotypes + ground truth specification.
Uses real genotype data from Neo4j (preserves actual LD structure).
"""

from __future__ import annotations

import numpy as np
from dataclasses import dataclass, field

from . import config as _cfg
from .db import GraphGWASConnection
from .genotype import (
    get_all_indices,
    variant_iterator,
    build_dosage,
    build_carrier_set,
)


# ===================================================================
# Ground truth data classes
# ===================================================================

@dataclass
class EpistaticPair:
    """Ground truth for a pairwise epistatic interaction."""
    variant_1: str
    variant_2: str
    beta_marginal_1: float
    beta_marginal_2: float
    beta_interaction: float
    pos_1: int
    pos_2: int
    af_1: float
    af_2: float


@dataclass
class CausalVariant:
    """Ground truth for a single causal variant."""
    variant_id: str
    chr: str
    pos: int
    af: float
    beta: float
    gene: str | None = None
    pathway: str | None = None


@dataclass
class SimulationResult:
    """Complete simulation output with phenotype + ground truth."""
    phenotype: np.ndarray
    sample_ids: np.ndarray
    scenario: str
    realized_h2: float
    causal_variants: list[CausalVariant] = field(default_factory=list)
    epistatic_pairs: list[EpistaticPair] = field(default_factory=list)
    noise_variance: float = 0.0
    config: dict = field(default_factory=dict)

    def ground_truth_variants(self) -> set[str]:
        """Return set of all variant IDs involved in the causal architecture."""
        ids = {cv.variant_id for cv in self.causal_variants}
        for pair in self.epistatic_pairs:
            ids.add(pair.variant_1)
            ids.add(pair.variant_2)
        return ids


# ===================================================================
# Helper: select variants with AF constraints
# ===================================================================

def _select_variants_by_af(
    conn: GraphGWASConnection,
    chr: str,
    af_min: float,
    af_max: float,
    n_needed: int,
    start: int | None = None,
    end: int | None = None,
    min_spacing: int = 100_000,
    rng: np.random.Generator | None = None,
) -> list[dict]:
    """Find variants in AF range, spaced apart to avoid LD.

    Args:
        min_spacing: minimum bp between selected variants (to reduce LD)
    """
    if rng is None:
        rng = np.random.default_rng(42)

    # Use chr+pos index with AF filter applied after. The chr+pos index is
    # much more selective than af_total alone.
    # If no region specified, sample sparsely along the chromosome.
    if start is None or end is None:
        # Need to cover the whole chromosome — sample at spacing intervals
        # by using cursor-based pagination with wider steps
        cursor = -1
        selected = []
        while len(selected) < n_needed:
            result = conn.execute_read(
                """
                MATCH (v:Variant)
                WHERE v.chr = $chr AND v.pos > $cursor
                  AND v.af_total >= $af_min AND v.af_total <= $af_max
                  AND v.gt_packed IS NOT NULL
                RETURN v.variantId AS variantId, v.chr AS chr, v.pos AS pos,
                       v.ref AS ref, v.alt AS alt, v.af_total AS af_total,
                       v.gt_packed AS gt_packed
                ORDER BY v.pos
                LIMIT 1
                """,
                {"chr": chr, "cursor": cursor,
                 "af_min": af_min, "af_max": af_max},
            )
            recs = list(result)
            if not recs:
                break
            v = dict(recs[0])
            selected.append(v)
            cursor = v["pos"] + min_spacing  # skip forward
    else:
        # Region-bounded: load all in region then filter by spacing
        result = conn.execute_read(
            """
            MATCH (v:Variant)
            WHERE v.chr = $chr AND v.pos >= $start AND v.pos <= $end
              AND v.af_total >= $af_min AND v.af_total <= $af_max
              AND v.gt_packed IS NOT NULL
            RETURN v.variantId AS variantId, v.chr AS chr, v.pos AS pos,
                   v.ref AS ref, v.alt AS alt, v.af_total AS af_total,
                   v.gt_packed AS gt_packed
            ORDER BY v.pos
            """,
            {"chr": chr, "start": start, "end": end,
             "af_min": af_min, "af_max": af_max},
        )
        selected = []
        last_pos = -min_spacing - 1
        for rec in result:
            v = dict(rec)
            if v["pos"] - last_pos < min_spacing:
                continue
            selected.append(v)
            last_pos = v["pos"]
            if len(selected) >= n_needed:
                break

    return selected


def _get_dosage_vectors(
    conn: GraphGWASConnection,
    variants: list[dict],
    all_idx: np.ndarray,
) -> dict[str, np.ndarray]:
    """Get dosage vector for each variant, mean-imputing missing."""
    result = {}
    for v in variants:
        dosage = build_dosage(v["gt_packed"], all_idx, _cfg.N_SAMPLES)
        m = np.nanmean(dosage)
        dosage = np.where(np.isnan(dosage), m, dosage)
        result[v["variantId"]] = dosage
    return result


# ===================================================================
# EPISTASIS SCENARIOS
# ===================================================================

def simulate_S1_pure_interaction(
    conn: GraphGWASConnection,
    chr: str,
    n_pairs: int = 3,
    beta_interaction: float = 1.5,
    h2_target: float = 0.3,
    af_range: tuple[float, float] = (0.1, 0.4),
    seed: int = 42,
) -> SimulationResult:
    """S1: Pure interaction, no marginal effects.

    Y = Σ β_int × (G_i - mean) × (G_j - mean) + ε

    This is the hardest scenario: only interaction term is real, no main effects.
    Tests whether a method can find interactions WITHOUT marginal signals.
    """
    rng = np.random.default_rng(seed)
    all_idx = get_all_indices(conn)
    n = len(all_idx)

    # Need 2 variants per pair, well-spaced to avoid LD
    n_variants_needed = n_pairs * 2
    variants = _select_variants_by_af(
        conn, chr, af_range[0], af_range[1], n_variants_needed,
        min_spacing=1_000_000,  # 1Mb apart to ensure no LD
        rng=rng,
    )

    if len(variants) < n_variants_needed:
        raise ValueError(
            f"Could not find {n_variants_needed} variants in AF range {af_range} on {chr}"
        )

    # Pair up the variants
    dosages = _get_dosage_vectors(conn, variants, all_idx)
    pairs = []
    genetic = np.zeros(n)

    for i in range(n_pairs):
        v1 = variants[2 * i]
        v2 = variants[2 * i + 1]
        d1 = dosages[v1["variantId"]]
        d2 = dosages[v2["variantId"]]
        # Center so that interaction has no marginal component
        d1_c = d1 - d1.mean()
        d2_c = d2 - d2.mean()
        contribution = beta_interaction * d1_c * d2_c
        genetic += contribution

        pairs.append(
            EpistaticPair(
                variant_1=v1["variantId"],
                variant_2=v2["variantId"],
                beta_marginal_1=0.0,
                beta_marginal_2=0.0,
                beta_interaction=beta_interaction,
                pos_1=v1["pos"],
                pos_2=v2["pos"],
                af_1=v1["af_total"],
                af_2=v2["af_total"],
            )
        )

    # Add noise to achieve target h²
    var_g = np.var(genetic)
    if var_g == 0:
        raise RuntimeError("Genetic variance is zero — check dosages")
    var_e = var_g * (1 / h2_target - 1)
    noise = rng.normal(0, np.sqrt(var_e), n)
    phenotype = genetic + noise

    realized_h2 = var_g / np.var(phenotype)

    return SimulationResult(
        phenotype=phenotype,
        sample_ids=all_idx,
        scenario="S1_pure_interaction",
        realized_h2=realized_h2,
        epistatic_pairs=pairs,
        noise_variance=var_e,
        config={"n_pairs": n_pairs, "beta_interaction": beta_interaction,
                "h2_target": h2_target, "chr": chr},
    )


def simulate_S2_additive_plus_interaction(
    conn: GraphGWASConnection,
    chr: str,
    n_pairs: int = 3,
    beta_marginal: float = 0.3,
    beta_interaction: float = 0.5,
    h2_target: float = 0.3,
    af_range: tuple[float, float] = (0.1, 0.4),
    seed: int = 42,
) -> SimulationResult:
    """S2: Additive + interaction effects.

    Y = Σ (β_A·G_A + β_B·G_B + β_AB·G_A·G_B) + ε

    More realistic: both marginal and interaction effects present.
    Marginal tests find the main effects; interaction test should find β_AB.
    """
    rng = np.random.default_rng(seed)
    all_idx = get_all_indices(conn)
    n = len(all_idx)

    n_variants_needed = n_pairs * 2
    variants = _select_variants_by_af(
        conn, chr, af_range[0], af_range[1], n_variants_needed,
        min_spacing=1_000_000, rng=rng,
    )
    if len(variants) < n_variants_needed:
        raise ValueError("Not enough variants found")

    dosages = _get_dosage_vectors(conn, variants, all_idx)
    pairs = []
    genetic = np.zeros(n)

    for i in range(n_pairs):
        v1 = variants[2 * i]
        v2 = variants[2 * i + 1]
        d1 = dosages[v1["variantId"]]
        d2 = dosages[v2["variantId"]]
        contribution = beta_marginal * d1 + beta_marginal * d2 + beta_interaction * d1 * d2
        genetic += contribution

        pairs.append(
            EpistaticPair(
                variant_1=v1["variantId"],
                variant_2=v2["variantId"],
                beta_marginal_1=beta_marginal,
                beta_marginal_2=beta_marginal,
                beta_interaction=beta_interaction,
                pos_1=v1["pos"],
                pos_2=v2["pos"],
                af_1=v1["af_total"],
                af_2=v2["af_total"],
            )
        )

    var_g = np.var(genetic)
    var_e = var_g * (1 / h2_target - 1)
    noise = rng.normal(0, np.sqrt(var_e), n)
    phenotype = genetic + noise

    return SimulationResult(
        phenotype=phenotype,
        sample_ids=all_idx,
        scenario="S2_additive_plus_interaction",
        realized_h2=var_g / np.var(phenotype),
        epistatic_pairs=pairs,
        noise_variance=var_e,
        config={"n_pairs": n_pairs, "beta_marginal": beta_marginal,
                "beta_interaction": beta_interaction, "h2_target": h2_target, "chr": chr},
    )


def simulate_S3_synthetic_lethality(
    conn: GraphGWASConnection,
    chr: str,
    beta_marginal: float = 0.3,
    h2_target: float = 0.2,
    af_range: tuple[float, float] = (0.2, 0.4),
    seed: int = 42,
) -> SimulationResult:
    """S3: Synthetic lethality — double mutant is removed from sample.

    Simulates: Y ~ G_A + G_B, but samples where BOTH G_A=2 AND G_B=2 are excluded.
    This creates a signature where the pair co-occurs less often than expected.

    Method M4 (dark matter) should detect this; other methods will miss it.
    """
    rng = np.random.default_rng(seed)
    all_idx = get_all_indices(conn)

    variants = _select_variants_by_af(
        conn, chr, af_range[0], af_range[1], n_needed=2,
        min_spacing=1_000_000, rng=rng,
    )
    if len(variants) < 2:
        raise ValueError("Could not find 2 variants")

    dosages = _get_dosage_vectors(conn, variants, all_idx)
    d1 = dosages[variants[0]["variantId"]]
    d2 = dosages[variants[1]["variantId"]]

    # Identify "lethal" samples (both hom-alt)
    lethal_mask = (d1 == 2) & (d2 == 2)
    # Keep samples not both hom-alt
    keep_mask = ~lethal_mask

    # Phenotype for surviving samples
    d1_k = d1[keep_mask]
    d2_k = d2[keep_mask]
    genetic = beta_marginal * d1_k + beta_marginal * d2_k
    var_g = np.var(genetic)
    var_e = var_g * (1 / h2_target - 1)
    noise = rng.normal(0, np.sqrt(var_e), keep_mask.sum())

    # Full phenotype vector with lethal samples as NaN
    phenotype = np.full(len(all_idx), np.nan)
    phenotype[keep_mask] = genetic + noise

    pairs = [
        EpistaticPair(
            variant_1=variants[0]["variantId"],
            variant_2=variants[1]["variantId"],
            beta_marginal_1=beta_marginal,
            beta_marginal_2=beta_marginal,
            beta_interaction=-np.inf,  # synthetic lethality = infinite negative interaction
            pos_1=variants[0]["pos"],
            pos_2=variants[1]["pos"],
            af_1=variants[0]["af_total"],
            af_2=variants[1]["af_total"],
        )
    ]

    return SimulationResult(
        phenotype=phenotype,
        sample_ids=all_idx,
        scenario="S3_synthetic_lethality",
        realized_h2=var_g / np.var(phenotype[keep_mask]),
        epistatic_pairs=pairs,
        noise_variance=var_e,
        config={"beta_marginal": beta_marginal, "h2_target": h2_target,
                "n_lethal_removed": int(lethal_mask.sum()), "chr": chr},
    )


def simulate_S5_multiway(
    conn: GraphGWASConnection,
    chr: str,
    n_way: int = 3,
    beta_multiway: float = 1.0,
    h2_target: float = 0.2,
    af_range: tuple[float, float] = (0.2, 0.4),
    seed: int = 42,
) -> SimulationResult:
    """S5: Higher-order interaction (3-way or 4-way).

    Y = β × G_A × G_B × G_C + ε  (for 3-way)

    Tests methods that claim to detect higher-order interactions (M5, M7).
    """
    rng = np.random.default_rng(seed)
    all_idx = get_all_indices(conn)
    n = len(all_idx)

    variants = _select_variants_by_af(
        conn, chr, af_range[0], af_range[1], n_needed=n_way,
        min_spacing=1_000_000, rng=rng,
    )
    if len(variants) < n_way:
        raise ValueError(f"Could not find {n_way} variants")

    dosages = _get_dosage_vectors(conn, variants, all_idx)
    vid_list = [v["variantId"] for v in variants]

    # Compute product of centered dosages
    genetic = np.ones(n) * beta_multiway
    for v in variants:
        d = dosages[v["variantId"]]
        d_c = d - d.mean()
        genetic = genetic * d_c

    var_g = np.var(genetic)
    if var_g == 0:
        raise RuntimeError("Genetic variance is zero")
    var_e = var_g * (1 / h2_target - 1)
    noise = rng.normal(0, np.sqrt(var_e), n)
    phenotype = genetic + noise

    # Store as single "pair" containing all k variants in variant_1 (comma-separated)
    # For multi-way we use causal_variants list instead
    causal = [
        CausalVariant(
            variant_id=v["variantId"],
            chr=chr,
            pos=v["pos"],
            af=v["af_total"],
            beta=beta_multiway,  # effect only in combination
        )
        for v in variants
    ]

    return SimulationResult(
        phenotype=phenotype,
        sample_ids=all_idx,
        scenario=f"S5_{n_way}way_interaction",
        realized_h2=var_g / np.var(phenotype),
        causal_variants=causal,
        noise_variance=var_e,
        config={"n_way": n_way, "beta_multiway": beta_multiway,
                "h2_target": h2_target, "chr": chr},
    )


# ===================================================================
# FINE-MAPPING SCENARIOS
# ===================================================================

def simulate_F1_single_causal_in_ld(
    conn: GraphGWASConnection,
    chr: str,
    locus_center: int,
    locus_window: int = 500_000,
    beta: float = 0.5,
    h2_target: float = 0.1,
    af_range: tuple[float, float] = (0.05, 0.5),
    seed: int = 42,
) -> SimulationResult:
    """F1: One causal variant in an LD block with ~50-200 correlated variants.

    Tests whether a fine-mapping method can distinguish the true causal variant
    from its many proxies.
    """
    rng = np.random.default_rng(seed)
    all_idx = get_all_indices(conn)
    n = len(all_idx)

    # Pick ONE causal variant in the center of the locus
    start = locus_center - 5000
    end = locus_center + 5000

    causal_variants = _select_variants_by_af(
        conn, chr, af_range[0], af_range[1], n_needed=1,
        start=start, end=end, min_spacing=0, rng=rng,
    )
    if not causal_variants:
        raise ValueError(f"No variant in AF range at {chr}:{locus_center}")

    causal = causal_variants[0]
    dosage = build_dosage(causal["gt_packed"], all_idx, _cfg.N_SAMPLES)
    dosage = np.where(np.isnan(dosage), np.nanmean(dosage), dosage)

    genetic = beta * dosage
    var_g = np.var(genetic)
    var_e = var_g * (1 / h2_target - 1)
    noise = rng.normal(0, np.sqrt(var_e), n)
    phenotype = genetic + noise

    return SimulationResult(
        phenotype=phenotype,
        sample_ids=all_idx,
        scenario="F1_single_causal_in_ld",
        realized_h2=var_g / np.var(phenotype),
        causal_variants=[
            CausalVariant(
                variant_id=causal["variantId"],
                chr=chr,
                pos=causal["pos"],
                af=causal["af_total"],
                beta=beta,
            )
        ],
        noise_variance=var_e,
        config={"locus_center": locus_center, "locus_window": locus_window,
                "beta": beta, "h2_target": h2_target, "chr": chr},
    )


def simulate_F2_two_independent_causal(
    conn: GraphGWASConnection,
    chr: str,
    locus_center: int,
    locus_window: int = 2_000_000,  # 2Mb — wide enough for 2 independent signals
    betas: tuple[float, float] = (0.5, 0.3),
    h2_target: float = 0.15,
    af_range: tuple[float, float] = (0.1, 0.4),
    seed: int = 42,
) -> SimulationResult:
    """F2: Two independent causal variants in same locus, not in LD.

    Tests whether fine-mapping can identify BOTH causal variants
    (multi-signal detection — where single-causal methods fail).
    """
    rng = np.random.default_rng(seed)
    all_idx = get_all_indices(conn)
    n = len(all_idx)

    start = locus_center - locus_window // 2
    end = locus_center + locus_window // 2

    causal_vars = _select_variants_by_af(
        conn, chr, af_range[0], af_range[1], n_needed=2,
        start=start, end=end,
        min_spacing=500_000,  # 500kb apart = likely independent
        rng=rng,
    )
    if len(causal_vars) < 2:
        raise ValueError("Could not find 2 independent causal variants in locus")

    dosages = _get_dosage_vectors(conn, causal_vars, all_idx)
    d1 = dosages[causal_vars[0]["variantId"]]
    d2 = dosages[causal_vars[1]["variantId"]]

    genetic = betas[0] * d1 + betas[1] * d2
    var_g = np.var(genetic)
    var_e = var_g * (1 / h2_target - 1)
    noise = rng.normal(0, np.sqrt(var_e), n)
    phenotype = genetic + noise

    causal_list = [
        CausalVariant(
            variant_id=causal_vars[i]["variantId"],
            chr=chr,
            pos=causal_vars[i]["pos"],
            af=causal_vars[i]["af_total"],
            beta=betas[i],
        )
        for i in range(2)
    ]

    return SimulationResult(
        phenotype=phenotype,
        sample_ids=all_idx,
        scenario="F2_two_independent_causal",
        realized_h2=var_g / np.var(phenotype),
        causal_variants=causal_list,
        noise_variance=var_e,
        config={"locus_center": locus_center, "locus_window": locus_window,
                "betas": list(betas), "h2_target": h2_target, "chr": chr},
    )


# ===================================================================
# Simulator registry
# ===================================================================

EPISTASIS_SCENARIOS = {
    "S1": simulate_S1_pure_interaction,
    "S2": simulate_S2_additive_plus_interaction,
    "S3": simulate_S3_synthetic_lethality,
    "S5": simulate_S5_multiway,
}

def simulate_F6_multi_locus_pathway(
    conn: GraphGWASConnection,
    chr: str,
    n_pathway_loci: int = 5,
    n_null_loci: int = 5,
    beta: float = 0.3,
    h2_target: float = 0.10,
    af_range: tuple[float, float] = (0.05, 0.5),
    seed: int = 42,
) -> SimulationResult:
    """F6: Multiple causal loci, some sharing a biological pathway.

    Selects causal variants in genes of a shared pathway (n_pathway_loci)
    plus causal variants NOT in any pathway (n_null_loci). Used to benchmark
    cross-locus fine-mapping (CLGF) which shares evidence via pathways.
    """
    rng = np.random.default_rng(seed)
    all_idx = get_all_indices(conn)
    n = len(all_idx)

    # Find genes in a pathway with many genes
    pathway_genes = conn.execute_read("""
        MATCH (g:Gene)-[:IN_PATHWAY]->(p:Pathway)
        WHERE g.symbol IS NOT NULL
        WITH p.name AS pathway, collect(DISTINCT g.symbol) AS genes
        WHERE size(genes) >= $min_genes
        RETURN pathway, genes
        ORDER BY size(genes) DESC
        LIMIT 5
    """, {"min_genes": n_pathway_loci})

    pathway_data = [(r["pathway"], r["genes"]) for r in pathway_genes]
    if not pathway_data:
        raise ValueError("No pathway with enough genes found")

    # Pick the largest pathway
    pw_name, pw_genes = pathway_data[0]
    rng.shuffle(pw_genes)
    selected_genes = pw_genes[:n_pathway_loci]

    # For each selected gene, find a variant with HAS_CONSEQUENCE in AF range
    causal_variants = []
    genetic = np.zeros(n)

    for gene in selected_genes:
        result = conn.execute_read("""
            MATCH (v:Variant)-[:HAS_CONSEQUENCE]->(g:Gene {symbol: $gene})
            WHERE v.af_total >= $af_min AND v.af_total <= $af_max
                  AND v.gt_packed IS NOT NULL
            RETURN v.variantId AS vid, v.pos AS pos, v.chr AS chr,
                   v.af_total AS af, v.gt_packed AS gtp
            ORDER BY abs(v.af_total - 0.2) ASC
            LIMIT 5
        """, {"gene": gene, "af_min": af_range[0], "af_max": af_range[1]})

        candidates = list(result)
        if not candidates:
            continue

        pick = candidates[rng.integers(len(candidates))]
        dosage = build_dosage(pick["gtp"], all_idx, _cfg.N_SAMPLES)
        dosage = np.where(np.isnan(dosage), np.nanmean(dosage), dosage)
        genetic += beta * dosage

        causal_variants.append(CausalVariant(
            variant_id=pick["vid"], chr=pick["chr"], pos=pick["pos"],
            af=pick["af"], beta=beta, gene=gene, pathway=pw_name,
        ))

    # Add null-pathway causal variants (far from annotated genes)
    null_variants = conn.execute_read("""
        MATCH (v:Variant)
        WHERE v.chr = $chr
              AND v.af_total >= $af_min AND v.af_total <= $af_max
              AND v.gt_packed IS NOT NULL
              AND NOT (v)-[:HAS_CONSEQUENCE]->(:Gene)
        RETURN v.variantId AS vid, v.pos AS pos, v.chr AS chr,
               v.af_total AS af, v.gt_packed AS gtp
        ORDER BY rand()
        LIMIT $n
    """, {"chr": chr, "af_min": af_range[0], "af_max": af_range[1],
          "n": n_null_loci * 3})

    null_picks = list(null_variants)
    rng.shuffle(null_picks)
    for pick in null_picks[:n_null_loci]:
        dosage = build_dosage(pick["gtp"], all_idx, _cfg.N_SAMPLES)
        dosage = np.where(np.isnan(dosage), np.nanmean(dosage), dosage)
        genetic += beta * dosage

        causal_variants.append(CausalVariant(
            variant_id=pick["vid"], chr=pick["chr"], pos=pick["pos"],
            af=pick["af"], beta=beta, gene=None, pathway=None,
        ))

    if len(causal_variants) == 0:
        raise ValueError("No causal variants could be selected")

    var_g = np.var(genetic)
    if var_g < 1e-10:
        raise ValueError("No genetic variance")
    var_e = var_g * (1 / h2_target - 1)
    noise = rng.normal(0, np.sqrt(var_e), n)
    phenotype = genetic + noise

    return SimulationResult(
        phenotype=phenotype,
        sample_ids=all_idx,
        scenario="F6_multi_locus_pathway",
        realized_h2=var_g / np.var(phenotype),
        causal_variants=causal_variants,
        noise_variance=var_e,
        config={
            "chr": chr, "n_pathway_loci": n_pathway_loci,
            "n_null_loci": n_null_loci, "beta": beta,
            "h2_target": h2_target, "pathway": pw_name,
            "pathway_genes": selected_genes,
        },
    )


FINEMAPPING_SCENARIOS = {
    "F1": simulate_F1_single_causal_in_ld,
    "F2": simulate_F2_two_independent_causal,
    "F6": simulate_F6_multi_locus_pathway,
}


def load_simulation_as_phenotype(
    conn: GraphGWASConnection,
    sim: SimulationResult,
    trait_name: str = "simulated",
) -> dict:
    """Load a simulated phenotype into Neo4j as gwas_value on Sample nodes."""
    valid = ~np.isnan(sim.phenotype)
    indices = sim.sample_ids[valid]
    values = sim.phenotype[valid]

    batch = [
        {"idx": int(indices[i]), "val": float(values[i])}
        for i in range(len(indices))
    ]

    # Set is_case/is_control to null (quantitative trait)
    for i in range(0, len(batch), 500):
        sub = batch[i:i + 500]
        conn.execute_write(
            """
            UNWIND $batch AS row
            MATCH (s:Sample {packed_index: row.idx})
            SET s.gwas_value = row.val,
                s.gwas_phenotype = $trait,
                s.is_case = null,
                s.is_control = null
            """,
            {"batch": sub, "trait": trait_name},
        )

    return {
        "n_samples": len(values),
        "scenario": sim.scenario,
        "realized_h2": sim.realized_h2,
        "trait_name": trait_name,
    }
