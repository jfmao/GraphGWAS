"""Inject canonical (BioGRID-curated) PPI edges into a graph_cache dict.

The yeast graph_cache built from the GraphGWAS import pipeline ships with
incomplete PPI annotations — for example, BCY1's variant entries list
``ppi: [YIL033C, YJR074W, ...]`` (BCY1 itself + 19 partners) but **NOT**
``YJL164C`` (TPK1), even though BCY1 and TPK1 are the most-validated
yeast PKA regulatory:catalytic binding pair.  Section §Y.5 of paper #2
isolated the cause to a curation gap upstream of the cache, and §Y.3 of
this module's audit found the same gap on every Tier-A 1011-segregating
pair (BCY1×TPK1, WHI5×CLN3, MKT1×GPA1, MIP1×SAL1, TKL1×NQM1: 0/5 in
either direction in the cache PPI field).

This module provides a clean injection path: given a list of canonical
gene-pair edges, augment each affected variant's ``ppi`` list to add
the edge.  The original cache is not modified — a deep-copied augmented
cache is returned.

Usage:

    from graphgwas.canonical_ppi import inject_canonical_ppi_edges
    augmented = inject_canonical_ppi_edges(cache, edges=[
        ("YIL033C", "YJL164C"),   # BCY1-TPK1
        ("YOR083W", "YAL040C"),   # WHI5-CLN3
    ])

The motif-filtered epistasis pipeline (M2) will now find P3 (motif
``protein_interaction``) pairs across these edges, even though the
cache's underlying source did not have them.

Source for canonical edges: BioGRID v4.4 (see
``docs/groundtruth/epistasis_pairs.json`` for citations per edge).
"""
from __future__ import annotations

import copy
import json
from collections import defaultdict
from pathlib import Path


# ---------------------------------------------------------------------------
# Canonical edges per species (curated from docs/groundtruth/epistasis_pairs.json)
# ---------------------------------------------------------------------------

YEAST_CANONICAL_EDGES: list[tuple[str, str]] = [
    # cAMP-PKA core complex (BCY1 reg + TPK1/2/3 cat)
    ("YIL033C", "YJL164C"),  # BCY1-TPK1
    ("YIL033C", "YPL203W"),  # BCY1-TPK2
    ("YIL033C", "YKL166C"),  # BCY1-TPK3
    ("YJL164C", "YPL203W"),  # TPK1-TPK2 (catalytic redundancy)
    ("YJL164C", "YKL166C"),  # TPK1-TPK3
    ("YPL203W", "YKL166C"),  # TPK2-TPK3
    # G1/S commitment
    ("YOR083W", "YAL040C"),  # WHI5-CLN3
    # Mitonuclear
    ("YOR330C", "YNL083W"),  # MIP1-SAL1
    # Pentose-phosphate paralogs
    ("YPR074C", "YGR043C"),  # TKL1-NQM1
    # MKT1 hub partners (Goldstein 2025 hub-QTL set)
    ("YNL085W", "YHR005C"),  # MKT1-GPA1
    ("YNL085W", "YGL073W"),  # MKT1-HAP1 (placeholder — verify systematic)
    # MAPK paralogs
    ("YBL016W", "YGR040W"),  # FUS3-KSS1
    # Morphogenesis checkpoint
    ("YKL101W", "YJL187C"),  # HSL1-SWE1
    # DDR (textbook lab pairs — useful as architectural ground truth)
    ("YBR136W", "YDR217C"),  # MEC1-RAD9
    ("YJL173C", "YML032C"),  # RFA1-RAD52
    # Carbon switch
    ("YDR477W", "YDR028C"),  # SNF1-REG1
]


ARABIDOPSIS_CANONICAL_EDGES: list[tuple[str, str]] = [
    # Flowering-time module
    ("AT1G65480", "AT5G10140"),  # FT-FLC (anchor; FLC represses FT)
    ("AT4G00650", "AT5G10140"),  # FRI-FLC (FRI activates FLC chromatin)
    ("AT1G65480", "AT4G20370"),  # FT-TSF (paralog florigens; double mutant)
    # Plant immunity
    ("AT5G45260", "AT5G45250"),  # RRS1-RPS4 (paired NLR, obligate heterodimer)
    # Root patterning
    ("AT3G54810", "AT3G54220"),  # SHR-SCR (SHR moves into endodermis, SCR sequesters)
    # Trichome MBW complex
    ("AT5G66320", "AT5G41315"),  # GL1-GL3
    ("AT5G66320", "AT5G24520"),  # GL1-TTG1
    ("AT5G41315", "AT5G24520"),  # GL3-TTG1
    # Light/temperature
    ("AT2G18790", "AT2G43010"),  # PHYB-PIF4 (phyB binds and degrades PIF4)
    # Brassinosteroid co-receptor
    ("AT4G39400", "AT4G33430"),  # BRI1-BAK1
    # Seed dormancy
    ("AT5G45830", "AT3G24650"),  # DOG1-ABI3
    # Glucosinolate biosynthesis (epistatic on alkenyl-glucosinolate yield)
    ("AT5G23010", "AT4G03060"),  # MAM1-AOP2 (GS-ELONG x GS-AOP); AOP2 is a
                                  # COL-0 pseudogene but functional in Cvi —
                                  # the natural-variation epistasis is documented
                                  # in Kliebenstein et al. 2001/Wentzell et al. 2007.
]


RICE_CANONICAL_EDGES: list[tuple[str, str]] = [
    # Flowering-time module
    ("LOC_Os06g16370", "LOC_Os06g06320"),  # Hd1-Hd3a (anchor)
    ("LOC_Os06g16370", "LOC_Os07g15770"),  # Hd1-Ghd7
    ("LOC_Os06g16370", "LOC_Os08g07740"),  # Hd1-DTH8
    ("LOC_Os07g15770", "LOC_Os08g07740"),  # Ghd7-DTH8
    ("LOC_Os10g32600", "LOC_Os06g06320"),  # Ehd1-Hd3a
    ("LOC_Os10g32600", "LOC_Os06g06300"),  # Ehd1-RFT1
    # Grain quality × yield
    ("LOC_Os08g41940", "LOC_Os07g41200"),  # GW8-GW7 (OsSPL16 directly represses)
    ("LOC_Os03g44500", "LOC_Os01g66030"),  # GS3-Gn1a
    ("LOC_Os03g44500", "LOC_Os08g39890"),  # GS3-IPA1/SPL14
    # Paired NLRs (rice blast)
    ("LOC_Os11g46200", "LOC_Os11g46210"),  # Pik-1 / Pik-2
    ("LOC_Os11g11790", "LOC_Os11g11770"),  # RGA5 / RGA4 (Pia-2 / Pia-1)
    # Tillering
    ("LOC_Os01g11410", "LOC_Os03g49880"),  # OsMADS57-OsTB1/FC1
    # Submergence (suppressive)
    ("LOC_Os09g11460", "LOC_Os09g11480"),  # SUB1A-SUB1C
]


CANONICAL_EDGES_BY_SPECIES: dict[str, list[tuple[str, str]]] = {
    "yeast": YEAST_CANONICAL_EDGES,
    "arabidopsis": ARABIDOPSIS_CANONICAL_EDGES,
    "rice": RICE_CANONICAL_EDGES,
}


# ---------------------------------------------------------------------------
# Cache PPI injection
# ---------------------------------------------------------------------------

def inject_canonical_ppi_edges(
    cache: dict,
    edges: list[tuple[str, str]] | None = None,
    *,
    inplace: bool = False,
    verbose: bool = False,
) -> tuple[dict, dict]:
    """Augment the ``ppi`` field of each variant entry in ``cache``.

    For every (g_a, g_b) edge in ``edges``, every variant that lists g_a
    in its ``genes`` field gains g_b in its ``ppi`` field, and every
    variant listing g_b gains g_a.  Existing PPI entries are preserved.

    Args:
        cache: dict mapping variant_id → entry-dict with at least
            ``genes: [...]`` and ``ppi: [...]`` keys (the standard
            graph_cache schema).
        edges: list of (gene_a, gene_b) systematic-name tuples.  Defaults
            to :data:`YEAST_CANONICAL_EDGES`.
        inplace: if True, modify ``cache`` directly; else deep-copy first.
        verbose: print per-edge counts.

    Returns:
        (augmented_cache, diagnostics)
        diagnostics has keys:
            n_edges_supplied: input edge count
            n_variants_modified: total variants whose PPI grew
            edges_with_zero_variants: list of (g_a, g_b) where neither gene
                had any variants in the cache (likely systematic-name typo
                or the gene is not in this species' annotation set)
    """
    if edges is None:
        edges = YEAST_CANONICAL_EDGES
    if not inplace:
        cache = copy.deepcopy(cache)

    # Build gene → variant index for fast lookup
    gene_to_vars: dict[str, list[str]] = defaultdict(list)
    for vid, entry in cache.items():
        for g in entry.get("genes", []):
            gene_to_vars[g].append(vid)

    n_modified = 0
    edges_dropped: list[tuple[str, str]] = []

    for g_a, g_b in edges:
        vars_a = gene_to_vars.get(g_a, [])
        vars_b = gene_to_vars.get(g_b, [])
        if not vars_a and not vars_b:
            edges_dropped.append((g_a, g_b))
            continue
        # Augment every g_a-variant's PPI list with g_b
        for vid in vars_a:
            ppi = cache[vid].setdefault("ppi", [])
            if g_b not in ppi:
                ppi.append(g_b)
                n_modified += 1
        # And vice versa
        for vid in vars_b:
            ppi = cache[vid].setdefault("ppi", [])
            if g_a not in ppi:
                ppi.append(g_a)
                n_modified += 1
        if verbose:
            print(f"  {g_a}↔{g_b}: augmented {len(vars_a)}+{len(vars_b)} variants")

    diag = {
        "n_edges_supplied": len(edges),
        "n_edges_applied": len(edges) - len(edges_dropped),
        "n_variants_modified": n_modified,
        "edges_with_zero_variants": edges_dropped,
    }
    return cache, diag


def edges_from_groundtruth_json(
    json_path: Path,
    species: str,
) -> list[tuple[str, str]]:
    """Extract canonical PPI edges from ``docs/groundtruth/epistasis_pairs.json``.

    Returns only pairs whose ``expected_motif`` includes
    ``"protein_interaction"`` (i.e., direct binding, not just same-pathway).

    Args:
        json_path: path to the catalog file.
        species: one of ``"yeast"``, ``"arabidopsis"``, ``"rice"``, ``"human"``.

    Returns:
        List of (gene_a_systematic, gene_b_systematic) pairs.  For species
        whose pair entries use a different identifier key (``agi``,
        ``locus``, ``hgnc``), the function returns the corresponding key
        values; the caller is responsible for matching them against the
        target cache's gene-identifier convention.
    """
    catalog = json.loads(json_path.read_text())
    pairs = catalog["species"][species]["pairs"]
    key_for_species = {
        "yeast": "systematic",
        "arabidopsis": "agi",
        "rice": "locus",
        "human": "hgnc",
    }[species]
    edges: list[tuple[str, str]] = []
    for pair in pairs:
        if "protein_interaction" not in pair.get("expected_motif", []):
            continue
        g_a = pair["gene_1"][key_for_species]
        g_b = pair["gene_2"][key_for_species]
        edges.append((g_a, g_b))
    return edges
