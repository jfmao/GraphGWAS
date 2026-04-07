"""Epistasis engine — variant co-occurrence network + community detection.

The first method impossible in any matrix-based GWAS tool.

Builds a sparse variant co-occurrence graph where edges represent
case-enriched co-carriage. Community detection (Louvain) identifies
epistatic modules — variant sets whose joint co-occurrence in cases
exceeds independent expectation. Permutation testing validates modules.

Complexity: O(V² × N/64) for V variants in a window using numpy bool AND.
Bounded by: max_variants cap, min_cocarriers threshold, window size ≤ 1 Mb.
"""

from __future__ import annotations

import numpy as np
import networkx as nx

from . import config as _cfg
from .db import GraphGWASConnection
from .genotype import (
    build_carrier_set,
    build_count_table,
    get_phenotype_indices,
    unpack_genotypes,
    variant_iterator,
)


def build_cooccurrence_graph(conn: GraphGWASConnection, chr: str,
                             start: int, end: int,
                             min_cocarriers: int = 3,
                             min_enrichment: float = 1.5,
                             max_variants: int = 2000,
                             verbose: bool = True) -> nx.Graph:
    """Build variant co-occurrence graph for a genomic window.

    Nodes = variants with case carrier count ≥ min_cocarriers.
    Edges = variant pairs with case co-carriers ≥ min_cocarriers
            AND case/control enrichment ≥ min_enrichment.
    Edge weight = enrichment score.

    Returns networkx Graph. Node attributes: variantId, pos, case_carriers, ctrl_carriers.
    Edge attributes: case_co, ctrl_co, enrichment.
    """
    case_idx, ctrl_idx = get_phenotype_indices(conn)
    n_case, n_ctrl = len(case_idx), len(ctrl_idx)

    if verbose:
        print(f"Building co-occurrence graph: {chr}:{start}-{end} "
              f"| {n_case} cases, {n_ctrl} controls")

    # Load variants and pre-filter by case carrier count
    variants = []
    for v in variant_iterator(conn, chr, start, end):
        gt_packed = v["gt_packed"]
        if gt_packed is None:
            continue
        case_carriers = build_carrier_set(gt_packed, case_idx, _cfg.N_SAMPLES)
        n_case_carriers = int(np.sum(case_carriers))
        if n_case_carriers < min_cocarriers:
            continue

        ctrl_carriers = build_carrier_set(gt_packed, ctrl_idx, _cfg.N_SAMPLES)
        n_ctrl_carriers = int(np.sum(ctrl_carriers))

        variants.append({
            "variantId": v["variantId"],
            "pos": v["pos"],
            "case_carriers": case_carriers,
            "ctrl_carriers": ctrl_carriers,
            "n_case_carriers": n_case_carriers,
            "n_ctrl_carriers": n_ctrl_carriers,
        })

        if len(variants) >= max_variants:
            break

    if verbose:
        print(f"  {len(variants)} variants after pre-filtering (case carriers ≥ {min_cocarriers})")

    if len(variants) < 2:
        return nx.Graph()

    # Build graph: all-pairs co-occurrence
    G = nx.Graph()
    eps = 0.5 / max(n_ctrl, 1)  # pseudocount for enrichment denominator

    for i, v in enumerate(variants):
        G.add_node(v["variantId"], pos=v["pos"],
                   case_carriers=v["n_case_carriers"],
                   ctrl_carriers=v["n_ctrl_carriers"])

    n_edges = 0
    for i in range(len(variants)):
        for j in range(i + 1, len(variants)):
            v1, v2 = variants[i], variants[j]

            # Case co-carriers: bitwise AND
            case_co = int(np.sum(v1["case_carriers"] & v2["case_carriers"]))
            if case_co < min_cocarriers:
                continue

            # Control co-carriers
            ctrl_co = int(np.sum(v1["ctrl_carriers"] & v2["ctrl_carriers"]))

            # Enrichment: (case_co / n_case) / (ctrl_co / n_ctrl + eps)
            case_rate = case_co / n_case
            ctrl_rate = ctrl_co / n_ctrl + eps
            enrichment = case_rate / ctrl_rate

            if enrichment >= min_enrichment:
                G.add_edge(v1["variantId"], v2["variantId"],
                           case_co=case_co, ctrl_co=ctrl_co,
                           enrichment=float(enrichment))
                n_edges += 1

    if verbose:
        print(f"  Co-occurrence graph: {G.number_of_nodes()} nodes, {n_edges} edges")

    return G


def detect_epistatic_modules(G: nx.Graph,
                             resolution: float = 1.0,
                             min_module_size: int = 2) -> list[list[str]]:
    """Run Louvain community detection on co-occurrence graph.

    Each community = a candidate epistatic module.

    Returns list of variant ID lists, sorted by module size (descending).
    """
    if G.number_of_nodes() < 2:
        return []

    communities = nx.community.louvain_communities(
        G, weight="enrichment", resolution=resolution, seed=42
    )

    # Filter by minimum size and sort
    modules = [sorted(c) for c in communities if len(c) >= min_module_size]
    modules.sort(key=len, reverse=True)
    return modules


def compute_module_score(case_carriers_list: list[np.ndarray],
                         ctrl_carriers_list: list[np.ndarray],
                         n_case: int, n_ctrl: int) -> float:
    """Compute mean pairwise case enrichment for a module.

    Higher score = stronger case-specific co-occurrence pattern.
    """
    k = len(case_carriers_list)
    if k < 2:
        return 0.0

    eps = 0.5 / max(n_ctrl, 1)
    enrichments = []

    for i in range(k):
        for j in range(i + 1, k):
            case_co = int(np.sum(case_carriers_list[i] & case_carriers_list[j]))
            ctrl_co = int(np.sum(ctrl_carriers_list[i] & ctrl_carriers_list[j]))
            case_rate = case_co / n_case
            ctrl_rate = ctrl_co / n_ctrl + eps
            enrichments.append(case_rate / ctrl_rate)

    return float(np.mean(enrichments))


def test_epistatic_module(conn: GraphGWASConnection,
                          module_variant_ids: list[str],
                          n_permutations: int = 1000,
                          verbose: bool = True) -> dict:
    """Permutation test for an epistatic module's significance.

    Shuffles case/control labels, recomputes module score under null.
    """
    case_idx, ctrl_idx = get_phenotype_indices(conn)
    n_case, n_ctrl = len(case_idx), len(ctrl_idx)
    all_idx = np.concatenate([case_idx, ctrl_idx])

    # Load carrier sets for module variants
    case_carriers_list = []
    ctrl_carriers_list = []
    all_carriers_list = []

    for vid in module_variant_ids:
        rec = conn.execute_read(
            "MATCH (v:Variant {variantId: $vid}) RETURN v.gt_packed AS gtp",
            {"vid": vid},
        ).single()
        if not rec or rec["gtp"] is None:
            continue
        gt_packed = rec["gtp"]
        case_carriers_list.append(build_carrier_set(gt_packed, case_idx, _cfg.N_SAMPLES))
        ctrl_carriers_list.append(build_carrier_set(gt_packed, ctrl_idx, _cfg.N_SAMPLES))
        all_carriers_list.append(build_carrier_set(gt_packed, all_idx, _cfg.N_SAMPLES))

    if len(case_carriers_list) < 2:
        return {"module_score": 0.0, "p_value": 1.0,
                "n_variants": len(module_variant_ids), "n_permutations": 0}

    # Observed module score
    observed_score = compute_module_score(case_carriers_list, ctrl_carriers_list, n_case, n_ctrl)

    # Permutation null
    rng = np.random.default_rng(42)
    null_scores = []
    n_all = n_case + n_ctrl

    for perm in range(n_permutations):
        # Shuffle case/control assignment
        perm_labels = rng.permutation(n_all)
        perm_case_mask = perm_labels < n_case  # first n_case get case label

        # Recompute carrier sets under permuted labels
        perm_case_carriers = [cs[perm_case_mask] for cs in all_carriers_list]
        perm_ctrl_carriers = [cs[~perm_case_mask] for cs in all_carriers_list]

        perm_score = compute_module_score(perm_case_carriers, perm_ctrl_carriers, n_case, n_ctrl)
        null_scores.append(perm_score)

    null_scores = np.array(null_scores)
    p_value = float((np.sum(null_scores >= observed_score) + 1) / (n_permutations + 1))

    if verbose:
        print(f"  Module ({len(module_variant_ids)} variants): "
              f"score={observed_score:.3f}, p={p_value:.4f}")

    return {
        "module_score": float(observed_score),
        "p_value": p_value,
        "n_variants": len(case_carriers_list),
        "n_permutations": n_permutations,
        "variant_ids": module_variant_ids,
    }


def epistasis_scan(conn: GraphGWASConnection, chr: str,
                   window_size: int = 1_000_000,
                   step_size: int = 500_000,
                   min_cocarriers: int = 3,
                   min_enrichment: float = 1.5,
                   n_permutations: int = 1000,
                   verbose: bool = True) -> list[dict]:
    """Sliding-window epistasis scan across a chromosome.

    For each window: build co-occurrence graph → detect communities →
    test significant modules.

    Returns list of significant epistatic modules with p-values.
    """
    # Get chromosome extent
    result = conn.execute_read(
        "MATCH (v:Variant) WHERE v.chr = $chr "
        "RETURN min(v.pos) AS min_pos, max(v.pos) AS max_pos",
        {"chr": chr},
    ).single()
    if not result:
        return []

    min_pos, max_pos = result["min_pos"], result["max_pos"]
    all_modules = []

    pos = min_pos
    window_idx = 0
    while pos < max_pos:
        w_start = pos
        w_end = min(pos + window_size, max_pos)
        window_idx += 1

        if verbose:
            print(f"\nWindow {window_idx}: {chr}:{w_start}-{w_end}")

        # Build co-occurrence graph
        G = build_cooccurrence_graph(
            conn, chr, w_start, w_end,
            min_cocarriers=min_cocarriers,
            min_enrichment=min_enrichment,
            verbose=verbose,
        )

        if G.number_of_nodes() < 2:
            pos += step_size
            continue

        # Detect communities
        modules = detect_epistatic_modules(G)

        if verbose:
            print(f"  Detected {len(modules)} modules")

        # Test each module
        for module in modules:
            result = test_epistatic_module(
                conn, module, n_permutations=n_permutations, verbose=verbose
            )
            result["chr"] = chr
            result["window_start"] = w_start
            result["window_end"] = w_end
            all_modules.append(result)

        pos += step_size

    # Sort by p-value
    all_modules.sort(key=lambda m: m["p_value"])

    if verbose:
        sig = sum(1 for m in all_modules if m["p_value"] < 0.05)
        print(f"\nEpistasis scan complete: {len(all_modules)} modules, {sig} significant (p < 0.05)")

    return all_modules
