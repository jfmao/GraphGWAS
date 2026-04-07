"""LD fine-mapping via graph centrality.

The causal variant is the hub through which all proxy associations flow
in the LD graph. Betweenness centrality replaces Bayesian posterior
probability computation — milliseconds instead of hours per locus.

The LD graph structure IS the prior that methods like SuSiE/FINEMAP
spend hours estimating from the LD matrix.
"""

from __future__ import annotations

import numpy as np
import networkx as nx

from . import config as _cfg
from .db import GraphGWASConnection
from .genotype import (
    build_dosage,
    get_all_indices,
    variant_iterator,
)


def build_ld_graph(conn: GraphGWASConnection, chr: str,
                   start: int, end: int,
                   r2_threshold: float = 0.3,
                   min_af: float = 0.001,
                   verbose: bool = True) -> nx.Graph:
    """Build LD graph for a genomic locus.

    Nodes = variants. Edges = r² > threshold.
    Edge weight = r² value.

    Uses Pearson correlation of dosage vectors for all variant pairs.
    """
    indices = get_all_indices(conn)

    # Load variants and dosage vectors
    variants = []
    dosages = []
    for v in variant_iterator(conn, chr, start, end):
        gt_packed = v["gt_packed"]
        if gt_packed is None:
            continue
        af = v.get("af_total", 0)
        if af < min_af or af > (1 - min_af):
            continue  # skip monomorphic

        d = build_dosage(gt_packed, indices, _cfg.N_SAMPLES)
        # Replace NaN with mean for correlation
        mean_d = np.nanmean(d)
        if np.isnan(mean_d):
            continue
        d = np.where(np.isnan(d), mean_d, d)

        variants.append({
            "variantId": v["variantId"],
            "pos": v["pos"],
            "af": af,
        })
        dosages.append(d)

    if verbose:
        print(f"LD graph: {len(variants)} variants in {chr}:{start}-{end}")

    if len(variants) < 2:
        return nx.Graph()

    # Compute pairwise correlation matrix
    dosage_matrix = np.array(dosages)  # (V, N)
    corr_matrix = np.corrcoef(dosage_matrix)
    r2_matrix = corr_matrix ** 2
    np.fill_diagonal(r2_matrix, 0)  # no self-edges

    # Build graph
    G = nx.Graph()
    for i, v in enumerate(variants):
        G.add_node(v["variantId"], pos=v["pos"], af=v["af"])

    n_edges = 0
    for i in range(len(variants)):
        for j in range(i + 1, len(variants)):
            r2 = r2_matrix[i, j]
            if r2 >= r2_threshold and not np.isnan(r2):
                G.add_edge(variants[i]["variantId"], variants[j]["variantId"],
                           r2=float(r2), weight=float(r2))
                n_edges += 1

    if verbose:
        print(f"  LD graph: {G.number_of_nodes()} nodes, {n_edges} edges (r² ≥ {r2_threshold})")

    return G


def centrality_finemapping(ld_graph: nx.Graph,
                           p_values: dict[str, float] | None = None,
                           method: str = "betweenness",
                           top_k: int = 10,
                           verbose: bool = True) -> list[dict]:
    """Fine-map a locus using LD graph centrality.

    The causal variant is the hub — highest centrality node.

    Args:
        ld_graph: LD graph from build_ld_graph.
        p_values: optional dict mapping variantId → p-value (for annotation).
        method: 'betweenness', 'degree', 'pagerank', 'eigenvector'.
        top_k: number of top variants to return.

    Returns list sorted by centrality score (descending).
    """
    if ld_graph.number_of_nodes() == 0:
        return []

    # Compute centrality
    if method == "betweenness":
        centrality = nx.betweenness_centrality(ld_graph, weight="r2")
    elif method == "degree":
        # Weighted degree: sum of r² to all neighbors
        centrality = {}
        for node in ld_graph.nodes():
            centrality[node] = sum(
                ld_graph[node][nbr].get("r2", 0) for nbr in ld_graph.neighbors(node)
            )
    elif method == "pagerank":
        centrality = nx.pagerank(ld_graph, weight="r2")
    elif method == "eigenvector":
        try:
            centrality = nx.eigenvector_centrality(ld_graph, weight="r2", max_iter=1000)
        except nx.PowerIterationFailedConvergence:
            centrality = nx.degree_centrality(ld_graph)
    else:
        raise ValueError(f"Unknown centrality method: {method}")

    # Build results
    results = []
    for vid, score in centrality.items():
        node_data = ld_graph.nodes[vid]
        result = {
            "variantId": vid,
            "centrality_score": float(score),
            "centrality_method": method,
            "pos": node_data.get("pos"),
            "af": node_data.get("af"),
            "n_ld_neighbors": ld_graph.degree(vid),
        }
        if p_values and vid in p_values:
            result["p_value"] = p_values[vid]
        results.append(result)

    results.sort(key=lambda r: r["centrality_score"], reverse=True)

    # Add rank
    for i, r in enumerate(results):
        r["rank"] = i + 1

    if verbose and results:
        print(f"\nFine-mapping credible set (top {min(top_k, len(results))}, method={method}):")
        for r in results[:top_k]:
            p_str = f"p={r['p_value']:.2e}" if "p_value" in r else ""
            print(f"  #{r['rank']} {r['variantId']:<50} "
                  f"centrality={r['centrality_score']:.4f} "
                  f"n_neighbors={r['n_ld_neighbors']} {p_str}")

    return results[:top_k]


def finemap_locus(conn: GraphGWASConnection, chr: str,
                  start: int, end: int,
                  r2_threshold: float = 0.3,
                  method: str = "betweenness",
                  run_id: str | None = None,
                  verbose: bool = True) -> list[dict]:
    """Complete fine-mapping pipeline for a locus.

    1. Build LD graph
    2. Optionally load p-values from stored AssociationResult nodes
    3. Compute centrality-based credible set
    """
    # Build LD graph
    ld_graph = build_ld_graph(conn, chr, start, end,
                              r2_threshold=r2_threshold, verbose=verbose)

    # Load p-values if a run_id is available
    p_values = None
    if run_id:
        result = conn.execute_read(
            """
            MATCH (ar:AssociationResult {run_id: $run_id})-[:FOR_VARIANT]->(v:Variant)
            WHERE v.chr = $chr AND v.pos >= $start AND v.pos <= $end
            RETURN v.variantId AS vid, ar.p_value AS p
            """,
            {"run_id": run_id, "chr": chr, "start": start, "end": end},
        )
        p_values = {r["vid"]: r["p"] for r in result}

    # Fine-map
    return centrality_finemapping(ld_graph, p_values=p_values,
                                  method=method, verbose=verbose)
