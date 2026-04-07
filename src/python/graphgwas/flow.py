"""Max-flow / min-cut from phenotype to biology.

Unifies variant testing, gene mapping, and pathway analysis into a single
graph computation. The min-cut IS the most parsimonious disease architecture.

Flow network: Source → Case Samples → Variants → Genes → Pathways → Sink
Edge capacities encode biological signal strength (1/AF for rare variants,
impact weight for consequences).

ONE statistical test (permutation of case labels) replaces the multi-step
pipeline of separate variant, gene, and pathway enrichment tests.
"""

from __future__ import annotations

import numpy as np
import networkx as nx

from . import config as _cfg
from .db import GraphGWASConnection
from .genotype import (
    build_carrier_set,
    get_phenotype_indices,
    variant_iterator,
)


def build_flow_network(conn: GraphGWASConnection,
                       chr: str | None = None,
                       af_threshold: float = 0.05,
                       consequence_filter: list[str] | None = None,
                       verbose: bool = True) -> nx.DiGraph:
    """Build directed flow network: Source → Sample → Variant → Gene → Pathway.

    Only includes:
    - Variants with af_total < af_threshold
    - Variants with qualifying HAS_CONSEQUENCE to Gene nodes
    - Case samples who carry qualifying variants (decoded from gt_packed)

    Edge capacities:
    - Source → Case Sample: 1.0
    - Case Sample → Variant: 1/AF (rare variants carry more signal)
    - Variant → Gene: impact_weight (HIGH=1.0, MODERATE=0.5, LOW=0.1)
    - Gene → Pathway: 1.0 / sqrt(pathway_gene_count)

    Returns: nx.DiGraph with 'capacity' on all edges.
    """
    if consequence_filter is None:
        consequence_filter = ["HIGH", "MODERATE"]

    case_idx, ctrl_idx = get_phenotype_indices(conn)
    n_case = len(case_idx)

    G = nx.DiGraph()
    G.add_node("SOURCE", type="source")

    impact_weights = {"HIGH": 1.0, "MODERATE": 0.5, "LOW": 0.1, "MODIFIER": 0.01}

    # Query variants with functional annotation
    consequence_list = ", ".join(f"'{c}'" for c in consequence_filter)
    chr_filter = f"AND v.chr = '{chr}'" if chr else ""

    result = conn.execute_read(
        f"""
        MATCH (v:Variant)-[r:HAS_CONSEQUENCE]->(g:Gene)
        WHERE v.af_total < $af_threshold
          AND r.impact IN [{consequence_list}]
          {chr_filter}
        OPTIONAL MATCH (g)-[:IN_PATHWAY]->(p:Pathway)
        RETURN v.variantId AS vid, v.af_total AS af, v.gt_packed AS gtp,
               r.impact AS impact, g.symbol AS gene,
               collect(DISTINCT p.name) AS pathways
        """,
        {"af_threshold": af_threshold},
    )

    variants_data = [dict(r) for r in result]
    if verbose:
        print(f"Flow network: {len(variants_data)} variant-gene annotations "
              f"(AF < {af_threshold}, consequence in {consequence_filter})")

    # Track which case samples are added
    case_samples_added = set()
    genes_added = set()
    pathways_added = set()

    for vdata in variants_data:
        gt_packed = vdata["gtp"]
        if gt_packed is None:
            continue

        vid = vdata["vid"]
        af = max(vdata["af"], 1e-6)
        impact = vdata["impact"]
        gene = vdata["gene"]
        pathways = [p for p in (vdata["pathways"] or []) if p]

        # Decode which case samples carry this variant
        case_carriers = build_carrier_set(gt_packed, case_idx, _cfg.N_SAMPLES)
        carrier_indices = np.where(case_carriers)[0]

        if len(carrier_indices) == 0:
            continue

        # Add variant node
        G.add_node(f"V:{vid}", type="variant", af=af, impact=impact)

        # Sample → Variant edges (for each case carrier)
        sample_variant_capacity = 1.0 / af  # rare variants carry more signal
        for ci in carrier_indices:
            sample_id = f"S:{int(case_idx[ci])}"
            if sample_id not in case_samples_added:
                G.add_node(sample_id, type="sample")
                G.add_edge("SOURCE", sample_id, capacity=1.0)
                case_samples_added.add(sample_id)
            G.add_edge(sample_id, f"V:{vid}", capacity=sample_variant_capacity)

        # Variant → Gene edge
        if gene:
            gene_node = f"G:{gene}"
            if gene not in genes_added:
                G.add_node(gene_node, type="gene")
                genes_added.add(gene)
            impact_w = impact_weights.get(impact, 0.1)
            # Use max capacity if multiple consequences exist
            if G.has_edge(f"V:{vid}", gene_node):
                G[f"V:{vid}"][gene_node]["capacity"] = max(
                    G[f"V:{vid}"][gene_node]["capacity"], impact_w)
            else:
                G.add_edge(f"V:{vid}", gene_node, capacity=impact_w)

            # Gene → Pathway edges
            for pathway in pathways:
                pw_node = f"P:{pathway}"
                if pathway not in pathways_added:
                    G.add_node(pw_node, type="pathway")
                    pathways_added.add(pathway)
                if not G.has_edge(gene_node, pw_node):
                    G.add_edge(gene_node, pw_node, capacity=1.0)

    if verbose:
        n_samples = len(case_samples_added)
        n_variants = sum(1 for _, d in G.nodes(data=True) if d.get("type") == "variant")
        print(f"  Network: {n_samples} case samples, {n_variants} variants, "
              f"{len(genes_added)} genes, {len(pathways_added)} pathways, "
              f"{G.number_of_edges()} edges")

    return G


def compute_pathway_flows(flow_network: nx.DiGraph,
                          verbose: bool = True) -> list[dict]:
    """Compute max-flow from Source to each pathway node.

    For each pathway:
    1. Add temporary sink connected to pathway
    2. Compute max-flow using Edmonds-Karp (networkx default)
    3. Extract min-cut: the bottleneck variant-gene set
    4. Remove temporary sink

    Returns list of dicts sorted by flow_value descending.
    """
    pathway_nodes = [n for n, d in flow_network.nodes(data=True)
                     if d.get("type") == "pathway"]

    if verbose:
        print(f"Computing flow to {len(pathway_nodes)} pathways...")

    results = []
    for pw_node in pathway_nodes:
        # Add temporary sink
        sink = "SINK_TEMP"
        flow_network.add_edge(pw_node, sink, capacity=float("inf"))

        try:
            flow_value, flow_dict = nx.maximum_flow(
                flow_network, "SOURCE", sink, capacity="capacity"
            )

            # Extract min-cut
            cut_value, (reachable, non_reachable) = nx.minimum_cut(
                flow_network, "SOURCE", sink, capacity="capacity"
            )

            # Identify cut edges (variants/genes in the cut)
            cut_variants = []
            cut_genes = []
            for u in reachable:
                for v in flow_network.successors(u):
                    if v in non_reachable:
                        node_data = flow_network.nodes[u]
                        if node_data.get("type") == "variant":
                            cut_variants.append(u.replace("V:", ""))
                        elif node_data.get("type") == "gene":
                            cut_genes.append(u.replace("G:", ""))

            results.append({
                "pathway": pw_node.replace("P:", ""),
                "flow_value": float(flow_value),
                "cut_value": float(cut_value),
                "min_cut_variants": cut_variants,
                "min_cut_genes": cut_genes,
                "n_cut_elements": len(cut_variants) + len(cut_genes),
            })

        except nx.NetworkXError:
            pass
        finally:
            flow_network.remove_node(sink)

    results.sort(key=lambda r: r["flow_value"], reverse=True)

    if verbose and results:
        print(f"  Top pathway: {results[0]['pathway']} (flow={results[0]['flow_value']:.2f})")

    return results


def test_pathway_flow(conn: GraphGWASConnection,
                      chr: str | None = None,
                      af_threshold: float = 0.05,
                      n_permutations: int = 100,
                      verbose: bool = True) -> list[dict]:
    """Permutation test for pathway flow significance.

    1. Build flow network with true case labels
    2. Compute observed flow to each pathway
    3. For each permutation: shuffle case/control labels, rebuild Source edges,
       recompute flows
    4. Empirical p-value per pathway

    Returns pathway results with p-values.
    """
    # Build the base flow network
    base_network = build_flow_network(conn, chr=chr, af_threshold=af_threshold, verbose=verbose)

    if base_network.number_of_nodes() < 5:
        if verbose:
            print("Flow network too small for testing.")
        return []

    # Observed flows
    observed_flows = compute_pathway_flows(base_network, verbose=verbose)
    observed_dict = {r["pathway"]: r["flow_value"] for r in observed_flows}

    if not observed_dict:
        return []

    # Permutation null
    case_idx, ctrl_idx = get_phenotype_indices(conn)
    all_sample_nodes = [n for n, d in base_network.nodes(data=True)
                        if d.get("type") == "sample"]
    n_case = len(case_idx)

    if verbose:
        print(f"Running {n_permutations} permutations...")

    rng = np.random.default_rng(42)
    null_flows = {pw: [] for pw in observed_dict}

    for perm in range(n_permutations):
        # Shuffle: randomly select n_case sample nodes as "cases"
        perm_case_nodes = set(rng.choice(all_sample_nodes, size=n_case, replace=False))

        # Rebuild Source edges: only connect Source to permuted case samples
        for sample_node in all_sample_nodes:
            if base_network.has_edge("SOURCE", sample_node):
                base_network.remove_edge("SOURCE", sample_node)
            if sample_node in perm_case_nodes:
                base_network.add_edge("SOURCE", sample_node, capacity=1.0)

        # Compute permuted flows
        perm_flows = compute_pathway_flows(base_network, verbose=False)
        perm_dict = {r["pathway"]: r["flow_value"] for r in perm_flows}

        for pw in null_flows:
            null_flows[pw].append(perm_dict.get(pw, 0.0))

    # Restore original Source edges
    for sample_node in all_sample_nodes:
        if base_network.has_edge("SOURCE", sample_node):
            base_network.remove_edge("SOURCE", sample_node)
        # Re-add original case connections
        base_network.add_edge("SOURCE", sample_node, capacity=1.0)

    # Compute p-values
    results = []
    for r in observed_flows:
        pw = r["pathway"]
        null_arr = np.array(null_flows.get(pw, [0.0]))
        p_value = float((np.sum(null_arr >= r["flow_value"]) + 1) / (n_permutations + 1))
        r["p_value"] = p_value
        r["n_permutations"] = n_permutations
        results.append(r)

    results.sort(key=lambda r: r["p_value"])

    if verbose:
        sig = sum(1 for r in results if r["p_value"] < 0.05)
        print(f"Flow analysis: {len(results)} pathways tested, {sig} significant (p < 0.05)")

    return results


def extract_disease_architecture(conn: GraphGWASConnection,
                                 chr: str | None = None,
                                 af_threshold: float = 0.05,
                                 p_threshold: float = 0.05,
                                 n_permutations: int = 100,
                                 verbose: bool = True) -> dict:
    """Extract the disease architecture: significant pathways + their min-cut.

    The union of all min-cuts = the disease architecture subgraph.
    """
    results = test_pathway_flow(conn, chr=chr, af_threshold=af_threshold,
                                n_permutations=n_permutations, verbose=verbose)

    sig_pathways = [r for r in results if r["p_value"] < p_threshold]
    all_variants = set()
    all_genes = set()

    for r in sig_pathways:
        all_variants.update(r.get("min_cut_variants", []))
        all_genes.update(r.get("min_cut_genes", []))

    architecture = {
        "significant_pathways": [r["pathway"] for r in sig_pathways],
        "disease_variants": sorted(all_variants),
        "disease_genes": sorted(all_genes),
        "n_pathways": len(sig_pathways),
        "n_variants": len(all_variants),
        "n_genes": len(all_genes),
        "pathway_details": sig_pathways,
    }

    if verbose:
        print(f"\nDisease architecture: {architecture['n_pathways']} pathways, "
              f"{architecture['n_variants']} variants, {architecture['n_genes']} genes")

    return architecture
