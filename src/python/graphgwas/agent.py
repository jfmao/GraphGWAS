"""LangGraph GWAS Agent — Natural Language Interface to GraphGWAS.

Phase 5: ReAct-style agent that translates natural language GWAS queries
into GraphGWAS operations. Combines structured graph retrieval (Cypher)
with statistical analysis tools and LLM interpretation.

The agent can:
- Run single-locus GWAS on specified regions
- Find epistatic modules in a region
- Extract disease architecture via max-flow
- Fine-map a locus using LD centrality
- Query stored results and interpret them
- Generate hypotheses about genetic architecture

Usage:
    from graphgwas.agent import create_agent, run_query
    agent = create_agent(conn)
    result = run_query(agent, "What variants are associated with this disease?")
"""

from __future__ import annotations

import json
from typing import Annotated

try:
    from langchain_core.messages import HumanMessage, SystemMessage, AIMessage
    from langchain_core.tools import tool
    from langchain_anthropic import ChatAnthropic
    from langgraph.prebuilt import create_react_agent
    HAS_LANGGRAPH = True
except ImportError:
    HAS_LANGGRAPH = False

from .db import GraphGWASConnection


# ---------------------------------------------------------------------------
# 5.1 Tool Definitions
# ---------------------------------------------------------------------------

def _make_tools(conn: GraphGWASConnection):
    """Create tool functions bound to a database connection."""

    @tool
    def run_association_scan(chromosome: str, start: int, end: int,
                            method: str = "auto") -> str:
        """Run single-locus GWAS on a genomic region.

        Args:
            chromosome: e.g. 'chr19'
            start: start position
            end: end position
            method: 'auto', 'chi2', 'fisher', 'logistic', 'firth'

        Returns top hits as formatted string.
        """
        from .assoc import single_locus_scan
        results = single_locus_scan(conn, chromosome, start, end,
                                     method=method, verbose=False)
        if not results:
            return "No variants found in region."

        # Sort by p-value and return top hits
        results.sort(key=lambda r: r["p_value"])
        top = results[:20]
        lines = [f"Scanned {len(results)} variants in {chromosome}:{start}-{end}"]
        sig = sum(1 for r in results if r["p_value"] < 5e-8)
        lines.append(f"Genome-wide significant: {sig}")
        lines.append("")
        for r in top:
            lines.append(f"  {r['variantId']} p={r['p_value']:.2e} "
                        f"beta={r.get('beta', 0):.3f} OR={r.get('odds_ratio', 1):.2f} "
                        f"AF={r.get('af', 0):.4f} method={r.get('method', 'unknown')}")
        return "\n".join(lines)

    @tool
    def find_epistatic_modules(chromosome: str, start: int, end: int,
                               min_cocarriers: int = 5) -> str:
        """Find epistatic variant modules in a region using co-occurrence network analysis.

        This is a graph-native method impossible in matrix-based GWAS tools.
        It builds a variant co-occurrence graph and detects communities.

        Args:
            chromosome: e.g. 'chr19'
            start: start position
            end: end position
            min_cocarriers: minimum case co-carriers for edge
        """
        from .epistasis import build_cooccurrence_graph, detect_epistatic_modules
        G = build_cooccurrence_graph(conn, chromosome, start, end,
                                      min_cocarriers=min_cocarriers, verbose=False)
        modules = detect_epistatic_modules(G)

        lines = [f"Co-occurrence graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges"]
        lines.append(f"Epistatic modules found: {len(modules)}")
        for i, m in enumerate(modules[:5]):
            lines.append(f"  Module {i+1}: {len(m)} variants — {', '.join(m[:5])}")
        return "\n".join(lines)

    @tool
    def extract_disease_architecture(chromosome: str = None,
                                     af_threshold: float = 0.05) -> str:
        """Extract disease architecture using max-flow/min-cut analysis.

        Finds the most parsimonious set of variants, genes, and pathways
        connecting phenotype to biology through the knowledge graph.

        Args:
            chromosome: optional chromosome filter (None = genome-wide)
            af_threshold: allele frequency threshold for rare variants
        """
        from .flow import build_flow_network, compute_pathway_flows
        G = build_flow_network(conn, chr=chromosome, af_threshold=af_threshold,
                               verbose=False)
        if G.number_of_nodes() < 5:
            return "Flow network too small. Try larger af_threshold or different region."

        flows = compute_pathway_flows(G, verbose=False)
        lines = [f"Flow network: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges"]
        lines.append(f"Top pathways by phenotype flow:")
        for f in flows[:10]:
            lines.append(f"  {f['pathway']}: flow={f['flow_value']:.2f}, "
                        f"min-cut={f['n_cut_elements']} elements")
            if f.get('min_cut_genes'):
                lines.append(f"    Key genes: {', '.join(f['min_cut_genes'][:5])}")
        return "\n".join(lines)

    @tool
    def finemap_locus(chromosome: str, start: int, end: int,
                      r2_threshold: float = 0.3) -> str:
        """Fine-map a locus using LD graph centrality.

        The causal variant is the hub with highest betweenness centrality
        in the LD graph. This replaces hours of Bayesian computation
        with millisecond graph centrality analysis.

        Args:
            chromosome: e.g. 'chr19'
            start: start position
            end: end position
            r2_threshold: LD threshold for graph edges
        """
        from .finemapping import finemap_locus as _finemap
        results = _finemap(conn, chromosome, start, end,
                          r2_threshold=r2_threshold, verbose=False)
        if not results:
            return "No variants found in locus."

        lines = ["Credible set (centrality-based fine-mapping):"]
        for r in results[:10]:
            p_str = f" p={r['p_value']:.2e}" if 'p_value' in r else ""
            lines.append(f"  #{r['rank']} {r['variantId']} "
                        f"centrality={r['centrality_score']:.4f} "
                        f"neighbors={r['n_ld_neighbors']}{p_str}")
        return "\n".join(lines)

    @tool
    def run_mpat_gene_test(chromosome: str = None,
                           genes: str = None) -> str:
        """Run MPAT (Message-Passing Association Test) on genes.

        Aggregates per-variant z-scores through gene→pathway topology
        with LD correction. Provides gene-level p-values.

        Args:
            chromosome: chromosome to scan (all genes on this chr)
            genes: comma-separated gene symbols (e.g. 'APOE,TOMM40,BCAM')
        """
        from .mpat import mpat_scan
        gene_list = genes.split(",") if genes else None
        results = mpat_scan(conn, chr=chromosome, genes=gene_list, verbose=False)

        lines = [f"MPAT gene-level results ({len(results)} genes tested):"]
        for r in results[:20]:
            lines.append(f"  {r['gene']}: Z={r['z_gene']:.3f} p={r['p_value']:.2e} "
                        f"({r['n_variants']} variants)")
        return "\n".join(lines)

    @tool
    def query_graph(cypher_query: str) -> str:
        """Execute a Cypher query against the Neo4j graph database.

        Use this for custom queries about variants, genes, pathways,
        or stored association results.

        Args:
            cypher_query: a READ-ONLY Cypher query
        """
        # Safety: only allow read queries
        query_lower = cypher_query.lower().strip()
        write_keywords = ['create', 'merge', 'delete', 'set ', 'remove', 'drop']
        if any(kw in query_lower for kw in write_keywords):
            return "ERROR: Only read-only queries are allowed."

        try:
            result = conn.execute_read(cypher_query)
            data = result.data()
            if not data:
                return "Query returned no results."
            # Format as JSON (limit output)
            return json.dumps(data[:50], indent=2, default=str)
        except Exception as e:
            return f"Query error: {e}"

    @tool
    def get_database_status() -> str:
        """Get current database status: node counts, phenotype info, GWAS runs."""
        from .schema import audit_schema
        report = audit_schema(conn)
        lines = ["Database Status:"]
        for k, v in report.items():
            lines.append(f"  {k}: {v:,}" if isinstance(v, int) else f"  {k}: {v}")
        return "\n".join(lines)

    return [
        run_association_scan,
        find_epistatic_modules,
        extract_disease_architecture,
        finemap_locus,
        run_mpat_gene_test,
        query_graph,
        get_database_status,
    ]


# ---------------------------------------------------------------------------
# 5.2 Agent Creation
# ---------------------------------------------------------------------------

SYSTEM_PROMPT = """You are GraphGWAS Agent — an AI assistant for graph-native genome-wide association studies.

You have access to a Neo4j database containing 70.7M genetic variants from 3,202 human samples (1000 Genomes Phase 3), with full functional annotation (91,973 genes, 2,279 pathways, 18,858 GO terms).

Your capabilities (unique to GraphGWAS — impossible in matrix tools like PLINK):

1. **Single-locus GWAS** — standard association testing (chi-squared, Fisher, Firth logistic)
2. **Epistasis detection** — variant co-occurrence networks + community detection
3. **Disease architecture** — max-flow/min-cut from phenotype through variant→gene→pathway
4. **LD fine-mapping** — graph centrality on LD topology (milliseconds, not hours)
5. **MPAT gene-level** — message-passing association test with LD-corrected analytical null
6. **Custom Cypher queries** — direct graph database access

When a user asks about genetic associations:
- Start with single-locus GWAS to identify significant regions
- If few hits: try MPAT gene-level or disease architecture (flow analysis)
- For significant loci: fine-map with LD centrality
- For rare disease with few cases: use epistasis and flow analysis
- Always report effect sizes, p-values, and biological context (genes, pathways)
- Interpret results in the context of known biology

Important: all p-values are from the simulated Alzheimer's-like phenotype based on the APOE locus.
The phenotype was designed with rs429358 (chr19:45411941) as the causal variant.
"""


def create_agent(conn: GraphGWASConnection,
                 model: str = "claude-sonnet-4-20250514"):
    """Create a LangGraph ReAct agent with GraphGWAS tools.

    Args:
        conn: active database connection.
        model: Anthropic model to use.

    Returns:
        LangGraph agent (compiled graph).
    """
    if not HAS_LANGGRAPH:
        raise ImportError(
            "LangGraph required. pip install langgraph langchain-core langchain-anthropic"
        )

    tools = _make_tools(conn)
    llm = ChatAnthropic(model=model, temperature=0)

    agent = create_react_agent(
        llm,
        tools,
        prompt=SystemMessage(content=SYSTEM_PROMPT),
    )
    return agent


def run_query(agent, query: str, verbose: bool = True) -> str:
    """Run a natural language query through the GWAS agent.

    Args:
        agent: LangGraph agent from create_agent.
        query: natural language question about genetics/GWAS.

    Returns:
        Agent's response as string.
    """
    if not HAS_LANGGRAPH:
        raise ImportError("LangGraph required")

    messages = [HumanMessage(content=query)]
    result = agent.invoke({"messages": messages})

    # Extract final response
    final_messages = result.get("messages", [])
    if final_messages:
        response = final_messages[-1].content
        if verbose:
            print(f"\nAgent Response:\n{response}")
        return response
    return "No response generated."


# ---------------------------------------------------------------------------
# 5.3 Standalone query function (no LLM required)
# ---------------------------------------------------------------------------

def interpret_results(conn: GraphGWASConnection, run_id: str,
                      verbose: bool = True) -> dict:
    """Structured interpretation of GWAS results without LLM.

    Provides:
    - Summary statistics (n_tested, n_significant, lambda_GC)
    - Top loci with gene annotations
    - Pathway enrichment via stored results
    - Recommendations for follow-up analyses

    This is the non-LLM fallback for environments without API keys.
    """
    from .results import query_top_hits, query_manhattan_data

    # Get results
    top_hits = query_top_hits(conn, run_id, p_threshold=5e-8, limit=50)
    all_data = query_manhattan_data(conn, run_id)

    n_tested = len(all_data)
    n_significant = len(top_hits)

    # Lambda GC
    import numpy as np
    from scipy import stats as sp_stats
    p_values = [d["log10p"] for d in all_data if d.get("log10p")]
    if p_values:
        chi2_obs = [sp_stats.chi2.isf(10 ** (-lp), 1) for lp in p_values if lp > 0]
        lambda_gc = float(np.median(chi2_obs) / 0.4549) if chi2_obs else 1.0
    else:
        lambda_gc = None

    # Unique loci (clump by 500kb windows)
    loci = []
    if top_hits:
        sorted_hits = sorted(top_hits, key=lambda h: h["p_value"])
        used_positions = []
        for hit in sorted_hits:
            pos = hit.get("pos", 0)
            chr_val = hit.get("chr", "")
            # Check if within 500kb of an existing locus
            is_new = True
            for used_chr, used_pos in used_positions:
                if chr_val == used_chr and abs(pos - used_pos) < 500000:
                    is_new = False
                    break
            if is_new:
                used_positions.append((chr_val, pos))
                loci.append(hit)

    interpretation = {
        "summary": {
            "n_variants_tested": n_tested,
            "n_genome_wide_significant": n_significant,
            "n_independent_loci": len(loci),
            "lambda_gc": lambda_gc,
        },
        "top_loci": loci[:20],
        "recommendations": [],
    }

    # Generate recommendations
    if n_significant == 0:
        interpretation["recommendations"].append(
            "No genome-wide significant hits. Consider: "
            "(1) MPAT gene-level test for aggregated signal, "
            "(2) Disease architecture via max-flow for pathway-level signal, "
            "(3) Epistasis scan for multi-variant interactions."
        )
    elif n_significant > 0 and n_significant < 10:
        interpretation["recommendations"].append(
            f"Found {n_significant} significant variants in {len(loci)} loci. "
            "Recommended follow-up: LD fine-mapping per locus, "
            "MPAT for gene-level confirmation, "
            "epistasis scan to check for multi-variant architecture."
        )
    else:
        interpretation["recommendations"].append(
            f"Strong polygenic signal: {n_significant} significant variants. "
            "Consider: GNN for multi-locus prediction, "
            "disease architecture for pathway-level summary."
        )

    if lambda_gc and lambda_gc > 1.1:
        interpretation["recommendations"].append(
            f"Genomic inflation λ_GC = {lambda_gc:.3f} > 1.1. "
            "Consider adding PCA covariates or spectral correction."
        )

    if verbose:
        print(f"\n=== GWAS Interpretation (run: {run_id}) ===")
        print(f"Variants tested: {n_tested:,}")
        print(f"Genome-wide significant: {n_significant}")
        print(f"Independent loci: {len(loci)}")
        if lambda_gc:
            print(f"Lambda GC: {lambda_gc:.3f}")
        print(f"\nTop loci:")
        for locus in loci[:10]:
            genes = ", ".join(locus.get("genes", [])[:3]) or "intergenic"
            print(f"  {locus.get('chr', '')}:{locus.get('pos', '')} "
                  f"p={locus.get('p_value', 1):.2e} genes={genes}")
        print(f"\nRecommendations:")
        for rec in interpretation["recommendations"]:
            print(f"  - {rec}")

    return interpretation
