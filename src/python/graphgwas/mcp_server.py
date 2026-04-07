"""GraphGWAS MCP Server — Model Context Protocol interface.

Exposes GraphGWAS capabilities as MCP tools, enabling any MCP-compatible client
(Claude Desktop, Claude Code, Cursor, etc.) to run GWAS analyses conversationally.

This replaces the need for the hardcoded LangGraph agent — the LLM client IS the
agent, and GraphGWAS just exposes tools.

Usage:
    graphgwas mcp                    # stdio transport (default)
    graphgwas mcp --transport sse    # SSE transport for web clients
"""

from __future__ import annotations

import json

try:
    from mcp.server.fastmcp import FastMCP
    HAS_MCP = True
except ImportError:
    HAS_MCP = False

from .db import GraphGWASConnection


# ---------------------------------------------------------------------------
# Server instance
# ---------------------------------------------------------------------------

if HAS_MCP:
    mcp = FastMCP(
        "GraphGWAS",
        instructions=(
            "Graph-native GWAS platform on Neo4j. Provides single-locus GWAS, "
            "epistasis detection, disease architecture extraction, LD fine-mapping, "
            "gene-level MPAT, spectral correction, GNN association, and direct "
            "Cypher queries over a 70.7M-variant human genomics graph."
        ),
    )
else:
    mcp = None

# Module-level connection — initialized by create_mcp_server()
_conn: GraphGWASConnection | None = None


def _get_conn() -> GraphGWASConnection:
    if _conn is None:
        raise RuntimeError(
            "MCP server not initialized. Call create_mcp_server() first."
        )
    return _conn


# ---------------------------------------------------------------------------
# Tool 1: Association scan
# ---------------------------------------------------------------------------

if HAS_MCP:

    @mcp.tool()
    def run_association_scan(chromosome: str, start: int, end: int,
                             method: str = "auto") -> str:
        """Run single-locus GWAS on a genomic region.

        Scans all variants in the region and returns top hits with p-values,
        effect sizes, and allele frequencies.

        Args:
            chromosome: Chromosome name (e.g. 'chr19')
            start: Start position (bp)
            end: End position (bp)
            method: Statistical method — 'auto' (recommended), 'chi2', 'fisher',
                    'logistic', 'firth', or 'linear'

        Returns:
            Formatted text with variant count, significant hits, and top results.
        """
        from .assoc import single_locus_scan
        conn = _get_conn()
        results = single_locus_scan(conn, chromosome, start, end,
                                    method=method, verbose=False)
        if not results:
            return "No variants found in region."

        results.sort(key=lambda r: r["p_value"])
        top = results[:20]
        sig = sum(1 for r in results if r["p_value"] < 5e-8)

        lines = [f"Scanned {len(results)} variants in {chromosome}:{start}-{end}",
                 f"Genome-wide significant: {sig}", ""]
        for r in top:
            lines.append(
                f"  {r['variantId']} p={r['p_value']:.2e} "
                f"beta={r.get('beta', 0):.3f} OR={r.get('odds_ratio', 1):.2f} "
                f"AF={r.get('af', 0):.4f} method={r.get('method', 'unknown')}"
            )
        return "\n".join(lines)


    # -----------------------------------------------------------------------
    # Tool 2: Epistasis
    # -----------------------------------------------------------------------

    @mcp.tool()
    def find_epistatic_modules(chromosome: str, start: int, end: int,
                                min_cocarriers: int = 5) -> str:
        """Find epistatic variant modules using co-occurrence network analysis.

        Graph-native method impossible in matrix-based GWAS tools. Builds a
        variant co-occurrence graph and detects communities via Louvain.

        Args:
            chromosome: Chromosome name (e.g. 'chr19')
            start: Start position (bp)
            end: End position (bp)
            min_cocarriers: Minimum case co-carriers for an edge

        Returns:
            Co-occurrence graph stats and detected epistatic modules.
        """
        from .epistasis import build_cooccurrence_graph, detect_epistatic_modules
        conn = _get_conn()
        G = build_cooccurrence_graph(conn, chromosome, start, end,
                                     min_cocarriers=min_cocarriers, verbose=False)
        modules = detect_epistatic_modules(G)

        lines = [
            f"Co-occurrence graph: {G.number_of_nodes()} nodes, "
            f"{G.number_of_edges()} edges",
            f"Epistatic modules found: {len(modules)}",
        ]
        for i, m in enumerate(modules[:5]):
            lines.append(f"  Module {i+1}: {len(m)} variants — {', '.join(m[:5])}")
        return "\n".join(lines)


    # -----------------------------------------------------------------------
    # Tool 3: Disease architecture
    # -----------------------------------------------------------------------

    @mcp.tool()
    def extract_disease_architecture(chromosome: str = "",
                                     af_threshold: float = 0.05) -> str:
        """Extract disease architecture using max-flow/min-cut analysis.

        Finds the most parsimonious set of variants, genes, and pathways
        connecting phenotype to biology through the knowledge graph.

        Args:
            chromosome: Chromosome filter (empty string = genome-wide)
            af_threshold: Allele frequency threshold for rare variants

        Returns:
            Top pathways by phenotype flow with key genes.
        """
        from .flow import build_flow_network, compute_pathway_flows
        conn = _get_conn()
        chr_arg = chromosome if chromosome else None
        G = build_flow_network(conn, chr=chr_arg, af_threshold=af_threshold,
                               verbose=False)
        if G.number_of_nodes() < 5:
            return "Flow network too small. Try larger af_threshold or different region."

        flows = compute_pathway_flows(G, verbose=False)
        lines = [
            f"Flow network: {G.number_of_nodes()} nodes, "
            f"{G.number_of_edges()} edges",
            "Top pathways by phenotype flow:",
        ]
        for f in flows[:10]:
            lines.append(
                f"  {f['pathway']}: flow={f['flow_value']:.2f}, "
                f"min-cut={f['n_cut_elements']} elements"
            )
            if f.get("min_cut_genes"):
                lines.append(f"    Key genes: {', '.join(f['min_cut_genes'][:5])}")
        return "\n".join(lines)


    # -----------------------------------------------------------------------
    # Tool 4: Fine-mapping
    # -----------------------------------------------------------------------

    @mcp.tool()
    def finemap_locus(chromosome: str, start: int, end: int,
                      r2_threshold: float = 0.3) -> str:
        """Fine-map a locus using LD graph centrality.

        The causal variant is the hub with highest betweenness centrality
        in the LD graph. Replaces hours of Bayesian computation with
        millisecond graph centrality analysis.

        Args:
            chromosome: Chromosome name (e.g. 'chr19')
            start: Start position (bp)
            end: End position (bp)
            r2_threshold: LD threshold for graph edges

        Returns:
            Credible set ranked by centrality score.
        """
        from .finemapping import finemap_locus as _finemap
        conn = _get_conn()
        results = _finemap(conn, chromosome, start, end,
                           r2_threshold=r2_threshold, verbose=False)
        if not results:
            return "No variants found in locus."

        lines = ["Credible set (centrality-based fine-mapping):"]
        for r in results[:10]:
            p_str = f" p={r['p_value']:.2e}" if "p_value" in r else ""
            lines.append(
                f"  #{r['rank']} {r['variantId']} "
                f"centrality={r['centrality_score']:.4f} "
                f"neighbors={r['n_ld_neighbors']}{p_str}"
            )
        return "\n".join(lines)


    # -----------------------------------------------------------------------
    # Tool 5: MPAT gene-level
    # -----------------------------------------------------------------------

    @mcp.tool()
    def run_mpat_gene_test(chromosome: str = "", genes: str = "") -> str:
        """Run MPAT (Message-Passing Association Test) on genes.

        Aggregates per-variant z-scores through gene-pathway topology
        with LD correction. Provides gene-level p-values.

        Args:
            chromosome: Chromosome to scan (empty = use genes arg)
            genes: Comma-separated gene symbols (e.g. 'APOE,TOMM40,BCAM')

        Returns:
            Gene-level results with z-scores and p-values.
        """
        from .mpat import mpat_scan
        conn = _get_conn()
        gene_list = genes.split(",") if genes else None
        chr_arg = chromosome if chromosome else None
        results = mpat_scan(conn, chr=chr_arg, genes=gene_list, verbose=False)

        lines = [f"MPAT gene-level results ({len(results)} genes tested):"]
        for r in results[:20]:
            lines.append(
                f"  {r['gene']}: Z={r['z_gene']:.3f} p={r['p_value']:.2e} "
                f"({r['n_variants']} variants)"
            )
        return "\n".join(lines)


    # -----------------------------------------------------------------------
    # Tool 6: Cypher query
    # -----------------------------------------------------------------------

    @mcp.tool()
    def query_graph(cypher_query: str) -> str:
        """Execute a read-only Cypher query against the Neo4j graph database.

        Use this for custom queries about variants, genes, pathways,
        or stored association results.

        Args:
            cypher_query: A READ-ONLY Cypher query

        Returns:
            Query results as JSON (limited to 50 rows).
        """
        conn = _get_conn()
        query_lower = cypher_query.lower().strip()
        write_keywords = ["create", "merge", "delete", "set ", "remove", "drop"]
        if any(kw in query_lower for kw in write_keywords):
            return "ERROR: Only read-only queries are allowed."

        try:
            result = conn.execute_read(cypher_query)
            data = result.data()
            if not data:
                return "Query returned no results."
            return json.dumps(data[:50], indent=2, default=str)
        except Exception as e:
            return f"Query error: {e}"


    # -----------------------------------------------------------------------
    # Tool 7: Database status
    # -----------------------------------------------------------------------

    @mcp.tool()
    def get_database_status() -> str:
        """Get current database status: node counts, phenotype info, GWAS runs."""
        from .schema import audit_schema
        conn = _get_conn()
        report = audit_schema(conn)
        lines = ["Database Status:"]
        for k, v in report.items():
            lines.append(f"  {k}: {v:,}" if isinstance(v, int) else f"  {k}: {v}")
        return "\n".join(lines)


    # -----------------------------------------------------------------------
    # Tool 8: Spectral correction
    # -----------------------------------------------------------------------

    @mcp.tool()
    def run_spectral_correction(n_components: int = 10,
                                af_threshold: float = 0.01,
                                chromosome: str = "") -> str:
        """Compute spectral phenotype correction from rare-variant similarity.

        Generalizes PCA-based population structure correction by decomposing
        the phenotype signal in the frequency domain of a rare-variant-based
        sample similarity graph.

        Args:
            n_components: Number of spectral components to compute
            af_threshold: Allele frequency threshold for rare variants
            chromosome: Single chromosome for fast testing (empty = all)

        Returns:
            Spectral decomposition summary with eigenvalue range and variance explained.
        """
        from .spectral import spectral_gwas_correction
        conn = _get_conn()
        chroms = [chromosome] if chromosome else None
        result = spectral_gwas_correction(
            conn, n_components=n_components,
            af_threshold=af_threshold, chromosomes=chroms, verbose=False,
        )
        lines = [
            f"Spectral correction computed:",
            f"  Components: {result['n_components']}",
            f"  Eigenvalue range: [{result['eigenvalues'][0]:.4f}, "
            f"{result['eigenvalues'][-1]:.4f}]",
        ]
        return "\n".join(lines)


    # -----------------------------------------------------------------------
    # Tool 9: GNN pipeline
    # -----------------------------------------------------------------------

    @mcp.tool()
    def run_gnn_pipeline(chromosome: str, start: int = 0, end: int = 0,
                         hidden_channels: int = 64, epochs: int = 100) -> str:
        """Run full GNN pipeline: export graph, train model, explain variants, import embeddings.

        Heterogeneous Graph Neural Network over the Variant-Sample-Gene graph
        learns multi-locus epistatic architecture through message passing.

        Args:
            chromosome: Chromosome name (e.g. 'chr19')
            start: Start position (0 = use full chromosome)
            end: End position (0 = use full chromosome)
            hidden_channels: GNN hidden dimension
            epochs: Training epochs

        Returns:
            Best AUROC, number of imported embeddings, top variant attributions.
        """
        from .gnn import run_gnn_pipeline as _run_gnn
        conn = _get_conn()
        start_pos = start if start > 0 else None
        end_pos = end if end > 0 else None
        result = _run_gnn(conn, chromosome, start_pos, end_pos,
                          hidden_channels=hidden_channels, epochs=epochs,
                          verbose=False)
        if not result:
            return "No variants in region. Aborting."

        lines = [
            f"GNN Pipeline Results:",
            f"  Best AUROC: {result['best_auroc']:.4f}",
            f"  Embeddings imported: {result['n_embeddings_imported']}",
            "",
            "Top variants by GNN importance:",
        ]
        for v in result.get("top_variants", [])[:10]:
            lines.append(
                f"  #{v['rank']} {v['variantId']} score={v['importance_score']:.6f}"
            )
        return "\n".join(lines)


    # -----------------------------------------------------------------------
    # Tool 10: Heritability estimation
    # -----------------------------------------------------------------------

    @mcp.tool()
    def estimate_heritability(method: str = "spectral",
                              chromosome: str = "",
                              af_threshold: float = 0.01,
                              permutations: int = 200) -> str:
        """Estimate heritability using graph-native methods.

        Novel multi-resolution heritability decomposition that no classical
        tool (GCTA, LDSC) can provide.  Decomposes h² by biological scale
        (variant → gene → pathway) using graph topology.

        Args:
            method: 'spectral' (eigenspectrum), 'conductance' (graph cut),
                    'flow' (pathway decomposition), 'multi' (all three scales),
                    'report' (run spectral + conductance + flow)
            chromosome: Optional chromosome filter (empty = genome-wide)
            af_threshold: Allele frequency threshold for rare variants
            permutations: Number of permutations (for conductance/flow)

        Returns:
            Heritability estimates with decomposition details.
        """
        conn = _get_conn()
        chr_arg = chromosome if chromosome else None
        chroms = [chromosome] if chromosome else None

        if method == "spectral":
            from .heritability import spectral_heritability
            r = spectral_heritability(conn, af_threshold=af_threshold,
                                      chromosomes=chroms, verbose=False)
            return (f"Spectral Heritability:\n"
                    f"  h²_spectral = {r['h2']:.4f}\n"
                    f"  Components: {r['optimal_cutoff']} / {r['n_components']}")

        elif method == "conductance":
            from .heritability import conductance_heritability
            r = conductance_heritability(conn, af_threshold=af_threshold,
                                         chromosomes=chroms,
                                         n_permutations=permutations,
                                         verbose=False)
            return (f"Conductance Heritability:\n"
                    f"  h²_conductance = {r['h2_conductance']:.4f}\n"
                    f"  p-value = {r['p_value']:.4f}\n"
                    f"  Conductance: {r['conductance']:.6f} "
                    f"(null: {r['phi_null_mean']:.6f})")

        elif method == "flow":
            from .heritability import flow_heritability
            r = flow_heritability(conn, chr=chr_arg,
                                  af_threshold=af_threshold,
                                  n_permutations=permutations,
                                  verbose=False)
            lines = [f"Flow Heritability Decomposition:",
                     f"  h²_flow = {r['h2_flow_total']:.4f} (p={r['p_value']:.4f})",
                     f"  Top pathway contributions:"]
            for c in r.get("pathway_contributions", [])[:5]:
                lines.append(f"    {c['pathway']}: h²={c['h2_contribution']:.4f} "
                             f"({c['flow_fraction']*100:.1f}%)")
            return "\n".join(lines)

        elif method == "multi":
            from .heritability import multiresolution_heritability
            r = multiresolution_heritability(conn, verbose=False)
            return (f"Multi-Resolution Heritability:\n"
                    f"  h²_variant  = {r['h2_variant']:.4f}\n"
                    f"  h²_gene     = {r['h2_gene']:.4f}\n"
                    f"  h²_pathway  = {r['h2_pathway']:.4f}\n"
                    f"  h²_intergenic        = {r['h2_intergenic']:.4f}\n"
                    f"  h²_non_pathway_genic = {r['h2_non_pathway_genic']:.4f}")

        elif method == "report":
            from .heritability import heritability_report
            r = heritability_report(conn, chromosomes=chroms, verbose=False)
            lines = ["Heritability Report:"]
            for name, h2 in r["summary"].items():
                lines.append(f"  h²_{name} = {h2:.4f}")
            return "\n".join(lines)

        else:
            return (f"Unknown method '{method}'. "
                    f"Use: spectral, conductance, flow, multi, report")


    # -----------------------------------------------------------------------
    # Tool 11: Multi-trait genetic correlation
    # -----------------------------------------------------------------------

    @mcp.tool()
    def genetic_correlation(trait1: str, trait2: str,
                            method: str = "spectral",
                            af_threshold: float = 0.01) -> str:
        """Estimate genetic correlation between two traits.

        Computes r_G, Cov_G, and bivariate heritability using graph-native
        spectral methods.  Can also decompose by biological scale
        (variant/gene/pathway) or compute spectral coherence.

        Args:
            trait1: First trait name (e.g. 'case_control')
            trait2: Second trait name (e.g. 'bmi')
            method: 'spectral' (r_G + Cov_G), 'coherence' (frequency-resolved),
                    'multi' (variant/gene/pathway decomposition)
            af_threshold: Allele frequency threshold for similarity graph

        Returns:
            Genetic correlation estimates with decomposition details.
        """
        conn = _get_conn()

        if method == "spectral":
            from .multivariate import spectral_genetic_covariance
            r = spectral_genetic_covariance(conn, trait1, trait2,
                                            af_threshold=af_threshold,
                                            verbose=False)
            if "error" in r:
                return f"Error: {r['error']}"
            return (f"Genetic Correlation: {trait1} × {trait2}\n"
                    f"  r_G    = {r['r_g']:.4f}\n"
                    f"  Cov_G  = {r['cov_g']:.6f}\n"
                    f"  h²_biv = {r['h2_bivariate']:.4f}\n"
                    f"  h²({trait1}) = {r['h2_trait1']:.4f}\n"
                    f"  h²({trait2}) = {r['h2_trait2']:.4f}")

        elif method == "coherence":
            from .multivariate import spectral_coherence
            r = spectral_coherence(conn, trait1, trait2,
                                   af_threshold=af_threshold, verbose=False)
            if "error" in r:
                return f"Error: {r['error']}"
            return (f"Spectral Coherence: {trait1} × {trait2}\n"
                    f"  Genetic coherence:     {r['mean_genetic_coherence']:.4f}\n"
                    f"  Environmental coherence: {r['mean_env_coherence']:.4f}\n"
                    f"  Genetic signed:        {r['mean_genetic_signed']:+.4f}\n"
                    f"  Environmental signed:  {r['mean_env_signed']:+.4f}")

        elif method == "multi":
            from .multivariate import multiresolution_correlation
            r = multiresolution_correlation(conn, trait1, trait2,
                                            af_threshold=af_threshold,
                                            verbose=False)
            if "error" in r:
                return f"Error: {r['error']}"
            return (f"Multi-Resolution r_G: {trait1} × {trait2}\n"
                    f"  r_G_variant  = {r['r_g_variant']:.4f}\n"
                    f"  r_G_gene     = {r['r_g_gene']:.4f}\n"
                    f"  r_G_pathway  = {r['r_g_pathway']:.4f}\n"
                    f"  LD-driven    = {r['r_g_ld_driven']:.4f}\n"
                    f"  Gene-specific = {r['r_g_gene_specific']:.4f}")

        return f"Unknown method '{method}'. Use: spectral, coherence, multi"


    # -----------------------------------------------------------------------
    # Tool 12: Gene Pleiotropy Strength Index
    # -----------------------------------------------------------------------

    @mcp.tool()
    def gene_pleiotropy_score(run_ids: str, p_threshold: float = 5e-8) -> str:
        """Compute Pleiotropy Strength Index (PSI) for every gene across traits.

        PSI = breadth × magnitude × mechanism_bonus, where:
        - breadth: fraction of traits where gene is significant
        - magnitude: mean -log10(p) across significant traits
        - mechanism_bonus: 1 + shared_variant_ratio (biological > mediated)

        High PSI = gene strongly and consistently affects many traits.

        Args:
            run_ids: Comma-separated GWAS run IDs (one per trait)
            p_threshold: Significance threshold per trait

        Returns:
            Top genes ranked by PSI with per-trait breakdown.
        """
        conn = _get_conn()
        runs = [r.strip() for r in run_ids.split(",")]
        from .multivariate import gene_pleiotropy_score as _psi
        result = _psi(conn, runs, p_threshold=p_threshold, verbose=False)

        lines = [f"Gene Pleiotropy Strength Index ({result['n_traits']} traits):",
                 f"  Genes scored: {result['n_genes_scored']}",
                 f"  Multi-trait (≥2): {result['n_multi_trait']}",
                 ""]
        for g in result["genes"][:15]:
            lines.append(
                f"  #{g['rank']} {g['gene']:<15} PSI={g['psi']:.3f} "
                f"breadth={g['breadth']:.2f} "
                f"mag={g['magnitude']:.1f} "
                f"({g['pleiotropy_type']})"
            )
        return "\n".join(lines)


    # -----------------------------------------------------------------------
    # Tool 13: Polygenic Risk Score
    # -----------------------------------------------------------------------

    @mcp.tool()
    def compute_prs(run_id: str, method: str = "report",
                    p_threshold: float = 5e-8) -> str:
        """Compute Polygenic Risk Score using graph-native methods.

        Graph-native PRS goes beyond weighted sums: LD-aware pruning via
        graph centrality, pathway-partitioned decomposition, and GNN-learned
        non-linear scores.

        Args:
            run_id: GWAS run ID for effect sizes
            method: 'classical' (baseline), 'graph_pruned' (LD-aware),
                    'pathway' (decomposed), 'report' (compare all)
            p_threshold: p-value threshold for variant inclusion

        Returns:
            PRS evaluation with AUROC, R², and method comparison.
        """
        conn = _get_conn()
        if method == "report":
            from .prs import prs_report
            r = prs_report(conn, run_id, p_threshold, verbose=False)
            lines = [f"PRS Report (run: {run_id}):"]
            for name, s in r.get("summary", {}).items():
                auroc = s.get("auroc", "N/A")
                if isinstance(auroc, float):
                    lines.append(f"  {name}: AUROC={auroc:.4f}, "
                                 f"n_variants={s.get('n_variants', '?')}")
            return "\n".join(lines)
        elif method == "classical":
            from .prs import classical_prs, prs_evaluation
            r = classical_prs(conn, run_id, p_threshold, verbose=False)
            if r.get("prs_scores"):
                ev = prs_evaluation(r["prs_scores"], r["n_cases"],
                                    r["n_controls"], verbose=False)
                return (f"Classical PRS: {r['n_variants']} variants\n"
                        f"  AUROC={ev['auroc']:.4f}, R²={ev['nagelkerke_r2']:.4f}")
            return f"No significant variants found"
        elif method == "pathway":
            from .prs import pathway_partitioned_prs
            r = pathway_partitioned_prs(conn, run_id, p_threshold, verbose=False)
            lines = [f"Pathway PRS: {r['n_pathways']} pathways"]
            for p in r.get("pathway_prs", [])[:5]:
                lines.append(f"  {p['pathway']}: d={p['cohens_d']:.3f}")
            return "\n".join(lines)
        return f"Unknown method '{method}'. Use: classical, graph_pruned, pathway, report"


    # -----------------------------------------------------------------------
    # Tool 14: Mendelian Randomization
    # -----------------------------------------------------------------------

    @mcp.tool()
    def mendelian_randomization(exposure_run_id: str, outcome_run_id: str,
                                p_threshold: float = 5e-8) -> str:
        """Graph-native Mendelian Randomization for causal inference.

        Uses graph centrality for instrument selection, detects pleiotropy
        from pathway structure, and decomposes causal effects by pathway.

        Args:
            exposure_run_id: GWAS run ID for exposure trait
            outcome_run_id: GWAS run ID for outcome trait
            p_threshold: instrument significance threshold

        Returns:
            MR estimates (IVW, Egger, weighted median) with sensitivity analyses.
        """
        conn = _get_conn()
        from .mr import mr_report
        r = mr_report(conn, exposure_run_id, outcome_run_id,
                      p_threshold=p_threshold, verbose=False)
        if "error" in r:
            return f"Error: {r['error']}"
        s = r["summary"]
        lines = [
            f"MR: {exposure_run_id} → {outcome_run_id} ({r['n_instruments']} instruments)",
            f"  IVW:    β={s['ivw_beta']:.4f} (p={s['ivw_p']:.2e})",
            f"  Egger:  β={s['egger_beta']:.4f} (p={s['egger_p']:.2e})",
            f"  W.Med:  β={s['wm_beta']:.4f} (p={s['wm_p']:.2e})",
            f"  Egger intercept p={s['egger_intercept_p']:.4f}",
            f"  Pleiotropic instruments: {s['n_pleiotropic']}",
        ]
        return "\n".join(lines)


    # -----------------------------------------------------------------------
    # Tool 15: Multi-Environment Trial analysis
    # -----------------------------------------------------------------------

    @mcp.tool()
    def analyze_met(trial_id: str, trait: str,
                    method: str = "report") -> str:
        """Analyze multi-environment trial data.

        Graph-native MET analysis: environment similarity, mega-environment
        detection, G×E decomposition, and design connectivity diagnostics.
        The experimental design IS the graph topology.

        Args:
            trial_id: Trial identifier
            trait: Trait name for analysis
            method: 'report' (full), 'similarity' (env correlations),
                    'mega-env' (community detection), 'gxe' (variance decomp),
                    'diagnostics' (design connectivity)

        Returns:
            Analysis results with design metrics and G×E decomposition.
        """
        conn = _get_conn()

        if method == "diagnostics":
            from .met import design_diagnostics
            r = design_diagnostics(conn, trial_id, verbose=False)
            if "error" in r:
                return f"Error: {r['error']}"
            return (f"Design Diagnostics ({trial_id}):\n"
                    f"  {r['n_genotypes']} genotypes × {r['n_environments']} environments\n"
                    f"  Balance: {r['balance_ratio']:.1%}\n"
                    f"  Connected components: {r['n_connected_components']}\n"
                    f"  Algebraic connectivity: {r['algebraic_connectivity']:.4f}")

        elif method == "similarity":
            from .met import environment_similarity
            r = environment_similarity(conn, trial_id, trait, verbose=False)
            pairs = sorted(r["pairwise"], key=lambda p: p["r_g"], reverse=True)
            lines = [f"Environment Similarity ({r['n_environments']} environments):"]
            for p in pairs[:10]:
                lines.append(f"  {p['env1']} × {p['env2']}: r_G={p['r_g']:.4f}")
            return "\n".join(lines)

        elif method == "mega-env":
            from .met import mega_environments
            r = mega_environments(conn, trial_id, trait, verbose=False)
            lines = [f"Mega-Environments: {r['n_clusters']} clusters"]
            for ci, info in r["clusters"].items():
                lines.append(f"  Cluster {ci}: {info['n']} envs, "
                             f"mean r_G={info['mean_within_r_g']:.4f}")
            return "\n".join(lines)

        elif method == "gxe":
            from .met import gxe_decomposition
            r = gxe_decomposition(conn, trial_id, trait, verbose=False)
            if "error" in r:
                return f"Error: {r['error']}"
            return (f"G×E Decomposition ({trial_id} / {trait}):\n"
                    f"  V_G={r['v_g']:.4f}, V_E={r['v_e']:.4f}, "
                    f"V_GxE={r['v_gxe']:.4f}\n"
                    f"  h²_broad={r['h2_broad']:.4f}\n"
                    f"  Type B r_G={r['type_b_correlation']:.4f}")

        elif method == "report":
            from .met import met_report
            r = met_report(conn, trial_id, trait, verbose=False)
            d = r.get("diagnostics", {})
            g = r.get("gxe", {})
            m = r.get("mega_environments", {})
            lines = [
                f"MET Report: {trial_id} / {trait}",
                f"  Design: {d.get('n_genotypes', '?')} × {d.get('n_environments', '?')} "
                f"({d.get('balance_ratio', 0):.1%} balanced)",
            ]
            if "h2_broad" in g:
                lines.append(f"  h²={g['h2_broad']:.4f}, G×E={g['gxe_ratio']:.4f}")
            lines.append(f"  Mega-environments: {m.get('n_clusters', '?')}")
            return "\n".join(lines)

        return f"Unknown method '{method}'. Use: report, similarity, mega-env, gxe, diagnostics"


    # -----------------------------------------------------------------------
    # Tool 14: Graph-based MET imputation
    # -----------------------------------------------------------------------

    @mcp.tool()
    def impute_missing_met(trial_id: str, trait: str,
                           genotype: str, environment: str,
                           k_neighbors: int = 10) -> str:
        """Predict missing genotype-environment phenotype via graph neighbors.

        Finds genetically similar genotypes tested in the target environment
        and computes similarity-weighted prediction.  No matrix inversion —
        works with any pattern of missingness.

        Args:
            trial_id: Trial identifier
            trait: Trait name
            genotype: Genotype ID to impute for
            environment: Environment ID to impute in
            k_neighbors: Number of nearest neighbors

        Returns:
            Imputed value with confidence and neighbor details.
        """
        conn = _get_conn()
        from .met import graph_impute_missing
        r = graph_impute_missing(conn, trial_id, trait,
                                 genotype, environment,
                                 k_neighbors=k_neighbors, verbose=False)
        if "error" in r:
            return f"Error: {r['error']}"

        lines = [
            f"Imputation: {genotype} in {environment}",
            f"  Predicted {trait} = {r['imputed_value']:.4f}",
            f"  Confidence: {r['confidence']:.4f}",
            f"  Neighbors: {r['n_neighbors']}",
        ]
        for n in r["neighbors"][:5]:
            lines.append(f"    {n['genotype']}: sim={n['similarity']:.4f}, "
                         f"val={n['observed_value']:.4f}")
        return "\n".join(lines)


# ---------------------------------------------------------------------------
# Server lifecycle
# ---------------------------------------------------------------------------

def create_mcp_server(uri: str, user: str, password: str,
                      database: str | None = None) -> "FastMCP":
    """Initialize the MCP server with a database connection.

    Args:
        uri: Neo4j bolt URI
        user: Neo4j username
        password: Neo4j password
        database: Neo4j database name (None = default)

    Returns:
        Configured FastMCP server instance.
    """
    if not HAS_MCP:
        raise ImportError(
            "MCP package required. pip install 'mcp>=1.0'"
        )

    global _conn
    _conn = GraphGWASConnection(uri, user, password, database)
    _conn.__enter__()
    return mcp


def run_server(transport: str = "stdio"):
    """Entry point for CLI: graphgwas mcp."""
    if not HAS_MCP:
        raise ImportError(
            "MCP package required. pip install 'mcp>=1.0'"
        )

    from .config import NEO4J_URI, NEO4J_USER, NEO4J_PASSWORD, NEO4J_DATABASE
    create_mcp_server(NEO4J_URI, NEO4J_USER, NEO4J_PASSWORD, NEO4J_DATABASE)
    mcp.run(transport=transport)
