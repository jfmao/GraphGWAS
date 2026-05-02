"""GraphGWAS CLI — Click-based command-line interface."""


import click

from .config import NEO4J_URI, NEO4J_USER, NEO4J_PASSWORD, NEO4J_DATABASE
from .db import GraphGWASConnection


class _AutoDetectConnection(GraphGWASConnection):
    """GraphGWASConnection that auto-detects N_SAMPLES on enter."""
    def __enter__(self):
        result = super().__enter__()
        from .config import detect_n_samples
        detect_n_samples(self)
        return result


def _connect(ctx):
    """Get connection from Click context. Auto-detects N_SAMPLES."""
    return _AutoDetectConnection(
        ctx.obj["uri"], ctx.obj["user"], ctx.obj["password"], ctx.obj["database"]
    )


@click.group()
@click.option("--neo4j-uri", default=NEO4J_URI, help="Neo4j URI")
@click.option("--neo4j-user", default=NEO4J_USER, help="Neo4j user")
@click.option("--neo4j-password", default=NEO4J_PASSWORD, help="Neo4j password")
@click.option("--database", default=NEO4J_DATABASE, help="Neo4j database")
@click.pass_context
def cli(ctx, neo4j_uri, neo4j_user, neo4j_password, database):
    """GraphGWAS — Graph-native GWAS platform."""
    ctx.ensure_object(dict)
    ctx.obj["uri"] = neo4j_uri
    ctx.obj["user"] = neo4j_user
    ctx.obj["password"] = neo4j_password
    ctx.obj["database"] = database


# ---------------------------------------------------------------------------
# Status
# ---------------------------------------------------------------------------

@cli.command()
@click.pass_context
def status(ctx):
    """Show database summary and GWAS status."""
    from .schema import audit_schema, verify_graphmana_compat

    with _connect(ctx) as conn:
        report = audit_schema(conn)
        issues = verify_graphmana_compat(conn)

    click.echo("=== GraphGWAS Database Status ===")
    for key, val in report.items():
        click.echo(f"  {key}: {val:,}" if isinstance(val, int) else f"  {key}: {val}")
    if issues:
        click.echo("\nCompatibility issues:")
        for issue in issues:
            click.echo(f"  WARNING: {issue}")
    else:
        click.echo("\nGraphMana compatibility: OK")


# ---------------------------------------------------------------------------
# Phenotype commands
# ---------------------------------------------------------------------------

@cli.group()
def phenotype():
    """Phenotype import and management."""
    pass


@phenotype.command("load")
@click.option("--phenotype-file", required=True, type=click.Path(exists=True))
@click.option("--sample-id-col", required=True)
@click.option("--traits", default=None, help="Comma-separated trait columns (default: all)")
@click.pass_context
def phenotype_load(ctx, phenotype_file, sample_id_col, traits):
    """Load phenotype CSV into Sample.phenotypes MAP."""
    from .phenotype import load_phenotypes

    trait_cols = traits.split(",") if traits else None
    with _connect(ctx) as conn:
        result = load_phenotypes(conn, phenotype_file, sample_id_col, trait_cols)

    click.echo(f"Loaded: {result['n_matched']} matched, {result['n_unmatched']} unmatched")
    click.echo(f"Traits: {', '.join(result['traits_loaded'])}")
    for trait, summary in result["trait_summary"].items():
        if summary["type"] == "categorical":
            click.echo(f"  {trait}: {summary['type']} — {summary['levels']}")
        else:
            click.echo(f"  {trait}: {summary['type']} — mean={summary['mean']:.2f}")


@phenotype.command("activate")
@click.option("--trait", required=True)
@click.option("--case-value", default=None, help="Value indicating case (binary trait)")
@click.option("--control-value", default=None, help="Value indicating control (binary trait)")
@click.pass_context
def phenotype_activate(ctx, trait, case_value, control_value):
    """Activate a trait for GWAS (set is_case/is_control flags)."""
    from .phenotype import activate_phenotype

    with _connect(ctx) as conn:
        result = activate_phenotype(conn, trait, case_value, control_value)

    click.echo(f"Activated: {result['trait']} ({result['trait_type']})")
    if result["trait_type"] == "binary":
        click.echo(f"  Cases: {result['n_cases']}, Controls: {result['n_controls']}, "
                    f"Missing: {result['n_missing']}")
    else:
        click.echo(f"  With value: {result['n_with_value']}, Missing: {result['n_missing']}")


# ---------------------------------------------------------------------------
# Schema commands
# ---------------------------------------------------------------------------

@cli.command("ensure-indexes")
@click.pass_context
def ensure_indexes(ctx):
    """Create all GWAS indexes."""
    from .schema import ensure_indexes as _ensure

    with _connect(ctx) as conn:
        results = _ensure(conn)
    for r in results:
        click.echo(r)


# ---------------------------------------------------------------------------
# QC commands
# ---------------------------------------------------------------------------

@cli.group()
def qc():
    """Quality control procedures."""
    pass


@qc.command("run")
@click.option("--chromosomes", default=None, help="Comma-separated chr list (default: all)")
@click.pass_context
def qc_run(ctx, chromosomes):
    """Run full variant and sample QC scan."""
    from .qc import run_variant_qc

    chrom_list = chromosomes.split(",") if chromosomes else None
    with _connect(ctx) as conn:
        result = run_variant_qc(conn, chromosomes=chrom_list)

    click.echo("\nQC Complete:")
    click.echo(f"  Variants scanned: {result['n_variants_scanned']:,}")
    click.echo(f"  HWE failures: {result['n_hwe_fail']:,}")
    click.echo(f"  MAC < 5 in cases: {result['n_mac_fail']:,}")


# ---------------------------------------------------------------------------
# Association commands
# ---------------------------------------------------------------------------

@cli.group()
def assoc():
    """Association analysis."""
    pass


@assoc.command("scan")
@click.option("--chr", "chromosome", required=True, help="Chromosome (e.g. chromosome4) or 'all'")
@click.option("--start", type=int, default=None)
@click.option("--end", type=int, default=None)
@click.option("--method", default="auto",
              type=click.Choice(["auto", "chi2", "fisher", "logistic", "firth", "linear"]))
@click.option("--covariates", default=None, help="Comma-separated covariate property names")
@click.option("--output", "-o", default=None, type=click.Path(), help="Output TSV file path")
@click.option("--store/--no-store", default=False, help="Also store results in Neo4j graph")
@click.option("--parallel/--no-parallel", default=True, help="Parallelize across chromosomes")
@click.option("--workers", default=None, type=int, help="Number of parallel workers")
@click.pass_context
def assoc_scan(ctx, chromosome, start, end, method, covariates, output, store, parallel, workers):
    """Run single-locus association scan on a region.

    Results are written to a TSV file (like PLINK/SAIGE). Use --store to
    also integrate results into the Neo4j graph as AssociationResult nodes.
    """
    import csv
    import time
    import numpy as np
    from .assoc import single_locus_scan, generate_run_id

    run_id = generate_run_id()
    cov_names = covariates.split(",") if covariates else None
    t0 = time.time()

    with _connect(ctx) as conn:
        # Get active phenotype
        rec = conn.execute_read(
            "MATCH (s:Sample) WHERE s.gwas_phenotype IS NOT NULL "
            "RETURN s.gwas_phenotype AS p LIMIT 1"
        ).single()
        phenotype_key = rec["p"] if rec else "unknown"

        # Default output path
        if output is None:
            output = f"gwas_{phenotype_key}_{method}_{run_id}.tsv"

        # Genome-wide mode
        if chromosome == "all":
            chroms_result = conn.execute_read(
                "MATCH (v:Variant) RETURN DISTINCT v.chr AS c ORDER BY c"
            )
            chromosomes = [r["c"] for r in chroms_result]
        else:
            chromosomes = [chromosome]

        # Run scan — parallel or sequential
        if parallel and len(chromosomes) > 1:
            from .parallel import parallel_genome_scan
            all_results = parallel_genome_scan(
                conn, chromosomes, method=method,
                covariate_names=cov_names, max_workers=workers,
            )
        else:
            all_results = []
            for chrom in chromosomes:
                results = single_locus_scan(
                    conn, chrom, start if len(chromosomes) == 1 else None,
                    end if len(chromosomes) == 1 else None,
                    method=method, covariate_names=cov_names,
                )
                all_results.extend(results)

        elapsed = time.time() - t0

        # Write TSV output (default behavior, like conventional GWAS tools)
        tsv_cols = ["variantId", "chr", "pos", "ref", "alt", "beta", "se",
                     "p_value", "p_value_log10", "r_squared", "af", "mac",
                     "method", "n_cases", "n_controls"]
        with open(output, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=tsv_cols, delimiter="\t",
                                    extrasaction="ignore")
            writer.writeheader()
            for r in all_results:
                writer.writerow(r)

        # Summary — compute lambda from common variants only (AF 0.05-0.95)
        # to avoid deflation from low-power rare variant tests
        from scipy import stats as sp_stats
        pvals_all = np.array([r["p_value"] for r in all_results if r["p_value"] > 0])
        pvals_common = np.array([r["p_value"] for r in all_results
                                 if r["p_value"] > 0 and r.get("af", 0) >= 0.05
                                 and r.get("af", 1) <= 0.95])
        n_gw = int((pvals_all < 5e-8).sum()) if len(pvals_all) > 0 else 0
        n_sug = int((pvals_all < 1e-5).sum()) if len(pvals_all) > 0 else 0
        if len(pvals_common) > 100:
            chi2 = sp_stats.chi2.isf(pvals_common, 1)
            lam = float(np.median(chi2) / 0.4549)
        else:
            lam = float(np.median(-2 * np.log(pvals_all)) / 1.386) if len(pvals_all) > 0 else 0

        click.echo(f"Done: {len(all_results):,} variants, {elapsed:.1f}s")
        click.echo(f"  Genome-wide significant (p<5e-8): {n_gw:,}")
        click.echo(f"  Suggestive (p<1e-5): {n_sug:,}")
        click.echo(f"  Lambda (common variants, AF>5%): {lam:.4f}")
        click.echo(f"  Output: {output}")

        # Optionally store in Neo4j
        if store:
            from .results import create_gwas_study, store_results, link_study
            region = f"{chromosome}:{start}-{end}" if start and chromosome != "all" else chromosome
            n_cases = all_results[0].get("n_cases", 0) if all_results else 0
            n_controls = all_results[0].get("n_controls", 0) if all_results else 0
            create_gwas_study(conn, run_id, phenotype_key, method,
                              n_cases, n_controls, len(all_results), region)
            n_stored = store_results(conn, all_results, run_id, phenotype_key)
            link_study(conn, run_id)
            click.echo(f"  Stored {n_stored:,} results in Neo4j (run_id: {run_id})")


@assoc.command("mpat")
@click.option("--chr", "chromosome", default=None, help="Chromosome to scan")
@click.option("--genes", default=None, help="Comma-separated gene symbols")
@click.option("--method", default="directed", type=click.Choice(["directed", "undirected"]))
@click.option("--min-variants", default=2, type=int)
@click.pass_context
def assoc_mpat(ctx, chromosome, genes, method, min_variants):
    """Run MPAT gene-level association test."""
    from .mpat import mpat_scan

    gene_list = genes.split(",") if genes else None

    with _connect(ctx) as conn:
        results = mpat_scan(conn, chr=chromosome, genes=gene_list,
                            method=method, min_variants=min_variants)

    click.echo(f"\n{'Gene':<15} {'Z_gene':<10} {'P-value':<12} {'N_variants':<12} {'Method'}")
    click.echo("-" * 65)
    for r in results[:50]:
        click.echo(f"{r['gene']:<15} {r['z_gene']:<10.3f} {r['p_value']:<12.2e} "
                    f"{r['n_variants']:<12} {r['method']}")


@assoc.command("top-hits")
@click.option("--run-id", required=True)
@click.option("--p-threshold", default=5e-8, type=float)
@click.pass_context
def assoc_top_hits(ctx, run_id, p_threshold):
    """Show top hits from a previous GWAS run."""
    from .results import query_top_hits

    with _connect(ctx) as conn:
        hits = query_top_hits(conn, run_id, p_threshold)

    if not hits:
        click.echo("No hits found.")
        return

    click.echo(f"{'Variant':<50} {'Chr':<6} {'Pos':<12} {'P-value':<12} {'Beta':<10} {'Genes'}")
    click.echo("-" * 110)
    for h in hits:
        genes = ",".join(h.get("genes", [])[:3])
        click.echo(f"{h['variant']:<50} {h['chr']:<6} {h['pos']:<12} "
                    f"{h['p_value']:.2e}  {h.get('beta', 'N/A'):<10} {genes}")


# ---------------------------------------------------------------------------
# Gene counts
# ---------------------------------------------------------------------------

@cli.command("gene-counts")
@click.pass_context
def gene_counts(ctx):
    """Pre-compute per-gene variant counts for burden test filtering."""
    with _connect(ctx) as conn:
        conn.execute_write(
            """
            MATCH (v:Variant)-[r:HAS_CONSEQUENCE]->(g:Gene)
            WITH g, count(v) AS total,
                 sum(CASE WHEN v.af_total < 0.01 AND r.impact IN ['HIGH','MODERATE'] THEN 1 ELSE 0 END) AS rare_func,
                 sum(CASE WHEN r.impact = 'HIGH' THEN 1 ELSE 0 END) AS high_impact,
                 sum(CASE WHEN r.impact = 'MODERATE' THEN 1 ELSE 0 END) AS mod_impact
            SET g.n_variants_total = total,
                g.n_rare_functional = rare_func,
                g.n_high_impact = high_impact,
                g.n_moderate_impact = mod_impact
            """
        )
    click.echo("Gene variant counts updated.")


# ---------------------------------------------------------------------------
# Plot commands
# ---------------------------------------------------------------------------

@cli.group()
def plot():
    """Visualization and plotting."""
    pass


@plot.command("manhattan")
@click.option("--run-id", required=True)
@click.option("--output", default="manhattan.png", help="Output file path")
@click.option("--title", default="", help="Plot title")
@click.pass_context
def plot_manhattan(ctx, run_id, output, title):
    """Generate Manhattan plot from stored results."""
    from .plots import manhattan_plot

    with _connect(ctx) as conn:
        manhattan_plot(conn, run_id, output, title=title)


@plot.command("qq")
@click.option("--run-id", required=True)
@click.option("--output", default="qq.png", help="Output file path")
@click.option("--title", default="", help="Plot title")
@click.pass_context
def plot_qq(ctx, run_id, output, title):
    """Generate QQ plot with genomic inflation factor."""
    from .plots import qq_plot

    with _connect(ctx) as conn:
        qq_plot(conn, run_id, output, title=title)


# ---------------------------------------------------------------------------
# Results export
# ---------------------------------------------------------------------------

@cli.group()
def results():
    """Result management and export."""
    pass


@results.command("export")
@click.option("--run-id", required=True)
@click.option("--output", required=True, help="Output TSV file path")
@click.pass_context
def results_export(ctx, run_id, output):
    """Export association results to TSV (GNExT/PheWeb compatible)."""
    from .plots import results_to_tsv

    with _connect(ctx) as conn:
        results_to_tsv(conn, run_id, output)


@results.command("delete")
@click.option("--run-id", required=True)
@click.confirmation_option(prompt="Delete all results for this run?")
@click.pass_context
def results_delete(ctx, run_id):
    """Delete all results for a GWAS run."""
    from .results import delete_run

    with _connect(ctx) as conn:
        n = delete_run(conn, run_id)
    click.echo(f"Deleted {n} results for run {run_id}")


# ---------------------------------------------------------------------------
# PCA import
# ---------------------------------------------------------------------------

@cli.command("import-pca")
@click.option("--eigenvec-file", required=True, type=click.Path(exists=True))
@click.option("--n-components", default=10, type=int)
@click.pass_context
def import_pca_cmd(ctx, eigenvec_file, n_components):
    """Import PCA coordinates from PLINK2 eigenvec file."""
    from .genotype import import_pca

    with _connect(ctx) as conn:
        n = import_pca(conn, eigenvec_file, n_components)
    click.echo(f"Updated {n} samples with {n_components} PC coordinates")


@cli.group()
def popstruct():
    """Population structure: graph-spectral PCs and correction."""
    pass


@popstruct.command("spectral-pcs")
@click.option("--n-components", default=20, type=int, help="Number of spectral PCs")
@click.option("--af-threshold", default=0.05, type=float, help="AF threshold for rare variants")
@click.pass_context
def popstruct_spectral_pcs(ctx, n_components, af_threshold):
    """Compute graph-spectral PCs from rare-variant similarity.

    Stores spectral_pc_1..spectral_pc_K on Sample nodes. Use as
    --covariates spectral_pc_1,...,spectral_pc_K in GWAS scans.
    This is a graph-native replacement for PCA that captures
    rare-variant-based population structure.
    """
    from .popstruct import compute_spectral_pcs

    with _connect(ctx) as conn:
        result = compute_spectral_pcs(conn, n_components=n_components,
                                       af_threshold=af_threshold)
    click.echo(f"\nSpectral PCs computed: {result['n_components']} components, "
               f"{result['n_variants_used']:,} variants, {result['n_samples']} samples")


@popstruct.command("grammar")
@click.option("--trait", required=True, help="Phenotype trait name (must be loaded)")
@click.option("--n-pcs", default=5, type=int, help="Number of PCs as fixed effects in GRAMMAR")
@click.pass_context
def popstruct_grammar(ctx, trait, n_pcs):
    """GRAMMAR mixed-model correction: residualize phenotype against GRM.

    Computes GRM from common variants in the graph, estimates h² via REML,
    and stores GRAMMAR residuals as gwas_value on Sample nodes. Then run
    'graphgwas assoc scan --method linear' on the residuals.

    Also reports GRAMMAR+ calibration factor for post-hoc p-value correction.
    """
    import time
    import numpy as np
    from .genotype import get_all_indices, get_phenotype_values, variant_iterator, build_dosage
    from .popstruct import compute_grm, grammar_residualize, grammar_plus_calibrate

    t0 = time.time()

    with _connect(ctx) as conn:
        # Activate trait
        from .phenotype import activate_phenotype
        info = activate_phenotype(conn, trait)
        click.echo(f"Trait: {trait} ({info.get('n_with_value', 0)} samples)")

        all_idx = get_all_indices(conn)
        pheno = get_phenotype_values(conn, all_idx)
        valid = ~np.isnan(pheno)
        all_idx_valid = all_idx[valid]
        pheno_valid = pheno[valid]
        n = len(pheno_valid)
        click.echo(f"Samples with phenotype: {n}")

        # Build genotype matrix (common variants)
        click.echo("Building genotype matrix from graph (common variants, AF 5-95%)...")
        chromosomes = [r["c"] for r in conn.execute_read(
            "MATCH (v:Variant) RETURN DISTINCT v.chr AS c ORDER BY c")]

        geno_rows = []
        for chrom in chromosomes:
            for v in variant_iterator(conn, chrom):
                af = v.get("af_total", 0)
                if af < 0.05 or af > 0.95:
                    continue
                gtp = v["gt_packed"]
                if gtp is None:
                    continue
                from .config import N_SAMPLES
                dosage = build_dosage(gtp, all_idx_valid, N_SAMPLES)
                m = np.nanmean(dosage)
                dosage = np.where(np.isnan(dosage), m, dosage)
                geno_rows.append(dosage)

        G = np.array(geno_rows, dtype=np.float64)
        click.echo(f"  {G.shape[0]:,} common variants × {G.shape[1]} samples")

        # GRM
        click.echo("Computing GRM...")
        K = compute_grm(G, verbose=True)

        # Get stored PCs as covariates (if available)
        pcs = None
        if n_pcs > 0:
            try:
                from .genotype import get_covariate_matrix
                pc_names = [f"grm_pc_{i+1}" for i in range(n_pcs)]
                pcs = get_covariate_matrix(conn, pc_names, all_idx_valid)
                if np.isnan(pcs).any():
                    click.echo(f"  Warning: NaN in PCs, using {n_pcs} PCs from GRM eigendecomposition")
                    from .popstruct import compute_pca
                    pca = compute_pca(K, n_components=n_pcs, verbose=False)
                    pcs = pca["pcs"][:, :n_pcs]
            except Exception:
                from .popstruct import compute_pca
                pca = compute_pca(K, n_components=n_pcs, verbose=False)
                pcs = pca["pcs"][:, :n_pcs]
                click.echo(f"  Using {n_pcs} PCs from GRM eigendecomposition")

        # GRAMMAR
        click.echo("GRAMMAR residualization...")
        grammar = grammar_residualize(pheno_valid, K, covariates=pcs, verbose=True)
        residuals = grammar["residuals"]
        click.echo(f"  h² = {grammar['h2']:.4f}")

        # GRAMMAR+ calibration
        click.echo("GRAMMAR+ calibration...")
        factor = grammar_plus_calibrate(G, residuals, verbose=True)

        # Store residuals on Sample nodes
        click.echo("Storing GRAMMAR residuals as gwas_value on Sample nodes...")
        batch = [{"idx": int(all_idx_valid[i]), "val": float(residuals[i])} for i in range(n)]
        for i in range(0, len(batch), 500):
            sub = batch[i:i+500]
            conn.execute_write("""
                UNWIND $batch AS row
                MATCH (s:Sample {packed_index: row.idx})
                SET s.gwas_value = row.val,
                    s.gwas_phenotype = $trait_name
            """, {"batch": sub, "trait_name": f"{trait}_grammar"})

        elapsed = time.time() - t0
        click.echo(f"\nGRAMMAR complete ({elapsed:.0f}s)")
        click.echo(f"  h² = {grammar['h2']:.4f}")
        click.echo(f"  GRAMMAR+ calibration factor = {factor:.4f}")
        click.echo(f"  Residuals stored as gwas_value (phenotype: {trait}_grammar)")
        click.echo("\nNext: run GWAS on residuals:")
        click.echo(f"  graphgwas assoc scan --chr all --method linear --parallel -o gwas_{trait}_grammar.tsv")
        click.echo(f"  Then apply GRAMMAR+ factor to p-values: p_cal = chi2.sf(chi2.isf(p, 1) * {factor:.4f}, 1)")


# ---------------------------------------------------------------------------
# Phase 3: Graph-native methods
# ---------------------------------------------------------------------------

@cli.group()
def epistasis():
    """Epistasis detection via co-occurrence networks."""
    pass


@epistasis.command("scan")
@click.option("--chr", "chromosome", required=True)
@click.option("--start", type=int, default=None)
@click.option("--end", type=int, default=None)
@click.option("--window", default=1000000, type=int, help="Window size in bp")
@click.option("--min-cocarriers", default=3, type=int)
@click.option("--min-enrichment", default=1.5, type=float)
@click.option("--permutations", default=100, type=int)
@click.option("--parallel/--no-parallel", default=True, help="Parallelize across windows")
@click.option("--workers", default=None, type=int, help="Number of parallel workers")
@click.pass_context
def epistasis_scan(ctx, chromosome, start, end, window, min_cocarriers,
                   min_enrichment, permutations, parallel, workers):
    """Scan for epistatic modules via co-occurrence networks."""
    from .epistasis import epistasis_scan as _scan, build_cooccurrence_graph, detect_epistatic_modules

    with _connect(ctx) as conn:
        if start is not None and end is not None:
            # Single window
            G = build_cooccurrence_graph(conn, chromosome, start, end,
                                         min_cocarriers=min_cocarriers,
                                         min_enrichment=min_enrichment)
            modules = detect_epistatic_modules(G)
            click.echo(f"Co-occurrence graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
            click.echo(f"Detected {len(modules)} epistatic modules")
            for i, m in enumerate(modules[:10]):
                click.echo(f"  Module {i+1}: {len(m)} variants — {', '.join(m[:5])}")
        elif parallel:
            # Parallel sliding window scan
            from .parallel import parallel_epistasis_scan as _par_scan
            # Build window list
            result = conn.execute_read(
                "MATCH (v:Variant) WHERE v.chr = $chr "
                "RETURN min(v.pos) AS mn, max(v.pos) AS mx",
                {"chr": chromosome},
            ).single()
            if not result:
                click.echo("No variants found.")
                return
            step = window // 2
            windows = []
            pos = result["mn"]
            while pos < result["mx"]:
                windows.append((pos, min(pos + window, result["mx"])))
                pos += step
            results = _par_scan(
                conn, chromosome, windows,
                min_cocarriers=min_cocarriers,
                min_enrichment=min_enrichment,
                n_permutations=permutations,
                max_workers=workers,
            )
            for r in results[:20]:
                click.echo(f"  Module: {r['n_variants']} variants, "
                           f"score={r['module_score']:.3f}, p={r['p_value']:.4f}")
        else:
            # Sequential sliding window scan
            results = _scan(conn, chromosome, window_size=window,
                            min_cocarriers=min_cocarriers,
                            min_enrichment=min_enrichment,
                            n_permutations=permutations)
            for r in results[:20]:
                click.echo(f"  Module: {r['n_variants']} variants, "
                           f"score={r['module_score']:.3f}, p={r['p_value']:.4f}")


@cli.group()
def flow():
    """Max-flow / min-cut disease architecture analysis."""
    pass


@flow.command("architecture")
@click.option("--chr", "chromosome", default=None)
@click.option("--af-threshold", default=0.05, type=float)
@click.option("--permutations", default=100, type=int)
@click.option("--p-threshold", default=0.05, type=float)
@click.pass_context
def flow_architecture(ctx, chromosome, af_threshold, permutations, p_threshold):
    """Extract disease architecture via max-flow analysis."""
    from .flow import extract_disease_architecture

    with _connect(ctx) as conn:
        arch = extract_disease_architecture(
            conn, chr=chromosome, af_threshold=af_threshold,
            p_threshold=p_threshold, n_permutations=permutations
        )

    click.echo("\nDisease Architecture:")
    click.echo(f"  Significant pathways: {arch['n_pathways']}")
    click.echo(f"  Disease variants: {arch['n_variants']}")
    click.echo(f"  Disease genes: {arch['n_genes']}")
    for pw in arch.get("pathway_details", [])[:10]:
        click.echo(f"    {pw['pathway']}: flow={pw['flow_value']:.2f} p={pw['p_value']:.4f}")


@cli.command("finemap")
@click.option("--chr", "chromosome", required=True)
@click.option("--start", required=True, type=int)
@click.option("--end", required=True, type=int)
@click.option("--r2-threshold", default=0.3, type=float)
@click.option("--method", default="betweenness",
              type=click.Choice(["betweenness", "degree", "pagerank"]))
@click.option("--run-id", default=None, help="GWAS run ID for p-value annotation")
@click.pass_context
def finemap_cmd(ctx, chromosome, start, end, r2_threshold, method, run_id):
    """Fine-map a locus using LD graph centrality."""
    from .finemapping import finemap_locus

    with _connect(ctx) as conn:
        results = finemap_locus(conn, chromosome, start, end,
                                r2_threshold=r2_threshold, method=method,
                                run_id=run_id)


@cli.group()
def spectral():
    """Graph spectral phenotype analysis."""
    pass


@spectral.command("correction")
@click.option("--n-components", default=10, type=int)
@click.option("--af-threshold", default=0.01, type=float)
@click.option("--chr", "chromosome", default=None, help="Single chr for fast testing")
@click.pass_context
def spectral_correction(ctx, n_components, af_threshold, chromosome):
    """Compute spectral phenotype correction from rare-variant similarity."""
    from .spectral import spectral_gwas_correction

    chroms = [chromosome] if chromosome else None
    with _connect(ctx) as conn:
        result = spectral_gwas_correction(
            conn, n_components=n_components,
            af_threshold=af_threshold, chromosomes=chroms
        )

    click.echo("\nSpectral correction computed:")
    click.echo(f"  Components: {result['n_components']}")
    click.echo(f"  Eigenvalue range: [{result['eigenvalues'][0]:.4f}, {result['eigenvalues'][-1]:.4f}]")


# ---------------------------------------------------------------------------
# Phase 4: GNN commands
# ---------------------------------------------------------------------------

@cli.group()
def gnn():
    """GNN association engine."""
    pass


@gnn.command("run")
@click.option("--chr", "chromosome", required=True)
@click.option("--start", type=int, default=None)
@click.option("--end", type=int, default=None)
@click.option("--hidden-channels", default=64, type=int)
@click.option("--epochs", default=100, type=int)
@click.option("--max-variants", default=5000, type=int)
@click.option("--device", default="auto", type=click.Choice(["auto", "cpu", "cuda"]),
              help="Compute device for GNN training")
@click.pass_context
def gnn_run(ctx, chromosome, start, end, hidden_channels, epochs, max_variants, device):
    """Run full GNN pipeline: export → train → explain → import embeddings."""
    from .gnn import run_gnn_pipeline

    with _connect(ctx) as conn:
        result = run_gnn_pipeline(
            conn, chromosome, start, end,
            hidden_channels=hidden_channels, epochs=epochs,
            max_variants=max_variants, device=device,
        )
        if result:
            click.echo("\nGNN Results:")
            click.echo(f"  Device: {result.get('device', 'cpu')}")
            click.echo(f"  Best AUROC: {result['best_auroc']:.4f}")
            click.echo(f"  Embeddings imported: {result['n_embeddings_imported']}")
            if result.get("top_variants"):
                click.echo(f"  Top variant: {result['top_variants'][0]['variantId']}")


@gnn.command("export")
@click.option("--chr", "chromosome", required=True)
@click.option("--start", type=int, default=None)
@click.option("--end", type=int, default=None)
@click.option("--output", default="graph_data.pt", help="Output .pt file")
@click.option("--max-variants", default=5000, type=int)
@click.pass_context
def gnn_export(ctx, chromosome, start, end, output, max_variants):
    """Export subgraph to PyTorch Geometric format."""
    from .gnn import export_to_pyg
    import torch

    with _connect(ctx) as conn:
        data = export_to_pyg(conn, chromosome, start, end, max_variants=max_variants)
        torch.save(data, output)
        click.echo(f"Saved HeteroData to {output}")


# ---------------------------------------------------------------------------
# Phase 5: Interpretation commands (LLM-free)
# ---------------------------------------------------------------------------

@cli.command("interpret")
@click.option("--run-id", required=True)
@click.pass_context
def interpret_cmd(ctx, run_id):
    """Interpret GWAS results (rule-based; no LLM required)."""
    from .interpret import interpret_results

    with _connect(ctx) as conn:
        interpret_results(conn, run_id)


# ---------------------------------------------------------------------------
# Heritability commands
# ---------------------------------------------------------------------------

@cli.group()
def heritability():
    """Graph-native heritability estimation."""
    pass


@heritability.command("spectral")
@click.option("--af-threshold", default=0.01, type=float)
@click.option("--n-components", default=20, type=int)
@click.option("--chr", "chromosome", default=None, help="Single chr for fast testing")
@click.pass_context
def h2_spectral(ctx, af_threshold, n_components, chromosome):
    """Estimate heritability from Laplacian eigenspectrum."""
    from .heritability import spectral_heritability
    chroms = [chromosome] if chromosome else None
    with _connect(ctx) as conn:
        result = spectral_heritability(conn, af_threshold=af_threshold,
                                       n_components=n_components,
                                       chromosomes=chroms)
    click.echo(f"\nh²_spectral = {result['h2']:.4f}")


@heritability.command("conductance")
@click.option("--af-threshold", default=0.01, type=float)
@click.option("--chr", "chromosome", default=None)
@click.option("--permutations", default=200, type=int)
@click.pass_context
def h2_conductance(ctx, af_threshold, chromosome, permutations):
    """Estimate heritability from graph conductance."""
    from .heritability import conductance_heritability
    chroms = [chromosome] if chromosome else None
    with _connect(ctx) as conn:
        result = conductance_heritability(conn, af_threshold=af_threshold,
                                          chromosomes=chroms,
                                          n_permutations=permutations)
    click.echo(f"\nh²_conductance = {result['h2_conductance']:.4f}  "
               f"(p = {result['p_value']:.4f})")


@heritability.command("flow")
@click.option("--chr", "chromosome", default=None)
@click.option("--af-threshold", default=0.05, type=float)
@click.option("--permutations", default=100, type=int)
@click.pass_context
def h2_flow(ctx, chromosome, af_threshold, permutations):
    """Decompose heritability by pathway using max-flow."""
    from .heritability import flow_heritability
    with _connect(ctx) as conn:
        result = flow_heritability(conn, chr=chromosome,
                                   af_threshold=af_threshold,
                                   n_permutations=permutations)
    click.echo(f"\nh²_flow = {result['h2_flow_total']:.4f}  "
               f"(p = {result['p_value']:.4f})")
    for c in result.get("pathway_contributions", [])[:10]:
        click.echo(f"  {c['pathway']}: h²={c['h2_contribution']:.4f} "
                   f"({c['flow_fraction']*100:.1f}%)")


@heritability.command("gnn")
@click.option("--chr", "chromosome", required=True)
@click.option("--start", type=int, default=None)
@click.option("--end", type=int, default=None)
@click.option("--epochs", default=100, type=int)
@click.option("--device", default="auto", type=click.Choice(["auto", "cpu", "cuda"]))
@click.pass_context
def h2_gnn(ctx, chromosome, start, end, epochs, device):
    """Estimate heritability from GNN prediction accuracy."""
    from .heritability import gnn_heritability
    with _connect(ctx) as conn:
        result = gnn_heritability(conn, chr=chromosome, start=start, end=end,
                                  epochs=epochs, device=device)
    click.echo(f"\nh²_GNN = {result.get('h2_gnn', 0):.4f}  "
               f"(AUROC = {result.get('auroc', 0):.4f})")
    for i, (label, h2, inc) in enumerate(zip(
            result.get("layer_labels", []),
            result.get("h2_per_layer", []),
            result.get("h2_layer_increments", []))):
        click.echo(f"  Layer {i+1} ({label}): h²={h2:.4f}, "
                   f"increment={inc:.4f}")


@heritability.command("multi")
@click.option("--n-components", default=20, type=int)
@click.pass_context
def h2_multi(ctx, n_components):
    """Multi-resolution heritability (variant / gene / pathway scales)."""
    from .heritability import multiresolution_heritability
    with _connect(ctx) as conn:
        result = multiresolution_heritability(conn, n_components=n_components)
    click.echo(f"\nh²_variant  = {result['h2_variant']:.4f}")
    click.echo(f"h²_gene     = {result['h2_gene']:.4f}")
    click.echo(f"h²_pathway  = {result['h2_pathway']:.4f}")
    click.echo(f"h²_intergenic        = {result['h2_intergenic']:.4f}")
    click.echo(f"h²_non_pathway_genic = {result['h2_non_pathway_genic']:.4f}")


@heritability.command("report")
@click.option("--chr", "chromosome", default=None)
@click.option("--include-multi", is_flag=True, help="Include multi-resolution (slow)")
@click.option("--include-gnn", is_flag=True, help="Include GNN (requires PyTorch)")
@click.pass_context
def h2_report(ctx, chromosome, include_multi, include_gnn):
    """Run all heritability estimators and produce unified report."""
    from .heritability import heritability_report
    chroms = [chromosome] if chromosome else None
    with _connect(ctx) as conn:
        result = heritability_report(conn, chromosomes=chroms,
                                     include_multiresolution=include_multi,
                                     include_gnn=include_gnn,
                                     gnn_chr=chromosome or "chr19")
    click.echo("\nSummary:")
    for name, h2 in result["summary"].items():
        click.echo(f"  h²_{name:<20s} = {h2:.4f}")


# ---------------------------------------------------------------------------
# Multivariate commands
# ---------------------------------------------------------------------------

@cli.group()
def multivariate():
    """Multi-trait genetic correlation and G-matrix."""
    pass


@multivariate.command("correlation")
@click.option("--trait1", required=True, help="First trait name")
@click.option("--trait2", required=True, help="Second trait name")
@click.option("--af-threshold", default=0.01, type=float)
@click.option("--n-components", default=20, type=int)
@click.pass_context
def mv_correlation(ctx, trait1, trait2, af_threshold, n_components):
    """Estimate genetic correlation between two traits."""
    from .multivariate import spectral_genetic_covariance
    with _connect(ctx) as conn:
        result = spectral_genetic_covariance(
            conn, trait1, trait2, af_threshold=af_threshold,
            n_components=n_components,
        )
    click.echo(f"\nr_G({trait1}, {trait2}) = {result.get('r_g', 0):.4f}")
    click.echo(f"Cov_G = {result.get('cov_g', 0):.6f}")
    click.echo(f"h²_biv = {result.get('h2_bivariate', 0):.4f}")


@multivariate.command("g-matrix")
@click.option("--traits", required=True, help="Comma-separated trait names")
@click.option("--af-threshold", default=0.01, type=float)
@click.option("--n-components", default=20, type=int)
@click.pass_context
def mv_gmatrix(ctx, traits, af_threshold, n_components):
    """Compute the T×T genetic variance-covariance matrix."""
    from .multivariate import g_matrix
    trait_list = [t.strip() for t in traits.split(",")]
    with _connect(ctx) as conn:
        result = g_matrix(conn, trait_list, af_threshold=af_threshold,
                          n_components=n_components)
    if "error" in result:
        click.echo(f"Error: {result['error']}")
        return
    click.echo("\nHeritabilities:")
    for t, h2 in zip(result["trait_names"], result["heritabilities"]):
        click.echo(f"  h²({t}) = {h2:.4f}")


@multivariate.command("multi-resolution")
@click.option("--trait1", required=True)
@click.option("--trait2", required=True)
@click.option("--af-threshold", default=0.01, type=float)
@click.pass_context
def mv_multiresolution(ctx, trait1, trait2, af_threshold):
    """Genetic correlation at variant/gene/pathway scales."""
    from .multivariate import multiresolution_correlation
    with _connect(ctx) as conn:
        result = multiresolution_correlation(conn, trait1, trait2,
                                             af_threshold=af_threshold)
    click.echo(f"\nr_G_variant  = {result['r_g_variant']:.4f}")
    click.echo(f"r_G_gene     = {result['r_g_gene']:.4f}")
    click.echo(f"r_G_pathway  = {result['r_g_pathway']:.4f}")


@multivariate.command("coherence")
@click.option("--trait1", required=True)
@click.option("--trait2", required=True)
@click.option("--af-threshold", default=0.01, type=float)
@click.option("--n-components", default=20, type=int)
@click.pass_context
def mv_coherence(ctx, trait1, trait2, af_threshold, n_components):
    """Frequency-resolved genetic vs environmental covariance."""
    from .multivariate import spectral_coherence
    with _connect(ctx) as conn:
        result = spectral_coherence(conn, trait1, trait2,
                                    af_threshold=af_threshold,
                                    n_components=n_components)
    click.echo(f"\nGenetic coherence:     {result['mean_genetic_coherence']:.4f}")
    click.echo(f"Environmental coherence: {result['mean_env_coherence']:.4f}")


@multivariate.command("pleiotropy")
@click.option("--run1", required=True, help="GWAS run ID for trait 1")
@click.option("--run2", required=True, help="GWAS run ID for trait 2")
@click.option("--p-threshold", default=5e-8, type=float)
@click.pass_context
def mv_pleiotropy(ctx, run1, run2, p_threshold):
    """Find pleiotropic genes shared between two GWAS runs."""
    from .multivariate import pleiotropic_genes
    with _connect(ctx) as conn:
        result = pleiotropic_genes(conn, run1, run2, p_threshold=p_threshold)
    n = len(result["pleiotropic_genes"])
    click.echo(f"\nPleiotropic genes: {n}")
    click.echo(f"  Biological: {result['n_biological_pleiotropy']}")
    click.echo(f"  Mediated:   {result['n_mediated_pleiotropy']}")
    for g in result["pleiotropic_genes"][:10]:
        click.echo(f"  {g['gene']} ({g['pleiotropy_type']}): "
                   f"{g['n_shared_variants']} shared variants")


@multivariate.command("psi")
@click.option("--run-ids", required=True, help="Comma-separated GWAS run IDs")
@click.option("--labels", default=None, help="Comma-separated trait labels")
@click.option("--p-threshold", default=5e-8, type=float)
@click.pass_context
def mv_psi(ctx, run_ids, labels, p_threshold):
    """Compute Gene Pleiotropy Strength Index (PSI) across traits."""
    from .multivariate import gene_pleiotropy_score
    runs = [r.strip() for r in run_ids.split(",")]
    lbls = [l.strip() for l in labels.split(",")] if labels else None
    with _connect(ctx) as conn:
        result = gene_pleiotropy_score(conn, runs, trait_labels=lbls,
                                       p_threshold=p_threshold)
    click.echo(f"\nGenes scored: {result['n_genes_scored']}")
    click.echo(f"Multi-trait (≥2): {result['n_multi_trait']}")
    for g in result["genes"][:20]:
        click.echo(f"  #{g['rank']} {g['gene']:<15} PSI={g['psi']:.3f} "
                   f"({g['n_traits_significant']}/{result['n_traits']} traits, "
                   f"{g['pleiotropy_type']})")


@multivariate.command("report")
@click.option("--trait1", required=True)
@click.option("--trait2", required=True)
@click.option("--af-threshold", default=0.01, type=float)
@click.pass_context
def mv_report(ctx, trait1, trait2, af_threshold):
    """Full multivariate report: covariance + coherence."""
    from .multivariate import multivariate_report
    with _connect(ctx) as conn:
        result = multivariate_report(conn, trait1, trait2,
                                     af_threshold=af_threshold)
    s = result["summary"]
    click.echo("\nSummary:")
    click.echo(f"  r_G           = {s['r_g']:.4f}")
    click.echo(f"  h²_bivariate  = {s['h2_bivariate']:.4f}")
    click.echo(f"  Genetic coh.  = {s['genetic_coherence']:.4f}")


# ---------------------------------------------------------------------------
# Multi-Environment Trial commands
# ---------------------------------------------------------------------------

@cli.group()
def met():
    """Multi-Environment Trial analysis."""
    pass


@met.command("load")
@click.option("--csv", "csv_path", required=True, type=click.Path(exists=True))
@click.option("--genotype-col", required=True)
@click.option("--env-col", required=True)
@click.option("--traits", required=True, help="Comma-separated trait columns")
@click.option("--env-covariates", default=None, help="Comma-separated env covariate columns")
@click.option("--rep-col", default=None)
@click.option("--block-col", default=None)
@click.option("--trial-id", default="default")
@click.pass_context
def met_load(ctx, csv_path, genotype_col, env_col, traits, env_covariates,
             rep_col, block_col, trial_id):
    """Import MET data from long-format CSV."""
    from .met import load_met_data
    trait_list = [t.strip() for t in traits.split(",")]
    env_cov_list = [c.strip() for c in env_covariates.split(",")] if env_covariates else None
    with _connect(ctx) as conn:
        result = load_met_data(conn, csv_path, genotype_col, env_col,
                               trait_list, env_covariate_cols=env_cov_list,
                               rep_col=rep_col, block_col=block_col,
                               trial_id=trial_id)
    click.echo(f"\nLoaded: {result['n_genotypes']} genotypes, "
               f"{result['n_environments']} environments, "
               f"{result['n_observations']} observations")


@met.command("status")
@click.option("--trial-id", default="default")
@click.pass_context
def met_status(ctx, trial_id):
    """Show MET trial summary."""
    with _connect(ctx) as conn:
        result = conn.execute_read(
            """
            MATCH (t:Trial {id: $tid})
            RETURN t.n_genotypes AS ng, t.n_environments AS ne, t.n_observations AS no
            """,
            {"tid": trial_id},
        ).single()
    if result:
        click.echo(f"Trial: {trial_id}")
        click.echo(f"  Genotypes: {result['ng']}, Environments: {result['ne']}, "
                   f"Observations: {result['no']}")
    else:
        click.echo(f"Trial '{trial_id}' not found.")


@met.command("env-similarity")
@click.option("--trial-id", default="default")
@click.option("--trait", required=True)
@click.pass_context
def met_env_sim(ctx, trial_id, trait):
    """Compute genetic correlation between environments."""
    from .met import environment_similarity
    with _connect(ctx) as conn:
        result = environment_similarity(conn, trial_id, trait)
    click.echo(f"\n{result['n_environments']} environments, "
               f"{len(result['pairwise'])} pairs computed")


@met.command("mega-env")
@click.option("--trial-id", default="default")
@click.option("--trait", required=True)
@click.pass_context
def met_mega_env(ctx, trial_id, trait):
    """Detect mega-environments via community detection."""
    from .met import mega_environments
    with _connect(ctx) as conn:
        result = mega_environments(conn, trial_id, trait)
    click.echo(f"\n{result['n_clusters']} mega-environments detected")


@met.command("gxe")
@click.option("--trial-id", default="default")
@click.option("--trait", required=True)
@click.pass_context
def met_gxe(ctx, trial_id, trait):
    """G×E variance decomposition."""
    from .met import gxe_decomposition
    with _connect(ctx) as conn:
        result = gxe_decomposition(conn, trial_id, trait)
    if "error" not in result:
        click.echo(f"\nh²_broad = {result['h2_broad']:.4f}, "
                   f"G×E ratio = {result['gxe_ratio']:.4f}, "
                   f"Type B r_G = {result['type_b_correlation']:.4f}")


@met.command("diagnostics")
@click.option("--trial-id", default="default")
@click.pass_context
def met_diagnostics(ctx, trial_id):
    """Graph connectivity analysis of experimental design."""
    from .met import design_diagnostics
    with _connect(ctx) as conn:
        result = design_diagnostics(conn, trial_id)
    if "error" not in result:
        click.echo(f"\nBalance: {result['balance_ratio']:.1%}, "
                   f"Components: {result['n_connected_components']}, "
                   f"Fiedler: {result['algebraic_connectivity']:.4f}")


@met.command("reaction-norm")
@click.option("--trial-id", default="default")
@click.option("--trait", required=True)
@click.option("--env-covariate", default=None, help="Environment property for regression")
@click.option("--top-n", default=100, type=int)
@click.pass_context
def met_reaction_norm(ctx, trial_id, trait, env_covariate, top_n):
    """Variant-level reaction norms (effect × environment)."""
    from .met import variant_reaction_norms
    with _connect(ctx) as conn:
        result = variant_reaction_norms(conn, trial_id, trait,
                                        env_covariate=env_covariate,
                                        top_n_variants=top_n)
    click.echo(f"\n{result['n_variants']} variants with multi-env effects")


@met.command("impute")
@click.option("--trial-id", default="default")
@click.option("--trait", required=True)
@click.option("--genotype", required=True)
@click.option("--environment", required=True)
@click.option("--k-neighbors", default=10, type=int)
@click.pass_context
def met_impute(ctx, trial_id, trait, genotype, environment, k_neighbors):
    """Impute missing G×E cell via graph neighbors."""
    from .met import graph_impute_missing
    with _connect(ctx) as conn:
        result = graph_impute_missing(conn, trial_id, trait,
                                      genotype, environment,
                                      k_neighbors=k_neighbors)
    if "error" in result:
        click.echo(f"Error: {result['error']}")
    else:
        click.echo(f"\nImputed {trait}({genotype}, {environment}) = "
                   f"{result['imputed_value']:.4f} "
                   f"(confidence: {result['confidence']:.4f})")


@met.command("report")
@click.option("--trial-id", default="default")
@click.option("--trait", required=True)
@click.pass_context
def met_report_cmd(ctx, trial_id, trait):
    """Comprehensive MET analysis report."""
    from .met import met_report
    with _connect(ctx) as conn:
        met_report(conn, trial_id, trait)


# ---------------------------------------------------------------------------
# PRS commands
# ---------------------------------------------------------------------------

@cli.group()
def prs():
    """Polygenic Risk Score computation."""
    pass


@prs.command("classical")
@click.option("--run-id", required=True)
@click.option("--p-threshold", default=5e-8, type=float)
@click.pass_context
def prs_classical(ctx, run_id, p_threshold):
    """Standard weighted-sum PRS."""
    from .prs import classical_prs, prs_evaluation
    with _connect(ctx) as conn:
        result = classical_prs(conn, run_id, p_threshold)
    if result.get("prs_scores"):
        ev = prs_evaluation(result["prs_scores"], result["n_cases"],
                            result["n_controls"], verbose=True)


@prs.command("graph-pruned")
@click.option("--run-id", required=True)
@click.option("--p-threshold", default=5e-8, type=float)
@click.option("--r2-threshold", default=0.2, type=float)
@click.option("--method", default="centrality",
              type=click.Choice(["centrality", "pagerank", "independent_set"]))
@click.pass_context
def prs_pruned(ctx, run_id, p_threshold, r2_threshold, method):
    """LD-aware PRS using graph centrality pruning."""
    from .prs import graph_pruned_prs, prs_evaluation
    with _connect(ctx) as conn:
        result = graph_pruned_prs(conn, run_id, p_threshold,
                                  r2_threshold=r2_threshold,
                                  pruning_method=method)
    if result.get("prs_scores"):
        prs_evaluation(result["prs_scores"], result["n_cases"],
                       result["n_controls"], verbose=True)


@prs.command("pathway")
@click.option("--run-id", required=True)
@click.option("--p-threshold", default=1e-5, type=float)
@click.pass_context
def prs_pathway(ctx, run_id, p_threshold):
    """Pathway-partitioned PRS."""
    from .prs import pathway_partitioned_prs
    with _connect(ctx) as conn:
        pathway_partitioned_prs(conn, run_id, p_threshold)


@prs.command("report")
@click.option("--run-id", required=True)
@click.option("--p-threshold", default=5e-8, type=float)
@click.option("--include-gnn", is_flag=True)
@click.pass_context
def prs_report_cmd(ctx, run_id, p_threshold, include_gnn):
    """Compare all PRS methods."""
    from .prs import prs_report
    with _connect(ctx) as conn:
        prs_report(conn, run_id, p_threshold, include_gnn=include_gnn)


# ---------------------------------------------------------------------------
# MR commands
# ---------------------------------------------------------------------------

@cli.group()
def mr():
    """Mendelian Randomization."""
    pass


@mr.command("report")
@click.option("--exposure-run", required=True, help="GWAS run ID for exposure")
@click.option("--outcome-run", required=True, help="GWAS run ID for outcome")
@click.option("--p-threshold", default=5e-8, type=float)
@click.pass_context
def mr_report_cmd(ctx, exposure_run, outcome_run, p_threshold):
    """Full MR analysis: IVW + Egger + weighted median + pleiotropy."""
    from .mr import mr_report
    with _connect(ctx) as conn:
        mr_report(conn, exposure_run, outcome_run, p_threshold=p_threshold)


# ---------------------------------------------------------------------------
# Phase 6: Server commands
# ---------------------------------------------------------------------------

@cli.command("serve")
@click.option("--host", default="0.0.0.0", help="Bind host")
@click.option("--port", default=8000, type=int, help="Bind port")
@click.pass_context
def serve(ctx, host, port):
    """Start GraphGWAS REST API server (FastAPI + Uvicorn)."""
    import os
    os.environ["GRAPHGWAS_NEO4J_URI"] = ctx.obj["uri"]
    os.environ["GRAPHGWAS_NEO4J_USER"] = ctx.obj["user"]
    os.environ["GRAPHGWAS_NEO4J_PASSWORD"] = ctx.obj["password"]
    if ctx.obj["database"]:
        os.environ["GRAPHGWAS_NEO4J_DATABASE"] = ctx.obj["database"]

    try:
        import uvicorn
    except ImportError:
        click.echo("FastAPI/Uvicorn required. pip install 'graphgwas[api]'")
        return

    click.echo(f"Starting GraphGWAS API on {host}:{port}")
    click.echo(f"  Docs: http://{host}:{port}/docs")
    uvicorn.run("graphgwas.api:app", host=host, port=port)


@cli.command("mcp")
@click.option("--transport", default="stdio",
              type=click.Choice(["stdio", "sse"]), help="MCP transport")
@click.pass_context
def mcp_cmd(ctx, transport):
    """Start GraphGWAS MCP server for AI assistant integration."""
    try:
        from .mcp_server import create_mcp_server
    except ImportError:
        click.echo("MCP package required. pip install 'graphgwas[mcp]'")
        return

    click.echo(f"Starting GraphGWAS MCP server (transport: {transport})")
    server = create_mcp_server(
        ctx.obj["uri"], ctx.obj["user"],
        ctx.obj["password"], ctx.obj["database"],
    )
    server.run(transport=transport)


if __name__ == "__main__":
    cli()
