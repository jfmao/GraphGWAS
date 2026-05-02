"""GraphGWAS REST API — FastAPI server for programmatic access.

Provides HTTP endpoints for all GraphGWAS operations. Long-running operations
(genome-wide scan, epistasis, flow analysis, GNN training) run as background
jobs with polling status endpoints.

Usage:
    graphgwas serve --host 0.0.0.0 --port 8000
    # or directly:
    uvicorn graphgwas.api:app --host 0.0.0.0 --port 8000
"""

from __future__ import annotations

import os
import threading
import uuid
from contextlib import asynccontextmanager
from datetime import datetime, timezone
from enum import Enum
from typing import Any

from fastapi import FastAPI, HTTPException, Request
from pydantic import BaseModel, Field


# ---------------------------------------------------------------------------
# Pydantic models
# ---------------------------------------------------------------------------

class JobStatusEnum(str, Enum):
    pending = "pending"
    running = "running"
    completed = "completed"
    failed = "failed"


class JobResponse(BaseModel):
    job_id: str
    status: JobStatusEnum
    created_at: str
    result: Any | None = None
    error: str | None = None
    progress: str | None = None


class ScanRequest(BaseModel):
    chromosome: str = Field("all", description="Chromosome (e.g. 'chr19') or 'all'")
    start: int | None = None
    end: int | None = None
    method: str = Field("auto", description="auto, chi2, fisher, logistic, firth, linear")
    covariates: list[str] | None = None
    parallel: bool = True
    workers: int | None = None
    store: bool = True


class EpistasisRequest(BaseModel):
    chromosome: str
    start: int | None = None
    end: int | None = None
    window_size: int = 1_000_000
    min_cocarriers: int = 3
    min_enrichment: float = 1.5
    permutations: int = 100
    parallel: bool = True
    workers: int | None = None


class FlowRequest(BaseModel):
    chromosome: str | None = None
    af_threshold: float = 0.05
    permutations: int = 100
    p_threshold: float = 0.05


class FinemapRequest(BaseModel):
    chromosome: str
    start: int
    end: int
    r2_threshold: float = 0.3
    method: str = "betweenness"
    run_id: str | None = None


class MpatRequest(BaseModel):
    chromosome: str | None = None
    genes: list[str] | None = None
    method: str = "directed"
    min_variants: int = 2


class GnnRequest(BaseModel):
    chromosome: str
    start: int | None = None
    end: int | None = None
    hidden_channels: int = 64
    epochs: int = 100
    max_variants: int = 5000
    device: str = "auto"


class PhenotypeLoadRequest(BaseModel):
    phenotype_file: str
    sample_id_col: str
    traits: list[str] | None = None


class PhenotypeActivateRequest(BaseModel):
    trait: str
    case_value: str | None = None
    control_value: str | None = None


class CypherRequest(BaseModel):
    query: str


# ---------------------------------------------------------------------------
# Job store (in-memory, single-server)
# ---------------------------------------------------------------------------

_jobs: dict[str, dict] = {}
_jobs_lock = threading.Lock()


def _create_job() -> str:
    job_id = str(uuid.uuid4())
    with _jobs_lock:
        _jobs[job_id] = {
            "status": JobStatusEnum.pending,
            "created_at": datetime.now(timezone.utc).isoformat(),
            "result": None,
            "error": None,
            "progress": None,
        }
    return job_id


def _update_job(job_id: str, **kwargs):
    with _jobs_lock:
        if job_id in _jobs:
            _jobs[job_id].update(kwargs)


def _get_job(job_id: str) -> dict | None:
    with _jobs_lock:
        return _jobs.get(job_id, {}).copy() if job_id in _jobs else None


# ---------------------------------------------------------------------------
# Background job runner
# ---------------------------------------------------------------------------

def _run_job_sync(job_id: str, func, *args, **kwargs):
    """Run a function synchronously in a thread, updating job state."""
    _update_job(job_id, status=JobStatusEnum.running)
    try:
        result = func(*args, **kwargs)
        # Ensure result is JSON-serializable
        if isinstance(result, list):
            result = {"items": result, "count": len(result)}
        elif not isinstance(result, dict):
            result = {"result": str(result)}
        _update_job(job_id, status=JobStatusEnum.completed, result=result)
    except Exception as e:
        _update_job(job_id, status=JobStatusEnum.failed, error=str(e))


def _submit_job(func, *args, **kwargs) -> str:
    """Create a job and run it in a background thread."""
    job_id = _create_job()
    t = threading.Thread(target=_run_job_sync, args=(job_id, func, *args),
                         kwargs=kwargs, daemon=True)
    t.start()
    return job_id


# ---------------------------------------------------------------------------
# App lifecycle
# ---------------------------------------------------------------------------

@asynccontextmanager
async def lifespan(app: FastAPI):
    """Manage Neo4j connection lifecycle."""
    from .db import GraphGWASConnection

    uri = os.environ.get("GRAPHGWAS_NEO4J_URI", "bolt://localhost:7688")
    user = os.environ.get("GRAPHGWAS_NEO4J_USER", "neo4j")
    password = os.environ.get("GRAPHGWAS_NEO4J_PASSWORD", "graphgwas")
    database = os.environ.get("GRAPHGWAS_NEO4J_DATABASE") or None

    conn = GraphGWASConnection(uri, user, password, database)
    conn.__enter__()
    app.state.conn = conn
    app.state.conn_params = {
        "uri": uri, "user": user, "password": password, "database": database,
    }
    yield
    conn.__exit__(None, None, None)


app = FastAPI(
    title="GraphGWAS API",
    description="Graph-native GWAS platform REST API",
    version="0.1.0",
    lifespan=lifespan,
)


def _conn(request: Request):
    return request.app.state.conn


# ---------------------------------------------------------------------------
# Status endpoint
# ---------------------------------------------------------------------------

@app.get("/api/v1/status")
async def status(request: Request):
    """Database summary and GWAS status."""
    from .schema import audit_schema
    return audit_schema(_conn(request))


# ---------------------------------------------------------------------------
# Phenotype endpoints
# ---------------------------------------------------------------------------

@app.post("/api/v1/phenotype/load")
async def phenotype_load(request: Request, body: PhenotypeLoadRequest):
    """Load phenotype CSV into Sample nodes."""
    from .phenotype import load_phenotypes
    return load_phenotypes(_conn(request), body.phenotype_file,
                           body.sample_id_col, body.traits)


@app.post("/api/v1/phenotype/activate")
async def phenotype_activate(request: Request, body: PhenotypeActivateRequest):
    """Activate a trait for GWAS (set is_case/is_control flags)."""
    from .phenotype import activate_phenotype
    return activate_phenotype(_conn(request), body.trait,
                              body.case_value, body.control_value)


# ---------------------------------------------------------------------------
# Association endpoints
# ---------------------------------------------------------------------------

@app.post("/api/v1/assoc/scan")
async def assoc_scan(request: Request, body: ScanRequest):
    """Submit association scan (async job for genome-wide, sync for single region)."""
    conn = _conn(request)

    def _run():
        from .assoc import single_locus_scan, generate_run_id
        from .results import create_gwas_study, store_results

        run_id = generate_run_id()

        # Determine chromosomes
        if body.chromosome == "all":
            chroms_result = conn.execute_read(
                "MATCH (v:Variant) RETURN DISTINCT v.chr AS c ORDER BY c"
            )
            chromosomes = [r["c"] for r in chroms_result]
        else:
            chromosomes = [body.chromosome]

        # Parallel or sequential
        if body.parallel and len(chromosomes) > 1:
            from .parallel import parallel_genome_scan
            all_results = parallel_genome_scan(
                conn, chromosomes, method=body.method,
                covariate_names=body.covariates,
                max_workers=body.workers, verbose=False,
            )
        else:
            all_results = []
            for chrom in chromosomes:
                results = single_locus_scan(
                    conn, chrom,
                    body.start if len(chromosomes) == 1 else None,
                    body.end if len(chromosomes) == 1 else None,
                    method=body.method,
                    covariate_names=body.covariates,
                    verbose=False,
                )
                all_results.extend(results)

        # Store results
        if body.store and all_results:
            rec = conn.execute_read(
                "MATCH (s:Sample) WHERE s.gwas_phenotype IS NOT NULL "
                "RETURN s.gwas_phenotype AS p LIMIT 1"
            ).single()
            phenotype_key = rec["p"] if rec else "unknown"
            n_cases = all_results[0].get("n_cases", 0)
            n_controls = all_results[0].get("n_controls", 0)

            create_gwas_study(conn, run_id, phenotype_key, body.method,
                              n_cases, n_controls, len(all_results),
                              body.chromosome)
            store_results(conn, all_results, run_id, phenotype_key)

        sig = sum(1 for r in all_results if r["p_value"] < 5e-8)
        return {
            "run_id": run_id,
            "n_variants": len(all_results),
            "n_significant": sig,
            "top_hits": sorted(all_results, key=lambda r: r["p_value"])[:20],
        }

    job_id = _submit_job(_run)
    return {"job_id": job_id}


@app.get("/api/v1/assoc/results/{run_id}")
async def assoc_results(request: Request, run_id: str,
                        p_threshold: float = 5e-8, limit: int = 100):
    """Query top hits from a stored GWAS run."""
    from .results import query_top_hits
    hits = query_top_hits(_conn(request), run_id, p_threshold, limit)
    return {"run_id": run_id, "hits": hits, "count": len(hits)}


# ---------------------------------------------------------------------------
# Epistasis endpoints
# ---------------------------------------------------------------------------

@app.post("/api/v1/epistasis/scan")
async def epistasis_scan(request: Request, body: EpistasisRequest):
    """Submit epistasis scan (async job)."""
    conn = _conn(request)

    def _run():
        if body.start is not None and body.end is not None:
            from .epistasis import build_cooccurrence_graph, detect_epistatic_modules
            G = build_cooccurrence_graph(
                conn, body.chromosome, body.start, body.end,
                min_cocarriers=body.min_cocarriers,
                min_enrichment=body.min_enrichment, verbose=False,
            )
            modules = detect_epistatic_modules(G)
            return {
                "nodes": G.number_of_nodes(),
                "edges": G.number_of_edges(),
                "modules": [{"variants": m, "size": len(m)} for m in modules[:20]],
            }
        else:
            from .epistasis import epistasis_scan as _scan
            results = _scan(conn, body.chromosome,
                            window_size=body.window_size,
                            min_cocarriers=body.min_cocarriers,
                            min_enrichment=body.min_enrichment,
                            n_permutations=body.permutations,
                            verbose=False)
            return {"modules": results}

    job_id = _submit_job(_run)
    return {"job_id": job_id}


# ---------------------------------------------------------------------------
# Flow endpoints
# ---------------------------------------------------------------------------

@app.post("/api/v1/flow/architecture")
async def flow_architecture(request: Request, body: FlowRequest):
    """Submit flow architecture analysis (async job)."""
    conn = _conn(request)

    def _run():
        from .flow import extract_disease_architecture
        return extract_disease_architecture(
            conn, chr=body.chromosome, af_threshold=body.af_threshold,
            p_threshold=body.p_threshold,
            n_permutations=body.permutations, verbose=False,
        )

    job_id = _submit_job(_run)
    return {"job_id": job_id}


# ---------------------------------------------------------------------------
# Fine-mapping endpoint
# ---------------------------------------------------------------------------

@app.post("/api/v1/finemap")
async def finemap(request: Request, body: FinemapRequest):
    """Fine-map a locus using LD graph centrality (synchronous)."""
    from .finemapping import finemap_locus
    results = finemap_locus(_conn(request), body.chromosome, body.start, body.end,
                            r2_threshold=body.r2_threshold, method=body.method,
                            run_id=body.run_id, verbose=False)
    return {"results": results}


# ---------------------------------------------------------------------------
# MPAT endpoint
# ---------------------------------------------------------------------------

@app.post("/api/v1/mpat")
async def mpat(request: Request, body: MpatRequest):
    """Run MPAT gene-level test (synchronous for small gene sets)."""
    from .mpat import mpat_scan
    results = mpat_scan(_conn(request), chr=body.chromosome,
                        genes=body.genes, method=body.method,
                        min_variants=body.min_variants, verbose=False)
    return {"results": results, "count": len(results)}


# ---------------------------------------------------------------------------
# GNN endpoints
# ---------------------------------------------------------------------------

@app.post("/api/v1/gnn/run")
async def gnn_run(request: Request, body: GnnRequest):
    """Submit GNN pipeline (async job)."""
    conn = _conn(request)

    def _run():
        from .gnn import run_gnn_pipeline
        result = run_gnn_pipeline(
            conn, body.chromosome, body.start, body.end,
            hidden_channels=body.hidden_channels, epochs=body.epochs,
            max_variants=body.max_variants, device=body.device,
            verbose=False,
        )
        if not result:
            return {"error": "No variants in region"}
        return {
            "best_auroc": result["best_auroc"],
            "n_embeddings_imported": result["n_embeddings_imported"],
            "top_variants": result.get("top_variants", [])[:20],
        }

    job_id = _submit_job(_run)
    return {"job_id": job_id}


# ---------------------------------------------------------------------------
# Job status endpoint
# ---------------------------------------------------------------------------

@app.get("/api/v1/jobs/{job_id}")
async def job_status(job_id: str):
    """Poll job status."""
    job = _get_job(job_id)
    if job is None:
        raise HTTPException(status_code=404, detail=f"Job {job_id} not found")
    return {"job_id": job_id, **job}


@app.get("/api/v1/jobs")
async def list_jobs(limit: int = 50):
    """List recent jobs."""
    with _jobs_lock:
        items = [
            {"job_id": jid, **jdata}
            for jid, jdata in sorted(_jobs.items(),
                                     key=lambda x: x[1]["created_at"],
                                     reverse=True)[:limit]
        ]
    return {"jobs": items, "count": len(items)}


# ---------------------------------------------------------------------------
# Cypher query endpoint
# ---------------------------------------------------------------------------

@app.post("/api/v1/cypher/query")
async def cypher_query(request: Request, body: CypherRequest):
    """Execute a read-only Cypher query."""
    query_lower = body.query.lower().strip()
    write_keywords = ["create", "merge", "delete", "set ", "remove", "drop"]
    if any(kw in query_lower for kw in write_keywords):
        raise HTTPException(status_code=400, detail="Only read-only queries allowed")

    try:
        result = _conn(request).execute_read(body.query)
        data = result.data()
        return {"data": data[:500], "count": len(data)}
    except Exception as e:
        raise HTTPException(status_code=400, detail=f"Query error: {e}")


# ---------------------------------------------------------------------------
# Interpret endpoint
# ---------------------------------------------------------------------------

@app.post("/api/v1/interpret/{run_id}")
async def interpret(request: Request, run_id: str):
    """Structured GWAS interpretation (no LLM required)."""
    from .interpret import interpret_results
    return interpret_results(_conn(request), run_id, verbose=False)


# ---------------------------------------------------------------------------
# Heritability endpoints
# ---------------------------------------------------------------------------

class HeritabilitySpectralRequest(BaseModel):
    af_threshold: float = 0.01
    n_components: int = 20
    chromosome: str | None = None


class HeritabilityConductanceRequest(BaseModel):
    af_threshold: float = 0.01
    chromosome: str | None = None
    permutations: int = 200


class HeritabilityFlowRequest(BaseModel):
    chromosome: str | None = None
    af_threshold: float = 0.05
    permutations: int = 100


class HeritabilityGnnRequest(BaseModel):
    chromosome: str = "chr19"
    start: int | None = None
    end: int | None = None
    epochs: int = 100
    device: str = "auto"


class HeritabilityReportRequest(BaseModel):
    chromosome: str | None = None
    include_multiresolution: bool = False
    include_gnn: bool = False


@app.post("/api/v1/heritability/spectral")
async def h2_spectral(request: Request, body: HeritabilitySpectralRequest):
    """Spectral heritability from Laplacian eigenspectrum (synchronous)."""
    from .heritability import spectral_heritability
    chroms = [body.chromosome] if body.chromosome else None
    return spectral_heritability(_conn(request), af_threshold=body.af_threshold,
                                 n_components=body.n_components,
                                 chromosomes=chroms, verbose=False)


@app.post("/api/v1/heritability/conductance")
async def h2_conductance(request: Request, body: HeritabilityConductanceRequest):
    """Conductance heritability (synchronous)."""
    from .heritability import conductance_heritability
    chroms = [body.chromosome] if body.chromosome else None
    return conductance_heritability(_conn(request), af_threshold=body.af_threshold,
                                    chromosomes=chroms,
                                    n_permutations=body.permutations,
                                    verbose=False)


@app.post("/api/v1/heritability/flow")
async def h2_flow(request: Request, body: HeritabilityFlowRequest):
    """Flow heritability decomposition (async job — permutations)."""
    conn = _conn(request)

    def _run():
        from .heritability import flow_heritability
        return flow_heritability(conn, chr=body.chromosome,
                                 af_threshold=body.af_threshold,
                                 n_permutations=body.permutations,
                                 verbose=False)

    job_id = _submit_job(_run)
    return {"job_id": job_id}


@app.post("/api/v1/heritability/gnn")
async def h2_gnn(request: Request, body: HeritabilityGnnRequest):
    """GNN heritability (async job — training)."""
    conn = _conn(request)

    def _run():
        from .heritability import gnn_heritability
        return gnn_heritability(conn, chr=body.chromosome,
                                start=body.start, end=body.end,
                                epochs=body.epochs, device=body.device,
                                verbose=False)

    job_id = _submit_job(_run)
    return {"job_id": job_id}


@app.post("/api/v1/heritability/report")
async def h2_report(request: Request, body: HeritabilityReportRequest):
    """Integrated heritability report (async job)."""
    conn = _conn(request)

    def _run():
        from .heritability import heritability_report
        chroms = [body.chromosome] if body.chromosome else None
        return heritability_report(conn, chromosomes=chroms,
                                   include_multiresolution=body.include_multiresolution,
                                   include_gnn=body.include_gnn,
                                   verbose=False)

    job_id = _submit_job(_run)
    return {"job_id": job_id}


# ---------------------------------------------------------------------------
# Multivariate endpoints
# ---------------------------------------------------------------------------

class GeneticCorrelationRequest(BaseModel):
    trait1: str
    trait2: str
    af_threshold: float = 0.01
    n_components: int = 20


class GMatrixRequest(BaseModel):
    trait_names: list[str]
    af_threshold: float = 0.01
    n_components: int = 20


class PleiotropyRequest(BaseModel):
    trait1_run_id: str
    trait2_run_id: str
    p_threshold: float = 5e-8


class PsiRequest(BaseModel):
    run_ids: list[str]
    trait_labels: list[str] | None = None
    p_threshold: float = 5e-8


@app.post("/api/v1/multivariate/correlation")
async def mv_correlation(request: Request, body: GeneticCorrelationRequest):
    """Spectral genetic correlation between two traits."""
    from .multivariate import spectral_genetic_covariance
    return spectral_genetic_covariance(_conn(request), body.trait1, body.trait2,
                                       af_threshold=body.af_threshold,
                                       n_components=body.n_components,
                                       verbose=False)


@app.post("/api/v1/multivariate/g-matrix")
async def mv_gmatrix(request: Request, body: GMatrixRequest):
    """Compute T×T G-matrix."""
    from .multivariate import g_matrix
    return g_matrix(_conn(request), body.trait_names,
                    af_threshold=body.af_threshold,
                    n_components=body.n_components, verbose=False)


@app.post("/api/v1/multivariate/coherence")
async def mv_coherence(request: Request, body: GeneticCorrelationRequest):
    """Spectral coherence between two traits."""
    from .multivariate import spectral_coherence
    return spectral_coherence(_conn(request), body.trait1, body.trait2,
                              af_threshold=body.af_threshold,
                              n_components=body.n_components, verbose=False)


@app.post("/api/v1/multivariate/pleiotropy")
async def mv_pleiotropy(request: Request, body: PleiotropyRequest):
    """Find pleiotropic genes between two GWAS runs."""
    from .multivariate import pleiotropic_genes
    return pleiotropic_genes(_conn(request), body.trait1_run_id,
                             body.trait2_run_id,
                             p_threshold=body.p_threshold, verbose=False)


@app.post("/api/v1/multivariate/psi")
async def mv_psi(request: Request, body: PsiRequest):
    """Compute Gene Pleiotropy Strength Index across multiple traits."""
    from .multivariate import gene_pleiotropy_score
    return gene_pleiotropy_score(_conn(request), body.run_ids,
                                 trait_labels=body.trait_labels,
                                 p_threshold=body.p_threshold,
                                 verbose=False)


# ---------------------------------------------------------------------------
# MET endpoints
# ---------------------------------------------------------------------------

class MetLoadRequest(BaseModel):
    csv_path: str
    genotype_col: str
    environment_col: str
    trait_cols: list[str]
    env_covariate_cols: list[str] | None = None
    rep_col: str | None = None
    block_col: str | None = None
    trial_id: str = "default"


class MetTraitRequest(BaseModel):
    trial_id: str = "default"
    trait: str


class MetImputeRequest(BaseModel):
    trial_id: str = "default"
    trait: str
    genotype_id: str
    environment_id: str
    k_neighbors: int = 10


@app.post("/api/v1/met/load")
async def met_load(request: Request, body: MetLoadRequest):
    """Import MET data from long-format CSV."""
    from .met import load_met_data
    return load_met_data(_conn(request), body.csv_path, body.genotype_col,
                         body.environment_col, body.trait_cols,
                         env_covariate_cols=body.env_covariate_cols,
                         rep_col=body.rep_col, block_col=body.block_col,
                         trial_id=body.trial_id, verbose=False)


@app.get("/api/v1/met/status/{trial_id}")
async def met_status(request: Request, trial_id: str):
    """Get MET trial summary."""
    result = _conn(request).execute_read(
        "MATCH (t:Trial {id: $tid}) RETURN t", {"tid": trial_id},
    ).single()
    if not result:
        raise HTTPException(404, f"Trial '{trial_id}' not found")
    return dict(result["t"])


@app.post("/api/v1/met/env-similarity")
async def met_env_similarity(request: Request, body: MetTraitRequest):
    """Genetic correlation between environments."""
    from .met import environment_similarity
    return environment_similarity(_conn(request), body.trial_id,
                                  body.trait, verbose=False)


@app.post("/api/v1/met/mega-environments")
async def met_mega_env(request: Request, body: MetTraitRequest):
    """Detect mega-environments."""
    from .met import mega_environments
    return mega_environments(_conn(request), body.trial_id,
                             body.trait, verbose=False)


@app.post("/api/v1/met/gxe")
async def met_gxe(request: Request, body: MetTraitRequest):
    """G×E variance decomposition."""
    from .met import gxe_decomposition
    return gxe_decomposition(_conn(request), body.trial_id,
                              body.trait, verbose=False)


@app.get("/api/v1/met/diagnostics/{trial_id}")
async def met_diagnostics(request: Request, trial_id: str):
    """Design connectivity diagnostics."""
    from .met import design_diagnostics
    return design_diagnostics(_conn(request), trial_id, verbose=False)


@app.post("/api/v1/met/impute")
async def met_impute(request: Request, body: MetImputeRequest):
    """Impute missing G×E cell."""
    from .met import graph_impute_missing
    return graph_impute_missing(_conn(request), body.trial_id, body.trait,
                                body.genotype_id, body.environment_id,
                                k_neighbors=body.k_neighbors, verbose=False)


@app.post("/api/v1/met/report")
async def met_report_api(request: Request, body: MetTraitRequest):
    """Comprehensive MET report (async job)."""
    conn = _conn(request)

    def _run():
        from .met import met_report
        return met_report(conn, body.trial_id, body.trait, verbose=False)

    job_id = _submit_job(_run)
    return {"job_id": job_id}


# ---------------------------------------------------------------------------
# PRS endpoints
# ---------------------------------------------------------------------------

class PrsRequest(BaseModel):
    run_id: str
    p_threshold: float = 5e-8


class PrsReportRequest(BaseModel):
    run_id: str
    p_threshold: float = 5e-8
    include_gnn: bool = False


@app.post("/api/v1/prs/classical")
async def prs_classical(request: Request, body: PrsRequest):
    """Standard weighted-sum PRS."""
    from .prs import classical_prs
    return classical_prs(_conn(request), body.run_id,
                          body.p_threshold, verbose=False)


@app.post("/api/v1/prs/graph-pruned")
async def prs_pruned(request: Request, body: PrsRequest):
    """LD-aware graph-pruned PRS."""
    from .prs import graph_pruned_prs
    return graph_pruned_prs(_conn(request), body.run_id,
                             body.p_threshold, verbose=False)


@app.post("/api/v1/prs/pathway")
async def prs_pathway(request: Request, body: PrsRequest):
    """Pathway-partitioned PRS."""
    from .prs import pathway_partitioned_prs
    return pathway_partitioned_prs(_conn(request), body.run_id,
                                    body.p_threshold, verbose=False)


@app.post("/api/v1/prs/report")
async def prs_report_api(request: Request, body: PrsReportRequest):
    """Compare all PRS methods (async job)."""
    conn = _conn(request)

    def _run():
        from .prs import prs_report
        return prs_report(conn, body.run_id, body.p_threshold,
                          include_gnn=body.include_gnn, verbose=False)

    job_id = _submit_job(_run)
    return {"job_id": job_id}


# ---------------------------------------------------------------------------
# MR endpoints
# ---------------------------------------------------------------------------

class MrRequest(BaseModel):
    exposure_run_id: str
    outcome_run_id: str
    p_threshold: float = 5e-8


@app.post("/api/v1/mr/report")
async def mr_report_api(request: Request, body: MrRequest):
    """Full MR analysis (async job)."""
    conn = _conn(request)

    def _run():
        from .mr import mr_report
        return mr_report(conn, body.exposure_run_id, body.outcome_run_id,
                         p_threshold=body.p_threshold, verbose=False)

    job_id = _submit_job(_run)
    return {"job_id": job_id}
