"""Association result storage and querying.

Stores results as (:AssociationResult) nodes connected to (:Variant) via FOR_VARIANT,
and to (:GWASStudy) via IN_STUDY. All results are queryable graph objects.
"""

from __future__ import annotations


from .db import GraphGWASConnection


def create_gwas_study(conn: GraphGWASConnection, run_id: str, phenotype_key: str,
                      method: str, n_cases: int, n_controls: int,
                      n_variants_tested: int, region: str = "genome-wide",
                      notes: str = "") -> str:
    """Create a GWASStudy audit node for this run.

    Returns the study node ID.
    """
    study_id = f"study_{run_id}"
    conn.execute_write(
        """
        MERGE (gs:GWASStudy {id: $id})
        ON CREATE SET
            gs.run_id = $run_id,
            gs.name = $name,
            gs.phenotype_key = $phenotype_key,
            gs.method = $method,
            gs.n_cases = $n_cases,
            gs.n_controls = $n_controls,
            gs.n_variants_tested = $n_variants_tested,
            gs.region = $region,
            gs.significance_threshold = 5e-8,
            gs.timestamp = datetime(),
            gs.notes = $notes,
            gs.software_version = 'graphgwas-0.1.0'
        """,
        {
            "id": study_id,
            "run_id": run_id,
            "name": f"GWAS {phenotype_key} ({method})",
            "phenotype_key": phenotype_key,
            "method": method,
            "n_cases": n_cases,
            "n_controls": n_controls,
            "n_variants_tested": n_variants_tested,
            "region": region,
            "notes": notes,
        },
    )
    return study_id


def store_results(conn: GraphGWASConnection, results: list[dict],
                  run_id: str, phenotype_key: str,
                  batch_size: int = 5000) -> int:
    """Store association results as AssociationResult nodes.

    Creates FOR_VARIANT edges to existing Variant nodes.
    IN_STUDY edges are linked in bulk after all batches via link_study().

    Returns number of results stored.
    """
    n_stored = 0

    for i in range(0, len(results), batch_size):
        batch = []
        for r in results[i : i + batch_size]:
            batch.append({
                "ar_id": f"gwas_{run_id}_{r.get('variantId', 'unknown')}",
                "run_id": run_id,
                "variant_id": r.get("variantId", ""),
                "phenotype_key": phenotype_key,
                "method": r.get("method", "unknown"),
                "beta": r.get("beta", None),
                "se": r.get("se", None),
                "p_value": r.get("p_value", 1.0),
                "p_value_log10": r.get("p_value_log10", 0.0),
                "odds_ratio": r.get("odds_ratio", None),
                "ci_lower": r.get("ci_lower", None),
                "ci_upper": r.get("ci_upper", None),
                "n_cases": r.get("n_cases", 0),
                "n_controls": r.get("n_controls", 0),
                "maf": r.get("af", None),
                "mac": r.get("mac", None),
            })

        result = conn.execute_write(
            """
            UNWIND $batch AS row
            CREATE (ar:AssociationResult {
                id: row.ar_id,
                run_id: row.run_id,
                variant_id: row.variant_id,
                phenotype_key: row.phenotype_key,
                method: row.method,
                beta: row.beta,
                se: row.se,
                p_value: row.p_value,
                p_value_log10: row.p_value_log10,
                odds_ratio: row.odds_ratio,
                ci_lower: row.ci_lower,
                ci_upper: row.ci_upper,
                n_cases: row.n_cases,
                n_controls: row.n_controls,
                maf: row.maf,
                mac: row.mac,
                timestamp: datetime()
            })
            WITH ar, row
            MATCH (v:Variant {variantId: row.variant_id})
            CREATE (ar)-[:FOR_VARIANT]->(v)
            RETURN count(ar) AS stored
            """,
            {"batch": batch},
        )
        rec = result.single()
        n_stored += rec["stored"] if rec else 0

    return n_stored


def link_study(conn: GraphGWASConnection, run_id: str) -> int:
    """Link all AssociationResult nodes for a run to their GWASStudy node.

    Called once after all results are stored, avoiding per-batch MATCH on GWASStudy.
    """
    study_id = f"study_{run_id}"
    result = conn.execute_write(
        """
        MATCH (ar:AssociationResult {run_id: $run_id})
        WITH ar
        MATCH (gs:GWASStudy {id: $study_id})
        CREATE (ar)-[:IN_STUDY]->(gs)
        RETURN count(ar) AS linked
        """,
        {"run_id": run_id, "study_id": study_id},
    )
    rec = result.single()
    return rec["linked"] if rec else 0


def query_top_hits(conn: GraphGWASConnection, run_id: str,
                   p_threshold: float = 5e-8, limit: int = 100) -> list[dict]:
    """Query significant association results for a run."""
    result = conn.execute_read(
        """
        MATCH (ar:AssociationResult {run_id: $run_id})-[:FOR_VARIANT]->(v:Variant)
        WHERE ar.p_value < $threshold
        OPTIONAL MATCH (v)-[:HAS_CONSEQUENCE]->(g:Gene)
        RETURN ar.variant_id AS variant, v.chr AS chr, v.pos AS pos,
               ar.p_value AS p_value, ar.p_value_log10 AS log10p,
               ar.beta AS beta, ar.odds_ratio AS or,
               ar.method AS method, ar.maf AS maf,
               collect(DISTINCT g.symbol) AS genes
        ORDER BY ar.p_value ASC
        LIMIT $limit
        """,
        {"run_id": run_id, "threshold": p_threshold, "limit": limit},
    )
    return result.data()


def query_manhattan_data(conn: GraphGWASConnection, run_id: str) -> list[dict]:
    """Query Manhattan plot data: chr, pos, -log10(p) for all results."""
    result = conn.execute_read(
        """
        MATCH (ar:AssociationResult {run_id: $run_id})-[:FOR_VARIANT]->(v:Variant)
        RETURN v.chr AS chr, v.pos AS pos, ar.p_value_log10 AS log10p
        ORDER BY v.chr, v.pos
        """,
        {"run_id": run_id},
    )
    return result.data()


def delete_run(conn: GraphGWASConnection, run_id: str) -> int:
    """Delete all results and study node for a GWAS run."""
    result = conn.execute_write(
        """
        MATCH (ar:AssociationResult {run_id: $run_id})
        DETACH DELETE ar
        WITH count(*) AS deleted_results
        OPTIONAL MATCH (gs:GWASStudy {run_id: $run_id})
        DETACH DELETE gs
        RETURN deleted_results
        """,
        {"run_id": run_id},
    )
    rec = result.single()
    return rec["deleted_results"] if rec else 0
