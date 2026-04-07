"""Schema management — index creation, schema audit, GraphMana compatibility check."""

from .db import GraphGWASConnection

GWAS_INDEXES = [
    "CREATE RANGE INDEX variant_chr_pos IF NOT EXISTS FOR (v:Variant) ON (v.chr, v.pos)",
    "CREATE INDEX sample_is_case IF NOT EXISTS FOR (s:Sample) ON (s.is_case)",
    "CREATE INDEX sample_is_control IF NOT EXISTS FOR (s:Sample) ON (s.is_control)",
    "CREATE RANGE INDEX variant_af_total IF NOT EXISTS FOR (v:Variant) ON (v.af_total)",
    "CREATE INDEX gene_symbol IF NOT EXISTS FOR (g:Gene) ON (g.symbol)",
    "CREATE INDEX pathway_name IF NOT EXISTS FOR (p:Pathway) ON (p.name)",
    "CREATE INDEX assoc_run_id IF NOT EXISTS FOR (ar:AssociationResult) ON (ar.run_id)",
    "CREATE RANGE INDEX assoc_pvalue IF NOT EXISTS FOR (ar:AssociationResult) ON (ar.p_value)",
    "CREATE INDEX assoc_phenotype IF NOT EXISTS FOR (ar:AssociationResult) ON (ar.phenotype_key)",
    "CREATE CONSTRAINT gwas_study_id IF NOT EXISTS FOR (gs:GWASStudy) REQUIRE gs.id IS UNIQUE",
]


def ensure_indexes(conn: GraphGWASConnection) -> list[str]:
    """Create all GWAS indexes. Returns list of created/existing index names."""
    results = []
    for ddl in GWAS_INDEXES:
        try:
            conn.execute_write(ddl)
            results.append(f"OK: {ddl[:60]}...")
        except Exception as e:
            results.append(f"SKIP: {ddl[:60]}... ({e})")
    return results


def audit_schema(conn: GraphGWASConnection) -> dict:
    """Return JSON-like report of database schema state."""
    report = {}

    # Node counts
    for label in ["Variant", "Sample", "Gene", "Pathway", "GOTerm",
                  "Population", "AssociationResult", "GWASStudy"]:
        result = conn.execute_read(f"MATCH (n:{label}) RETURN count(n) AS c")
        rec = result.single()
        report[f"node_{label}"] = rec["c"] if rec else 0

    # Relationship counts
    for rel in ["HAS_CONSEQUENCE", "IN_PATHWAY", "HAS_GO_TERM",
                "ON_CHROMOSOME", "IN_POPULATION", "NEXT",
                "FOR_VARIANT", "IN_STUDY"]:
        result = conn.execute_read(
            f"MATCH ()-[r:{rel}]->() RETURN count(r) AS c"
        )
        rec = result.single()
        report[f"rel_{rel}"] = rec["c"] if rec else 0

    # Phenotype status
    for prop, label in [("is_case", "cases"), ("is_control", "controls")]:
        result = conn.execute_read(
            f"MATCH (s:Sample) WHERE s.{prop} = true RETURN count(s) AS c"
        )
        rec = result.single()
        report[f"n_{label}"] = rec["c"] if rec else 0

    result = conn.execute_read(
        "MATCH (s:Sample) WHERE s.phenotypes_loaded = true RETURN count(s) AS c"
    )
    rec = result.single()
    report["n_phenotyped"] = rec["c"] if rec else 0

    return report


def verify_graphmana_compat(conn: GraphGWASConnection) -> list[str]:
    """Check that expected GraphMana schema exists. Returns list of issues."""
    issues = []

    # Check required node labels
    result = conn.execute_read("CALL db.labels() YIELD label RETURN collect(label) AS labels")
    labels = set(result.single()["labels"])
    for required in ["Variant", "Sample", "Gene", "Pathway"]:
        if required not in labels:
            issues.append(f"Missing required node label: {required}")

    # Check gt_packed exists on Variants
    result = conn.execute_read(
        "MATCH (v:Variant) WHERE v.gt_packed IS NOT NULL RETURN count(v) AS c LIMIT 1"
    )
    rec = result.single()
    if not rec or rec["c"] == 0:
        issues.append("No Variant nodes have gt_packed property")

    # Check packed_index on Samples
    result = conn.execute_read(
        "MATCH (s:Sample) WHERE s.packed_index IS NOT NULL RETURN count(s) AS c"
    )
    rec = result.single()
    if not rec or rec["c"] == 0:
        issues.append("No Sample nodes have packed_index property")

    # Check HAS_CONSEQUENCE edges
    result = conn.execute_read(
        "MATCH ()-[r:HAS_CONSEQUENCE]->() RETURN count(r) AS c LIMIT 1"
    )
    rec = result.single()
    if not rec or rec["c"] == 0:
        issues.append("No HAS_CONSEQUENCE edges found")

    return issues
