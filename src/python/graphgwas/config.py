"""Connection defaults and GWAS constants."""

import os

# Neo4j connection defaults (shared with GraphMana)
NEO4J_URI = os.environ.get("GRAPHGWAS_NEO4J_URI", "bolt://localhost:7688")
NEO4J_USER = os.environ.get("GRAPHGWAS_NEO4J_USER", "neo4j")
NEO4J_PASSWORD = os.environ.get("GRAPHGWAS_NEO4J_PASSWORD", "graphgwas")
NEO4J_DATABASE = os.environ.get("GRAPHGWAS_NEO4J_DATABASE", None)

# GWAS significance thresholds
P_GENOME_WIDE = 5e-8
P_SUGGESTIVE = 1e-5

# QC defaults
HWE_P_THRESHOLD = 1e-10
CALL_RATE_MIN = 0.95
MAC_MIN = 5
SAMPLE_MISSINGNESS_MAX = 0.05
HET_RATE_SD = 3.0  # exclude samples outside mean ± 3 SD

# Genotype encoding constants (must match GraphMana PackedGenotypeReader)
GT_HOM_REF = 0
GT_HET = 1
GT_HOM_ALT = 2
GT_MISSING = 3

# Number of samples in current database
# Auto-detected on first use via detect_n_samples(); falls back to this default.
N_SAMPLES = 1011
_N_SAMPLES_DETECTED = False


def detect_n_samples(conn) -> int:
    """Auto-detect N_SAMPLES from gt_packed byte array size.

    N_SAMPLES = size(gt_packed) * 4, because each byte holds 4 genotypes
    (2 bits each). This is the authoritative sample count for genotype decoding.
    """
    global N_SAMPLES, _N_SAMPLES_DETECTED
    if _N_SAMPLES_DETECTED:
        return N_SAMPLES
    try:
        result = conn.execute_read(
            "MATCH (v:Variant) WHERE v.gt_packed IS NOT NULL "
            "RETURN size(v.gt_packed) AS sz LIMIT 1"
        )
        rec = result.single()
        if rec and rec["sz"]:
            N_SAMPLES = rec["sz"] * 4
            _N_SAMPLES_DETECTED = True
    except Exception:
        pass
    return N_SAMPLES
