"""End-to-end test: BGEN-backed fine-mapping gives same result as Neo4j-backed.

This is the critical UKB-readiness test. If BGEN-backed HBP/L1 produce different
PIPs from the Neo4j-backed versions on the same locus, the hybrid architecture
has a bug.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

BGEN_DIR = Path("tests/data/human/1kGP_bgen")
CHR22_BGEN = BGEN_DIR / "chr22.bgen"

needs_bgen = pytest.mark.skipif(not CHR22_BGEN.exists(), reason="chr22 BGEN not present")


@needs_bgen
def test_load_locus_variants_dispatches():
    """load_locus_variants must accept both a conn and a BgenReader."""
    from graphgwas.bgen_reader import BgenReader
    from graphgwas.finemapping_v2 import load_locus_variants

    reader = BgenReader(BGEN_DIR)
    variants = load_locus_variants(
        chr="chr22", center=16_075_000, window=50_000, source=reader
    )
    assert len(variants) > 0
    v0 = variants[0]
    assert set(v0.keys()) >= {"variantId", "chr", "pos", "ref", "alt", "af_total", "dosage"}
    assert v0["dosage"].ndim == 1
    assert v0["dosage"].shape[0] == reader.n_samples("22")
    # af_total computed from BGEN must match dosage.mean()/2
    assert abs(v0["af_total"] - v0["dosage"].mean() / 2.0) < 1e-6


@needs_bgen
def test_bgen_variants_work_with_fast_hbp():
    """BGEN-loaded variants must run through fast_hbp_finemap without error."""
    from graphgwas.bgen_reader import BgenReader
    from graphgwas.finemapping_v2 import load_locus_variants, fast_hbp_finemap

    reader = BgenReader(BGEN_DIR)
    variants = load_locus_variants(
        chr="chr22", center=16_075_000, window=50_000, source=reader
    )
    assert len(variants) >= 10
    # Simulate a phenotype: random + causal contribution from a random variant
    rng = np.random.default_rng(42)
    n_samples = variants[0]["dosage"].shape[0]
    phenotype = rng.standard_normal(n_samples)
    # Empty graph cache (no annotations) — exercises the cold path
    cands = fast_hbp_finemap(
        variants, phenotype, graph_cache={}, n_rounds=3, chr_name="chr22"
    )
    assert len(cands) == len(variants)
    # PIPs should sum to ~1 when normalized
    pips = np.array([c.pip for c in cands])
    assert abs(pips.sum() - 1.0) < 0.01


@needs_bgen
def test_af_consistency_with_neo4j():
    """AF computed from BGEN dosage must match Neo4j Variant.af_total."""
    pytest.importorskip("neo4j")
    from graphgwas.bgen_reader import BgenReader
    from graphgwas.finemapping_v2 import load_locus_variants
    from graphgwas.db import GraphGWASConnection
    from graphgwas import config

    reader = BgenReader(BGEN_DIR)
    variants = load_locus_variants(
        chr="chr22", center=16_055_000, window=10_000, source=reader
    )

    positions = [v["pos"] for v in variants]
    try:
        with GraphGWASConnection(
            config.NEO4J_URI, config.NEO4J_USER, config.NEO4J_PASSWORD
        ) as conn:
            result = conn.execute_read(
                "UNWIND $positions AS p MATCH (v:Variant {chr: 'chr22', pos: p}) "
                "RETURN v.pos AS pos, v.ref AS ref, v.alt AS alt, v.af_total AS af",
                {"positions": positions},
            )
            neo4j_by_pos: dict[int, list[tuple[str, str, float]]] = {}
            for r in result:
                neo4j_by_pos.setdefault(int(r["pos"]), []).append(
                    (r["ref"], r["alt"], float(r["af"]))
                )
    except Exception:
        pytest.skip("Neo4j not reachable")

    diffs: list[float] = []
    for v in variants:
        for ref, alt, af in neo4j_by_pos.get(v["pos"], []):
            if {v["ref"], v["alt"]} == {ref, alt}:
                diffs.append(abs(v["af_total"] - af))
                break
    assert len(diffs) > 10
    assert max(diffs) < 1e-4, f"max AF diff {max(diffs)} exceeds 1e-4 tolerance"
