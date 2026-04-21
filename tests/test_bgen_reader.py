"""Unit tests for graphgwas.bgen_reader — validates hybrid architecture.

Requires the 1KG chr22 BGEN at tests/data/human/1kGP_bgen/chr22.bgen and
a live Neo4j DB at the default port with the human 1KG graph loaded.

Run with: pytest tests/test_bgen_reader.py -v
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

BGEN_DIR = Path("tests/data/human/1kGP_bgen")
CHR22_BGEN = BGEN_DIR / "chr22.bgen"


needs_bgen = pytest.mark.skipif(
    not CHR22_BGEN.exists(), reason="chr22 BGEN not present — skip"
)


@needs_bgen
def test_load_locus_shape():
    from graphgwas.bgen_reader import BgenReader
    reader = BgenReader(BGEN_DIR)
    variants, dosage = reader.load_locus("chr22", 16_050_000, 16_060_000)
    assert len(variants) > 0
    assert dosage.shape == (reader.n_samples("22"), len(variants))
    assert dosage.dtype == np.float32
    assert variants.columns.tolist() == ["chr", "pos", "a1", "a2"]
    assert (variants["pos"] >= 16_050_000).all()
    assert (variants["pos"] <= 16_060_000).all()


@needs_bgen
def test_dosage_range():
    from graphgwas.bgen_reader import BgenReader
    reader = BgenReader(BGEN_DIR)
    _, dosage = reader.load_locus("chr22", 16_050_000, 16_055_000)
    assert dosage.min() >= 0.0
    assert dosage.max() <= 2.0


@needs_bgen
def test_hardcall_rounding():
    from graphgwas.bgen_reader import BgenReader
    reader = BgenReader(BGEN_DIR)
    _, hc = reader.load_locus("chr22", 16_050_000, 16_051_000, format="hardcall")
    assert hc.dtype == np.int8
    assert set(np.unique(hc).tolist()).issubset({0, 1, 2})


@needs_bgen
def test_empty_locus():
    from graphgwas.bgen_reader import BgenReader
    reader = BgenReader(BGEN_DIR)
    variants, dosage = reader.load_locus("chr22", 10_000, 10_001)
    assert len(variants) == 0
    assert dosage.shape == (reader.n_samples("22"), 0)


@needs_bgen
def test_chr_prefix_normalization():
    from graphgwas.bgen_reader import BgenReader
    reader = BgenReader(BGEN_DIR)
    v1, _ = reader.load_locus("chr22", 16_050_000, 16_051_000)
    v2, _ = reader.load_locus("22", 16_050_000, 16_051_000)
    assert len(v1) == len(v2)
    assert (v1["pos"].values == v2["pos"].values).all()


@needs_bgen
def test_af_matches_neo4j():
    """BGEN minor_allele_dosage should equal Neo4j Variant.af_total exactly."""
    pytest.importorskip("neo4j")
    from graphgwas.bgen_reader import BgenReader
    from graphgwas.db import GraphGWASConnection
    from graphgwas import config

    reader = BgenReader(BGEN_DIR)
    variants, dosage = reader.load_locus("chr22", 16_050_000, 16_052_000)
    bgen_af = dosage.mean(axis=0) / 2.0

    positions = variants["pos"].tolist()
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
    for i, row in variants.iterrows():
        pos = int(row["pos"])
        for ref, alt, af in neo4j_by_pos.get(pos, []):
            if {row["a1"], row["a2"]} == {ref, alt}:
                diffs.append(abs(bgen_af[i] - af))
                break
    assert len(diffs) > 10, f"matched too few variants: {len(diffs)}"
    diffs_arr = np.array(diffs)
    # Allele frequencies should match exactly (same samples, same genotypes)
    assert diffs_arr.max() < 1e-4, f"max AF diff {diffs_arr.max()} exceeds tolerance"
