"""Unit tests for graphgwas.summary_import.

Parses synthetic PLINK2/regenie files — no Neo4j required.
"""
from __future__ import annotations

from pathlib import Path

import pytest


def test_parse_plink2_linear(tmp_path: Path):
    from graphgwas.summary_import import parse_plink2_glm, _build_variant_id
    p = tmp_path / "height.glm.linear"
    p.write_text(
        "#CHROM\tPOS\tID\tREF\tALT\tA1\tTEST\tOBS_CT\tBETA\tSE\tT_STAT\tP\n"
        "22\t16050005\t.\tC\tT\tT\tADD\t3000\t0.05\t0.01\t5.0\t5.7e-7\n"
        "22\t16050012\t.\tA\tG\tG\tADD\t3000\t-0.02\t0.015\t-1.33\t0.18\n"
    )
    rows = list(parse_plink2_glm(p, p_threshold=1e-3))
    # Only the significant row survives
    assert len(rows) == 1
    r = rows[0]
    assert r["variantId"] == "chr22:16050005:C:T"
    assert r["method"] == "plink2_linear"
    assert r["beta"] == pytest.approx(0.05)
    assert r["p_value"] == pytest.approx(5.7e-7)
    assert r["p_value_log10"] > 6.0


def test_parse_plink2_logistic_flips_for_ref_a1(tmp_path: Path):
    from graphgwas.summary_import import parse_plink2_glm
    p = tmp_path / "bmi.glm.logistic"
    p.write_text(
        "#CHROM\tPOS\tID\tREF\tALT\tA1\tTEST\tOBS_CT\tOR\tSE\tZ_STAT\tP\n"
        "1\t12345\t.\tA\tG\tA\tADD\t5000\t2.0\t0.1\t7.0\t1e-10\n"
    )
    rows = list(parse_plink2_glm(p))
    # A1 is REF — effect should be flipped so OR becomes 1/2.0 = 0.5
    assert rows[0]["odds_ratio"] == pytest.approx(0.5)


def test_parse_regenie(tmp_path: Path):
    from graphgwas.summary_import import parse_regenie
    p = tmp_path / "trait.regenie"
    p.write_text(
        "CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ INFO N TEST BETA SE CHISQ LOG10P\n"
        "22 16050005 . C T 0.01 1.0 3000 ADD-WGR-LR 0.05 0.01 25 7.3\n"
        "22 16050012 . A G 0.2 1.0 3000 ADD-WGR-LR -0.02 0.015 1.77 0.25\n"
    )
    rows = list(parse_regenie(p, p_threshold=1e-3))
    assert len(rows) == 1
    r = rows[0]
    assert r["variantId"] == "chr22:16050005:C:T"
    assert r["p_value_log10"] == pytest.approx(7.3)
    assert r["beta"] == pytest.approx(0.05)
