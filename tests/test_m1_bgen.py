"""Unit test: M1 LD-pruned co-occurrence runs end-to-end on BGEN input.

Uses chr22 1KG BGEN with a synthetic quantitative phenotype to exercise
both the source-agnostic core and the BGEN convenience wrapper.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

BGEN_DIR = Path("tests/data/human/1kGP_bgen")
CHR22_BGEN = BGEN_DIR / "chr22.bgen"

needs_bgen = pytest.mark.skipif(
    not CHR22_BGEN.exists(), reason="chr22 BGEN not present"
)


@needs_bgen
def test_m1_from_data_runs():
    """M1 source-agnostic core runs on small synthetic case."""
    from graphgwas.epistasis_v2 import ld_pruned_cooccurrence_from_data

    rng = np.random.default_rng(0)
    n = 200
    # Simulate 20 variants with two interacting pair
    vs = []
    ds = []
    for i in range(20):
        af = rng.uniform(0.1, 0.45)
        d = rng.binomial(2, af, size=n).astype(float)
        ds.append(d)
        vs.append({"variantId": f"v{i}", "pos": 1000 + i * 200_000,
                   "chr": "chr22", "af_total": float(d.mean() / 2)})
    # Build pheno influenced by variant-0 × variant-10
    pheno = 1.5 * (ds[0] - ds[0].mean()) * (ds[10] - ds[10].mean())
    pheno += rng.standard_normal(n) * 0.5

    results = ld_pruned_cooccurrence_from_data(
        variants=vs, dosage_list=ds, phenotype=pheno,
        r2_prune=0.5, min_cocarriers=3,
        min_distance_bp=100_000, max_variants=50,
        verbose=False,
    )
    # Should test at least some pairs without errors
    assert isinstance(results, list)


@needs_bgen
def test_m1_bgen_wrapper_smoke():
    """The BGEN wrapper should return without error on a real locus."""
    from graphgwas.bgen_reader import BgenReader
    from graphgwas.epistasis_v2 import ld_pruned_cooccurrence_bgen

    reader = BgenReader(BGEN_DIR)
    samples = [str(s) for s in reader.samples("22")]
    rng = np.random.default_rng(42)
    # Random quantitative phenotype per BGEN sample
    pheno_by_sample = {s: float(rng.standard_normal()) for s in samples}
    results = ld_pruned_cooccurrence_bgen(
        reader, chr="22", start=16_050_000, end=16_150_000,
        phenotype_by_sample=pheno_by_sample,
        r2_prune=0.5, max_variants=100,
        min_cocarriers=3, min_distance_bp=10_000,
        verbose=False,
    )
    assert isinstance(results, list)


@needs_bgen
def test_m1_bgen_wrapper_reuses_core():
    """Wrapper should produce identical results as calling core directly
    with the same loaded data."""
    from graphgwas.bgen_reader import BgenReader
    from graphgwas.epistasis_v2 import (
        ld_pruned_cooccurrence_from_data, ld_pruned_cooccurrence_bgen,
    )

    reader = BgenReader(BGEN_DIR)
    samples = [str(s) for s in reader.samples("22")]
    rng = np.random.default_rng(123)
    pheno_by_sample = {s: float(rng.standard_normal()) for s in samples}

    # Path A: wrapper
    res_a = ld_pruned_cooccurrence_bgen(
        reader, chr="22", start=16_050_000, end=16_100_000,
        phenotype_by_sample=pheno_by_sample,
        r2_prune=0.5, max_variants=50,
        min_cocarriers=3, min_distance_bp=10_000,
        verbose=False,
    )
    # Path B: manual load + core
    vdf, dos = reader.load_locus("22", 16_050_000, 16_100_000, format="dosage")
    vs, ds = [], []
    pheno = np.array([pheno_by_sample[s] for s in samples], dtype=float)
    for i in range(len(vdf)):
        d = dos[:, i]
        m = np.nanmean(d)
        if not np.isfinite(m):
            continue
        af = float(m / 2)
        if af < 0.01 or af > 0.99:
            continue
        di = np.where(np.isnan(d), m, d)
        mac = min(int(di.sum()), int(2 * len(di) - di.sum()))
        if mac < 10:
            continue
        row = vdf.iloc[i]
        vs.append({"variantId": f"{row['chr']}:{row['pos']}:{row['a1']}:{row['a2']}",
                   "pos": int(row["pos"]), "chr": row["chr"], "af_total": af})
        ds.append(di)
        if len(vs) >= 50 * 3:
            break
    res_b = ld_pruned_cooccurrence_from_data(
        variants=vs, dosage_list=ds, phenotype=pheno,
        r2_prune=0.5, max_variants=50,
        min_cocarriers=3, min_distance_bp=10_000,
        verbose=False,
    )
    # Result shape should be identical (same random pheno → same pairs)
    assert len(res_a) == len(res_b)
    # Top result (if any) should have the same variant pair
    if res_a:
        assert res_a[0].variant_1 == res_b[0].variant_1
        assert res_a[0].variant_2 == res_b[0].variant_2
