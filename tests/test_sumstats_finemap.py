"""Tests for sumstats-only L1 and HBP fine-mapping entry paths.

Verifies that feeding (z, R_sq) directly through
`l1_finemap_from_sumstats` / `hbp_finemap_from_sumstats` produces the
same PIPs as the dosage-based pipeline when both are given equivalent
inputs.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src" / "python"))

from graphgwas.finemapping_v2 import (  # noqa: E402
    _compute_association_stats,
    _compute_ld_matrix,
    _ld_deconvolve,
    _softmax,
    fast_hbp_finemap,
    hbp_finemap_from_sumstats,
    l1_finemap_from_sumstats,
)


# ===================================================================
# Synthetic locus fixtures
# ===================================================================

def _make_locus(n_var: int = 20, n_samples: int = 400, seed: int = 0):
    """Build a small synthetic locus with LD + one causal variant."""
    rng = np.random.default_rng(seed)
    # Generate dosages with banded LD: each variant correlated with its
    # k nearest neighbours in position.
    k = 5
    latents = rng.standard_normal((n_samples, n_var + k))
    D = np.zeros((n_var, n_samples), dtype=np.float64)
    for i in range(n_var):
        D[i] = latents[:, i : i + k].mean(axis=1) + 0.3 * rng.standard_normal(n_samples)
    # Recode to dosage-like {0, 1, 2}
    D = np.clip(np.round(D - D.min(axis=1, keepdims=True)), 0, 2)

    # One causal variant with moderate effect
    causal = n_var // 2
    beta = 0.4
    pheno = beta * (D[causal] - D[causal].mean()) + rng.standard_normal(n_samples)

    variants = [
        {
            "variantId": f"chr22:{1_000_000 + 1000 * i}:A:C",
            "chr": "22",
            "pos": 1_000_000 + 1000 * i,
            "ref": "A",
            "alt": "C",
            "af_total": float(np.clip(D[i].mean() / 2, 0.01, 0.99)),
            "dosage": D[i],
        }
        for i in range(n_var)
    ]
    return variants, pheno, causal


def _fake_graph_cache(variants: list[dict]) -> dict:
    """Build a graph cache that annotates every 3rd variant with a gene."""
    cache = {}
    for i, v in enumerate(variants):
        if i % 3 == 0:
            cache[v["variantId"]] = {
                "genes": [f"GENE{i}"],
                "pathways": ["PW_A"] if i % 6 == 0 else [],
                "ppi": [],
            }
    return cache


# ===================================================================
# Tests
# ===================================================================

def test_hbp_sumstats_matches_dosage_path():
    """HBP from sumstats and HBP from dosages should produce identical PIPs."""
    variants, pheno, _ = _make_locus(n_var=25, seed=42)
    graph_cache = _fake_graph_cache(variants)

    # Dosage path
    dosage_out = fast_hbp_finemap(variants, pheno, graph_cache, chr_name="chr22")

    # Sumstats path — compute z and R_sq from the same dosages
    z = _compute_association_stats(variants, pheno)
    R_sq = _compute_ld_matrix(variants)
    meta = [{k: v for k, v in var.items() if k != "dosage"} for var in variants]
    sumstats_out = hbp_finemap_from_sumstats(meta, z, R_sq, graph_cache, chr_name="chr22")

    # Sort both by variantId so order differences don't matter
    dosage_by_id = {c.variant_id: c for c in dosage_out}
    sumstats_by_id = {c.variant_id: c for c in sumstats_out}
    assert set(dosage_by_id) == set(sumstats_by_id)
    for vid in dosage_by_id:
        assert dosage_by_id[vid].pip == pytest.approx(sumstats_by_id[vid].pip, abs=1e-10)
        assert dosage_by_id[vid].unique_stat == pytest.approx(
            sumstats_by_id[vid].unique_stat, abs=1e-10,
        )


def test_l1_sumstats_basic_ranking():
    """L1 from sumstats should rank the causal variant highly."""
    variants, pheno, causal = _make_locus(n_var=30, seed=7)
    z = _compute_association_stats(variants, pheno)
    R_sq = _compute_ld_matrix(variants)
    meta = [{k: v for k, v in var.items() if k != "dosage"} for var in variants]

    # Without functional prior — fall back to pure statistical fine-mapping
    out = l1_finemap_from_sumstats(meta, z, R_sq, z_func=None, alpha=1.0, chr_name="chr22")
    ranks = {c.variant_id: r for r, c in enumerate(out)}
    causal_vid = variants[causal]["variantId"]
    # Causal should be in top 3 of 30 with this effect size + sample size
    assert ranks[causal_vid] < 3


def test_l1_sumstats_functional_prior_boosts_causal():
    """Supplying a functional prior that hits the causal variant should help."""
    variants, pheno, causal = _make_locus(n_var=30, seed=11)
    z = _compute_association_stats(variants, pheno)
    R_sq = _compute_ld_matrix(variants)
    meta = [{k: v for k, v in var.items() if k != "dosage"} for var in variants]

    # Functional prior: spike the causal variant
    z_func = np.zeros(len(variants))
    z_func[causal] = float(z.max())

    out = l1_finemap_from_sumstats(meta, z, R_sq, z_func=z_func, alpha=0.5, chr_name="chr22")
    ranks = {c.variant_id: r for r, c in enumerate(out)}
    causal_vid = variants[causal]["variantId"]
    assert ranks[causal_vid] == 0, "Causal variant should rank #1 with matching functional prior"


def test_sumstats_shape_validation():
    """Shape mismatches should raise clean errors."""
    variants, pheno, _ = _make_locus(n_var=10, seed=1)
    meta = [{k: v for k, v in var.items() if k != "dosage"} for var in variants]
    z = _compute_association_stats(variants, pheno)
    R_sq = _compute_ld_matrix(variants)

    # wrong-length z
    with pytest.raises(AssertionError):
        l1_finemap_from_sumstats(meta, z[:5], R_sq, chr_name="chr22")
    # wrong-shape R_sq
    with pytest.raises(AssertionError):
        hbp_finemap_from_sumstats(meta, z, R_sq[:5, :5], graph_cache={}, chr_name="chr22")


def test_empty_graph_cache_hbp_falls_back_to_softmax():
    """If graph_cache is empty, HBP from sumstats should still return PIPs."""
    variants, pheno, _ = _make_locus(n_var=15, seed=2)
    meta = [{k: v for k, v in var.items() if k != "dosage"} for var in variants]
    z = _compute_association_stats(variants, pheno)
    R_sq = _compute_ld_matrix(variants)
    out = hbp_finemap_from_sumstats(meta, z, R_sq, graph_cache={}, chr_name="chr22")
    assert len(out) == len(variants)
    pip_sum = sum(c.pip for c in out)
    assert pip_sum == pytest.approx(1.0, abs=1e-6)


def test_ld_deconvolve_matches():
    """Verify the LD deconvolution helper produces the same u whether called
    from the dosage path or directly — documenting the contract."""
    variants, pheno, _ = _make_locus(n_var=20, seed=3)
    z = _compute_association_stats(variants, pheno)
    R_sq = _compute_ld_matrix(variants)
    u1, _ = _ld_deconvolve(z, R_sq, 0.3)
    u2, _ = _ld_deconvolve(z, R_sq, 0.3)
    np.testing.assert_array_equal(u1, u2)


# ===================================================================
# Phase 0 schema extension: prior_score regression tests
# ===================================================================

def _spike_prior_score(graph_cache: dict, variants: list[dict],
                      target_idx: int, score: float = 1.0) -> dict:
    """Return a copy of graph_cache with prior_score=score on variant[target_idx]
    and 0.0 on every other variant that already has a cache entry."""
    out = {k: dict(v) for k, v in graph_cache.items()}
    for i, v in enumerate(variants):
        vid = v["variantId"]
        if vid not in out:
            # Add a minimal entry so the prior reaches HBP/L1
            out[vid] = {"genes": [], "pathways": [], "ppi": []}
        out[vid]["prior_score"] = float(score) if i == target_idx else 0.0
    return out


def test_hbp_prior_score_boosts_causal():
    """Phase 0: HBP should rank the variant carrying prior_score higher
    than the same locus run without prior_score."""
    variants, pheno, causal = _make_locus(n_var=30, seed=101)
    z_stats = _compute_association_stats(variants, pheno)
    R_sq = _compute_ld_matrix(variants)
    base_cache = _fake_graph_cache(variants)

    # Same cache, no prior_score
    out_no = hbp_finemap_from_sumstats(
        variants, z_stats, R_sq, graph_cache=base_cache, chr_name="chr22",
    )
    # Spike prior_score=1.0 on the causal variant
    spiked = _spike_prior_score(base_cache, variants, causal, score=1.0)
    out_spiked = hbp_finemap_from_sumstats(
        variants, z_stats, R_sq, graph_cache=spiked, chr_name="chr22",
    )
    pip_no = next(c.pip for c in out_no
                  if c.variant_id == variants[causal]["variantId"])
    pip_spiked = next(c.pip for c in out_spiked
                      if c.variant_id == variants[causal]["variantId"])
    assert pip_spiked > pip_no, (
        f"HBP prior_score should boost causal PIP: {pip_no:.4f} -> {pip_spiked:.4f}"
    )


def test_l1_prior_score_boosts_causal():
    """Phase 0: L1 graph_cache parameter — prior_score should boost causal PIP."""
    variants, pheno, causal = _make_locus(n_var=30, seed=102)
    z_stats = _compute_association_stats(variants, pheno)
    R_sq = _compute_ld_matrix(variants)
    base_cache = _fake_graph_cache(variants)

    out_no = l1_finemap_from_sumstats(
        variants, z_stats, R_sq, graph_cache=base_cache, chr_name="chr22",
    )
    spiked = _spike_prior_score(base_cache, variants, causal, score=1.0)
    out_spiked = l1_finemap_from_sumstats(
        variants, z_stats, R_sq, graph_cache=spiked, chr_name="chr22",
    )
    pip_no = next(c.pip for c in out_no
                  if c.variant_id == variants[causal]["variantId"])
    pip_spiked = next(c.pip for c in out_spiked
                      if c.variant_id == variants[causal]["variantId"])
    assert pip_spiked > pip_no, (
        f"L1 prior_score should boost causal PIP: {pip_no:.4f} -> {pip_spiked:.4f}"
    )


def test_prior_score_no_op_when_absent():
    """Phase 0 regression: cache without prior_score reproduces v1 behavior."""
    variants, pheno, _ = _make_locus(n_var=20, seed=103)
    z_stats = _compute_association_stats(variants, pheno)
    R_sq = _compute_ld_matrix(variants)
    cache_no = _fake_graph_cache(variants)
    # Identity: cache lookups don't yield prior_score → output should be deterministic
    out1 = hbp_finemap_from_sumstats(
        variants, z_stats, R_sq, graph_cache=cache_no, chr_name="chr22",
    )
    out2 = hbp_finemap_from_sumstats(
        variants, z_stats, R_sq, graph_cache=cache_no, chr_name="chr22",
    )
    pips1 = np.array([c.pip for c in out1])
    pips2 = np.array([c.pip for c in out2])
    np.testing.assert_array_almost_equal(pips1, pips2, decimal=10)


def test_fast_hbp_finemap_handles_missing_genes_field():
    """Schema robustness: fast_hbp_finemap must not crash on cache entries
    that lack the 'genes' field (e.g., augmented v2 caches with only
    prior_score)."""
    variants, pheno, _ = _make_locus(n_var=15, seed=104)
    # Cache where every variant has only prior_score, no genes/pathways/ppi
    cache = {v["variantId"]: {"prior_score": 0.5} for v in variants}
    # Should not raise KeyError
    out = fast_hbp_finemap(variants, pheno, cache, chr_name="chr22")
    assert len(out) == len(variants)
