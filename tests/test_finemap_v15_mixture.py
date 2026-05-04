"""Unit tests for v0.1.5 GAFM/HBP enhancements:
λ_GC deflation, SBayesRC-style mixture-prior posterior, credible sets.

Run: pytest tests/test_finemap_v15_mixture.py -v
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src" / "python"))

from graphgwas.finemapping_v2 import (  # noqa: E402
    DEFAULT_MIXTURE_GAMMA,
    DEFAULT_MIXTURE_PI,
    apply_mixture_posterior,
    credible_set_from_pips,
    deflate_z_for_lambda_gc,
)


def test_deflation_basic():
    z = np.array([10.0, 5.0, 1.0])
    assert np.allclose(deflate_z_for_lambda_gc(z, 4.0), [5.0, 2.5, 0.5])


def test_deflation_identity_when_lambda_one():
    z = np.array([3.0, -2.0, 0.5])
    assert np.allclose(deflate_z_for_lambda_gc(z, 1.0), z)


def test_deflation_passthrough_on_invalid():
    z = np.array([3.0, -2.0])
    # Negative or zero λ_GC: passthrough (defensive)
    assert np.allclose(deflate_z_for_lambda_gc(z, 0.0), z)
    assert np.allclose(deflate_z_for_lambda_gc(z, -1.0), z)


def test_mixture_preserves_sum():
    rng = np.random.default_rng(0)
    n = 100
    pips = rng.dirichlet(np.ones(n))
    z = rng.standard_normal(n) * 0.5
    z[10] = 8.0
    new = apply_mixture_posterior(pips, z, n_samples=2400)
    assert abs(new.sum() - pips.sum()) < 1e-6


def test_mixture_sharpens_at_large_z():
    """Variant with z=8 should attract most of the mass after reweighting."""
    rng = np.random.default_rng(0)
    n = 100
    pips = rng.dirichlet(np.ones(n))  # roughly uniform
    z = rng.standard_normal(n) * 0.5
    z[10] = 8.0
    new = apply_mixture_posterior(pips, z, n_samples=2400)
    assert int(np.argmax(new)) == 10
    assert new[10] > 0.9, f"top PIP only {new[10]:.3f}"


def test_credible_set_single_variant_at_high_pip():
    pips = np.array([0.0, 0.99, 0.005, 0.005])
    cs = credible_set_from_pips(pips, coverage=0.95)
    assert cs == [1]


def test_credible_set_handles_all_zero():
    pips = np.zeros(5)
    cs = credible_set_from_pips(pips, coverage=0.95)
    # all-zero edge case: return single argmax
    assert len(cs) == 1


def test_default_mixture_components_match_sbayesrc():
    """SBayesRC's startPi[1:5] and gamma[1:5]; verifies we're not drifting."""
    assert DEFAULT_MIXTURE_PI == (0.005, 0.003, 0.001, 0.001)
    assert DEFAULT_MIXTURE_GAMMA == (0.001, 0.01, 0.1, 1.0)


def test_mixture_dampens_at_noise_z():
    """A variant with z near 0 should LOSE mass after reweighting."""
    rng = np.random.default_rng(1)
    n = 50
    pips = np.full(n, 1.0 / n)
    z = rng.standard_normal(n) * 0.3  # all noise
    z[5] = 5.0  # one mid-strength signal
    new = apply_mixture_posterior(pips, z, n_samples=2400)
    # the noise variants should each have lower PIP than they started with
    noise_idx = [i for i in range(n) if i != 5]
    assert (new[noise_idx] < pips[noise_idx]).sum() > n / 2, (
        f"only {(new[noise_idx] < pips[noise_idx]).sum()} of {n-1} noise "
        f"variants got dampened"
    )


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
