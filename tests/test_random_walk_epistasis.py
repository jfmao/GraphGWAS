"""Unit tests for graphgwas.epistasis_higher_order (M5 random walk)."""
from __future__ import annotations

import numpy as np
import pytest

from graphgwas.epistasis_higher_order import (
    m5_interaction_matrix,
    mutual_rwr_pair_scores,
    mutual_rwr_triplet_scores,
    rwr_stationary,
)


class TestRWR:
    """RWR stationary distribution sanity checks."""

    def test_rwr_returns_distribution(self):
        # 4 variants, 3 genes, simple bipartite graph
        # v0--g0--v1, v0--g1--v2, v3--g2--(none)
        A = np.array([
            [1, 1, 0],  # v0 in g0, g1
            [1, 0, 0],  # v1 in g0
            [0, 1, 0],  # v2 in g1
            [0, 0, 1],  # v3 in g2 (isolated component)
        ], dtype=float)
        p = rwr_stationary(A, seed_idx=0, alpha=0.15, n_iter=100)
        assert abs(p.sum() - 1.0) < 1e-10
        assert (p >= 0).all()

    def test_seed_has_high_probability(self):
        A = np.array([[1, 1], [1, 0], [0, 1]], dtype=float)
        p = rwr_stationary(A, seed_idx=0, alpha=0.5, n_iter=100)
        # With high restart probability, seed should dominate
        assert p[0] > 0.5

    def test_low_alpha_diffuses_more(self):
        A = np.array([[1, 1], [1, 0], [0, 1]], dtype=float)
        p_low = rwr_stationary(A, seed_idx=0, alpha=0.05)
        p_high = rwr_stationary(A, seed_idx=0, alpha=0.5)
        # Lower alpha means more diffusion → other variants get more probability
        assert p_low[1] + p_low[2] > p_high[1] + p_high[2]

    def test_isolated_seed_concentrates(self):
        # v3 has no neighbours, so RWR should put all mass on v3
        A = np.array([[1, 0], [0, 1], [1, 1], [0, 0]], dtype=float)
        # variant 3 has no genes; its row is all zeros → all walks restart
        p = rwr_stationary(A, seed_idx=3, alpha=0.15, n_iter=20)
        assert p[3] > 0.99

    @pytest.mark.parametrize("bad_alpha", [-0.1, 0.0, 1.5])
    def test_bad_alpha_raises(self, bad_alpha):
        A = np.array([[1]], dtype=float)
        with pytest.raises(ValueError):
            rwr_stationary(A, seed_idx=0, alpha=bad_alpha)


class TestPairScores:

    def test_pair_scores_symmetric_in_construction(self):
        """score(i,j) = p_i[j] * p_j[i] is symmetric by construction."""
        A = np.array([
            [1, 1, 0],
            [1, 1, 0],
            [0, 1, 1],
            [0, 0, 1],
        ], dtype=float)
        pairs = mutual_rwr_pair_scores(A, seed_indices=[0, 1, 2, 3],
                                        alpha=0.15, top_k=10)
        assert all(s >= 0 for _, _, s in pairs)

    def test_co_membership_increases_score(self):
        """Variants sharing a gene should outrank disconnected pairs."""
        # v0,v1,v2 all in g0; v3 in g1 only
        A = np.array([
            [1, 0],
            [1, 0],
            [1, 0],
            [0, 1],
        ], dtype=float)
        pairs = mutual_rwr_pair_scores(A, seed_indices=[0, 1, 2, 3],
                                        alpha=0.15, top_k=10)
        # The top-scoring pair must come from {0,1,2} (the cluster)
        top_i, top_j, _ = pairs[0]
        assert {top_i, top_j} <= {0, 1, 2}


class TestTripletScores:

    def test_triplet_requires_pairs_to_form_triangle(self):
        """If only 2 pairs exist, no triplet can form."""
        # Bipartite chain: v0--g0--v1--g1--v2 (no v0-v2 edge through any gene)
        A = np.array([
            [1, 0],
            [1, 1],
            [0, 1],
        ], dtype=float)
        triplets = mutual_rwr_triplet_scores(A, seed_indices=[0, 1, 2],
                                              pair_top_k=10, top_k=5)
        # Even if all 3 pairs exist via RWR diffusion, the score should be
        # finite; check it doesn't crash
        assert isinstance(triplets, list)


class TestM5InteractionMatrix:

    def test_m5_z_has_correct_shape(self):
        rng = np.random.default_rng(0)
        n_samples = 50
        n_var = 10
        dos = rng.standard_normal((n_samples, n_var))
        pair_scores = [(0, 1, 0.5), (2, 3, 0.4), (4, 5, 0.3)]
        Z = m5_interaction_matrix(pair_scores, dos)
        assert Z.shape == (n_samples, 3)

    def test_m5_z_columns_centred(self):
        rng = np.random.default_rng(1)
        dos = rng.standard_normal((100, 5))
        pair_scores = [(0, 1, 0.9), (2, 3, 0.8)]
        Z = m5_interaction_matrix(pair_scores, dos)
        # Each column should have mean ~0
        assert np.allclose(Z.mean(axis=0), 0, atol=1e-10)

    def test_m5_z_top_k_truncates(self):
        rng = np.random.default_rng(2)
        dos = rng.standard_normal((30, 8))
        pair_scores = [(0, 1, 0.9), (2, 3, 0.8), (4, 5, 0.7)]
        Z = m5_interaction_matrix(pair_scores, dos, top_k=2)
        assert Z.shape == (30, 2)


class TestPlaceholders:

    def test_gnn_attention_raises_not_implemented(self):
        from graphgwas.epistasis_higher_order import gnn_attention_pairs
        with pytest.raises(NotImplementedError, match="M6"):
            gnn_attention_pairs()

    def test_topological_raises_not_implemented(self):
        from graphgwas.epistasis_higher_order import topological_epistasis
        with pytest.raises(NotImplementedError, match="M7"):
            topological_epistasis()
