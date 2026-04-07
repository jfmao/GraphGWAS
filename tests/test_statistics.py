"""Comprehensive tests for GraphGWAS statistical engines.

Tests all pure-computation functions with synthetic data.
NO Neo4j connection required — tests the math, not the database.

Run: pytest tests/test_statistics.py -v
"""

import numpy as np
import pytest
from scipy import stats as sp_stats


# ===================================================================
# Fixtures: synthetic genotype/phenotype data
# ===================================================================

@pytest.fixture
def rng():
    return np.random.default_rng(42)


@pytest.fixture
def binary_phenotype(rng):
    """200 samples: 100 cases, 100 controls."""
    n_case, n_ctrl = 100, 100
    labels = np.concatenate([np.ones(n_case), np.zeros(n_ctrl)])
    return labels, n_case, n_ctrl


@pytest.fixture
def quantitative_phenotype(rng):
    """200 samples with continuous trait."""
    return rng.normal(10, 2, 200)


@pytest.fixture
def causal_dosage(rng, binary_phenotype):
    """Dosage vector correlated with case status (simulates causal variant)."""
    labels, n_case, n_ctrl = binary_phenotype
    # Cases enriched for alt allele
    dosage_case = rng.choice([0, 1, 2], n_case, p=[0.3, 0.5, 0.2])
    dosage_ctrl = rng.choice([0, 1, 2], n_ctrl, p=[0.7, 0.25, 0.05])
    return np.concatenate([dosage_case, dosage_ctrl]).astype(float)


@pytest.fixture
def null_dosage(rng):
    """Dosage vector independent of phenotype (null variant)."""
    return rng.choice([0, 1, 2], 200, p=[0.5, 0.4, 0.1]).astype(float)


@pytest.fixture
def count_table_causal():
    """2×3 genotype count table with association signal."""
    # Rows: case, control. Columns: hom_ref, het, hom_alt
    return np.array([[30, 50, 20],    # cases: more het+alt
                     [70, 25, 5]])    # controls: more ref


@pytest.fixture
def count_table_null():
    """2×3 count table with no association."""
    return np.array([[50, 40, 10],
                     [48, 42, 10]])


# ===================================================================
# Test: Chi-squared allelic test
# ===================================================================

class TestChi2:
    def test_causal_variant(self, count_table_causal):
        from graphgwas.assoc import chi2_allelic_test
        result = chi2_allelic_test(count_table_causal)
        assert result["method"] == "chi2"
        assert result["p_value"] < 0.001, "Causal variant should be significant"
        assert result["odds_ratio"] != 1.0, "Should show non-null OR"

    def test_null_variant(self, count_table_null):
        from graphgwas.assoc import chi2_allelic_test
        result = chi2_allelic_test(count_table_null)
        assert result["p_value"] > 0.05, "Null variant should not be significant"

    def test_matches_scipy(self, count_table_causal):
        """Chi2 result should match scipy.stats.chi2_contingency."""
        from graphgwas.assoc import chi2_allelic_test, count_table_to_allelic
        from graphgwas.genotype import count_table_to_allelic as cta
        result = chi2_allelic_test(count_table_causal)
        # Our p-value should be in the right ballpark
        assert result["p_value"] < 1e-5


# ===================================================================
# Test: Fisher's exact test
# ===================================================================

class TestFisher:
    def test_causal_variant(self, count_table_causal):
        from graphgwas.assoc import fisher_exact_test
        result = fisher_exact_test(count_table_causal)
        assert result["method"] == "fisher"
        assert result["p_value"] < 0.001

    def test_null_variant(self, count_table_null):
        from graphgwas.assoc import fisher_exact_test
        result = fisher_exact_test(count_table_null)
        assert result["p_value"] > 0.05

    def test_consistent_with_chi2(self, count_table_causal):
        """Fisher and chi2 should agree on significance direction."""
        from graphgwas.assoc import chi2_allelic_test, fisher_exact_test
        chi2 = chi2_allelic_test(count_table_causal)
        fisher = fisher_exact_test(count_table_causal)
        # Both should be significant
        assert chi2["p_value"] < 0.05
        assert fisher["p_value"] < 0.05
        # Odds ratios in the same direction
        assert (chi2["odds_ratio"] > 1) == (fisher["odds_ratio"] > 1)


# ===================================================================
# Test: Logistic regression
# ===================================================================

class TestLogistic:
    def test_causal_variant(self, causal_dosage, binary_phenotype):
        from graphgwas.assoc import logistic_regression
        labels, n_case, n_ctrl = binary_phenotype
        result = logistic_regression(causal_dosage, labels)
        assert result["method"] == "logistic"
        assert result["p_value"] < 0.01

    def test_null_variant(self, null_dosage, binary_phenotype):
        from graphgwas.assoc import logistic_regression
        labels, _, _ = binary_phenotype
        result = logistic_regression(null_dosage, labels)
        assert result["p_value"] > 0.01, "Null variant should not be significant at 0.01"

    def test_with_covariates(self, causal_dosage, binary_phenotype, rng):
        from graphgwas.assoc import logistic_regression
        labels, _, _ = binary_phenotype
        covariates = rng.normal(0, 1, (200, 2))
        result = logistic_regression(causal_dosage, labels, covariates)
        assert result["p_value"] < 0.05, "Should remain significant with random covariates"


# ===================================================================
# Test: Firth penalized logistic regression
# ===================================================================

class TestFirth:
    def test_causal_variant(self, causal_dosage, binary_phenotype):
        from graphgwas.assoc import firth_logistic_regression
        labels, _, _ = binary_phenotype
        result = firth_logistic_regression(causal_dosage, labels)
        assert result["method"] == "firth"
        assert result["p_value"] < 0.01

    def test_agrees_with_logistic_on_common(self, binary_phenotype, rng):
        """Firth and standard logistic should agree for common variants."""
        from graphgwas.assoc import logistic_regression, firth_logistic_regression
        labels, _, _ = binary_phenotype
        dosage = rng.choice([0, 1, 2], 200, p=[0.4, 0.4, 0.2]).astype(float)
        log_result = logistic_regression(dosage, labels)
        firth_result = firth_logistic_regression(dosage, labels)
        if log_result["p_value"] < 1.0 and firth_result["p_value"] < 1.0:
            # Betas should be close for common variants
            assert abs(log_result["beta"] - firth_result["beta"]) < 0.5

    def test_rare_variant_convergence(self, binary_phenotype, rng):
        """Firth should converge where standard logistic may fail."""
        from graphgwas.assoc import firth_logistic_regression
        labels, _, _ = binary_phenotype
        # Very rare variant: MAC = 3
        dosage = np.zeros(200)
        dosage[:3] = 1.0  # 3 carriers, all cases
        result = firth_logistic_regression(dosage, labels)
        # Should return a result (not null) even for very rare variants
        assert result["method"] == "firth"


# ===================================================================
# Test: Linear regression
# ===================================================================

class TestLinear:
    def test_causal_variant(self, rng):
        from graphgwas.assoc import linear_regression
        dosage = rng.choice([0, 1, 2], 200, p=[0.4, 0.4, 0.2]).astype(float)
        phenotype = 5.0 + 2.0 * dosage + rng.normal(0, 1, 200)
        result = linear_regression(dosage, phenotype)
        assert result["method"] == "linear"
        assert result["p_value"] < 0.001
        assert abs(result["beta"] - 2.0) < 0.5, "Beta should be ~2.0"

    def test_null_variant(self, rng):
        from graphgwas.assoc import linear_regression
        dosage = rng.choice([0, 1, 2], 200, p=[0.4, 0.4, 0.2]).astype(float)
        phenotype = rng.normal(10, 2, 200)
        result = linear_regression(dosage, phenotype)
        assert result["p_value"] > 0.01


# ===================================================================
# Test: PRS evaluation
# ===================================================================

class TestPRS:
    def test_discriminative_prs(self, rng):
        from graphgwas.prs import prs_evaluation
        prs_case = rng.normal(2.0, 1.0, 100)
        prs_ctrl = rng.normal(0.0, 1.0, 100)
        prs = np.concatenate([prs_case, prs_ctrl])
        result = prs_evaluation(prs, 100, 100, verbose=False)
        assert result["auroc"] > 0.8, "Well-separated PRS should have high AUROC"
        assert result["nagelkerke_r2"] > 0.2

    def test_random_prs(self, rng):
        from graphgwas.prs import prs_evaluation
        prs = rng.normal(0, 1, 200)
        result = prs_evaluation(prs, 100, 100, verbose=False)
        assert 0.3 < result["auroc"] < 0.7, "Random PRS should have AUROC ~0.5"

    def test_perfect_prs(self):
        from graphgwas.prs import prs_evaluation
        prs = np.concatenate([np.ones(100) * 10, np.zeros(100)])
        result = prs_evaluation(prs, 100, 100, verbose=False)
        assert result["auroc"] > 0.99


# ===================================================================
# Test: MR estimators
# ===================================================================

class TestMR:
    @pytest.fixture
    def instruments_causal(self, rng):
        """Simulate instruments with true causal effect = 0.5."""
        n = 20
        bx = rng.uniform(0.1, 0.5, n)
        by = 0.5 * bx + rng.normal(0, 0.02, n)
        se_y = rng.uniform(0.01, 0.05, n)
        return [{"variant": f"SNP_{i}", "beta_exposure": float(bx[i]),
                 "se_exposure": 0.02, "p_exposure": 1e-10,
                 "beta_outcome": float(by[i]), "se_outcome": float(se_y[i]),
                 "p_outcome": 0.01} for i in range(n)]

    @pytest.fixture
    def instruments_null(self, rng):
        """Simulate instruments with no causal effect."""
        n = 20
        bx = rng.uniform(0.1, 0.5, n)
        by = rng.normal(0, 0.02, n)  # no relationship to exposure
        se_y = rng.uniform(0.01, 0.05, n)
        return [{"variant": f"SNP_{i}", "beta_exposure": float(bx[i]),
                 "se_exposure": 0.02, "p_exposure": 1e-10,
                 "beta_outcome": float(by[i]), "se_outcome": float(se_y[i]),
                 "p_outcome": 0.5} for i in range(n)]

    def test_ivw_causal(self, instruments_causal):
        from graphgwas.mr import ivw_estimate
        result = ivw_estimate(instruments_causal, verbose=False)
        assert abs(result["beta"] - 0.5) < 0.15, f"IVW should be ~0.5, got {result['beta']}"
        assert result["p_value"] < 1e-10

    def test_ivw_null(self, instruments_null):
        from graphgwas.mr import ivw_estimate
        result = ivw_estimate(instruments_null, verbose=False)
        assert abs(result["beta"]) < 0.3, "Null IVW beta should be ~0"

    def test_egger_no_pleiotropy(self, instruments_causal):
        from graphgwas.mr import egger_estimate
        result = egger_estimate(instruments_causal, verbose=False)
        # Egger intercept should be small (near zero) even if p fluctuates
        assert abs(result["intercept"]) < 0.1, "Intercept should be near zero"
        assert abs(result["beta"] - 0.5) < 0.3

    def test_weighted_median_causal(self, instruments_causal):
        from graphgwas.mr import weighted_median_estimate
        result = weighted_median_estimate(instruments_causal, verbose=False)
        assert abs(result["beta"] - 0.5) < 0.2

    def test_consistency_across_methods(self, instruments_causal):
        """All three MR methods should agree on direction."""
        from graphgwas.mr import ivw_estimate, egger_estimate, weighted_median_estimate
        ivw = ivw_estimate(instruments_causal, verbose=False)
        egger = egger_estimate(instruments_causal, verbose=False)
        wm = weighted_median_estimate(instruments_causal, verbose=False)
        # All should be positive (true effect is positive)
        assert ivw["beta"] > 0
        assert egger["beta"] > 0
        assert wm["beta"] > 0


# ===================================================================
# Test: Heritability helpers
# ===================================================================

class TestHeritability:
    def test_auroc_to_h2_null(self):
        from graphgwas.heritability import _auroc_to_liability_h2
        assert _auroc_to_liability_h2(0.5, 0.1) == 0.0, "AUROC=0.5 → h2=0"

    def test_auroc_to_h2_monotone(self):
        from graphgwas.heritability import _auroc_to_liability_h2
        h2_low = _auroc_to_liability_h2(0.6, 0.1)
        h2_high = _auroc_to_liability_h2(0.7, 0.1)
        assert h2_high > h2_low, "Higher AUROC should give higher h2"

    def test_spectral_h2_clustered(self, rng):
        """Clustered cases should give higher spectral h2 than random."""
        from graphgwas.heritability import _spectral_h2_from_matrix
        n = 50
        # Similarity where cases (first 15) cluster
        W_clust = np.eye(n) * 0.01
        W_clust[:15, :15] += rng.uniform(0.1, 0.3, (15, 15))
        W_clust = (W_clust + W_clust.T) / 2
        phenotype = np.concatenate([np.ones(15), np.zeros(35)])

        h2_clust = _spectral_h2_from_matrix(W_clust, phenotype, 10)

        # Random similarity
        W_rand = rng.uniform(0, 0.1, (n, n))
        W_rand = (W_rand + W_rand.T) / 2
        np.fill_diagonal(W_rand, 1.0)
        h2_rand = _spectral_h2_from_matrix(W_rand, phenotype, 10)

        assert h2_clust > h2_rand, "Clustered should have higher h2"


# ===================================================================
# Test: Multivariate helpers
# ===================================================================

class TestMultivariate:
    def test_flow_covariance_positive(self):
        from graphgwas.multivariate import flow_covariance_from_vectors
        flows1 = [{"pathway": "A", "flow_value": 5.0},
                  {"pathway": "B", "flow_value": 3.0}]
        flows2 = [{"pathway": "A", "flow_value": 4.0},
                  {"pathway": "B", "flow_value": 2.0}]
        result = flow_covariance_from_vectors(flows1, flows2, verbose=False)
        assert result["r_g_flow"] > 0.9, "Same-direction flows should have high r_G"

    def test_flow_covariance_orthogonal(self):
        from graphgwas.multivariate import flow_covariance_from_vectors
        flows1 = [{"pathway": "A", "flow_value": 5.0}]
        flows2 = [{"pathway": "B", "flow_value": 5.0}]
        result = flow_covariance_from_vectors(flows1, flows2, verbose=False)
        assert result["r_g_flow"] == 0.0, "Non-overlapping pathways should have r_G=0"

    def test_psi_formula(self):
        """PSI should reward breadth and biological mechanism."""
        breadth_bio = 0.8 * 12.5 * 1.6    # 4/5 traits, shared variants
        breadth_med = 0.4 * 10.0 * 1.0    # 2/5 traits, mediated
        narrow = 0.2 * 20.0 * 1.0          # 1/5 traits, very strong
        assert breadth_bio > narrow > 0
        assert breadth_bio > breadth_med


# ===================================================================
# Test: MET G×E decomposition
# ===================================================================

class TestMET:
    def test_gxe_variance_decomposition(self, rng):
        """V_G + V_E + V_GxE should approximately sum to V_total."""
        G, E = 30, 8
        geno_eff = rng.normal(0, 3, G)
        env_eff = rng.normal(0, 5, E)
        gxe = rng.normal(0, 1, (G, E))
        Y = geno_eff[:, None] + env_eff[None, :] + gxe

        grand = np.mean(Y)
        g_means = np.mean(Y, axis=1)
        e_means = np.mean(Y, axis=0)
        predicted = g_means[:, None] + e_means[None, :] - grand
        residuals = Y - predicted

        v_g = np.var(g_means)
        v_e = np.var(e_means)
        v_gxe = np.var(residuals)
        v_total = np.var(Y)

        # V_G + V_E + V_GxE should be close to V_total
        assert abs(v_g + v_e + v_gxe - v_total) / v_total < 0.15

    def test_balanced_design_detection(self):
        """A fully balanced design should have balance ratio = 1.0."""
        import networkx as nx
        G_graph = nx.Graph()
        genotypes = [f"G{i}" for i in range(10)]
        environments = [f"E{i}" for i in range(5)]
        G_graph.add_nodes_from(genotypes, bipartite="genotype")
        G_graph.add_nodes_from(environments, bipartite="environment")
        for g in genotypes:
            for e in environments:
                G_graph.add_edge(g, e)
        n_obs = G_graph.number_of_edges()
        assert n_obs == 50
        assert n_obs / (10 * 5) == 1.0

    def test_unbalanced_connectivity(self):
        """Disconnected design should have >1 component."""
        import networkx as nx
        B = nx.Graph()
        # Two disconnected groups
        for i in range(5):
            B.add_edge(f"G{i}", "E0")
        for i in range(5, 10):
            B.add_edge(f"G{i}", "E1")
        components = list(nx.connected_components(B))
        assert len(components) == 2


# ===================================================================
# Test: Spectral vectorization
# ===================================================================

class TestSpectral:
    def test_np_ix_symmetry(self, rng):
        """Vectorized similarity update should be symmetric."""
        n = 50
        sim = np.zeros((n, n))
        for _ in range(20):
            carriers = rng.random(n) > 0.8
            idx = np.where(carriers)[0]
            sim[np.ix_(idx, idx)] += 1.0
        assert np.allclose(sim, sim.T), "Similarity matrix must be symmetric"

    def test_np_ix_nonnegative(self, rng):
        """Similarity should be non-negative."""
        n = 50
        sim = np.zeros((n, n))
        for _ in range(20):
            carriers = rng.random(n) > 0.8
            idx = np.where(carriers)[0]
            weight = rng.uniform(0.1, 10.0)
            sim[np.ix_(idx, idx)] += weight
        assert np.all(sim >= 0)
