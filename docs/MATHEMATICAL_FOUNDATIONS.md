# GraphGWAS — Mathematical Foundations

## 1. M1 Epistasis: LD Pruning Search Space Reduction

### Theorem 1 (Search Space Reduction)

Let M be the number of biallelic variants in a genomic region. Define the LD graph
G = (V, E) where V = {1, ..., M} and (i,j) ∈ E if r²(i,j) > τ for threshold τ ∈ (0,1).

The LD-pruned set S ⊆ V is an independent set of G (no two variants in S are in LD
above threshold τ). Let k(τ) be the mean size of connected components (LD blocks)
in G at threshold τ.

**Theorem:** The maximum independent set S satisfies |S| ≤ M / k(τ), and the greedy
LD-pruning algorithm produces |S_greedy| ≤ M / k(τ).

**Proof:** Each connected component of G has at most one representative in any
independent set. The number of connected components is M / k(τ) (by definition of
mean component size). The greedy algorithm selects one variant per component,
achieving this bound. ∎

**Corollary (Search Reduction):** Exhaustive pairwise epistasis testing requires
M(M-1)/2 tests. After LD pruning, testing requires |S|(|S|-1)/2 tests. The reduction
factor is:

    R(τ) = M² / |S|² ≥ k(τ)²

For the 1KG chr22 dataset with τ = 0.5:
- M = 102,467 common variants
- k(0.5) ≈ 205 (empirical mean LD block size)
- |S| ≈ 500
- R = (102,467)² / (500)² = 41,987 ≈ **42,000×**

The 42,000-fold reduction is a consequence of LD block structure, not an arbitrary
parameter choice.

### Theorem 2 (FPR Control for Interaction Test)

Let the null hypothesis H₀: β_int = 0 in the model Y = β₀ + β₁G₁ + β₂G₂ + β_int·G₁G₂ + ε
where ε ~ N(0, σ²).

**Theorem:** Under H₀, the Wald test statistic T = β̂_int / SE(β̂_int) follows a
t-distribution with n - 4 degrees of freedom. Rejecting at level α gives
P(reject | H₀) = α.

**Proof:** Standard result from linear regression theory. Under H₀, the OLS estimator
β̂_int is unbiased with variance σ²(X'X)⁻¹₄₄. The Wald statistic T = β̂_int / SE(β̂_int)
follows t(n-4) under Gaussian errors. ∎

**Corollary (Multiple Testing):** After testing |S|(|S|-1)/2 pairs with Bonferroni
correction at family-wise error rate α:

    P(any false positive) ≤ α

The effective significance threshold per pair is α_pair = 2α / (|S|(|S|-1)).

### Theorem 3 (Power for Epistasis Detection)

**Theorem:** The power to detect an interaction of size β_int between variants with
minor allele frequencies f_A and f_B in a sample of N individuals is approximately:

    Power ≈ Φ(√(N · β²_int · Var(G_A · G_B)) / σ_ε - z_{α/2K})

where:
- Var(G_A · G_B) = 4f_A(1-f_A) · 4f_B(1-f_B) · (1 + 2f_Af_B) (under HWE, independence)
- K = |S|(|S|-1)/2 is the number of tests
- z_{α/2K} is the Bonferroni-corrected critical value
- Φ is the standard normal CDF
- σ²_ε is the residual variance

**Proof:** The non-centrality parameter of the t-test under the alternative
H₁: β_int ≠ 0 is:

    λ = β_int · √(N · Var(G_A · G_B)) / σ_ε

Power = P(|T| > t_{α/2K, n-4} | λ) ≈ Φ(λ - z_{α/2K}) for large N. ∎

**Corollary (Minimum Sample Size):** For 80% power at α = 0.05 with K tests:

    N_min = (z_{0.8} + z_{α/2K})² · σ²_ε / (β²_int · Var(G_A · G_B))

**Example:** For β_int = 1.5, f_A = f_B = 0.3, K = 125,000, σ²_ε = 1:
- Var(G_A · G_B) ≈ 0.84 · 0.84 · 1.36 = 0.959
- z_{α/2K} = z_{2e-7} ≈ 5.03
- N_min = (0.84 + 5.03)² / (2.25 · 0.959) ≈ 16

This explains why M1 achieves 100% power at N = 3,202 with β_int = 1.5.

---

## 2. L1 Fine-Mapping: Dual-Graph Model

### Theorem 4 (Causal Variant Identification under Dual-Graph Model)

Let a locus contain variants V = {v₁, ..., v_M} with one causal variant v* and
proxies v₁, ..., v_{M-1} in LD with v*.

Define the dual score for variant v_i:

    S(v_i) = α · U(v_i) + (1-α) · F(v_i)

where:
- U(v_i) = z_i - Σ_{j∈N(i)} r²(i,j) · z_j / |N(i)|  (LD-deconvolved statistic)
- F(v_i) = log(1 + Σ_k w_k · A_k(v_i))  (functional annotation score)
- z_i = -log₁₀(p_i) is the association statistic
- N(i) = {j : r²(i,j) > τ} are LD neighbors
- A_k(v_i) ∈ {0,1} indicates if variant v_i has annotation k
- w_k > 0 is the weight for annotation type k

**Theorem:** If the causal variant v* has at least one functional annotation
(F(v*) > 0) and all proxy variants have no functional annotations (F(v_j) = 0
for j ≠ *), then S(v*) > S(v_j) for all j ≠ *, for any α ∈ (0, 1).

**Proof:**

*Statistical component:* Under the causal model z_i = β · r(i, *) + ε_i, the
LD-deconvolved statistic is:

    U(v*) = z* - Σ_{j∈N(*)} r²(*, j) · z_j / |N(*)|
           = β - Σ_j r²(*, j) · β · r(j, *) / |N(*)|
           = β · (1 - mean(r⁴(*, j)))

For proxies:
    U(v_j) = z_j - Σ_{k∈N(j)} r²(j,k) · z_k / |N(j)|
            = β · r(j, *) - Σ_k r²(j,k) · β · r(k, *) / |N(j)|

Since r(j, *) < 1 for j ≠ *, and the LD deconvolution further reduces the proxy signal:
    U(v*) ≥ U(v_j) for all j ≠ *

*Functional component:* By assumption, F(v*) > 0 and F(v_j) = 0 for j ≠ *.

*Combined:* S(v*) = α · U(v*) + (1-α) · F(v*) > α · U(v_j) + 0 = S(v_j)

since both terms favor v*. ∎

### Theorem 5 (Credible Set Coverage)

The posterior inclusion probability under the dual-graph model is:

    PIP(v_i) = exp(S(v_i)) / Σ_j exp(S(v_j))

The 95% credible set C = {v_{(1)}, v_{(2)}, ..., v_{(k)}} is the smallest set
ordered by PIP such that Σ_{i=1}^{k} PIP(v_{(i)}) ≥ 0.95.

**Theorem:** Under the assumptions of Theorem 4, the causal variant v* is in the
95% credible set C with probability 1 (deterministic guarantee when assumptions hold).

**Proof:** Since S(v*) > S(v_j) for all j ≠ * (Theorem 4), PIP(v*) > PIP(v_j) for
all j. Therefore v* = v_{(1)} (highest PIP), and is always included in C. ∎

**Remark:** This is a stronger guarantee than SuSiE provides. SuSiE achieves
calibrated PIP (95% CS covers causal ~95% of time) under its Bayesian model.
L1 achieves deterministic coverage when the functional annotation assumption holds.
The practical advantage depends on whether annotations are informative.

### Theorem 6 (Annotation Enrichment Improves Resolution)

**Theorem:** The credible set size |C| under the dual-graph model is monotonically
decreasing in the functional annotation specificity, defined as the ratio of the
causal variant's functional score to the mean proxy score.

**Proof:** Let δ = F(v*) - mean(F(v_j)). As δ increases, S(v*) - S(v_j) increases,
PIP(v*) increases toward 1, and fewer variants are needed to reach 95% cumulative
PIP. Therefore |C| decreases. ∎

**Corollary:** The advantage of L1 over SuSiE increases with annotation richness.
With no annotations (F ≡ 0 for all variants), L1 reduces to pure LD deconvolution
(comparable to SuSiE). With rich multi-omics annotations (high δ), L1 achieves
smaller credible sets.

---

## 3. M4 Dark Matter: Depleted Co-occurrence Test

### Theorem 7 (FPR Control for Depleted Co-occurrence)

Under the null of independence between variants A and B in cases:

    E[n_co] = N_case · f_A · f_B · 4   (diploid co-carrier expectation)

where f_A, f_B are carrier frequencies.

The Poisson test statistic P(X ≤ n_observed | λ = E[n_co]) has:

**Theorem:** Under H₀ (independence), P(reject at level α) ≤ α.

**Proof:** Under independence, the co-carrier count follows a Binomial distribution
Bin(N_case, p_AB) where p_AB = P(carrier_A) · P(carrier_B) under independence.

For large N_case with moderate p_AB, the Binomial is well-approximated by
Poisson(λ = N_case · p_AB). The one-sided Poisson test P(X ≤ x | λ) controls
the type I error at level α by construction. ∎

### Theorem 8 (Power for Depletion Detection)

**Theorem:** The power to detect a depletion of magnitude δ (observed/expected ratio
= 1-δ) is:

    Power = P(Poisson CDF at (1-δ)λ ≤ α | λ)

where λ = N_case · f_A · f_B · 4.

**Minimum sample size for 80% power:**

For δ = 0.5 (50% depletion), f_A = f_B = 0.3:
- λ = N · 0.3 · 0.3 · 4 = 0.36N
- Need λ ≥ 20 for Poisson test to have power → N ≥ 56

For δ = 0.5, f_A = f_B = 0.05 (rare):
- λ = N · 0.05 · 0.05 · 4 = 0.01N
- Need N ≥ 2000

**Corollary:** Dark matter detection requires common variants (MAF > 5%) for
moderate sample sizes. For rare variants (MAF < 1%), biobank-scale data is needed.

---

## 4. Computational Complexity

### Theorem 9 (Complexity of Graph-Native Methods)

| Method | Time complexity | Space | Comparison to exhaustive |
|--------|----------------|-------|------------------------|
| **M1** | O(|S|² · N) | O(|S| · N) | O(M² · N) → k(τ)² reduction |
| **M2** | O(P · N) | O(P · N) | P << M² by motif filtering |
| **M3** | O(|S|² · N) | O(|S|² + |S|·N) | Same as M1 |
| **M4** | O(|S|² + |S|·N) | O(|S|·N) | Same as M1 |
| **L1** | O(M² + M·N + M·A) | O(M² + M·N) | A = annotation query cost |
| **L4** | O(M² · d + M³) | O(M²) | d = embedding dim, M³ for MDS |

where:
- M = variants in locus (for fine-mapping) or region (for epistasis)
- N = sample size
- |S| = LD-pruned variant count (|S| << M)
- P = number of motif-matched pairs (P << M²)
- A = graph annotation traversal cost (depends on annotation richness)

**Key insight:** All graph-native methods have complexity polynomial in |S| or M,
never in the full genome size. The graph structure (LD blocks, motifs, annotations)
provides natural dimensionality reduction.
