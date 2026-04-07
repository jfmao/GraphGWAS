This is exactly the right question to ask before implementation begins, and the honest answer is: **the statistical engine of GraphGWAS has not been fully specified yet**. Let me think through this rigorously.Good. The project knowledge confirms the gap — statistical tests are discussed at a high level but the actual computation inside each method has not been pinned down. Let me now give you the clearest possible answer to this question, structured from first principles.

---

## The Honest Starting Point

When you ask "what statistical test or algorithm actually does the GWAS?" the truthful answer has three layers that you must distinguish carefully:

1. **What classical GWAS uses** — well-established, mathematically precise
2. **What GraphGWAS's Layer 2 (single-locus) uses** — essentially the same statistics, but *sourced differently* from the graph
3. **What GraphGWAS's Layers 3–5 use** — genuinely new computations that have no direct classical equivalent

Conflating these three leads to confusion. Let me work through each.

---

## Layer 2: Single-Locus Association — The Statistics Are Classical, the Data Source Is Not

The single-locus test in GraphGWAS uses exactly the same statistical machinery as PLINK2 or SAIGE. The graph does not change the mathematics of association testing at this level. What it changes is **how the input data is assembled** and **what covariates are naturally available**.

### For Binary Traits (Case-Control)

**The standard test:** Logistic regression.

For each variant *v*, fit the model:

```
logit(P(case)) = β₀ + β₁·G_v + β₂·PC1 + β₃·PC2 + ... + βₖ·covariates
```

where `G_v` is the genotype dosage (0, 1, or 2 alt alleles). The test statistic is the Wald test on `β₁`, asymptotically chi-squared with 1 df. The p-value tests H₀: β₁ = 0, i.e., no association.

**The Firth correction** modifies the likelihood to reduce bias when:
- Minor allele count is small (MAC < 20 in cases)
- Case-control ratio is severely unbalanced
- Any cell in the 2×2 genotype × phenotype table is zero or near-zero

The penalized log-likelihood becomes:
```
L*_F(β) = L(β) + ½ log|I(β)|
```
where `I(β)` is the Fisher information matrix. This adds a data-dependent prior that prevents infinite effect size estimates when a minor allele happens to appear only in cases.

**What the graph contributes here:** Instead of reading dosage vectors from a flat matrix, GraphGWAS queries:
```cypher
MATCH (s:Sample)-[:CARRIES {gt: g}]->(v:Variant {id: $variant_id})
RETURN s.id, g, s.phenotypes.trait, s.pca_coords
```
This single traversal simultaneously retrieves the genotype, the phenotype, AND the PC covariates — all natively. In PLINK2, these come from three separate files (`.bed`, `.pheno`, `.eigenvec`) that must be joined. The math is identical; the data plumbing is graph-native.

### For Quantitative Traits

**The standard test:** Linear regression (ordinary least squares).

```
Y = β₀ + β₁·G_v + β₂·PC1 + ... + ε,   ε ~ N(0, σ²)
```

Test statistic: t-test on β₁, or equivalently F-test. Assumes normality of residuals — in practice, inverse-normal transformation of the phenotype is standard.

### The Population Structure Problem and Its Graph-Native Solution

Both logistic and linear regression above assume that after conditioning on covariates, residual confounding by ancestry is negligible. This is why PCA covariates (PC1–PC10) are included. The PCs are computed from the genotype matrix and represent directions of maximal ancestry variation.

In GraphGWAS, the PC coordinates come from `graphpop.structure.pca()` — graph-native GRM construction from CARRIES edges, then eigendecomposition. The math (LAPACK eigendecomposition) is classical. The data source (graph-computed GRM) is native.

**But GraphGWAS can go further.** Beyond classical PCA covariates, it can use:
- **Community membership** (Leiden/Louvain labels) from the Sample–Sample kinship graph as a discrete ancestry covariate
- **Node2vec or GraphSAGE embeddings** of the kinship graph as a continuous ancestry representation richer than PCA

This is where the graph starts to genuinely augment the statistics — not replace them, but enrich the covariate structure.

### Linear Mixed Model (LMM) — The SAIGE/BOLT-LMM Approach

For large biobank-scale data with population structure and relatedness, the current state-of-the-art is not plain logistic regression but a **generalized linear mixed model** (GLMM):

```
logit(P(case)) = Xβ + u,    u ~ MVN(0, σ²ₐ·K)
```

where `K` is the genetic relationship matrix (GRM) — a kinship matrix that simultaneously accounts for population structure AND cryptic relatedness between samples. SAIGE and BOLT-LMM are the dominant implementations.

**In GraphGWAS**, the GRM `K` can be computed from KINSHIP edges natively:
```cypher
MATCH (s1:Sample)-[k:KINSHIP]->(s2:Sample)
RETURN s1.id, s2.id, k.coefficient
```
This sparse kinship matrix IS the GRM. The GLMM fitting (MCMC or variational approximation) is expensive but mathematically identical to SAIGE's approach — GraphGWAS just assembles the kinship matrix from the graph rather than recomputing it from scratch.

This is important to be clear about: **GraphGWAS does not replace LMM with something else for single-locus tests.** It makes LMM *easier to set up* because the kinship matrix already exists as graph edges.

---

## Layer 3: Multi-Locus Methods — This Is Where the Statistics Genuinely Diverge

Here GraphGWAS departs from classical methods. Let me be precise about each one.

### 3A. Variant Co-Occurrence Enrichment — Fisher's Exact Test

The co-occurrence test asks: is variant pair (v1, v2) carried together more often in cases than expected by chance?

The statistical test is **Fisher's exact test** on a 2×2 contingency table:

```
                  Carries both v1+v2    Does not carry both
Cases:                    a                     b
Controls:                 c                     d
```

The odds ratio is `(a·d)/(b·c)`. Fisher's exact test gives the exact p-value under the hypergeometric null. This is not a new statistic — it is the classical epistasis detection test, just now computable at scale because the graph makes it fast to retrieve co-carrier counts.

For multiple testing correction across all pairs: Benjamini-Hochberg FDR, or permutation-based FDR (shuffle case/control labels on Sample nodes, recompute all pair enrichments, build empirical null).

**The graph's contribution:** In classical GWAS, testing all M×M pairs is O(M²) and computationally intractable for M = 10M variants. The graph imposes a natural sparsity constraint — you only build CO_OCCURS edges where co-carriage count exceeds a threshold — reducing the search space from M² to a manageable sparse graph, without changing the statistical test itself.

### 3B. Burden Test — Logistic Regression on Aggregate Score

BURDEN aggregates rare variants within a gene into a single score per sample:

```
burden_score(sample, gene) = Σ w(v) · G(sample, v)    for all v in gene with AF < threshold
```

where `w(v)` is a variant weight (e.g., `w = 1/AF` to upweight ultra-rare variants, or CADD score). Then test:

```
logit(P(case)) = β₀ + β₁·burden_score + covariates
```

This is again logistic regression — but now the predictor is a graph traversal aggregate rather than a single genotype. In Cypher:

```cypher
MATCH (s:Sample)-[:CARRIES]->(v:Variant)-[:HAS_CONSEQUENCE]->(g:Gene {symbol: $gene})
WHERE v.af_total < $threshold
RETURN s.id, sum(1.0 / v.af_total) AS burden_score, s.is_case
```

That `sum(1.0 / v.af_total)` — computed in a single graph traversal — is the burden score that goes into the logistic regression. Classical tools (SKAT-O) require parsing VCF records, subsetting by gene coordinates, building a separate burden vector. GraphGWAS computes it as a Cypher aggregate.

**SKAT (Sequence Kernel Association Test)** is the variance-component alternative to BURDEN. Instead of the sum of weighted genotypes, it tests whether the variance of the distribution of rare variant effects is non-zero:

```
Q_SKAT = (y - μ)ᵀ · W · Kₛ · W · (y - μ)
```

where `Kₛ` is the genetic kernel matrix (weighted genotype covariance within the gene) and `μ` is the fitted phenotype mean under the null. The test statistic Q follows a mixture of chi-squared distributions. SKAT-O combines BURDEN and SKAT optimally. This is genuinely complex to implement but mathematically well-specified. GraphGWAS's contribution is assembling the kernel matrix `Kₛ` from CARRIES traversal rather than from a VCF slice.

### 3C. Phenotype Label Diffusion — New Statistical Territory

This is where GraphGWAS enters genuinely novel statistical ground, and it must be stated honestly: **label propagation as a GWAS method has not been formally established in the GWAS literature.** It is conceptually motivated and computationally attractive, but the null distribution and the correct statistical test for its outputs are not yet specified.

The algorithm:

```
Initialize: score(s) = 1.0 for cases, 0.0 for controls
Iterate:
  score(v) ← Σ [score(s) · w(s→v)] for all samples s carrying v
  score(g) ← Σ [score(v) · w(v→g)] for all variants v in gene g  
  score(p) ← Σ [score(g) · w(g→p)] for all genes g in pathway p
  score(s) ← α · initial_score(s) + (1-α) · Σ [score(v) · w(v→s)]
Until convergence
```

The diffusion score on a Pathway node tells you how much phenotype signal flowed through that pathway. But what is its null distribution? The correct answer is: **permutation testing**. Shuffle the `is_case` labels on Sample nodes 1000 times, rerun diffusion each time, and build the empirical distribution of pathway scores under the null. The observed score's position in that null distribution gives the permutation p-value.

This is computationally expensive (1000 × diffusion iterations) but statistically sound. It is the appropriate approach when the null distribution is analytically intractable — which it is here, because the graph topology introduces complex correlations between variant scores.

**The key statistical claim** that needs rigorous validation: does phenotype diffusion have power to detect pathway-level enrichment when no individual variant is significant? This must be demonstrated via simulation before the method can be published. The design is: simulate a cohort where 20 rare variants in one pathway each cause 5% of cases to carry them — no individual variant is significant at 5e-8, but the pathway accumulates 20× the expected diffusion signal. Power analysis against permuted null.

### 3D. The Disease Subgraph — A Graph Optimization Problem, Not a Traditional Test

The disease subgraph is not a statistical test in the classical sense. It is a **constrained graph optimization** — find the minimum-weight connected subgraph (spanning variant, gene, and pathway nodes) that maximally separates case-carrier topology from control-carrier topology.

The closest classical analogue is **network-based enrichment analysis** (e.g., HotNet2 in cancer genomics, which finds connected subgraphs of mutated genes). The statistical validation approach is the same: permute sample labels, extract subgraphs under the null, compare the observed subgraph's enrichment score to the null distribution.

**In GraphGWAS**, the enrichment score per node:
```
enrichment(v) = (case_carrier_count / n_cases) / (control_carrier_count / n_controls + ε)
```

Aggregate to gene: `enrichment(g) = mean(enrichment(v) for v in gene g)`
Aggregate to pathway: `enrichment(p) = sum(enrichment(g) for g in pathway p) / sqrt(|pathway|)`

The sqrt(|pathway|) term penalizes large pathways that would accumulate signal trivially. This size correction is standard in gene-set enrichment analysis (GSEA).

---

## Layer 4: GNN — A Machine Learning Model, Not a Hypothesis Test

The GNN (Heterogeneous Graph Transformer or R-GCN) operates in a completely different statistical paradigm — it is **supervised machine learning**, not hypothesis testing. This distinction is critical.

**The GNN does not produce p-values.** It produces:
- A phenotype prediction score per sample (probability of being a case)
- Feature attribution scores per variant/edge (via GNNExplainer or SHAP)

**What the GNN tests:** whether the graph topology — the joint structure of which samples carry which variants, connected through LD and functional annotation edges — contains signal predictive of phenotype, beyond what any individual variant contributes.

**Evaluation metric:** AUROC (area under ROC curve) on held-out samples. Compare GNN AUROC to logistic regression AUROC on same data. If GNN > logistic regression, the multi-locus topology contains signal the single-locus model missed.

**The causal attribution problem:** High SHAP values on a variant from the GNN do not constitute a statistical association in the frequentist sense. They are feature importance scores from a model that may be capturing LD or population structure confounders. This is a real limitation — GNN outputs require post-hoc validation with classical tests, not replacement of them.

**What GNN adds to GraphGWAS specifically:** it can detect non-linear, high-order interactions that would require exponentially many classical tests to find. If three variants jointly predict disease (a three-way interaction) but no pair does, a 3-layer GNN can learn this. Classical GWAS testing all triplets is O(M³) — infeasible. The GNN encodes it implicitly through message-passing.

---

## The Complete Picture: What GraphGWAS Actually Computes

Let me make this fully concrete as a table:

| GraphGWAS Method | Statistical Engine | Test / Criterion | p-value Source | Novel vs Classical |
|---|---|---|---|---|
| Single-locus (common variants) | Logistic/linear regression | Wald test on β₁ | Chi-squared, 1df | Classical, graph-sourced inputs |
| Single-locus (rare variants) | Firth penalized logistic | Penalized likelihood ratio | Chi-squared, 1df | Classical, graph-sourced inputs |
| Single-locus at biobank scale | Generalized LMM (SAIGE-like) | Score test with saddle-point approx | Empirical | Classical, kinship from graph |
| Gene burden test | Weighted logistic regression | Wald/likelihood ratio test on burden β | Chi-squared | Classical aggregation, graph traversal |
| SKAT/SKAT-O | Variance component test | Mixture of chi-squared (Davies method) | Analytically approximated | Classical, kernel from graph |
| Pairwise epistasis | Fisher's exact test | Hypergeometric | Exact | Classical test, graph-enabled scale |
| Higher-order epistatic modules | Community detection enrichment | Permutation FDR | Empirical permutation | Novel — no classical equivalent |
| Phenotype label diffusion | Iterative graph diffusion | Permutation test on pathway score | Empirical permutation | **Novel — new to GWAS** |
| Disease subgraph extraction | Subgraph enrichment score | Permutation test on subgraph score | Empirical permutation | **Novel — adapted from cancer genomics** |
| GNN phenotype prediction | Message-passing neural network | AUROC on held-out set | None (ML evaluation) | **Novel — no classical equivalent** |
| GNN feature attribution | SHAP / GNNExplainer | Feature importance score | None (ML evaluation) | **Novel — requires classical validation** |

---

## The Layered Logic: Why You Need All of Them

Here is the scientific logic of why GraphGWAS needs all these methods, not just one:

**Layer 2 (classical logistic/Firth)** gives you p-values interpretable by the entire GWAS community. It is your baseline, your credibility anchor, and your validation against PLINK2. You cannot publish without it.

**Layer 3 burden/SKAT** gives you power for rare variants that Layer 2 misses. These are established methods with a publication track record. They extend your coverage without requiring new statistical theory.

**Layer 3 co-occurrence + Fisher's exact** gives you pairwise epistasis. The test itself is classical; the graph makes it scalable. This is incrementally novel — you can describe it as "classical epistasis testing, made scalable by graph topology."

**Layer 3 diffusion and disease subgraph** are genuinely novel methods. They need rigorous simulation-based validation before publication. The burden of proof is on you to show: (a) the permutation-based p-values are calibrated, (b) the method has power above the classical alternatives, (c) the results are interpretable biologically.

**Layer 4 GNN** is an exploratory, hypothesis-generating tool — not a confirmatory test. Its outputs should be validated by Layer 2 or Layer 3 tests. Treat it as "the GNN highlighted this pathway; now Layer 3 diffusion confirms the signal at p = 0.002 by permutation."

The synthesis: **the graph is the data structure; the statistics are a layered stack from proven to novel, with the novel methods requiring simulation validation before they can be claimed as contributions.**