# Graph-Native Epistasis Detection v2 — Design Document

**Related document:** [LD_FINEMAPPING_V2_DESIGN.md](LD_FINEMAPPING_V2_DESIGN.md) — parallel design for LD fine-mapping methods.

## Overview

GraphGWAS v1 detects "epistatic modules" via co-occurrence network + Louvain community detection. **v2 is a fundamental rethink:** seven independent detection methods, each targeting a distinct epistasis signature, benchmarked against ground truth, and exposed as separate user options.

**Core principle:** Don't force one algorithm — expose multiple graph-native detection strategies, let the user select based on their biological question.

**Note on scope:** This document covers epistasis detection. The related LD fine-mapping v2 design covers methods for identifying causal variants within LD blocks. M1 below uses LD pruning as a preprocessing step; the LD fine-mapping methods (L1-L4) are separate and complementary.

---

## Architecture

```
┌─────────────────────────────────────────────────────────────────────┐
│  Input: genotypes + phenotype + annotations (graph)                 │
└─────────────────────────────────────────────────────────────────────┘
                            │
     ┌────────┬────────┬────┴────┬────────┬────────┬────────┐
     ▼        ▼        ▼         ▼        ▼        ▼        ▼
┌────────┐┌────────┐┌────────┐┌────────┐┌────────┐┌────────┐┌────────┐
│ M1:    ││ M2:    ││ M3:    ││ M4:    ││ M5:    ││ M6:    ││ M7:    │
│ Co-    ││ Motif  ││ Diff.  ││ Dark   ││ Random ││ GNN    ││ Topo-  │
│ Module ││ Pairs  ││ Sub-   ││ Matter ││ Walk   ││ Inter- ││ logical│
│ (LD-   ││ (path- ││ graph  ││ (synth.││ (diff- ││ action ││ (persis│
│ pruned)││ way)   ││ (Δedge)││ lethal)││ usion) ││ extract││ homol.)│
└───┬────┘└───┬────┘└───┬────┘└───┬────┘└───┬────┘└───┬────┘└───┬────┘
    └────────┴────────┴────┬────┴────────┴────────┴────────┘
                           ▼
                ┌──────────────────────┐
                │ Unified Ranker       │
                │ (Borda / top-K union)│
                └──────────────────────┘
                           ▼
                ┌──────────────────────┐
                │ Output: ranked pairs │
                │ + p-values + scores  │
                └──────────────────────┘

Pairwise methods: M1, M2, M3, M4 (detect pairs)
Higher-order methods: M5, M6, M7 (detect k-way interactions)
```

---

## Method 1: LD-Aware Co-occurrence with Proper Statistics

**Addresses:** the broken v1 algorithm. Adds LD pruning and interaction testing.

### Algorithm

```
Input: variants V, cases C, controls K, LD graph G_LD, AF_min, r²_max
Output: epistatic modules with p-values

1. LD pruning (greedy set-cover):
   S = []
   rank variants by single-locus p-value (ascending)
   for v in ranked V:
       if no w in S with r²(v, w) > r²_max:
           S.append(v)

2. Build co-occurrence graph on S only:
   for each pair (v_i, v_j) in S with pos_j - pos_i > 1Mb OR different chromosomes:
       n_co_case = |carriers(v_i) ∩ carriers(v_j) ∩ C|
       n_co_ctrl = |carriers(v_i) ∩ carriers(v_j) ∩ K|
       if n_co_case < min_cocarriers: skip
       enrichment = (n_co_case / |C|) / (n_co_ctrl / |K| + pseudo)
       if enrichment > threshold: add edge(v_i, v_j, weight=log(enrichment))

3. Community detection on co-occurrence graph:
   modules = Louvain(co_occurrence_graph)

4. Test each module for interaction (NOT just enrichment):
   for module m:
       for pair (v_i, v_j) in m:
           fit Y ~ β_i·G_i + β_j·G_j + β_ij·(G_i × G_j) + ε
           test H_0: β_ij = 0 via Wald test
           collect p-values per pair
       module_p = Fisher's combined p-value

5. Permutation testing:
   repeat N=1000 times:
       shuffle case/control labels
       recompute module_p under null
   empirical p-value = (# null ≥ observed + 1) / (N + 1)
```

### Key changes from v1

| v1 (broken) | v2 Method 1 |
|------------|-------------|
| No LD pruning | LD pruning r² < 0.5 |
| Tests all pairs in window | Tests distant pairs (>1Mb or different chr) |
| Enrichment score only | Enrichment + interaction term test |
| No p-values | Fisher's combined + permutation |
| Modules 100-700 variants | Modules 2-10 variants (tight clusters) |

### Complexity

- LD pruning: O(M²) worst case, O(M log M) with LD blocks
- Pairwise co-occurrence: O(S²) where S << M after pruning
- Interaction testing: O(modules × pairs_per_module) regressions
- Permutation: O(N × scan_time)

### Output

```python
{
    "module_id": int,
    "variants": [str],  # 2-10 variants
    "enrichment": float,
    "interaction_p_combined": float,
    "permutation_p": float,
    "top_pair": (v1, v2, interaction_beta, interaction_p),
}
```

---

## Method 2: Motif-Filtered Pairwise Testing

**Addresses:** combinatorial explosion and multiple testing burden via biological priors.

### Supported Motifs

```cypher
Motif P1 (Same Pathway):
  MATCH (v1:Variant)-[:HAS_CONSEQUENCE]->(g1:Gene)-[:IN_PATHWAY]->(p:Pathway)
        <-[:IN_PATHWAY]-(g2:Gene)<-[:HAS_CONSEQUENCE]-(v2:Variant)
  WHERE g1.geneId < g2.geneId

Motif P2 (Protein-Protein Interaction):
  MATCH (v1:Variant)-[:HAS_CONSEQUENCE]->(g1:Gene)-[:INTERACTS_WITH]->(g2:Gene)
        <-[:HAS_CONSEQUENCE]-(v2:Variant)

Motif P3 (Regulatory):
  MATCH (v1:Variant)-[:IN_REGULATORY_REGION]->(g:Gene)
        <-[:HAS_CONSEQUENCE]-(v2:Variant)
  WHERE v1.variantId < v2.variantId

Motif P4 (Cis-Trans Regulator):
  MATCH (v1:Variant)-[:IN_REGULATORY_REGION]->(g:Gene),
        (v2:Variant)-[:HAS_CONSEQUENCE]->(tf:Gene)-[:REGULATES]->(g)
```

### Algorithm

```
Input: motif queries, genotypes, phenotype
Output: ranked candidate pairs with interaction p-values

1. Motif enumeration:
   candidate_pairs = []
   for motif in [P1, P2, P3, P4]:
       candidate_pairs.extend(query_motif(motif))
   deduplicate

2. Filter by MAC:
   keep pairs where both variants have MAC ≥ 10 in cases and controls

3. Interaction term test for each candidate:
   for (v1, v2) in candidate_pairs:
       fit: Y ~ β1·G1 + β2·G2 + β12·(G1×G2) + covariates + ε
       return p-value for β12

4. Multiple testing correction:
   Bonferroni over candidate_pairs (typically 10³-10⁵ pairs, not 10¹²)
   OR FDR via Benjamini-Hochberg
```

### Why this is graph-native

Matrix tools can't express "pairs of variants whose genes share a pathway" without pre-computing and storing the pairs externally. The Cypher motif queries do this in one database operation.

### Complexity

- Motif queries: O(edges) — essentially free in Neo4j
- Typical pair counts: 10³-10⁶ vs 10¹² exhaustive
- Testing: O(pairs × N) per regression
- **Testing burden reduction: 10⁶-10⁹ fold**

### Output

```python
{
    "variant_1": str,
    "variant_2": str,
    "motif": str,  # "P1", "P2", etc.
    "shared_entity": str,  # gene, pathway, TF
    "beta_marginal_1": float,
    "beta_marginal_2": float,
    "beta_interaction": float,
    "interaction_p": float,
    "interaction_p_corrected": float,
}
```

---

## Method 3: Differential Subgraph Analysis

**Addresses:** finding interactions that differ in structure (not just magnitude) between cases and controls.

### Conceptual basis

Build two co-occurrence graphs. Compare them as sets of edges.
- **Case-unique edges** (high case co-occurrence, low control): enrichment
- **Control-unique edges** (high control co-occurrence, low case): synthetic lethality
- **Weight-differential edges** (present in both but very different): frequency shift

### Algorithm

```
Input: variants S (LD-pruned), cases C, controls K
Output: case-unique, control-unique, and differential edges

1. Build separate co-occurrence graphs:
   G_case = {(v_i, v_j): case_cocarriers(v_i, v_j) ≥ k_case}
   G_ctrl = {(v_i, v_j): ctrl_cocarriers(v_i, v_j) ≥ k_ctrl}

2. Edge classification:
   case_unique = G_case.edges - G_ctrl.edges
   ctrl_unique = G_ctrl.edges - G_case.edges
   shared = G_case.edges ∩ G_ctrl.edges

3. For each shared edge, compute log2 frequency ratio:
   lfr(e) = log2((cocarriers_case/|C|) / (cocarriers_ctrl/|K|))
   differential = shared where |lfr| > threshold

4. Statistical significance:
   For case_unique edges, compute Fisher's exact test:
     contingency = [[case_co, case_solo], [ctrl_co, ctrl_solo]]
   For ctrl_unique edges, same test with reversed expectation
```

### Biological interpretations

| Edge type | Biological meaning | Example |
|-----------|-------------------|---------|
| case_unique | Positive epistasis / disease pathway combination | Two risk variants co-amplifying disease |
| ctrl_unique | Synthetic lethality / incompatibility | Two variants fatal together, never co-occur in cases |
| differential (lfr < 0) | Purifying selection in cases | Control samples tolerate combo; cases don't |
| differential (lfr > 0) | Enrichment (matches M1) | Both variants needed for disease |

### Complexity

- Graph construction: O(S²) but sparse
- Set operations: O(|E|)
- Fisher tests: O(|edges|)

### Output

```python
{
    "edge": (v1, v2),
    "type": "case_unique" | "ctrl_unique" | "differential",
    "case_cocarriers": int,
    "ctrl_cocarriers": int,
    "log2_frequency_ratio": float,
    "fisher_p": float,
}
```

---

## Method 4: Dark Matter (Negative Co-occurrence / Synthetic Incompatibility)

**Addresses:** epistasis that manifests as ABSENCE of combinations (never tested explicitly in current GWAS).

### Conceptual basis

Under independence, expected co-occurrence of variants A, B in cases:
```
E[n_co_case] = |C| × AF_A × AF_B × 4  (for diploid)
```

If observed co-occurrence << expected, the combination is under negative selection in cases. This detects:
- Synthetic lethality (never co-occur because organism dies)
- Compensatory mutations (wild-type + compensated mutant are fine, double-mutant is not)
- Dobzhansky-Muller incompatibilities

### Algorithm

```
Input: variants (with AF), cases C, controls K
Output: depleted pairs

1. For each pair (v_i, v_j), compute:
   observed_case = count(carriers_i ∩ carriers_j ∩ C)
   expected_case = |C| × 4 × AF_i × AF_j  # random independence
   observed_ctrl = count(carriers_i ∩ carriers_j ∩ K)
   expected_ctrl = |K| × 4 × AF_i × AF_j

2. Depletion score:
   z_case = (observed_case - expected_case) / sqrt(expected_case + 1)
   z_ctrl = (observed_ctrl - expected_ctrl) / sqrt(expected_ctrl + 1)
   depletion_differential = z_ctrl - z_case  # positive = depleted in cases only

3. Require:
   - Both variants have MAF > threshold (to expect some co-occurrences)
   - expected_case > 5 (enough power)
   - observed_case < expected_case * 0.5 (50% depletion or more)

4. Significance test:
   Poisson test: P(X ≤ observed_case | λ = expected_case)
   Adjust for testing burden (FDR)
```

### Why this is novel

No existing tool explicitly tests for **depleted** co-occurrence. Tests like SKAT, burden, PLINK --epistasis all look for enrichment or deviation from additive. This method tests the opposite hypothesis: pairs that should co-occur but don't.

**Biological precedent:** Dobzhansky-Muller incompatibilities are well-documented in evolutionary biology but rarely tested in human GWAS.

### Complexity

- Pair enumeration: O(M²) naive, but prefilter to MAF > 5% reduces to ~0.1M × 0.1M = 10¹⁰ still too many
- **Solution:** use motif prefilter (Method 2) to reduce candidates
- Poisson test per pair: O(1)

### Output

```python
{
    "variant_1": str,
    "variant_2": str,
    "maf_1": float,
    "maf_2": float,
    "expected_co_case": float,
    "observed_co_case": int,
    "expected_co_ctrl": float,
    "observed_co_ctrl": int,
    "depletion_differential": float,
    "poisson_p_case": float,
}
```

---

## Method 5: Random Walk with Restart on Bipartite Graph

**Addresses:** higher-order interactions (3-way, 4-way) via diffusion.

### Conceptual basis

Build a bipartite graph: Samples ↔ Variants (edge = CARRIES). Start random walks from case samples, with restart probability α. Variants reached often = candidates. **Joint hitting probability** of pairs encodes interaction.

### Algorithm

```
Input: bipartite graph B = (Samples, Variants, CARRIES), cases C
Output: ranked pairs by joint hitting probability

1. Initialize probability vector:
   p_0 = uniform over cases C, zero elsewhere

2. Random walk with restart (power iteration):
   repeat until convergence:
       p = (1-α) × T × p + α × p_0
   where T is the transition matrix of B

3. Extract variant hitting probabilities:
   h_v = p at variant v  # probability of being at variant v

4. For each pair (v_i, v_j):
   # Conditional hitting: P(reach v_j starting from carriers of v_i)
   start v_i-carriers subset of samples
   compute reachability to v_j via 2-step walk
   joint_hit(v_i, v_j) = observed 2-step transition probability

5. Compare to null:
   null_joint = E[joint_hit | independence]
           = h_i × h_j / Σ h_k  # approximate
   score(v_i, v_j) = joint_hit(v_i, v_j) / null_joint

6. Rank pairs by score, test significance via label permutation
```

### Multi-step extension (higher-order)

Same framework with k-step walks:
- k=2: pairs
- k=3: triads
- k=4: 4-way interactions

Each step adds O(|V|) computation but reveals higher-order structure.

### Complexity

- Power iteration: O(|E| × iterations) where |E| = total CARRIES edges
- For 3,202 samples × 70M variants: |E| ~ 10¹⁰ edges, iterations ~ 20-50
- **Needs GPU** (sparse matrix operations) to be practical
- Per-pair scoring: O(M² × k) for k-step

### Output

```python
{
    "variants": [v1, v2, ...],
    "walk_length": int,  # k
    "joint_hitting_prob": float,
    "null_expectation": float,
    "score": float,  # observed / null
    "permutation_p": float,
}
```

---

---

## Method 6: GNN Interaction Extraction

**Addresses:** learning interaction patterns from data without specifying candidate pairs a priori.

### Conceptual basis

Represent each sample as a **subgraph** of their carried variants. Train a Graph Neural Network to classify cases vs controls. The GNN learns which *combinations* of variants are predictive. Use explainability methods (GNNExplainer, integrated gradients) to extract the most important subgraphs — these are the epistatic patterns.

GraphGWAS already has `gnn.py` with PyTorch Geometric infrastructure. This method leverages that.

### Algorithm

```
Input: bipartite graph (Samples, Variants, CARRIES), cases C, controls K
Output: important variant subgraphs with scores

1. Build per-sample subgraphs:
   for each sample s in C ∪ K:
       S_s = subgraph induced by variants s carries
       add edges between variants in S_s that co-occur above frequency threshold
       OR: use genotype-weighted edges (sample carries G_i × G_j for all pairs)

2. Feature engineering:
   node features = [variant MAF, impact score, chr position encoded]
   edge features = [r², physical distance, shared gene indicator]

3. Train HeteroGNN (existing GraphGWASModel):
   input: per-sample subgraphs
   target: case/control label
   loss: weighted BCE
   split: 70% train / 15% val / 15% test

4. Extract interactions via explainability:
   for each correctly-classified case:
       compute GNNExplainer mask: which edges/nodes drove the prediction?
   aggregate masks across cases:
       edge_importance(v_i, v_j) = mean(mask[v_i, v_j] across all cases)

5. Rank pairs by edge importance:
   top pairs = highest aggregated importance
   filter: require minimum co-carrier count in cases

6. Validate with interaction term test:
   for top pair (v_i, v_j):
       fit: Y ~ β_i·G_i + β_j·G_j + β_ij·(G_i × G_j)
       test β_ij ≠ 0
```

### Key differences from M1-M4

- **Data-driven, not hypothesis-driven:** Doesn't pre-specify candidate pairs
- **Handles arbitrary interactions:** k-way, non-linear, conditional
- **Captures context:** An interaction may be important only in certain samples
- **Requires training data:** Needs cases/controls with balanced representation

### Complexity

- Subgraph construction: O(|samples| × mean_carriers²)
- GNN training: O(epochs × |edges| × hidden_dim)
- Explainability: O(|samples| × explanation_steps)
- **GPU required** for practical runtime

### When this shines

- Complex architectures where marginal + pairwise models miss the signal
- When biological annotations are sparse (no pathway info) so motif methods fail
- Exploratory analysis before hypothesis testing

### Output

```python
{
    "variant_1": str,
    "variant_2": str,
    "gnn_importance": float,  # 0-1, aggregated over cases
    "n_cases_important": int,  # how many cases relied on this edge
    "marginal_beta_1": float,
    "marginal_beta_2": float,
    "interaction_beta": float,
    "interaction_p": float,
    "model_auroc": float,  # overall model AUROC
}
```

---

## Method 7: Persistent Homology / Topological Epistasis

**Addresses:** higher-order interactions (3-way, 4-way, k-way) via topological features.

### Conceptual basis

The case co-occurrence graph has topological structure beyond just edges:
- **0-simplices:** variants (nodes)
- **1-simplices:** pairs co-occurring above threshold k (edges)
- **2-simplices:** triples co-occurring above threshold k (filled triangles)
- **k-simplices:** (k+1)-variant combinations

**Persistent homology** tracks how these features appear and disappear as k varies from high to low. Features that **persist** across many thresholds are stable, reliable interactions.

### Intuition

- At high k: only very strong co-occurrences exist → sparse graph
- At low k: weak co-occurrences appear → dense graph with many triangles
- A 3-way interaction A×B×C manifests as a **persistent triangle** (A-B, B-C, A-C edges that all appear at similar k and form a triangle)
- Higher simplices represent higher-order interactions

### Algorithm

```
Input: variants S (LD-pruned), cases C, controls K
Output: persistent k-simplices (k-way interactions)

1. Build filtered simplicial complex:
   for each threshold t in [k_max, k_max-1, ..., k_min]:
       include all edges (v_i, v_j) with case_cocarriers >= t
       include all triangles (v_i, v_j, v_k) with triple_cocarriers >= t
       ... (up to dimension d, typically d=3 or 4)

2. Compute persistent homology:
   use Gudhi or Giotto-TDA
   returns: birth/death thresholds for each topological feature

3. Extract persistent features:
   persistence(feature) = birth - death
   filter: persistence > threshold (e.g., spanning 50% of threshold range)

4. Test significance via permutation:
   shuffle case labels 1000 times
   recompute simplicial complex + persistence diagram
   empirical p-value: compare observed persistence to null distribution

5. Convert persistent simplices to interactions:
   persistent 1-simplex (edge) = pairwise interaction
   persistent 2-simplex (triangle) = 3-way interaction
   persistent 3-simplex = 4-way interaction
```

### Why this is uniquely graph-native

- Matrix representations of higher-order interactions require tensors of increasing dimension
- Computational cost scales as O(M^k) for k-way tensor
- Topological features are **intrinsically graph-theoretic** — Betti numbers, homology groups
- Persistence provides natural multi-scale analysis (not just one threshold)

### Complexity

- Simplicial complex construction: O(M^d) where d is max dimension
- **Practical: d=3 (triangles) is tractable, d≥4 requires strict pre-filtering**
- Persistent homology computation: O(|simplices|^3) with standard algorithms
- **Recommendation:** pre-filter with motif (M2) or LD pruning before building complex

### Libraries

- `gudhi` (recommended for this application)
- `giotto-tda` (scikit-learn style API)
- `ripser` (fast Rips complexes)

### Output

```python
{
    "simplex": [str, str, str],  # 2-3 variants for pairs/triples
    "dimension": int,  # 1=pair, 2=triple, 3=quad
    "birth_threshold": int,
    "death_threshold": int,
    "persistence": int,
    "permutation_p": float,
    "interpretation": str,  # "stable 3-way interaction"
}
```

### When this shines

- Discovering novel higher-order interactions (3+way)
- Quantifying interaction stability across thresholds
- Detecting interactions invisible to pairwise methods

---

## Unified Scoring Framework

Each method produces ranked interactions (pairs or higher-order) with its own scoring function. Combine via:

### Option A: Conservative union

Take top K from each method, output the union. Each pair gets flagged with which method(s) detected it.

### Option B: Rank aggregation (Borda count)

```
For each pair p:
    rank_m(p) = rank of p in method m's output (∞ if not in top K)
    borda_score(p) = Σ_m (K - rank_m(p))
```

Pairs detected by multiple methods rank higher.

### Option C: Weighted ensemble

Train weights on ground truth simulations:
```
score(p) = w_1·score_M1(p) + w_2·score_M2(p) + ... + w_5·score_M5(p)
```

Weights learned from synthetic data where we know the answer.

### Option D: Method selection per biological question

| Biological question | Recommended method |
|---------------------|-------------------|
| Any epistatic pair | M1 (co-occurrence with LD pruning) |
| Pairs in known pathway | M2 (motif-filtered) |
| Synthetic lethality | M4 (dark matter) |
| Structural differences cases vs controls | M3 (differential subgraph) |
| 3-way interactions (via diffusion) | M5 (random walk) |
| Complex/unknown patterns (data-driven) | M6 (GNN extraction) |
| k-way interactions, multi-scale | M7 (topological) |

**Recommendation:** expose all 7 methods as CLI options, default to M1+M2 (conservative union for pairwise). For higher-order, default to M7 (if pre-filtered) or M5.

---

## Simulation Framework for Ground Truth

Before running on real data, need phenotype simulations with known epistasis.

### Simulation scenarios

**S1: Pure interaction (no marginals)**
```
Y = β × (G_A - 0.5) × (G_B - 0.5) + ε
```
Marginal effects of A and B alone are zero.

**S2: Additive + interaction**
```
Y = β_A·G_A + β_B·G_B + β_AB·G_A·G_B + ε
```

**S3: Synthetic lethality (for Method 4)**
```
Remove samples where both G_A=1 AND G_B=1
```
Simulates selection against the double mutant.

**S4: Pathway-based (for Method 2)**
```
Pick 5 genes in glycolysis pathway
Pick 2 variants per gene
Y = Σ interactions between variants in different glycolysis genes
```

**S5: Multi-way (for Method 5)**
```
Y = β × G_A × G_B × G_C + ε
```

### Configuration

- N_samples: 1000, 2000, 5000 (use 1KG 3,202)
- N_variants: from real genotypes (yeast 1.9M or human chr22 1M)
- Effect sizes: small, medium, large (OR = 1.5, 3, 10)
- MAF: common (>5%), low-frequency (1-5%), rare (<1%)
- Noise: variance_explained = 5%, 10%, 20%
- LD structure: use real 1KG LD (not simulated)
- Number of causal pairs: 1, 5, 10, 50

### Replicates

10 replicates per configuration × 5 scenarios × 9 configs = 450 simulations.
Tractable: each simulation takes seconds.

---

## Comparison Baselines

### External tools

| Tool | Epistasis method | Strengths | Weaknesses |
|------|------------------|-----------|-----------|
| **PLINK --epistasis** | Exhaustive pairwise χ² | Gold standard, fast | No biological prior, LD blocks |
| **BOOST** | Boolean operations | Very fast, boolean encoding | Only binary traits |
| **MDR** | Classification tree | Captures higher-order | Prone to overfitting |
| **GCTA REML** | Variance decomposition | Estimates epistatic h² | No specific pairs |
| **Random Forest** | Tree-based interactions | Flexible, interpretable | Not statistical test |

### Internal (for ablation)

- v1 (current) — baseline to beat
- Each v2 method alone
- All v2 methods combined

---

## Metrics

### Per-method metrics

**Sensitivity (True Positive Rate):**
```
sensitivity = |detected ∩ true_pairs| / |true_pairs|
```

**False Positive Rate:**
```
FPR = |detected - true_pairs| / |total_tested_pairs|
```

**Precision at top K:**
```
precision@K = |top_K_detected ∩ true_pairs| / K
```

**Area under PR curve (over K):**
```
AUPRC
```

**Runtime per 1M variant pairs**

**Module size distribution** (should be 2-10, not 100+)

### Meta-metrics (across methods)

**Method uniqueness:** What fraction of detected pairs are unique to each method?

**Consensus reliability:** If multiple methods detect the same pair, how often is it a true positive?

**Coverage:** Does any combination of methods achieve >90% sensitivity?

---

## Implementation Plan

### New module: `src/python/graphgwas/epistasis_v2.py`

```python
# Main entry points
def ld_aware_cooccurrence(conn, phenotype, af_min=0.05, r2_max=0.5,
                           min_cocarriers=5, n_permutations=1000)
def motif_filtered_pairs(conn, phenotype, motifs=["P1","P2","P3"],
                          mac_min=10, correction="BH")
def differential_subgraph(conn, phenotype, min_cocarriers=5, lfr_threshold=2.0)
def dark_matter(conn, phenotype, maf_min=0.05, depletion_threshold=0.5)
def random_walk_pairs(conn, phenotype, k_steps=2, alpha=0.15, n_iterations=50)

# Unified interface
def detect_epistasis(conn, phenotype, method="auto", **kwargs)

# Simulation
def simulate_epistatic_phenotype(conn, scenario, config)

# Benchmarking
def benchmark_all_methods(conn, simulation, ground_truth)
```

### New module: `src/python/graphgwas/epistasis_simulate.py`

```python
def simulate_pure_interaction(genotypes, causal_pairs, beta, noise)
def simulate_additive_plus_interaction(...)
def simulate_synthetic_lethality(...)
def simulate_pathway_epistasis(...)
def simulate_multiway(...)
```

### New CLI commands

```bash
# Run specific method
graphgwas epistasis detect --method motif --trait ethanol -o pairs.tsv

# Run all methods and compare
graphgwas epistasis detect --method all --trait ethanol -o combined.tsv

# Benchmark against ground truth
graphgwas epistasis benchmark --simulation S1 --n-reps 10 -o benchmark.tsv
```

### Dependencies

- Reuse: `genotype.py` (dosage, carrier sets), `popstruct.py` (LD calculation)
- New: interaction term regression (statsmodels Logit/OLS with interaction)
- New: Bipartite random walk implementation (scipy sparse)

---

## Runtime Estimates

Based on yeast (1.9M variants, 1011 samples) × human (70M × 3,202) extrapolation:

| Method | Yeast | Human chr22 (1M) | Human full genome |
|--------|-------|-----------------|-------------------|
| M1 Co-occurrence (LD-pruned) | ~5 min | ~30 min | ~10 hr |
| M2 Motif-filtered pairwise | ~2 min | ~10 min | ~1 hr |
| M3 Differential subgraph | ~5 min | ~30 min | ~10 hr |
| M4 Dark matter (with motif prefilter) | ~3 min | ~15 min | ~2 hr |
| M5 Random walk (GPU needed) | ~10 min | ~1 hr | ~24 hr (without GPU) |
| M6 GNN extraction (GPU required) | ~15 min | ~2 hr | ~8 hr |
| M7 Topological (dim 2, pre-filtered) | ~10 min | ~1 hr | ~4 hr |

**M2 is fastest and most actionable** — the biological prefilter does most of the work.
**M6 and M7 require GPU or pre-filtering to be practical on full genome.**

---

## Success Criteria

### Minimum success (proof of concept)

- [ ] At least one method achieves >80% sensitivity on pure-interaction simulation (S1)
- [ ] At least one method achieves >60% precision at top-10
- [ ] All methods have calibrated FPR (<5% on null simulations)
- [ ] At least one method detects synthetic lethality (S3) that PLINK misses

### Strong success (publishable)

- [ ] Union of methods achieves >90% sensitivity
- [ ] Motif-filtered (M2) detects pairs PLINK misses due to MTC burden
- [ ] Dark matter (M4) detects validated synthetic lethal pairs
- [ ] Random walk (M5) detects 3-way interactions in S5 simulation
- [ ] Combined method ranks causal pairs in top-10 with >70% accuracy

### Stretch success (paradigm-shifting)

- [ ] Discover novel epistatic interactions in real trait (yeast or human)
- [ ] Reproduce findings across multiple traits
- [ ] Biological validation of novel interactions (literature / experiments)
- [ ] M7 detects 3-way interaction in S5 simulation (topological advantage)
- [ ] M6 (GNN) finds interactions that M1-M4 miss (learned patterns)

---

## Open Research Questions

1. **What LD pruning threshold is optimal?** r²=0.2 (strict) vs 0.5 (lenient) vs 0.8 (loose). Too strict loses power; too loose confounds with LD.

2. **How to handle rare variants in epistasis?** If both variants are AF<1%, you might have 0 co-carriers purely by chance. Need aggregation (burden × burden) or different null.

3. **What's the proper null distribution for differential subgraph (M3)?** Label permutation is standard, but the edges themselves are computed from case/control split — need case-control resampling.

4. **Random walk convergence on 70M-node graphs?** Needs GPU or sparse approximations. May need subsampling samples.

5. **How to correct for multiple testing across methods?** If we run 5 methods, each with its own FDR, the overall FDR is inflated. Need a meta-procedure.

6. **Should we condition on marginal effects?** Test interaction term AFTER conditioning on both variants' marginal effects (more conservative but correct).

---

## Proposed Timeline

**Week 1: Simulation framework**
- Implement S1-S5 scenarios
- Validate simulations produce detectable signal

**Week 2: Methods M1, M2 implementation + benchmarks**
- LD-aware co-occurrence with interaction test
- Motif-filtered pairwise
- Run against S1-S2 simulations

**Week 3: Methods M3, M4 implementation + benchmarks**
- Differential subgraph analysis
- Dark matter / synthetic lethality
- Run against S3 simulation

**Week 4: Method M5 (random walk) — research mode**
- Bipartite random walk implementation
- Higher-order detection
- GPU acceleration if needed
- Run against S5 simulation

**Week 5: Unified framework + comparison**
- Implement Borda rank aggregation
- External baseline comparisons (PLINK, BOOST)
- Full ablation study

**Week 6: Real data application**
- Apply to yeast (validated traits)
- Apply to human chr22 (simulated + real traits)
- Write up results

---

## Risks and Mitigations

| Risk | Likelihood | Mitigation |
|------|-----------|------------|
| No method beats PLINK --epistasis | Medium | Still valuable: 5 approaches, biological interpretation |
| Motif-filtered misses interactions between unannotated genes | Low | M1 (unbiased) complements M2 |
| Random walk doesn't converge on full genome | High | Fall back to subsampled walks or pathway-local walks |
| Differential subgraph has too many false positives | Medium | Strict thresholds + permutation testing |
| Dark matter: no truly depleted pairs exist at current sample size | Medium | Report even if negative — informative about limits |
| Overall: all methods find the same pairs (redundancy) | Low | That would actually be validation that they work |

---

## Final Deliverable

A paper-ready benchmark showing:
1. **Current v1 is broken** (module sizes confirm LD confounding)
2. **Each v2 method has distinct detection signature** (captures different biology)
3. **Union of methods outperforms any single method** (synergy)
4. **Method 2 (motif-filtered) outperforms PLINK** at same FPR (biological prior helps)
5. **Method 4 (dark matter) finds synthetic lethality PLINK can't** (novel signal)

This gives users a principled choice: "I want to find pathway epistasis → M2. I want synthetic lethality → M4. I want everything → combined."

---

## Ready for Implementation?

**Before coding:**
- [ ] Review this design with user
- [ ] Decide which methods to implement first (I recommend M1 + M2 + ground truth simulator)
- [ ] Decide on LD pruning threshold default
- [ ] Decide on multiple testing correction strategy
- [ ] Decide whether to implement M5 (random walk) initially or defer

**Initial implementation scope proposal:**
- Phase 1 (2 weeks): Simulator + M1 + M2 + basic benchmarks
- Phase 2 (2 weeks): M3 + M4 + PLINK comparison
- Phase 3 (2 weeks): M5 + unified framework
- Phase 4 (3 weeks): M6 (GNN) + M7 (Topological) — higher-order methods
- Phase 5 (1 week): Paper-ready results + cross-method comparison
