# Graph-Native LD Fine-Mapping v2 — Design Document

**Related document:** [EPISTASIS_V2_DESIGN.md](EPISTASIS_V2_DESIGN.md) — parallel design for epistasis detection methods.

## Overview

LD fine-mapping is the problem of identifying the causal variant(s) within a locus where many correlated variants all show association. Current GraphGWAS has a simple centrality-based fine-mapping method (`finemapping.py`). **v2 proposes four graph-native methods**, each exploiting different aspects of the genomic graph structure beyond pure LD.

**The fundamental insight:** LD is not noise to be filtered — it's *structure* that encodes evolutionary history, haplotype inheritance, and recombination architecture. Combined with functional annotations (gene/pathway edges), graph representations can resolve causality in ways that pure LD-based Bayesian methods cannot.

**Core principle:** Each method combines LD structure with a *different* type of auxiliary information (functional, regulatory, evolutionary) to distinguish causal from proxy variants.

---

## Background: The LD Fine-Mapping Problem

After GWAS identifies a significant locus, typically 20-200 variants show p-values < 5e-8 due to LD. Only 1-5 are causal. Fine-mapping aims to rank them.

### Current approaches

| Tool | Principle | Strengths | Weaknesses |
|------|-----------|-----------|-----------|
| **SuSiE** | Sum of single effects, Bayesian | Handles multiple causal | Minutes-hours per locus |
| **FINEMAP** | Shotgun stochastic search | Multiple causal, calibrated | Slow, needs LD reference |
| **PAINTOR** | Integrates functional annotations | Uses external annotations | Complex setup |
| **fastENLOC** | Colocalization with eQTL | Combines GWAS+eQTL | Requires eQTL data |
| **GraphGWAS v1 (centrality)** | Betweenness centrality on LD graph | Seconds per locus | Pure topology, no biology |

### The graph-native advantage

Graph databases naturally represent:
- LD relationships (Variant-Variant edges weighted by r²)
- Functional annotations (Variant→Gene→Pathway edges)
- Regulatory architecture (Variant→RegulatoryRegion→Gene edges)
- Tissue-specific effects (Variant→Gene→Tissue edges)
- Protein interactions (Gene-Gene interaction edges)

**No matrix tool can integrate these in a single, unified algorithm.** All existing multi-omic fine-mappers (PAINTOR, fastENLOC) pre-compute features and feed them to a separate statistical model. Graph-native methods can traverse the integrated structure directly.

---

## Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│  Input: locus (lead SNP + window), GWAS summary stats          │
└─────────────────────────────────────────────────────────────────┘
                                │
        ┌───────────────────────┼────────────────────────┬────────────────────┐
        ▼                       ▼                        ▼                    ▼
┌─────────────────┐  ┌─────────────────┐  ┌─────────────────┐  ┌─────────────────┐
│ L1: Dual-Graph  │  │ L2: LD-Regul-   │  │ L3: Haplotype   │  │ L4: Recombin-   │
│ Fine-Mapping    │  │ arized Assoc.   │  │ Graph Assoc.    │  │ ation-Aware     │
│ (LD + func)     │  │ (LD as negative │  │ (haplotype-level│  │ Embedding       │
│                 │  │  edges)         │  │  statistics)    │  │                 │
└────────┬────────┘  └────────┬────────┘  └────────┬────────┘  └────────┬────────┘
         │                    │                    │                    │
         └────────────────────┴────────┬───────────┴────────────────────┘
                                       ▼
                            ┌──────────────────────┐
                            │ Unified Credible Set │
                            │ (consensus ranking)  │
                            └──────────────────────┘
                                       ▼
                            ┌──────────────────────┐
                            │ Output: ranked       │
                            │ causal candidates    │
                            │ + posterior probs    │
                            └──────────────────────┘
```

---

## Method L1: Dual-Graph Fine-Mapping (LD + Functional)

**Addresses:** distinguishing causal from proxy variants by integrating LD with functional annotations in a single graph computation.

### Conceptual basis

A causal variant should exhibit BOTH:
- Strong statistical signal (low p-value from GWAS)
- Functional relevance (traverses to a gene/pathway consistent with the phenotype)

A proxy variant typically has:
- Strong statistical signal (by LD)
- NO direct functional relevance (just happens to be nearby on the chromosome)

**The dual-graph score combines both signals.** Variants with functional context get their score amplified; variants with only LD-driven signal get suppressed via smoothing.

### Algorithm

```
Input: locus variants V with GWAS p-values, LD graph G_LD, annotation graph G_annot
Output: ranked causal candidates with posterior probabilities

1. Initialize statistical score:
   for v in V:
       z_stat(v) = -log10(p_value(v))  # or chi² statistic

2. Compute functional score via graph traversal:
   for v in V:
       functional_hits = []
       # traverse Variant→Gene→Pathway edges
       for path in graph.traverse(v, max_depth=3):
           if path ends at a pathway/gene relevant to phenotype:
               functional_hits.append(path.score)  # weighted by edge types
       z_func(v) = log(1 + sum(functional_hits))

3. Apply LD smoothing (variants in LD should have similar scores):
   for v in V:
       neighbors = LD_graph.neighbors(v)  # r² > 0.3
       ld_weights = [r²(v, u) for u in neighbors]
       smoothed_stat(v) = weighted_mean(z_stat(u) for u in [v] + neighbors, weights)

4. Apply LD differentiation (variants should be distinguished):
   # compute variant's unique contribution after removing LD neighbors
   for v in V:
       unique_score(v) = z_stat(v) - sum(r²(v,u) * z_stat(u) for u in neighbors)

5. Combined score:
   final_score(v) = α * unique_score(v) + (1-α) * z_func(v)
   # α balances statistical vs functional evidence

6. Convert to posterior probability (PIP-like):
   PIP(v) = exp(final_score(v)) / sum(exp(final_score(u)) for u in locus)

7. Credible set:
   sort variants by PIP, take cumulative sum to 95%
```

### How to define "phenotype-relevant" pathways/genes

Three strategies:
- **User-specified:** User provides a list of candidate genes or pathways
- **GWAS-derived:** Use genes harboring other GW-significant loci as "relevant"
- **Pathway database:** Use KEGG/GO terms enriched in GWAS hits (MAGMA-style)

### Why this beats matrix-based multi-omic fine-mappers

| Aspect | PAINTOR/fastENLOC | L1 dual-graph |
|--------|-------------------|---------------|
| Annotation integration | Pre-computed features | Direct graph traversal |
| Extensibility | Requires new model per annotation type | Just add new edge types |
| Queryability | None (black box model) | Cypher queries explain rankings |
| Performance on new data types | Retrain model | Zero-config |

### Complexity

- Per-locus: O(|locus| + |edges_traversed|)
- For 200-variant locus with pathway annotations: seconds
- Scales linearly with annotation richness

### Output

```python
{
    "variant": str,
    "z_stat": float,
    "z_functional": float,
    "smoothed_stat": float,
    "unique_stat": float,
    "final_score": float,
    "PIP": float,  # posterior inclusion probability
    "in_credible_set": bool,
    "annotations": [str],  # which pathways/genes gave functional boost
}
```

---

## Method L2: LD-Regularized Association

**Addresses:** decomposing GWAS signal into independent contributions by treating LD as a graph regularizer.

### Conceptual basis

Standard GWAS tests each variant independently, producing correlated p-values across LD blocks. Fine-mapping tries to "deconvolve" this post-hoc. Instead: **treat LD as part of the regression model from the start**, using graph Laplacian regularization.

Mathematically: if the signal vector z (one element per variant in a locus) contains correlated signals, we want to find the sparse representation in LD-eigenvector space where each coefficient represents an independent signal.

### Algorithm

```
Input: locus variants V, LD matrix R (r² between all pairs), GWAS z-scores
Output: sparse signal decomposition

1. Compute graph Laplacian of LD:
   A = R (adjacency matrix, r² values)
   D = diag(sum(A, axis=1))  # degree matrix
   L = D - A  # graph Laplacian

2. Eigendecompose Laplacian:
   eigenvalues, eigenvectors = eig(L)
   # low eigenvalues = smooth signals (shared across LD block)
   # high eigenvalues = sharp signals (variant-specific)

3. Project z-scores onto eigenvectors:
   coefficients = eigenvectors^T @ z
   # coefficients show signal's spectral distribution

4. Sparse reconstruction (LASSO in spectral domain):
   # most LD blocks have 1-2 causal variants
   # enforce sparsity: few nonzero coefficients
   solve: min ||z - eigenvectors @ c||² + λ||c||_1

5. Transform back to variant space:
   refined_z = eigenvectors @ c_sparse

6. Rank variants by refined_z magnitude:
   variants with large |refined_z| after sparse decomposition = causal

7. Convert to posterior via softmax:
   PIP(v) = exp(|refined_z(v)|) / sum(exp(|refined_z(u)|))
```

### Intuition

If variants A, B, C are all in LD and all show p=1e-10:
- **Standard fine-mapping:** all three get similar PIP (≈0.33 each)
- **L2 spectral decomposition:** shows that ONE spectral component explains the signal, and only one variant (the causal one) has high loading on that component

### Why this is distinct from existing methods

- SuSiE does sum of single effects in variant space — we do it in **spectral space**
- FINEMAP does stochastic search over configurations — we do **convex optimization**
- Runtime: O(|locus|³) for eigendecomposition, but each locus is small (100-500 variants)
- No LD reference panel needed if using genotypes directly

### Connection to existing methods

L2 is mathematically related to:
- **Graph signal processing** (emerging field)
- **Spectral clustering** (eigendecomposition of graph Laplacian)
- **L1-regularized regression** (LASSO)

But the combination is novel for GWAS fine-mapping.

### Complexity

- Laplacian eigendecomposition: O(|locus|³)
- LASSO: O(|locus|² × iterations)
- Fast enough for loci up to 1000 variants

### Output

```python
{
    "variant": str,
    "z_score_raw": float,
    "z_score_refined": float,
    "spectral_coefficient": float,
    "sparsity_rank": int,  # which sparse component
    "PIP": float,
    "in_credible_set": bool,
}
```

---

## Method L3: Haplotype Graph Association

**Addresses:** working with natural inheritance units instead of individual variants, side-stepping the LD problem entirely.

### Conceptual basis

LD exists because variants are inherited as **haplotypes** — contiguous chromosomal segments that travel together. Instead of testing variants, test **haplotypes directly**.

- A haplotype is a specific combination of alleles across variants in a region
- In an LD block of K variants, there are 2^K possible haplotypes but usually only 10-50 observed
- Causal variants can be identified by comparing "disease haplotypes" to "healthy haplotypes"

### Algorithm

```
Input: phased genotypes in a locus window, phenotype
Output: causal variants necessary to make a haplotype "disease-associated"

1. Build haplotype graph:
   nodes = distinct observed haplotypes in the locus
   edges = recombination events (differ by 1 crossover)
   node attributes:
     - variant signature (binary string)
     - case frequency, control frequency
     - case/control log ratio

2. Identify disease haplotypes:
   for each haplotype h:
       case_freq(h) = count of h in cases / total case haplotypes
       ctrl_freq(h) = count of h in controls / total ctrl haplotypes
       OR(h) = case_freq(h) / ctrl_freq(h)
   disease_haps = {h : OR(h) > threshold and case_freq(h) > min_freq}

3. Find causal variants via haplotype comparison:
   # variants that differentiate disease from non-disease haplotypes
   for v in variants:
       disease_alleles = [h.variant[v] for h in disease_haps]
       healthy_alleles = [h.variant[v] for h in other_haps]
       if disease_alleles are consistently different from healthy_alleles:
           v is a candidate causal variant

4. Formal test (Fisher's exact on haplotype-variant table):
   for v:
       table = [[disease_haps_carrying_alt, disease_haps_carrying_ref],
                [healthy_haps_carrying_alt, healthy_haps_carrying_ref]]
       p_value = fisher_exact(table)

5. Rank by p-value, apply multiple testing correction over locus
```

### Why this works for fine-mapping

- LD confusion disappears: we're comparing *actual inherited units*, not correlated variant patterns
- The causal variant(s) are the ones whose alleles differ systematically between disease and healthy haplotypes
- Variants that are just "riding along" don't consistently differ

### Computational considerations

- Locus size: typically 100-500 variants (phase-informative window)
- Distinct haplotypes: usually 10-100 per locus
- Computation: per-locus, fast (seconds)
- **Requires phased genotypes** (1KG phase 3 is phased, most modern cohorts are too)

### When this shines

- LD blocks with many variants but few distinct haplotypes (typical in non-recombinant regions)
- Cases where one haplotype dominates disease risk (classic founder effects)
- When phasing is reliable (statistical phasing or long-read sequencing)

### Output

```python
{
    "variant": str,
    "haplotype_test_p": float,
    "disease_haps_containing_alt": int,
    "disease_haps_containing_ref": int,
    "healthy_haps_containing_alt": int,
    "healthy_haps_containing_ref": int,
    "differentiation_score": float,
    "in_credible_set": bool,
}
```

---

## Method L4: Recombination-Aware Embedding

**Addresses:** identifying independent causal signals within a locus via embedding in recombination space.

### Conceptual basis

Physical distance does not equal genetic distance. Two adjacent variants in a recombination coldspot have r²≈1; two adjacent variants across a hotspot have r²≈0.

**Embed variants in a metric space where distance reflects recombination probability, not base pairs.** Then:
- Variants in LD cluster tightly (small recombination distance)
- Recombination hotspots create gaps
- Number of clusters = number of independent causal signals
- Each cluster's "center" (via some criterion) = causal candidate

### Algorithm

```
Input: locus variants V, LD matrix R, GWAS z-scores
Output: causal candidates per independent signal

1. Define distance metric:
   d(v_i, v_j) = -log(r²(v_i, v_j) + ε)
   # perfectly linked: d=0
   # independent: d=infinity

2. Embed variants using MDS or UMAP:
   coords = MDS(distance_matrix=d, n_components=2)
   # 2D representation preserving recombination distances

3. Identify clusters (independent signals):
   use DBSCAN or hierarchical clustering
   # each cluster = one independent signal in the locus

4. For each cluster, find the causal candidate:
   # variants with highest signal in the center of the cluster
   for cluster in clusters:
       # candidate 1: highest z-score
       top_z = max(cluster, key=z_score)
       # candidate 2: most central (smallest mean distance to others)
       center = min(cluster, key=lambda v: mean(d(v, u) for u in cluster))
       # candidate 3: highest functional score (if annotations available)
       top_func = max(cluster, key=functional_score)

5. Rank candidates per cluster:
   primary candidate: top-z intersected with top-center
   secondary: functional candidate if different
```

### Connection to biology

- Recombination hotspots are well-documented (PRDM9 binding sites)
- LD blocks are defined by hotspot boundaries
- This method automatically discovers block structure from data

### Why this is graph-native

The "distance" is computed from the LD graph, not from physical coordinates. The embedding preserves the **intrinsic geometry** of the locus — how it recombines, not just where variants are located.

### Complexity

- Distance matrix: O(|locus|²)
- MDS/UMAP: O(|locus|² × iterations)
- Clustering: O(|locus|²)
- Per-locus: seconds

### Output

```python
{
    "variant": str,
    "cluster_id": int,  # which independent signal
    "embedding_x": float,
    "embedding_y": float,
    "distance_to_cluster_center": float,
    "z_score": float,
    "is_primary_candidate": bool,
    "PIP": float,
}
```

### When this shines

- Loci with multiple independent causal signals (not apparent from raw p-values)
- Loci spanning recombination hotspots
- Visual interpretation: you can plot the embedding and see the signal structure

---

## Unified Credible Set

Each method produces its own ranking. Combine via:

### Option A: Consensus credible set

For each variant, count how many methods include it in their top-K. Variants appearing in ≥3/4 methods form the consensus set.

### Option B: Borda rank aggregation

```
rank_aggregate(v) = sum(rank_method_i(v) for i in methods)
```

Variants with lowest aggregate rank are top candidates.

### Option C: Per-method, let user choose

Expose all four methods with CLI flags. Users select based on their question:
- "I have functional annotations" → L1
- "I have many correlated variants" → L2
- "I have phased genotypes" → L3
- "I want to visualize multi-signal structure" → L4

---

## Simulation Framework for Ground Truth

### Simulation scenarios

**S1: Single causal variant in LD block**
- Pick one variant as causal
- Simulate phenotype: Y = β × G_causal + noise
- Add 50 variants in strong LD (r² > 0.8) as proxies
- Test: does each method rank causal variant #1?

**S2: Two causal variants in same locus (independent)**
- Two variants in same 1Mb window
- r² between them: 0.0 - 0.3 (independent)
- Both contribute: Y = β1 × G_1 + β2 × G_2 + noise
- Test: do methods identify both signals?

**S3: Causal variant + proxy (different effect sizes)**
- Main causal variant (β=0.5)
- Weakly correlated "causal-looking" variant in LD
- Test: does the method avoid false positives?

**S4: Causal variant with functional annotation**
- Causal variant in an annotated gene
- Matched LD partners in unannotated regions
- Test: does L1 (dual-graph) boost functional signal?

**S5: Multi-signal locus (3+ independent signals)**
- 3-4 independent causal variants in same 5Mb region
- Test: can L4 (embedding) detect multiple clusters?

### Real-world validation

- Apply to known fine-mapped loci in literature (BMI, height, T2D)
- Compare to SuSiE credible sets on real data
- Measure agreement and complementarity

---

## Comparison Baselines

| Tool | Type | Use for comparison |
|------|------|-------------------|
| **SuSiE** | Bayesian, sum of single effects | Gold standard |
| **FINEMAP** | Stochastic search | Multiple causal variants |
| **PAINTOR** | Functional annotation-aware | Direct competitor to L1 |
| **GraphGWAS v1 (centrality)** | LD topology only | Internal baseline |
| **Top p-value** | Naive baseline | Worst case |

---

## Metrics

### Per-method metrics

**Causal variant rank:** rank of true causal variant in method's ranking (ideal: 1)

**Credible set size:** number of variants in 95% credible set (smaller = more precise)

**Causal inclusion rate:** fraction of simulations where causal variant is in 95% CS

**Posterior calibration:** does PIP=0.9 correspond to 90% of variants being truly causal?

**Runtime per locus**

### Comparative metrics

**Method complementarity:** do methods identify different causal candidates?

**Consensus accuracy:** when ≥3/4 methods agree, is the consensus correct?

**Functional enrichment:** do top candidates enrich for known functional annotations more than LD-based methods alone?

---

## Implementation Plan

### New module: `src/python/graphgwas/finemapping_v2.py`

```python
# Method entry points
def dual_graph_finemap(conn, lead_snp, window_bp=500_000,
                        annotation_weights={"pathway": 1.0, "gene": 0.8},
                        alpha=0.5)
def ld_regularized_finemap(conn, lead_snp, window_bp=500_000,
                            lambda_sparsity=0.1)
def haplotype_graph_finemap(conn, lead_snp, window_bp=200_000,
                             min_hap_freq=0.01)
def recombination_embedding_finemap(conn, lead_snp, window_bp=500_000,
                                      eps_cluster=0.5)

# Unified interface
def finemap(conn, lead_snp, method="auto", **kwargs)

# Simulation
def simulate_locus_with_causal(genotypes, n_causal=1, r2_range=(0.8, 0.99))

# Benchmarking
def benchmark_finemapping(conn, simulations, methods, baselines=["SuSiE"])
```

### New CLI commands

```bash
# Specific method
graphgwas finemap-v2 dual-graph --lead-snp chr22:17000061 --window 500000

# Run all methods on a locus
graphgwas finemap-v2 all --lead-snp chr22:17000061 -o credible_set.tsv

# Benchmark
graphgwas finemap-v2 benchmark --simulation S1 --n-reps 50 -o benchmark.tsv
```

### Dependencies

- Reuse: `finemapping.py` (existing centrality-based), `genotype.py` (LD computation)
- New: `scikit-learn` for UMAP/MDS/DBSCAN (already installed)
- New: hapix / ShapeIt output parsing for phased haplotypes (if implementing L3 fully)

---

## Runtime Estimates

Per-locus (500-variant window):

| Method | Yeast | Human chr22 | Human genome-wide (10K loci) |
|--------|-------|-------------|------------------------------|
| L1 dual-graph | 2s | 2s | 5.5 hr |
| L2 LD-regularized | 5s | 5s | 14 hr |
| L3 haplotype graph | 10s | 10s | 28 hr |
| L4 embedding | 3s | 3s | 8 hr |

All methods operate per-locus; parallelize over loci for full-genome runs.

---

## Success Criteria

### Minimum (proof of concept)

- [ ] At least one method achieves causal variant in top-5 for >70% of single-causal simulations
- [ ] Methods have calibrated credible sets (90% CS covers true causal ~90% of time)
- [ ] At least one method matches or beats SuSiE on S1 simulation

### Strong (publishable)

- [ ] L1 (dual-graph) outperforms SuSiE when functional annotations are informative
- [ ] L2 (LD-regularized) matches SuSiE on standard single-signal loci
- [ ] L3 (haplotype) identifies correct causal variant where SuSiE is confused
- [ ] L4 (embedding) correctly identifies multiple independent signals in S5

### Stretch

- [ ] Combined ranking (4 methods) outperforms any single method
- [ ] Methods discover causal variants missed by SuSiE on real data
- [ ] Method selection rule (L1 vs L2 vs L3 vs L4) is learnable

---

## Open Research Questions

1. **Annotation weighting:** In L1, how to weight different edge types (gene vs pathway vs regulatory)? Should weights be learned from data?

2. **Sparsity tuning (L2):** How to choose λ in LASSO? Cross-validation? BIC?

3. **Haplotype block definition (L3):** What window size? How to handle recombination?

4. **Multi-signal detection (L4):** How to choose cluster number automatically?

5. **Meta-calibration:** If four methods give PIPs, how to combine into a consensus PIP with proper calibration?

6. **LD reference:** Do we compute LD on the fly from genotypes, or use reference panels (1KG)?

---

## Proposed Timeline

**Week 1: Simulation framework for fine-mapping**
- S1-S5 scenarios with real 1KG LD structure
- Implementation of known-causal simulation

**Week 2: L1 + L2 implementation**
- Dual-graph fine-mapping (L1)
- LD-regularized association (L2)
- Benchmark against SuSiE on S1-S3

**Week 3: L3 + L4 implementation**
- Haplotype graph association (L3)
- Recombination-aware embedding (L4)
- Benchmark on S2-S5

**Week 4: Unified framework + real-data comparison**
- Consensus credible set construction
- Apply to known loci (BMI, height)
- Compare with SuSiE and FINEMAP

---

## Risks and Mitigations

| Risk | Mitigation |
|------|-----------|
| None of the methods beat SuSiE | Still valuable: faster, more interpretable, graph-queryable |
| L1 functional boost only works with good annotations | Make annotation usage optional; fall back to pure LD |
| L3 requires phased data we don't have | Use 1KG phased data as reference; statistical phasing if needed |
| L2 sparsity too aggressive (missing secondary signals) | Tune λ via cross-validation |
| L4 clusters are unstable with small loci | Use as secondary signal detector, not primary |

---

## Final Deliverable

A benchmark showing:
1. **Each method has distinct strengths** (L1 for annotated regions, L2 for simple loci, L3 for phased data, L4 for multi-signal)
2. **Methods complement SuSiE** rather than directly replacing it
3. **Combined approach improves credible set quality** over any single method
4. **Graph-native integration of annotations** outperforms matrix-based multi-omic fine-mappers (PAINTOR)

This gives users principled choices: "I have functional data → L1; I need speed → L4; I want complete picture → run all, compare."

---

## Cross-References to Epistasis v2

Some synergies between this document and EPISTASIS_V2_DESIGN.md:

1. **LD pruning in M1 (epistasis)** can use L2's spectral decomposition to identify truly independent variants (not just r² threshold)

2. **Motif queries in M2 (epistasis)** can be combined with L1's functional scoring to prioritize epistatic pairs in annotated contexts

3. **Haplotype graphs (L3)** could enable haplotype-level epistasis testing — a novel direction

4. **Embedding (L4)** could identify independent loci for epistasis testing (ensure pairs are not in LD)

---

## Ready for Implementation?

**Before coding:**
- [ ] Review this design with user
- [ ] Decide which methods to implement first (recommend: L1 + L4 first — most novel signals)
- [ ] Decide on annotation strategy for L1 (user-specified vs GWAS-derived)
- [ ] Decide on sparsity tuning approach for L2
- [ ] Confirm we have phased genotypes for L3 (yes, 1KG is phased)
