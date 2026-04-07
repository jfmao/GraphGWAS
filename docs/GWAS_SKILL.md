# GWAS_SKILL.md — Graph-Native GWAS Domain Knowledge

## Purpose
This file is the domain knowledge reference for the GraphGWAS implementation.
Read this before implementing any GWAS-related procedure. It explains WHY each
design decision was made, not just WHAT to build.

---

## 1. Why Classical GWAS Fails — The Structural Problem

Standard GWAS is a **massively parallel univariate regression**. For each of ~10M
variants, it fits: `phenotype ~ genotype_i + covariates`

This is done 10 million times, independently. The result is a collection of *marginal
associations*, not a *joint model of the genome*. This creates three structural failures:

### 1.1 The Epistasis Blindspot
If variant A only has an effect when variant B is also present (genetic interaction),
neither A nor B will reach significance in isolation. The signal is split and attenuated.
**Standard GWAS structurally cannot detect this.**

### 1.2 The LD Proxy Problem
What GWAS "finds" is rarely the causal variant — it finds LD neighbors. Fine-mapping
tries to recover from this post-hoc. The graph has the LD structure natively; GWAS
should use it during discovery, not after.

### 1.3 Missing Heritability
Significant single-locus associations explain only a fraction of trait heritability for
most complex traits. The rest is distributed across millions of variants, their
interactions, and their context-dependence on genetic background — none of which
single-locus tests can capture.

**The root cause:** The matrix representation (samples × variants) has no natural way
to represent relationships *between* variants, the functional context of variants, or
the network topology of genetic effects. The graph has all of this natively.

---

## 2. The Imbalanced Data Problem (Rare Disease Context)

### 2.1 Why Imbalance Breaks Classical GWAS

With a rare disease (e.g., 200 cases, 100,000 controls):
- Every test operates on a response vector that is 99.8% zeros
- Standard logistic regression loses power; asymptotic approximations break down
- Firth regression and saddle-point corrections are patches, not solutions
- Rare causal variants (AF < 0.01) create genotype vectors that are ALSO nearly all zeros
- The correlation of two near-zero vectors is statistically undetectable

### 2.2 Why Graph Architecture Is Structurally Better

**Key insight:** In a matrix, sparsity is the enemy. In GraphGWAS's CARRIES graph,
sparsity is the native representation.

- A rare disease with 200 cases = 200 Sample nodes labelled `case`
- A rare causal variant with 200 carriers = 200 CARRIES edges
- The question "do case samples disproportionately share CARRIES edges to the same
  rare variants?" is a **graph topology question** — exactly what the graph answers efficiently

### 2.3 The Five Graph-Native Solutions to Imbalance

**Solution 1 — Graph-Native BURDEN Test**
Instead of testing rare variants individually (no power), count rare, high-impact
variants reachable from each case sample. Compare case burden to control burden.
Can be conditioned on pathway membership in a single Cypher query — impossible in PLINK.

**Solution 2 — Case-Specific Variant Subgraph**
Build the variant co-occurrence subgraph *only* for case samples. Even with 200 cases,
the topology of this subgraph (which variants co-cluster, which are hubs) carries
association signal that no individual variant test would reveal.

**Solution 3 — Kinship-Aware Cosegregation**
For rare diseases, cases are often related (familial forms, founder populations).
Classical GWAS treats this as a confounder to correct away. GraphGWAS traverses
KINSHIP edges to find variants co-carried by related cases — cosegregation analysis
at population scale. Related cases sharing a rare variant are far more likely to
implicate a causal locus than unrelated cases sharing by chance.

```cypher
MATCH (s1:Sample {phenotype: 'case'})-[:KINSHIP {coefficient: coeff}]-(s2:Sample {phenotype: 'case'})
WHERE coeff > 0.125
MATCH (s1)-[:CARRIES]->(v:Variant)<-[:CARRIES]-(s2)
WHERE v.af_total < 0.01
RETURN v.id, count(DISTINCT s1) AS family_case_carriers, v.gene_symbol
ORDER BY family_case_carriers DESC
```

**Solution 4 — Pathway-Level Label Propagation (Ultra-Rare Diseases)**
For 20–50 cases, no single-variant or gene-level test has power.
Phenotype diffusion through the functional graph aggregates signal across many rare
variants simultaneously:
1. Label case Sample nodes with phenotype score 1.0
2. Propagate through CARRIES edges to Variant nodes (weighted by 1/AF to upweight rare variants)
3. Propagate through HAS_CONSEQUENCE → Gene → IN_PATHWAY edges
4. Pathways accumulating the most signal are candidate disease pathways

With 20 cases each carrying different rare variants in the same pathway, the pathway
node accumulates 20× the propagated signal even though no individual variant is significant.

**Solution 5 — Cross-Phenotype Borrowing**
Rare diseases rarely exist in isolation. Patients often have comorbidities with shared
genetic architecture. Multi-label Sample nodes allow borrowing statistical strength
from phenotypically adjacent conditions in a single query.

---

## 3. Multi-Locus GWAS — The Graph-Native Architecture

### 3.1 Variant Co-Occurrence as Native Graph Object

The first computation that PLINK cannot do natively: **variant co-occurrence structure**.

```cypher
// Find pairs of variants co-carried above expected frequency
MATCH (s:Sample)-[:CARRIES]->(v1:Variant),
      (s)-[:CARRIES]->(v2:Variant)
WHERE v1.pos < v2.pos AND v1.chr = v2.chr
WITH v1, v2, count(s) AS co_count
WHERE co_count > expectedCoOccurrence(v1, v2)
MERGE (v1)-[:CO_OCCURS {count: co_count, excess: co_count - expectedCoOccurrence(v1, v2)}]->(v2)
```

This builds a Variant–Variant co-occurrence graph conditioned on actual sample carriers.
This is LD, but richer — it is the *joint distribution* of alleles across samples.

### 3.2 Phenotype-Conditioned Subgraph Traversal

The key query no matrix tool can express in one step:

```cypher
// Find variant pairs where co-carriage is enriched in cases vs controls
MATCH (s:Sample {phenotype: 'case'})-[:CARRIES]->(v1:Variant),
      (s)-[:CARRIES]->(v2:Variant)
WHERE v1 <> v2
WITH v1, v2, count(DISTINCT s) AS case_co_count

MATCH (s2:Sample {phenotype: 'control'})-[:CARRIES]->(v1),
      (s2)-[:CARRIES]->(v2)
WITH v1, v2, case_co_count, count(DISTINCT s2) AS ctrl_co_count
WHERE case_co_count * 1.0 / (ctrl_co_count + 1) > enrichment_threshold
RETURN v1.id, v2.id, case_co_count, ctrl_co_count
```

This is a **graph-native epistasis screen** — identifying variant pairs whose *joint*
carriage is phenotype-associated, independent of their marginal effects.

### 3.3 Higher-Order Modules via Community Detection

Beyond pairs, variant *communities* in the case-phenotype co-occurrence graph are
*putative epistatic modules* — sets of variants that tend to co-occur in the same
genetic backgrounds and are jointly associated with the phenotype.

Use Neo4j GDS Louvain/Leiden on the case-conditioned co-occurrence graph projection.
Each community is a candidate multi-locus association module.

### 3.4 Phenotype Flow: The Theoretical Framing

Complex trait association is fundamentally a **flow problem on a multilayer graph**,
not a regression problem on independent variants:

- Layer 1: Genomic (variant ↔ variant via LD, haplotype structure)
- Layer 2: Functional (variant → gene → pathway → biological process)
- Layer 3: Population (sample ↔ sample via kinship, IBD, ancestry)
- Layer 4: Phenotypic (phenotype ↔ phenotype via genetic correlation)

A causal variant is not a "significant SNP" — it is a **node with high betweenness
centrality in phenotype-conditioned flow** across all four layers simultaneously.
A classical GWAS p-value is a projection of this 4D flow onto a single dimension.

---

## 4. Population Stratification in Graph GWAS

Population stratification is the main confounder in GWAS. Classical correction uses
PCA covariates derived from the genotype matrix. GraphGWAS can do this natively:

- **PCA from graph-derived GRM:** `graphpop.structure.pca()` gives PC coordinates
  that can be directly used as covariates. Store them on Sample nodes as `pca_coords`.
- **Community-based stratification:** Leiden/Louvain community labels from the
  Sample–Sample kinship graph are a graph-native alternative to discrete ancestry labels.
- **Embedding-based stratification:** Node2vec or GraphSAGE embeddings of the
  Sample–Sample kinship/IBD graph capture continuous ancestry gradients.

**Novel capability:** GraphGWAS can condition association tests on *graph-derived*
population structure — not just classical PCA, but community membership, embedding
coordinates, and local ancestry estimates, all computable natively without external tools.

---

## 5. The Disease Subgraph Concept

This is the deepest conceptual innovation of GraphGWAS.

**Classical GWAS asks:** for each variant i, does P(carrier | case) > P(carrier | control)?

**GraphGWAS asks:** what is the **minimal connected subgraph** of the
variant–gene–pathway graph that is maximally enriched in case carriers,
conditioned on population structure?

The answer is not a list of p-values. It is a **subgraph**: a network of variants,
genes, and pathways that together constitute the genetic architecture of the disease.

### Why This Is Powerful for Rare Diseases
- Causal architecture may be spread across 10–50 rare variants in 3–5 genes in 1–2 pathways
- No individual component reaches significance alone
- But the connected subgraph — all variants connected through gene and pathway nodes —
  is massively enriched in case samples
- The subgraph IS the association signal

### Subgraph Extraction Algorithm (Conceptual)
1. Label all case-carrying variants with enrichment score = (case_freq / control_freq)
2. Label Gene nodes with aggregate enrichment from connected variants
3. Label Pathway nodes with aggregate enrichment from connected genes
4. Extract the minimum spanning subgraph above enrichment threshold
5. Apply statistical correction for subgraph connectivity (permutation-based)

---

## 6. GNN Association: Message-Passing as GWAS

### Why GNNs Replace Regression for Multi-Locus Association

A heterogeneous GNN over the Variant–Sample–Gene–Pathway graph with:
- Node types: Variant, Sample, Gene, Pathway
- Edge types: CARRIES, LD, HAS_CONSEQUENCE, IN_PATHWAY, KINSHIP

Trained to predict phenotype from the *graph neighborhood* of each sample —
not from any single variant, but from the entire local structure of variants,
their LD relationships, their functional annotations, and their pathway memberships.

The GNN learns the multi-locus epistatic architecture implicitly because
**message-passing is exactly the mathematical operation that aggregates neighborhood
information**. It does naturally what epistasis testing tries to do combinatorially.

### Implementation Notes
- Use PyTorch Geometric (PyG) with HeteroData for heterogeneous graphs
- Export from Neo4j GDS using graph projections to PyG format
- Architecture: R-GCN or HGT (Heterogeneous Graph Transformer) for multi-relational graphs
- Training objective: phenotype prediction (binary cross-entropy for case-control,
  MSE for quantitative traits)
- Feature engineering: Variant features = [af, consequence_encoded, cadd_score, ld_degree];
  Sample features = [pca_coords, ancestry_proportions, sex_encoded]

---

## 7. Using GraphPop Statistics as GWAS Covariates

**This is a unique capability of the combined GraphPop + GraphGWAS platform.**

GraphPop computes evolutionary statistics that can serve as covariates in GWAS:
- `Variant.tajima_d` (from GenomicWindow): regions under balancing selection have
  elevated Tajima's D; conditioning on this separates selection signal from association signal
- `Variant.ihs_score`: haplotype-based selection statistics as variant-level covariate
- `Variant.fst` (between populations): variants with high F_ST may reflect local adaptation
  rather than phenotype association; condition out
- `Sample.embedding` (from GraphPop GNN): population structure embedding richer than PCA

**Example:** A variant with high F_ST in an admixed cohort may appear as a GWAS hit
simply because it tags ancestry, not because it is causal. GraphGWAS can condition on
graph-derived F_ST in a single procedure call, eliminating this class of false positive.

---

## 8. Statistical Considerations

### 8.1 Multiple Testing Correction
- Standard: Bonferroni at 5×10⁻⁸ for single-locus (10M independent tests)
- For multi-locus co-occurrence: permutation-based FDR (shuffle case/control labels,
  recompute co-occurrence enrichment, empirical null distribution)
- For pathway diffusion: permutation of Sample labels (not variant labels) to preserve
  LD structure in the null

### 8.2 Effect Size Representation
- Single-locus: odds ratio (binary) or beta (quantitative) with SE and 95% CI
- Multi-locus/epistasis: interaction term delta (joint effect minus sum of marginals)
- Pathway diffusion: normalized diffusion score (z-score against permuted null)
- GNN: SHAP values for variant-level feature attribution

### 8.3 Quality Control Filters (Apply Before Any Analysis)
- Hardy-Weinberg equilibrium: exclude variants with HWE p < 1e-10 in controls
- Call rate: exclude variants with call_rate < 0.95
- Minor allele count: exclude variants with MAC < 5 in cases (for case-control)
- Sample missingness: exclude samples with genotype missingness > 5%
- Relatedness: for single-locus GWAS, exclude one of each pair with kinship > 0.125
  (for rare disease cosegregation, retain and exploit relatedness instead)

### 8.4 Firth Logistic Regression for Rare Variants
Standard logistic regression is biased when variants are rare or when case counts are small.
Use Firth penalized likelihood correction for:
- Any variant with MAC < 20 in cases
- Any dataset with fewer than 500 cases
- Implementation: `graphgwas.assoc.firth_correction()` procedure

---

## 9. Comparison to Classical Tools

| Capability | PLINK2 | SAIGE | Hail | **GraphGWAS** |
|---|---|---|---|---|
| Single-locus marginal association | ✅ Fast | ✅ Mixed model | ✅ Fast | ✅ |
| Firth correction for rare variants | ✅ | ✅ | ⚠️ Partial | ✅ |
| Pairwise epistasis | ❌ O(M²) | ❌ | ❌ | ✅ Co-occurrence graph |
| Functional conditioning in one step | ❌ Multi-step | ❌ | ❌ | ✅ Single Cypher |
| Phenotype propagation through LD graph | ❌ | ❌ | ❌ | ✅ Native label propagation |
| Rare disease subgraph extraction | ❌ | ❌ | ❌ | ✅ Novel |
| Kinship-aware cosegregation | ❌ | ⚠️ Partial | ❌ | ✅ Graph traversal |
| Cross-phenotype borrowing | ❌ Manual | ❌ | ⚠️ Manual | ✅ Multi-label nodes |
| Incremental sample addition | ❌ Re-run | ❌ Re-run | ❌ Re-run | ✅ Add edges + update |
| Association results as queryable graph objects | ❌ | ❌ | ❌ | ✅ AssociationResult nodes |
| Evolutionary stats as covariates | ❌ | ❌ | ❌ | ✅ GraphPop integration |

---

## 10. Key References and Concepts

- **GWAS review:** Visscher et al. (2017) Nature Genetics — 10 years of GWAS
- **Epistasis in GWAS:** Cordell (2009) Nature Reviews Genetics — detecting gene-gene interactions
- **Rare variant methods:** Lee et al. (2012) AJHG — SKAT-O optimal unified test
- **Firth regression:** Firth (1993) Biometrika; Heinze & Schemper (2002)
- **Label propagation:** Zhu et al. (2003) — semi-supervised learning with label propagation
- **Heterogeneous GNN:** Yun et al. (2019) — Graph Transformer Networks; Schlichtkrull et al. (2018) — R-GCN
- **Missing heritability:** Manolio et al. (2009) Nature — finding the missing heritability
- **Methodological opportunities:** See `Methodological_opportunities.pdf` in project knowledge
  (Nature Reviews Genetics 2025) — current GWAS methodological gaps, especially population diversity
