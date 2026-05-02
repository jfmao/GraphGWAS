# Paper #2 §Y.4 — Yeast real-data validation: multi-method comparison

**Question:** does each of the five implemented epistasis methods (M1
LD-pruned, M2 motif-filtered, M3 differential subgraph, M4 dark matter,
M5 random-walk) recover the literature-canonical BCY1 × TPK1 cAMP-PKA
epistatic pair from yeast 1011 Genomes genotypes?

**Script:** `tests/validate_m2_yeast.py`
**Output:** `results/paper2_epistasis/yeast_validation.json`

## Coverage

| Method | Source-agnostic implementation? | Tested in this validation? |
|--------|---|---|
| **M1** | ✓ `ld_pruned_cooccurrence_from_data` (epistasis_v2.py:829) | ✓ |
| **M2** | ✓ `motif_filtered_epistasis_from_data` (epistasis_v2.py:454) | ✓ |
| **M3** | ✗ `differential_subgraph` is Neo4j-only; needs lifting | not yet — see Stage 2 |
| **M4** | ✗ `dark_matter_epistasis` is Neo4j-only; needs lifting | not yet — see Stage 2 |
| **M5** | ✓ `mutual_rwr_pair_scores` (epistasis_higher_order.py) | ✓ |

M3 and M4 also target *different signal types* (case/control
co-occurrence enrichment for M3; synthetic-lethal depletion for M4),
which requires different phenotype simulators (Scenarios B and C in the
plan).  Both are deferred to a follow-up session.

## Setup (shared across all three tested methods)

| Item | Value |
|------|-------|
| Genotype panel | 1011 Yeast Genomes, MAF ≥ 0.05 |
| Loaded variants | 6,625 SNPs across chromosomes 9, 10, 12 |
| Samples | 1,011 |
| Annotation source | `data/yeast/yeast_graph_cache_v2.json` |
| Cache hit rate | 4,572 / 6,625 = 69 % |
| Ground-truth pair | BCY1 (YIL033C) × TPK1 (YJL164C) — yeast PKA regulatory + catalytic subunits |
| Recovery matching rule | **gene-level** — *any* (BCY1-variant × TPK1-variant) pair in the result list counts as recovery; the simulated representative pair may be LD-equivalent to another pair in the same gene |
| Phenotype DGP | y = β · G_BCY1 · G_TPK1 + ε (β = 3, causal R² ≈ 10 %) |

## Results: side-by-side multi-method ranks

| Method | n_pairs_tested | sig@q<0.05 | **GT rank** | q-value | β̂ | runtime |
|--------|---:|---:|:---:|---:|---:|---:|
| M1 LD-pruned (Bonferroni) | 133,837 | 346 | **NF** | — | — | 55 s |
| M2 motif-filtered (BH-FDR) | 94,440 | 5,033 | **1,148** (top 1.2 %) | 1.66 × 10⁻⁵ | +2.32 | 12 s |
| M5 RWR + interaction (BH-FDR) | 1,977 | 89 | **NF** | — | — | 12 s |

M2 is the only method that recovers the BCY1 × TPK1 pair on this panel.
The other two methods miss it for **method-specific, paper-worthy reasons**:

### Why M1 misses

`ld_pruned_cooccurrence_from_data()` greedy-prunes variants at r² ≥ 0.5
within 10 kb, keeping the first variant by genomic position.  For our
panel:

- BCY1 representative (chromosome9:290536:C:T) is **pruned by
  chromosome9:282798:A:G** (r² = 0.51).  But chr9:282798 has **no gene
  annotation** in the cache — it's intergenic.  M1's "BCY1 LD block" is
  thereby represented by an unannotated variant.
- TPK1 representative (chromosome10:110637:A:G) is **pruned by
  chromosome10:105170:C:T** (r² = 0.89), which is annotated as **YJL167W**
  — a different gene.  M1's "TPK1 LD block" is represented by a
  YJL167W-tagged variant.

→ Even with gene-level recovery matching, M1's result list doesn't
contain *any* (BCY1-variant × TPK1-variant) pair — because the LD-pruning
step removed all BCY1 and TPK1 representatives from the pool.

**Method-improvement implication for paper #2:** M1 needs *gene-aware*
LD pruning (don't prune across gene boundaries) to be competitive on
gene-level epistasis tests.  This is a low-cost fix and will land in a
follow-up commit.

### Why M5 misses

M5 scores variant pairs by mutual random-walk-with-restart probabilities
on the bipartite variant–gene graph.  The yeast cache subset has
4,572 variants × 887 genes; only **2 BCY1 variants** and **5 TPK1
variants** survive into the cache ∩ dosage intersection (most BCY1/TPK1
cache entries are below the panel's MAF≥0.05 threshold).

With only 7 BCY1+TPK1 seeds and 4,565 random others, the cross-gene
mutual-RWR scores between BCY1 and TPK1 don't make the top 2,000 pair
list.  The 89 BH-FDR-significant pairs M5 returned are dominated by
high-RWR-density genes (mating-type, ribosome-biogenesis), not BCY1/TPK1.

**Method-improvement implication for paper #2:** M5 needs either (i) a
gene-pair seeding strategy that prioritises gene-pair-level mutual RWR
above variant-pair-level (paper #2 §Y.5 higher-order extension), or
(ii) lower MAF cutoff to include more BCY1/TPK1 variants in the pool.

### Why M2 wins

M2 enumerates pairs by *gene/pathway/PPI motif* — not by position-LD or
graph density.  Every (BCY1-variant × TPK1-variant) pair where both
variants are MAC ≥ 10 enters the testing pool, regardless of LD or
graph-walk reachability.  The motif-typed enumeration **preserves
gene-level context**, which is the unit of biological hypothesis.

This is the central paper-#2 claim: **biology-typed motif filtering is
the right inductive bias when the hypothesis space is gene-level
epistasis.**  M1 (position-LD-typed) and M5 (graph-density-typed) bring
their own valuable inductive biases for *other* signal types, but on a
gene-level test M2 dominates.

## Honest limitations + future work

- **Only 3 of 5 methods tested**: M3 and M4 require (a) lifting from
  Neo4j to source-agnostic and (b) Scenarios B/C with appropriate
  phenotype simulators.  Stage 2 + Stage 3 of the §Y.4 plan.
- **Single ground-truth pair**: BCY1 × TPK1 only.  Paper #2 should
  ideally also test on (HSP104 × Sup35), (DPY1/DAL5 × NPR1), and other
  known yeast epistatic gene pairs from Bloom et al. 2015.
- **Single phenotype DGP**: positive interaction y = β·g₁·g₂ + ε.
  M3/M4 require qualitatively different DGPs (case/control enrichment,
  synthetic-lethal depletion).
- **Yeast 1011 panel has predominantly rare variants in BCY1 and TPK1**
  (best representative MAFs 0.06 and 0.07).  Larger biobank-scale panels
  with higher-MAF variants in PKA-pathway genes would give stronger
  detectability across all methods.

## Reproduction

```bash
# Multi-method run with simulated BCY1×TPK1 interaction
python tests/validate_m2_yeast.py --beta 3.0 --seed 2026

# Null-FPR control
python tests/validate_m2_yeast.py --beta 0.0 --seed 2026
```

Each run takes ~80 seconds (M1: 55 s, M2 + M5: 12 s each).
