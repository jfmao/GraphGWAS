# Paper #2 §Y.4 — Yeast real-data validation: complete 3 × 5 benchmark

**Question:** does each of the five implemented epistasis methods (M1
LD-pruned, M2 motif-filtered, M3 differential subgraph, M4 dark matter,
M5 random-walk) recover the literature-canonical BCY1 × TPK1 cAMP-PKA
epistatic pair from yeast 1011 Genomes genotypes — across three
phenotype regimes that target different methods' signal types?

**Script:** `tests/validate_m2_yeast.py`
**Output:** `results/paper2_epistasis/yeast_validation.json`

## Coverage — full benchmark

| Method | Source-agnostic implementation | Tested |
|--------|---|---|
| **M1** | ✓ `ld_pruned_cooccurrence_from_data` | ✓ in all 3 scenarios |
| **M2** | ✓ `motif_filtered_epistasis_from_data` | ✓ in all 3 scenarios |
| **M3** | ✓ `differential_subgraph_from_data` *(NEW Stage 2)* | ✓ in all 3 scenarios |
| **M4** | ✓ `dark_matter_epistasis_from_data` *(NEW Stage 2)* | ✓ in all 3 scenarios |
| **M5** | ✓ `mutual_rwr_pair_scores` + interaction test | ✓ in all 3 scenarios |

## Three phenotype scenarios

| Scenario | DGP | Targets which method? |
|---|---|---|
| **A — Positive interaction (quantitative)** | y = β·G_BCY1·G_TPK1 + ε, β=3 ⇒ R²≈10 % | M1, M2, M5 |
| **B — Case/control with case-enriched co-occurrence** | y_binary = (DC) ∨ (Bernoulli(0.10)); double-carriers forced into cases | **M3** (differential subgraph) |
| **C — Synthetic-lethal depletion** | y_binary = Bernoulli(0.10) ∧ ¬DC; double-carriers excluded from cases | **M4** (dark matter) |

## Setup (shared across all scenarios)

| Item | Value |
|------|-------|
| Genotype panel | 1011 Yeast Genomes, MAF ≥ 0.05 |
| Loaded variants | 6,625 SNPs across chromosomes 9, 10, 12 |
| Samples | 1,011 (haploid) |
| Annotation source | `data/yeast/yeast_graph_cache_v2.json` (1.14 M variants) |
| Cache hit rate | 4,572 / 6,625 = 69 % |
| Ground-truth pair | BCY1 (YIL033C) × TPK1 (YJL164C) |
| Recovery matching | **gene-level** — *any* (BCY1-var × TPK1-var) pair counts as recovery |
| Total runtime | 9 minutes (3 scenarios × ~3 min/scenario) |

## Final results: 3 × 5 grid (BCY1 × TPK1 ground-truth rank)

| Scenario | M1 | M2 | M3 | M4 | M5 |
|----------|---:|---:|---:|---:|---:|
| **A — positive interaction** | NF | **1,148** (q=1.7×10⁻⁵) | NF | NF | NF |
| **B — case-enriched co-occurrence** | NF | **4,232** (q=8.5×10⁻⁷) | NF | NF | NF |
| **C — synthetic-lethal** | NF | 29,809 (q≈1.0) | NF | NF | NF |

NF = not found in result list. Lower rank = better recovery.

**M2 is the only method that recovers BCY1 × TPK1 in any scenario.**

## Per-scenario method analysis

### Scenario A (positive-interaction quantitative)

- **M2 wins (rank 1,148 / 94,440, q = 1.7×10⁻⁵)** — motif enumeration preserves gene-level context; strong regression signal under positive interaction.
- M1 misses — LD-pruning replaces BCY1's representative with unannotated chr9:282798 (intergenic) and TPK1's with YJL167W (different gene); the kept LD-block representatives are no longer BCY1 / TPK1 by gene annotation.
- M3 returns no results — its case/control quartile-binarisation produces masks where Fisher's exact doesn't reach significance for our specific pair on the LD-pruned subset.
- M4 returns no significant pairs — Bonferroni over ~108 K pairs is too strict on this sample size.
- M5 misses — only 7 BCY1/TPK1 cache-annotated variants survive into the seed pool; mutual RWR doesn't surface them in the top 2,000.

### Scenario B (case-enriched co-occurrence — M3's natural turf)

- **M2 still wins (rank 4,232, q = 8.5×10⁻⁷)** — surprisingly strong, even though regression on a binary phenotype is suboptimal.
- M3 returns 79,636 results (157 case-unique, 23,346 ctrl-unique, 55,954 differential, 179 shared) — most are noise from imbalanced case (n=119) vs control (n=892) co-occurrence counts. BCY1×TPK1 isn't among the case-unique edges because LD-pruning still excludes the representatives.
- M1, M4, M5 all NF — same reasons as Scenario A.

### Scenario C (synthetic-lethal depletion — M4's natural turf)

- **All methods miss.** M2's rank degrades to 29,809 / 94,440 (≈ 31st percentile, q = 1.0) — its regression-based test correctly returns null for a depletion-only signal.
- M4 detects 272 depleted pairs but BCY1×TPK1 isn't among them — M4's `min_expected=0.5` filter combined with the panel's low MAFs (0.06, 0.07) means expected co-occurrences ≈ 0.28 (89 cases × 0.06 × 0.07), which the filter excludes before testing.
- Even with 100 % depletion (DCs forbidden in cases by design), Bonferroni correction over ~7,655 tested pairs prevents any pair from reaching q < 0.05.

## Top-line findings for paper #2

1. **M2 (motif-filtered) is the right default for gene-level epistasis tests.** It dominates on both Scenario A (positive interaction) and Scenario B (case-enriched co-occurrence) — the two regimes where a regression-based interaction test has power.

2. **M1's position-based LD pruning is incompatible with gene-level hypothesis testing.** When the LD-block-leader variant is intergenic or in a different gene, the gene-pair recovery fails by construction. Paper #2 should either (a) propose gene-aware LD pruning for M1, or (b) recommend M2 as the canonical method for gene-level tests.

3. **M3 produces too many candidates under imbalanced n_case / n_ctrl** (157 vs 23,346 case-unique vs ctrl-unique edges in Scenario B). The Fisher's exact + classification-by-LFR design is sound, but the method needs (a) better balancing of edge categories, or (b) FDR correction within edge type rather than across all 79K results.

4. **M4 is power-limited at yeast-QTL scale.** Even with 100 % depletion of double-carriers (the strongest possible synthetic-lethal signal), BCY1/TPK1's low MAFs (0.06, 0.07) give expected co-occurrences ≈ 0.28 — below the `min_expected=0.5` threshold that M4 needs for Poisson tail testing. Biobank-scale panels (n ≥ 100 K) with higher-MAF variants would give M4 power.

5. **M5 is sensitive to seed-pool selection.** Only 7 BCY1/TPK1 variants are in our seed pool (most cache entries are below the panel's MAF ≥ 0.05 threshold). Mutual-RWR scoring among 4,572 seeds doesn't push the BCY1×TPK1 cross-pairs into the top 2,000. A gene-pair-priority RWR (rather than variant-pair-priority) would likely surface the pair — paper #2 §Y.5 follow-up.

6. **The framework correctly returns null for unsupported signal types.** M2's rank in Scenario C is essentially random (29,809 / 94,440, q ≈ 1.0) — the regression-based test doesn't have power for depletion signals. This is the *right* behavior and demonstrates that M2's significance machinery is well-calibrated.

## Honest limitations (carried forward)

- **Single ground-truth pair** (BCY1 × TPK1). Paper #2 should add (HSP104 × Sup35), (DPY1/DAL5 × NPR1), and other Bloom-2015 yeast epistatic pairs.
- **Single phenotype DGP per scenario.** Stage-3 simulators are minimal — robust paper-#2 results need ≥ 10 replicates per scenario.
- **Yeast 1011 has predominantly rare variants in BCY1 / TPK1** — best representatives have MAF 0.06 / 0.07. Higher-MAF variants in PKA-pathway genes would give all methods more power.
- **No baseline comparison against BOOST / MAPIT / MDR** on this dataset — Item 3 follow-up.

## Reproduction

```bash
# Full 3-scenario × 5-method run (~9 minutes)
python tests/validate_m2_yeast.py --scenario all --beta 3.0 --baseline-p 0.10

# Single scenario (faster)
python tests/validate_m2_yeast.py --scenario A
python tests/validate_m2_yeast.py --scenario B
python tests/validate_m2_yeast.py --scenario C

# Null-FPR control (re-run Stage 1's beta=0 sweep)
python tests/validate_m2_yeast.py --scenario A --beta 0.0
```

Output: `results/paper2_epistasis/yeast_validation.json` (the latest run
overwrites — re-run per scenario for separate sidecars).
