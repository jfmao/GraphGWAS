# Deferred Experiments — Final Report

Three experiments queued during the simulated review-revision loop have now been
executed end-to-end on the GraphGWAS chr22 testbed.

| Experiment | Output | Headline |
|---|---|---|
| 1. Null-overlap permutation | `results/null_overlap/null_overlap.{json,tsv}` | Rice + Arabidopsis: observed >> null (p<0.001); yeast + human need coordinate catalogues |
| 2. Hyperparameter sensitivity | `results/hyperparam_sensitivity/hyperparam_sensitivity.{json,tsv}` | HBP robust at 70% across all perturbations; GAFM highly α-sensitive — Methods/code inconsistency caught |
| 3. Placebo-prior control | `results/placebo_prior/placebo_prior.json` | In simulation: informative cCRE = placebo permuted = 13%; baseline = 10%. No biological-structure gain over heterogeneity per se |

---

## 1. Null-overlap permutation experiment

**Method.** For each species, generate 1,000 random "leads" matched to the chromosome distribution of the observed lead set; compute the 250 kb-window overlap rate against the species' validation catalogue; report the empirical p-value of the observed validation rate vs the null distribution.

**Results.**

| Species | Observed rate | Null mean | Null p95 | p-empirical | Enrichment |
|---|---:|---:|---:|---:|---:|
| Rice (3kRG) | **33.3%** (24/72) | 0.9% | 2.8% | <0.001 | **38.4×** |
| *Arabidopsis* (1001G) | 5.6% (3/54) | 0.0% (sparse catalogue) | 0.0% | <0.001 | (catalogue too sparse to estimate) |
| Yeast (1011) | 18.0% (44/245) | (placeholder) | (placeholder) | n/a | n/a |
| Human (Pan-UKB GW) | 6.5% (21/321) | (placeholder) | (placeholder) | n/a | n/a |

**Interpretation.** Rice and *Arabidopsis* validation rates **are significantly above what would be expected by chance**. The rice 38.4× enrichment is the strongest result. Yeast and human require a coordinate-level catalogue file (not currently in the repo) for proper null calibration; their permuted-flag null trivially equals the observed rate by construction.

**Action for paper.** Add the rice and *Arabidopsis* p-values + enrichments to Methods §"Cross-species fine-mapping protocol" or to Sup Note S4 (validation criterion). Add a footnote acknowledging that yeast and human null calibration is forthcoming.

---

## 2. Hyperparameter sensitivity

**Method.** 10 simulated F1 loci on chr22 (β=0.5, h²=0.10). For each hyperparameter, hold others at default and sweep ±33–50% around the default value.

### HBP (defaults: α=0.6, damping=0.5, n_rounds=5, r2_smooth=0.3)

| Parameter | Sweep | Rank-1 rate range |
|---|---|---|
| α | 0.4 – 0.8 | 70% (constant) |
| damping (λ) | 0.3 – 0.7 | 70% (constant) |
| n_rounds | 3 – 7 | 70% (constant) |
| r2_smooth (τ_L) | 0.2 – 0.4 | 60–70% (single dip at 0.35) |

**HBP is hyperparameter-robust** at 70% across all perturbations within ±50% of default values.

### GAFM (defaults per Methods: α=0.5, r2_smooth=0.3)

| α | Rank-1 rate |
|---:|---:|
| 0.3 | 0/10 (0%) |
| 0.4 | 0/10 (0%) |
| 0.5 (Methods default) | **1/10 (10%)** |
| 0.6 | 5/10 (50%) |
| 0.7 | 6/10 (60%) |
| 0.8 | 6/10 (60%) |
| **0.9 (actual benchmark default)** | **7/10 (70%)** |
| 1.0 (no eQTL prior, pure z-score) | 7/10 (70%) |

| r2_smooth | Rank-1 rate |
|---:|---:|
| 0.2 – 0.4 | 0–10% (highly degraded at α=0.5) |

**GAFM is highly α-sensitive.** Critically, the **Methods text states "α = 0.5 by default"** but the **actual benchmark code uses α = 0.9** (`tests/benchmark_hbp_vs_susie_finemap.py:45`). At α = 0.5, GAFM gets 10% rank-1; at α = 0.9, GAFM gets 70%, consistent with Sup Table S2 strong-signal headline of 73%. Different benchmark scripts use different α values:

| Script | α used |
|---|---:|
| `benchmark_hbp_vs_susie_finemap.py` (head-to-head main) | 0.9 |
| `benchmark_graph_native_finemapping.py` | 0.9 |
| `benchmark_null_simulations.py` | 0.9 |
| `benchmark_pip_calibration.py` | 0.9 |
| `benchmark_functional_causal.py` | 0.5 |
| `benchmark_tool_comparison_fast.py` | 0.5 |
| `benchmark_weak_signal.py` (the 27–2 result) | 0.5 |
| `benchmark_panukb_finemap.py` (sumstats-only) | 1.0 (no eQTL prior) |

**Action for paper.** Methods §GAFM must disclose the scenario-dependent α: α=0.9 for strong-signal head-to-head + null calibration; α=0.5 for weak-signal tissue-specific eQTL benchmark; α=1.0 for sumstats-only Pan-UKB. The Methods text claim "α = 0.5 by default" is inaccurate as a single value.

---

## 3. Placebo-prior control

**Method.** 30 simulated F1 loci on chr22 (β=0.5, h²=0.10). For each lead, run HBP and GAFM (α=0.9 to match benchmark) under three cache configurations: (a) baseline = no prior_score; (b) informative = real cCRE-class prior_score from the v2 chr22 cache; (c) placebo = same v2 cache but prior_score values **permuted** across variants (preserving the marginal distribution, randomising assignment). The placebo is averaged over 5 independent permutations.

**Results.**

| Cache configuration | HBP-tighter rate (CS smaller than GAFM) |
|---|---:|
| Baseline (no prior) | **3/30 (10%)** |
| Informative (real cCRE prior) | **4/30 (13%)** |
| Placebo (5 perms, permuted prior) | **4/30 (13%) ± 0** |

**Interpretation.**

- **In this simulation regime, the gain from adding a non-uniform prior is +3pp (10% → 13%), not the +88pp seen on real Pan-UKB lead loci.**
- **Critically, the gain is entirely captured by the placebo control** — informative and placebo are identical at 13%. This means the +3pp gain from adding any non-uniform prior is **fully explained by heterogeneity per se**, with **zero additional contribution from biological structure** in this regime.
- **This does not refute the original 0% → 88% genome-wide Pan-UKB result** but it bounds the interpretation: in a setting where the causal variant is randomly placed (most simulated F1 loci do not fall in regulatory elements), biological-structure priors and random non-uniform priors perform identically. The 0% → 88% result must therefore depend on enrichment of real disease-associated leads in cCRE-class regions — which is itself a known phenomenon, but is a different mechanism from "the graph prior adds biological information about the causal variant".

**Action for paper.** Discussion §Limitations must report this finding explicitly. Results §"A non-uniform per-variant prior...lets the graph break LD ties" must be updated to acknowledge that the simulation placebo control shows zero biological-structure gain over heterogeneity per se, while the Pan-UKB genome-wide result remains robust under a different interpretation (real-disease-lead enrichment in regulatory elements). The reframing the Round-1 review demanded — "non-uniform vs uniform priors" rather than "heterogeneity vs density" — is reinforced by this finding, with the additional clarification that the gain in simulation is from heterogeneity per se.

---

## Consolidated paper-side actions

1. **Methods §GAFM** — disclose scenario-dependent α (0.9 strong / 0.5 weak / 1.0 sumstats); remove the misleading "α = 0.5 by default" single-value claim.
2. **Methods §"Cross-species fine-mapping protocol"** — add rice and *Arabidopsis* null-overlap p-values and enrichment factors.
3. **Discussion §Limitations** — add the placebo-control finding: simulation-regime gain is fully captured by heterogeneity per se; the Pan-UKB genome-wide result depends on real-disease-lead enrichment in cCRE-class regions, not on the graph prior providing variant-specific biological information.
4. **Results §non-uniform-prior section closing paragraph** — replace "necessary-condition claim" with the simulation-vs-genome-wide reframing reflecting the placebo finding.
5. **Sup Note S2** — add hyperparameter-sensitivity numerical results as a Supplementary Figure (or paragraph) and add the null-overlap permutation results.
