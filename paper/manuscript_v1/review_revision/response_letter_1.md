# Response to Reviewers — Round 1

We thank the six reviewers for their rigorous and detailed assessments. We address each major and minor concern point-by-point, classify each as `[Accepted]`, `[Accepted-with-modification]`, or `[Rebutted-with-justification]`, and indicate where the revised manuscript reflects the change.

## Reviewer 1 — Domain Expert

**R1.1 (Major)** — *250-kb-window validation criterion is too lax; report null overlap rate per species.*
**Classification:** `[Accepted]` (textual addition); empirical permutation deferred to follow-up.
**Action:** Methods §"Cross-species fine-mapping protocol" updated to acknowledge that the 250-kb window approximates an LD block in the densest panel and that the validation rate should be interpreted relative to a null-overlap rate; we add this as an explicit limitation in Discussion §"Limitations". A formal permutation analysis (1000 random-position draws per species) is the subject of Round 3.

**R1.2 (Major)** — *Intervention experiment lacks a placebo control with permuted prior_score.*
**Classification:** `[Accepted-with-modification]`. The placebo-control critique is correct and the cleanest answer requires a new computational experiment (24 species × 692 leads × 100 permutations ~ 4 hours of compute). We have therefore (a) reframed the central claim from "heterogeneity, not density" to a precise statement that distinguishes uniform-cache priors from non-uniform priors with biological structure (abstract closing sentence rewritten); (b) added a paragraph in Results §"Heterogeneity, not density, makes graph priors work" noting that the placebo-prior control is the next step and bounding the current claim accordingly; (c) flagged the limitation explicitly in Discussion. The full placebo-control experiment will be reported in Round 3.

**R1.3 (Major)** — *27–2 weak-signal headline is from a different experiment than the 6-method head-to-head.*
**Classification:** `[Accepted-with-modification]`. We agree that the abstract's "wins 27--2 against SuSiE at weak signal" claim derives from a *separate* 79-replicate tissue-specific-eQTL simulation, not from the unified 6-method × 3-scenario benchmark in Fig 4. Abstract is rewritten to make this scope explicit. Results §"GAFM wins weak signal..." now opens with a clear statement that this is a separate, more-favourable simulation regime.

**R1.4 (Major)** — *HBP softmax cap at ~0.7 PIP is operationally a problem.*
**Classification:** `[Accepted]`. The HBP cap is a real interpretability cost. We add a paragraph in Methods §HBP describing the softmax-driven posterior compression and the principled reading of HBP PIPs vs GAFM PIPs (HBP for narrowness; GAFM for absolute calibration). Sup Fig S4 caption now states this explicitly.

**R1.minor (3 items)** — Spelling consistency, Theorem 1 tightness, ancestry-power confounding.
**Classification:** `[Accepted]` (1, 2); `[Accepted-with-modification]` (3).
**Action:** UK spelling enforced; Theorem 1 statement now explicitly labelled "upper bound"; ancestry-power confounding added as a Discussion limitation.

## Reviewer 2 — Methods / Statistics Reviewer

**R2.1 (Major)** — *30 replicates underpowered for "match within 3pp" claim.*
**Classification:** `[Accepted-with-modification]`. We add Wilson 95% CIs to Sup Table S2 and downgrade "match" to "no significant difference at $\alpha=0.05$" in Results / Fig 4 caption. Increasing $N_{\text{reps}}$ to 100 per cell is computationally feasible and queued for Round 3.

**R2.2 (Major)** — *PIP calibration error bars too wide; report Brier / ECE with bootstrap CIs.*
**Classification:** `[Accepted]` for ECE in Sup Fig S4 caption. Reformulating Sup Fig S4 with bootstrap CIs is queued for Round 3.

**R2.3 (Major)** — *No hyperparameter sensitivity analysis for HBP / GAFM.*
**Classification:** `[Accepted]`. Hyperparameter sensitivity Supplementary Figure (one-at-a-time $\pm 50\%$ perturbations) is queued as a new Sup Fig S7 in Round 3. Round-2 manuscript notes the gap explicitly in Methods.

**R2.4 (Major)** — *Restrict rice claim to 14/18 traits with $\lambda_{GC} \in [0.85, 1.20]$.*
**Classification:** `[Accepted]`. Results §"Rice (3,000 Rice Genomes...)" reformulated: primary claim restricted to 14 well-calibrated traits; the 4 high-$\lambda$ traits are now reported as a "secondary subset with residual subpopulation confounding". Sup Table S5 unchanged (all 18 traits remain) but main-text language tightened.

**R2.5 (Major)** — *Theorem 3 assumption statement informal.*
**Classification:** `[Accepted]`. Sup Note S1 Theorem 3 statement now explicitly conditions on $\max_{i \neq j, j \neq c} r^2_{ij} \leq \rho^*$ with $\rho^* = 1 - \alpha (1-\delta)$ for the eQTL-prior strength $\delta$.

**R2.minor (Wilcoxon vs sign-test, GraphGWAS-vs-PLINK2 within-species)** — `[Accepted]` (Wilcoxon now used in Methods §evaluation); `[Accepted-with-modification]` (within-species cross-validation already reported in Sup Table S4 footer for yeast at $r=1.0000$; we extend a similar note to Arabidopsis and rice as a Round-3 deliverable).

## Reviewer 3 — Software / Systems Reviewer

**R3.1 (Major)** — *17 GB graph dump is too large; provide streaming alternative.*
**Classification:** `[Accepted-with-modification]`. The Zenodo deposit will include both the full dump and a per-chromosome JSON shard set (chr1 alone is ~700 MB) so reviewers can verify functionality on a single chromosome before downloading the full 17 GB. Documented in `zenodo/graph_dumps/README.md` (Round 2 deliverable).

**R3.2 (Major)** — *Un-validated method entry points need explicit warnings.*
**Classification:** `[Accepted]`. CLGF, L4 and LPCE entry points will raise `UserWarning` at first invocation, citing the benchmark-status disclosure. Code change for Round 2.

**R3.3 (Major)** — *Single-command Fig 4 reproducibility.*
**Classification:** `[Accepted]`. Adding `tests/reproduce_fig4.sh` as a one-command driver. Round 2 deliverable.

**R3.4 (Major)** — *Pan-UKB local-cache mode.*
**Classification:** `[Accepted]`. `panukb.cache_locus()` documented in `docs/PANUKB_INTEGRATION_PLAN.md`; we add a fixture for BMI/chr22 to the test suite. Round 2 deliverable (code-side; not visible in manuscript text).

**R3.minor (3 items)** — Bibtex `fey2019fast`, README, lint sweep.
**Classification:** `[Accepted]`. Bibtex entry given a `pages = {n/a}` field; README, lint sweep, and `LICENSE`/`CITATION.cff` will be committed in the next push.

## Reviewer 4 — Adversarial Skeptic

**R4.1 (Critical)** — *"Heterogeneity vs density" is rhetorical; cCRE-flag is itself a form of density.*
**Classification:** `[Accepted]`. This is the cleanest and most damaging review point. We rewrite the abstract closing sentence and the Discussion paragraph to drop the "density" contrast and reframe in terms of *uniform vs non-uniform* priors. New abstract closing sentence:

> *"These results recast multi-omics fine-mapping as a non-uniform-prior-curation problem rather than a uniform-coverage problem, and reframe post-GWAS analysis as message passing over biological structure rather than weighted regression on flattened annotations."*

**R4.2 (Major)** — *"Same code" understates species-specific tuning.*
**Classification:** `[Accepted]`. Discussion §"Non-model species" now says: *"Within a fixed HBP/GAFM core, the same algorithm runs unchanged; the surrounding annotation graph, LD windows, validation catalogues and GWAS pipelines are tailored per species."*

**R4.3 (Major)** — *5–40× speed vs single-threaded baselines.*
**Classification:** `[Accepted-with-modification]`. SuSiE and FINEMAP were run with their default single-thread configuration (the configuration most users invoke). We add an explicit note in Methods §evaluation and Sup Table S2 footnote stating this. Multi-threaded baseline benchmarking is a Round-3 deliverable.

**R4.4 (Major)** — *LPCE in decision tree (Sup Fig S6) is a category error.*
**Classification:** `[Accepted]`. The LPCE node in Sup Fig S6 will be replaced with a "(forthcoming work)" marker in the figure regeneration. Code change for Round 2.

**R4.5 (Major)** — *5–40× excludes graph-cache + LD-load cost.*
**Classification:** `[Accepted-with-modification]`. We add an explicit "amortised vs first-locus" runtime breakdown to Sup Table S2 footnote. Round 3 will include a startup-cost row in the table.

**R4.minor (3 items)** — Lipid recovery cache configuration, 692-leads aggregation, title qualifier.
**Classification:** `[Accepted]` (1, 2); `[Rebutted-with-justification]` (3 — "under weak signal" is the reviewer-most-likely-to-resist claim and is the central methodological pitch; we keep the title).
**Action:** Abstract clarifies that *LDLR / APOE / LPL / GCKR / ANGPTL3* recoveries come from the genome-wide GAFM run (not the intervention cache); the 692-leads aggregation is now broken out by species in the abstract.

## Reviewer 5 — Reader-Advocate

**R5.1 (Major)** — *Acronym density too high.*
**Classification:** `[Accepted-with-modification]`. We add a one-paragraph "Abbreviations" box at the end of the Introduction listing the 12 most frequent acronyms with their expansions. Round 2 deliverable.

**R5.2 (Major)** — *Decision tree demoted to supplementary; readers will miss it.*
**Classification:** `[Rebutted-with-justification]`. The decision tree was deliberately moved to supplementary as part of the prior reorganisation that focuses Results on empirical findings only. We retain a short "method-choice" paragraph in Discussion §"When SuSiE and FINEMAP remain the right tool" that names which scenario maps to which method, so readers who skip the supplementary still see the recommendation.

**R5.3 (Major)** — *Table 1 caption + column headers not self-contained.*
**Classification:** `[Accepted]`. Table 1 reformulated: column headers expanded ("HBP narrows tighter than GAFM" instead of "HBP-tighter"); caption now defines "baseline cache" vs "intervention cache" inline.

**R5.4 (Major)** — *Inconsistent panel-count across main figures.*
**Classification:** `[Rebutted-with-justification]`. Main figures are sized to their content; forcing a uniform panel count would either pad or truncate. We retain the current sizing.

**R5.minor (3 items)** — Caption sub-panel formatting, orphans, long `\nameref{}` titles.
**Classification:** `[Accepted]` (1, 2); `[Accepted-with-modification]` (3 — Sup Note S2 title shortened).

## Reviewer 6 — Editor-in-Chief

We acknowledge the consolidated Major Revision verdict and have applied all in-text changes that do not require new computational experiments in Round 2. Three families of revisions are deferred to Round 3 with explicit caveats:

1. Placebo-prior control (R1.2 / R4.1) — full permutation experiment.
2. Hyperparameter sensitivity Sup Fig S7 (R2.3).
3. Null-overlap rate per species for the 250-kb validation criterion (R1.1).

For Round 2 advancement, we ask the reviewer panel to evaluate whether the textual reframing (heterogeneity → non-uniform priors), the abstract scope clarification (separate-experiment 27–2), the rice subset restriction (14 well-calibrated traits), the LPCE decision-tree removal, and the explicit limitations section are sufficient given that the three deferred experiments are explicitly named with target-of-completion dates.

— Author Team, on behalf of all five authors
