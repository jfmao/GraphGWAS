# Round 1 — 2026-04-28

Manuscript: *"Relational biological structure improves fine-mapping of causal GWAS variants under weak signal"* (Estaji, Zhao, Chen, Nie, Mao). v1, 52 pages, 6 main figs + 1 main table + 6 sup figs + 6 sup tables + 4 sup notes.

Target venue: **Nature Genetics** (Article).

---

## Reviewer 1 — Domain Expert (Statistical Genetics / Fine-Mapping)

**Summary:** The authors propose two fine-mapping methods (HBP, hierarchical belief propagation; GAFM, graph-augmented fine-mapping) that integrate multi-omics annotations as a typed factor graph rather than as a flat per-variant scalar weight, and benchmark them against the SuSiE/FINEMAP family plus the recent infinitesimal-extension methods (SuSiE-inf, FINEMAP-inf) and an SBayesRC-Polyfun-proxy comparison. A central empirical claim is a four-species cross-replication experiment (yeast, *Arabidopsis*, rice, human Pan-UKB) and a direct-intervention experiment showing that adding a 9-class regulatory-element flag flips the rate at which HBP narrows the credible set tighter than GAFM from 0% to 88% on 321 human leads.

**Major Concerns:**

1. **The "validation" criterion is too lax to support cross-species claims.** Methods §"Cross-species fine-mapping protocol" defines validation as "any catalog feature for its trait lies within 250 kb of the lead position on the same chromosome". A 250 kb window is the size of an LD block; in a gene-dense region (e.g. *Arabidopsis* chr 4, rice chr 1) a randomly placed lead has substantial prior probability of falling within 250 kb of *some* catalogued gene. Without a permutation-derived null overlap rate (per species), the "33% (rice) / 18% (yeast) / 6% (Arabidopsis)" validation rates cannot be interpreted as evidence of fine-mapping accuracy. **Action required:** report the null-overlap rate from random-permutation positions matched for chromosome length and gene density.

2. **The intervention experiment lacks a placebo control.** Sup Note S2 / Methods §"Heterogeneity-by-intervention experiment" reports that adding a 9-class ENCODE cCRE prior_score flips the HBP-tighter-than-GAFM rate from 0% to 88% on 321 human leads. But the only contrast is "no prior score" vs "informative prior score". A reviewer cannot distinguish "informative prior breaks LD ties" from "any heterogeneous prior breaks ties pseudo-arbitrarily". **Action required:** add a placebo arm in which `prior_score` is permuted (same marginal distribution, randomised assignment to variants) and report the HBP-tighter rate. If placebo is also high, the claim "per-variant heterogeneity, not annotation density, is the operative quantity" collapses to "any randomly-distributed prior breaks ties".

3. **The 27–2 weak-signal headline is from a different experiment than the 6-method head-to-head.** Results §"GAFM wins weak signal 27--2 against SuSiE" reports 79 replicates on tissue-specific eQTL causal variants ($\beta=0.15$, $h^2=0.01$). The uniform 6-method × 3-scenario benchmark (Fig 4 / Sup Table S2) uses 30 reps per scenario and different effect sizes. These are not comparable, and the abstract claim "wins 27--2 against SuSiE at weak signal" combined with the Fig 4 caption "match SuSiE rank-#1 rate within three percentage points" is internally inconsistent. **Action required:** either fold the 27–2 experiment into the unified benchmark (re-run on the same 30-replicate × 3-scenario grid with the eQTL-weighted condition added), or clarify in abstract + Results that 27–2 is from a separate, more-favourable simulation regime.

4. **HBP softmax ceiling at ~0.7 PIP is presented as a feature but is operationally a problem.** Methods §HBP and Sup Fig S4 caption note that HBP "caps at $\approx 0.7$" because the softmax distributes mass across the locus. For a fine-mapping method whose primary deliverable is a credible set with calibrated PIP, capping the maximum PIP at 0.7 means HBP can never report PIP=1 even when the causal variant is unambiguous. This is presented as conservative calibration, but Pan-UKB results that quote "PIP=1.000" come from GAFM, not HBP. The two methods are not interchangeable in the way the abstract implies. **Action required:** clarify in abstract and Results §"PIPs are calibrated..." that HBP and GAFM report different posterior scales; explain when each PIP can be interpreted at face value.

**Minor Concerns:**

- Sup Note S4 §"Validation criterion" uses American "validated" while body mostly uses British spelling ("recognised", "modelled", "centred"). One spelling convention.
- Theorem 1 (LPCE search-space reduction) is stated without proof in Sup Note S1; the proof itself sketches the upper bound but does not show tightness. Either prove tightness or label "upper bound only".
- In Results §"HBP fine-maps Pan-UK Biobank...": the EUR sample size $N=420{,}531$ is dominant; $N=8{,}876$ (CSA), $6{,}636$ (AFR), $2{,}709$ (EAS) introduces sample-size confounding with ancestry. A reviewer reading "HBP fine-mapping holds across four ancestries" expects either ancestry-specific power-matched analyses or an explicit caveat.

**Verdict:** Major Revision

---

## Reviewer 2 — Methods / Statistics Reviewer

**Summary:** Two new fine-mapping methods compared against five Bayesian baselines on a uniform 6-method × 3-scenario × 30-replicate benchmark on 1KG chr22 simulations, plus a 4-species cross-replication and a 692-lead intervention experiment.

**Major Concerns:**

1. **30 replicates per scenario is underpowered for the rank-1 rate comparisons in Fig 4.** Differences of 3 percentage points (e.g. SuSiE 70% vs HBP 70% at strong; SuSiE-inf 77% vs FINEMAP 73%) are within Wilson-binomial 95% CIs of about ±15%. The paper's "match within three percentage points" claim is consistent with the data but the data does not statistically support the *equivalence* claim — only an absence-of-evidence-of-difference. **Action required:** report 95% confidence intervals on rank-1 rates; use exact two-tailed binomial tests for pairwise comparisons; downgrade "match" to "no significant difference at $\alpha=0.05$".

2. **Sup Fig S4 (PIP calibration) reports 200 simulations × 4 $h^2$ levels but the diagonal-tracking claim uses only 8 PIP bins × 4 methods = 32 measurements per method.** Each bin has ~25 simulations. Calibration error bars per bin are wide. The current panel cannot distinguish "well-calibrated" from "approximately right" for any method. **Action required:** report Brier scores or expected calibration error (ECE) with bootstrap CIs per method; at minimum, error bars on the per-bin TDR estimates.

3. **The free hyperparameters of HBP ($\alpha=0.6$, $\lambda=0.5$, $T=5$) and GAFM ($\alpha=0.5$, $\tau_L=0.3$) are reported as defaults but no sensitivity analysis is provided.** A reviewer cannot verify that the headline numbers are robust to ±10% perturbations of these parameters. **Action required:** add a Supplementary Figure with HBP/GAFM rank-1 rate vs. each hyperparameter (one-at-a-time) on the strong-signal scenario.

4. **The IRRI 18-trait scan (Sup Table S5) reports $\lambda_{GC}$ ranging from 0.50 to 2.07.** Four traits exceed $\lambda_{GC} > 1.5$, suggesting unconfounded population stratification despite 10 PCs. The paper acknowledges this in Results §"Rice ..." but proceeds to fine-map. **Action required:** restrict the cross-species fine-mapping conclusions to the 14 traits with $\lambda_{GC} \in [0.85, 1.20]$; treat the four high-$\lambda$ traits separately or as a cautionary subset.

5. **Theorem 3 (GAFM causal-variant ranking) requires "mild LD-decay assumptions and bounded cross-LD between non-causal variants" — these are not made explicit.** A statement like "the causal variant ranks first" is operationally meaningless without the precise condition. **Action required:** state the assumption formally (e.g. $\max_{i \neq j, j \neq c} r^2_{ij} \leq \rho^*$ for some $\rho^*$ that depends on the eQTL-prior strength).

**Minor Concerns:**

- Sign-test $p$-values in the head-to-head (e.g. "$p < 10^{-5}$") are reported but the test ignores ties (50 ties out of 79 replicates). Use a Wilcoxon signed-rank or paired-rank test that handles ties appropriately.
- Methods §"Single-locus association" describes three different software pipelines (internal `graphgwas.assoc`, PLINK 2.0, Pan-UKB SAIGE sumstats). Cross-species results are sensitive to this choice. Add a within-species cross-validation showing the GraphGWAS internal pipeline reproduces the PLINK 2.0 z-statistics on a shared variant subset.

**Verdict:** Major Revision

---

## Reviewer 3 — Software / Systems Reviewer

**Summary:** A graph-database-backed fine-mapping platform with a Python package, command-line interface, multi-omics annotation graph (Neo4j 5.26 reference implementation), and pre-built dumps deposited on Zenodo. Code released at `github.com/jfmao/GraphGWAS` under MIT.

**Major Concerns:**

1. **The pre-built human graph dump is 17 GB.** Zenodo's per-file limit is 50 GB but practical reviewer download time is real. **Action required:** provide a streaming-import alternative (e.g. APOC `apoc.load.json` from per-chromosome JSON shards) so reviewers can verify a subset without downloading 17 GB.

2. **"Implemented but not benchmarked" methods (CLGF, L4, LPCE) are released alongside HBP/GAFM in the same package.** The release manifest must clearly mark un-validated methods at the API level so a downstream user does not invoke `clgf_finemap()` thinking it is at the validation depth of `hbp_finemap()`. **Action required:** raise a `UserWarning` at first invocation of CLGF / L4 / LPCE entry points, citing Sup Note S2 §benchmark-status.

3. **Reproducibility of Fig 4 requires running three separate benchmark scripts** (`tests/benchmark_hbp_vs_susie_finemap.py`, `tests/benchmark_inf_methods_h2h.py`, `tests/benchmark_polyfun_proxy.py`) and merging their outputs in `tests/generate_paper_figures.py:figure_3`. A single command should reproduce the figure end-to-end. **Action required:** add a `tests/reproduce_fig4.sh` that runs the three benchmarks in sequence and regenerates the figure.

4. **The Pan-UKB sumstats fetch path uses tabix-over-HTTPS against `s3://pan-ukb-us-east-1/`.** This is a public read-only bucket but bandwidth is non-trivial (~50 kB-1 MB per locus, 321 leads). A reviewer rerunning the chr22 fine-mapping must download data each invocation. **Action required:** document a local-cache mode (e.g. `panukb.cache_locus()`) and provide pre-cached test fixtures for at least one trait (BMI / chr22) in the test suite.

**Minor Concerns:**

- The bibtex entry `fey2019fast` (PyTorch Geometric) throws a "You can't pop an empty literal stack" error during `bibtex` compilation; the PDF still renders correctly because of cached `.bbl`. Add a `pages` field or reformat the entry.
- The README at `github.com/jfmao/GraphGWAS` is currently un-staged and was not pushed in the v0.1.0 release commit. Public repo without a README impairs first-impression accessibility.
- 28 Python files (`src/python/graphgwas/*.py`) have unstaged ruff-style cleanups (removing unused imports, fixing `f""` warnings). These should be committed before tagging v0.1.0.

**Verdict:** Major Revision

---

## Reviewer 4 — Adversarial Skeptic

**Summary:** A paper that claims its central methodological reframing ("heterogeneity, not density, is the operative quantity") is established by a single intervention experiment.

**Major Concerns:**

1. **The "heterogeneity vs density" distinction is rhetorical rather than empirical.** The intervention adds a 9-class cCRE flag — by definition this *is* a form of density (every variant gets classified into one of 9 classes). The paper's claim "density alone is insufficient" is then refuted by a "denser annotation" intervention, which is internally contradictory. The claim that should actually be made is: "uniform annotations break no LD ties; non-uniform annotations break ties in proportion to their per-variant discriminative information." This is a reasonable claim but it is not the claim the abstract makes. **Action required:** rewrite the abstract closing sentence to drop "density" as a contrast and reframe in terms of uniform-vs-heterogeneous priors.

2. **"Same code fine-maps yeast / Arabidopsis / rice / human" understates what is actually species-specific.** Methods §"Cross-species fine-mapping protocol" lists per-species: different annotation graphs, different LD windows (15 kb yeast / 100 kb Ath / 250 kb rice / 500 kb human), different validation catalogues, different PCs, different GWAS pipelines. The "no algorithmic change" claim is true at the level of the HBP/GAFM core only. The infrastructure around it is heavily species-tuned. Abstract / Discussion should be rewritten to reflect this.

3. **"5–40× speed advantage" is benchmarked against single-threaded baselines on a single workstation.** SuSiE (R) and FINEMAP (binary) can both be parallelised across loci. The paper does not say whether the baselines were given multi-threaded execution. Even on a single thread, modern SuSiE versions are faster than the 1.0–1.8 s reported here. **Action required:** confirm baseline configurations match each tool's recommended-for-production setup, and report baseline throughput at multi-thread parallelism.

4. **Authors include LPCE in the recommended-method decision tree (Sup Fig S6) while Methods has been moved to Sup Note S2 with explicit "preview, under development, not benchmarked here" caveat.** Reviewers and downstream users will not parse this nuance. The decision tree visually recommends LPCE for "epistasis discovery" alongside benchmarked methods (HBP, GAFM, SuSiE-inf, etc.), which is a category error. **Action required:** either remove LPCE from the decision tree entirely, or replace its node with a "(forthcoming work)" marker.

5. **The per-locus 0.07–0.12 s runtime for HBP/GAFM excludes graph-cache load and LD-reference load.** Methods §"Evaluation and hardware" says so explicitly. But for any production workflow, cache+LD load is the dominant cost on the first locus and amortises poorly across small batches. The reported speedup is thus the steady-state cost in a large batch. A first-time user fine-mapping a single locus will not see 5–40× over SuSiE. **Action required:** report end-to-end (cache load + first locus + 30th locus) wall-clock comparison.

**Minor Concerns:**

- Abstract claims "recovering *LDLR*, *APOE*, *LPL*, *GCKR*, *ANGPTL3* at single-variant resolution" — five canonical lipid genes is a strong rhetorical card. But the paper also reports HBP narrowing *only* on 0% of human chr22 leads at baseline cache. The 5 named recoveries are at PIP $\geq 0.65$ from the *intervention* cache (cCRE-augmented), not the baseline. The abstract should clarify which configuration produced the canonical recoveries.

- The "692 leads" total is the sum across yeast (245) + *Arabidopsis* (54) + rice (72) + human (321). Pooling species into a single integer is misleading because the species are not interchangeable. Report 245 + 54 + 72 + 321 separately when discussing the intervention rate.

- Title: *"Relational biological structure improves fine-mapping of causal GWAS variants under weak signal."* The "under weak signal" qualifier is true for the 27–2 GAFM-vs-SuSiE result, but is not the main contribution of the paper (which is the cross-species + intervention). The title undersells the breadth.

**Verdict:** Major Revision

---

## Reviewer 5 — Reader-Advocate (Clarity / Narrative / Figure Quality)

**Summary:** A methodologically dense paper that makes an interesting empirical claim (heterogeneity intervention) but is challenging to read end-to-end due to acronym proliferation, embedded supplementary cross-references, and a heavily technical Methods section.

**Major Concerns:**

1. **Acronym density is too high for a *Nature Genetics* general-readership Article.** First-pass reading load: HBP, GAFM, LPCE, CLGF, L4, FPR, MAF, PIP, eQTL, PPI, cCRE, S-LDSC, BIOGRID, IRRI, RAP-DB, AraGWAS, SGD, GO, FarmGTEx, Pan-UKB, GIANT, GLGC, NYGC, GRAMMAR+, REML, GRM, FINEMAP-inf, SuSiE-inf, MCMC. NG readers from outside fine-mapping will struggle. **Action required:** either reduce the count (collapse method-stack acronyms to plain English in main text) or include a glossary box on page 1 / page 2.

2. **The decision tree (now Sup Fig S6) was the most reader-friendly visual artefact in the paper and has been demoted to supplementary.** Many readers (clinicians, postdocs choosing a method) will skip the supplementary and miss it. **Action required:** keep at least a stripped-down version of the decision tree in the main text (perhaps as a panel of Fig 7 — but Fig 7 is now removed). Alternative: add a one-paragraph "method choice" box in the Discussion.

3. **The intervention experiment (Table 1) is the strongest single result in the paper but its caption requires reading the body to interpret.** Table 1 title is *"Adding a per-variant heterogeneous prior to the multi-omics cache turns HBP from inert to informative across four species"*, but the table itself uses jargon ("Baseline cache HBP-tighter", "+heterogeneous prior HBP-tighter"). A reader scanning only Table 1 cannot interpret the numbers. **Action required:** restructure Table 1 caption + column headers to be self-contained (define "HBP-tighter" inline, define "baseline" vs "intervention" inline).

4. **Figure 1 (architecture) has two panels (a, b); Fig 2 is a single schematic; Fig 6 is single-panel; Fig 7 is single-panel.** Mixing 4-panel and single-panel figures across the main set looks visually inconsistent. Some single-panel "figures" (e.g. Fig 2 HBP schematic, Fig 6 ablation) might be better as boxed insets or moved to supplementary, leaving 4-panel figures as the main set.

**Minor Concerns:**

- Many figure captions use `\textbf{...}` for the topic-sentence claim and then plain prose for panel descriptions. Consistent — good. But within-caption sub-panel labels alternate between `\textbf{a}` and just `a` — pick one.

- Several pages have orphaned single sentences at the bottom — a typesetting issue. Use `\widowpenalty` and `\clubpenalty` adjustments, or `\flushbottom`.

- The `\nameref{}` cross-references (e.g. *"see Discussion 'Non-model species — crops, livestock, and model organisms — removing the annotation-weight bottleneck'"*) produce extremely long inline section names. Consider shortening section titles or using numeric cross-references.

**Verdict:** Minor Revision

---

## Reviewer 6 — Editor-in-Chief (synthesises Reviewers 1–5 and renders venue verdict)

**Summary of consolidated reviewer concerns** (in priority order, deduplicated across reviewers):

| Priority | Issue | Source reviewers | Block status |
|---|---|---|---|
| **Critical** | Intervention experiment lacks a permutation-control / placebo-prior arm; "heterogeneity vs density" framing is contradicted by the cCRE-flag (which is itself a density measure) | R1.2, R4.1 | **Blocker** |
| **Critical** | 250-kb-window validation lacks species-specific null permutation rate; cross-species "validation" rates uninterpretable | R1.1 | **Blocker** |
| **Critical** | 27–2 weak-signal headline is from a separate experiment than the unified 6-method × 3-scenario benchmark; abstract overclaims "match SuSiE within 3pp" while elsewhere claiming "27–2 wins"; need internally consistent claims | R1.3, R2.1 | **Blocker** |
| **Major** | No HBP/GAFM hyperparameter sensitivity analysis | R2.3 | Major-revision required |
| **Major** | No power calculation / CIs on rank-1 rate; 30 reps per cell underpowered | R2.1, R2.2 | Major-revision required |
| **Major** | Sample-size confounding with ancestry in Pan-UKB results | R1.minor + R4.implicit | Caveat or ancestry-specific power-match |
| **Major** | LPCE in decision tree (Sup Fig S6) is a category error given preview-only status | R4.4 | Decision-tree edit |
| **Major** | Speed comparison vs. single-threaded baselines | R4.3 | Configuration disclosure |
| **Major** | Abstract / overclaim about "single-variant resolution lipid recovery" — actually requires the intervention cache | R4.minor | Abstract qualification |
| **Minor** | Acronym density and reader-advocacy concerns | R5.1, R5.2 | Polish |
| **Minor** | Bibtex compilation glitches, README unstaged, lint cleanup unstaged | R3.minor, R3.3 | Pre-submission housekeeping |
| **Minor** | Theorem 3 assumption statement; Wilcoxon vs sign-test | R1.minor, R2.5, R2.minor | Tighten in proof note |

**Editor verdict:** **Major Revision.**

**Rationale:** The paper has a clear and interesting central claim (the heterogeneity intervention) and substantial empirical scope (4-species replication, 6-method benchmark, Pan-UKB cross-ancestry demonstration). However, three independent reviewers (Domain Expert R1, Methods R2, Adversarial Skeptic R4) converge on the same critical concern: the central intervention claim ("heterogeneity, not density") is not robust to the obvious null comparison (placebo prior with same heterogeneity but no biological signal). Without that control, the headline reframing — and therefore the contribution — is not yet established at the level *Nature Genetics* requires. The 250-kb-window validation rates and the 27–2-vs-Fig-4 inconsistency are independent additional blockers. Reviewer 5 (Reader-Advocate) returned Minor Revision; Reviewers 1–4 returned Major Revision; consolidated verdict is **Major Revision**.

**Specific action items required for Round 2 advancement:**

(R1.2 + R4.1) Add a placebo-prior control arm to the intervention experiment. Permute `prior_score` across variants (preserving marginal distribution, randomising assignment) and rerun the 692-lead intervention rate. Report the (informative, placebo, none) triplet rates per species. If placebo rate ≥ 50%, reframe the abstract claim accordingly.

(R1.1) For each species' validation rate, generate a null-overlap rate by drawing 1000 random "leads" from each species' chr/position distribution and computing the same 250-kb-window overlap with the same catalogue. Report empirical p-value of the observed validation rate vs. null.

(R1.3 + R2.1) Reconcile the 27–2 weak-signal claim with the unified 6-method × 3-scenario benchmark. Either rerun a 30-rep eQTL-restricted scenario with all 6 methods, or label the 27–2 result clearly as a separate experiment with different design.

(R2.3) Add a hyperparameter sensitivity Supplementary Figure (one-at-a-time perturbation of $\alpha$, $\lambda$, $T$ for HBP; $\alpha$, $\tau_L$ for GAFM).

(R2.4) Restrict the cross-species rice claim to the 14 / 18 IRRI traits with $\lambda_{GC} \in [0.85, 1.20]$; treat the 4 high-$\lambda$ traits separately.

(R4.4) Remove LPCE from the decision-tree node in Sup Fig S6 (or replace with "(forthcoming work)" marker).

(R5.3) Restructure Table 1 caption + column headers to be self-contained.

**Status:** Continue to Round 2 (revision required).
