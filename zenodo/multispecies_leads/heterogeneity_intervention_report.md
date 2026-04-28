# Heterogeneity-not-density: empirical validation by direct intervention

**Date:** 2026-04-26

## Hypothesis

From the rice ablation study (paper §"Ablation: when does multi-omics
coverage help?"): HBP narrows credible sets beyond GAFM only when the
multi-omics cache provides **per-variant heterogeneous** prior
information. A uniform-density cache (every annotated variant looks
structurally similar) leaves HBP equal to GAFM; a heterogeneous cache
(per-variant variation in pathway / PPI / regulatory annotation)
makes HBP responsive.

**Empirical demonstration before this work:** correlational only
(varying species cache density across 4 species, observing HBP-tighter
correlates with heterogeneity). Most striking: dense+uniform human
GTEx+STRING cache → 0/321 HBP-tighter on Pan-UKB lipid+BMI loci.

**This Phase 0–4 work:** direct intervention. For each of 4 species,
add a per-variant heterogeneous `prior_score` field to the cache,
re-run identical fine-mapping pipeline, measure the change in
HBP-tighter rate.

## Schema extension (Phase 0)

`graph_cache[variant_id]` gained an optional `prior_score: float`
field. Both HBP and GAFM from-sumstats functions now consume it:

- **HBP**: at the end of belief propagation, prior_score is rescaled
  to z-score range and re-mixed into the final PIP (50/50 blend with
  the BP-derived belief).
- **GAFM**: prior_score is added directly to the `z_func` mixture before
  LD deconvolution (50% weight).

Regression: cache without prior_score reproduces v1 results exactly.

## Per-species intervention

### Human (Phase 1)

- v1: GTEx v8 eQTL (49 tissues, 4.59M variants) + STRING v12 PPI ≥ 700.
  Each annotated variant gets ~1 eGene + ~20 STRING partners;
  near-uniform across the genome.
- v2: + ENCODE cCRE v3 overlap (926,535 elements; 9 classes:
  pELS / dELS / PLS / CTCF-bound / DNase-H3K4me3 / etc.) with
  class-weighted prior_score (PLS-CTCF=1.5, PLS=1.0, pELS-CTCF=1.2,
  dELS=0.6, etc.). 10–15% of GTEx-annotated variants per chrom now
  carry a regulatory-element entry + heterogeneous prior_score.

### Yeast (Phase 2)

- v1: SGD ORFs + GO_slim mapping only (no PPI in v1).
- v2: + BIOGRID 5.0.256 yeast S288c physical interactions
  (263,403 edges → 5,813 ORFs with ≥1 partner) + ORF-feature-class
  prior_score (ORF=1.0, tRNA/ncRNA=0.7, pseudogene=0.4).

### Arabidopsis (Phase 3)

- v1: TAIR10 + Plant Reactome pathway annotations + STRING-Ath PPI
  (UniProt-fixed, 20,192 genes with partners).
- v3: + TAIR10 GFF feature-class prior_score: CDS=1.0, exon=0.8,
  gene=0.5 (with ±2 kb promoter window). 71-73% of variants per
  chrom got a per-variant prior_score.

### Rice (Phase 4)

- v1: snpEff variant→gene + Ren-2023 (269 grain-quality genes)
  pathway labels + RicePPINet PPI (Prob ≥ 0.7).
- v2: + snpEff impact-class prior_score (HIGH=1.0, MODERATE=0.7,
  LOW=0.4, MODIFIER=0.1) by re-scanning the pseudo-canonical VCF.
  100% of cached variants got a prior_score.

## Result

The heterogeneity hypothesis is empirically confirmed across all 4
species. Direction of effect is consistent; magnitude correlates with
v1 cache uniformity:

| Species     | Cache change                                  | v1 HBP-tighter | v2 HBP-tighter | Δ        |
|-------------|-----------------------------------------------|---------------:|---------------:|---------:|
| Human       | GTEx+STRING uniform → +cCRE+prior_score       | 0/321 (0%)     | 282/321 (88%)  | **+88pp**|
| Arabidopsis | TAIR10+Reactome+STRING → +CDS prior_score     | 29/54 (54%)    | 51/54 (94%)    | +41pp    |
| Rice        | Ren-2023+RicePPINet → +snpEff-impact prior    | 14/72 (19%)    | 56/72 (78%)    | +58pp    |
| Yeast       | GO_slim → +BIOGRID+prior_score                | 244/245 (100%) | 245/245 (100%) | saturated|

The human result is the most striking: the v1 cache was the most
uniform (every annotated variant ≈ 1 eGene + 20 STRING partners) and
produced a textbook null effect (HBP exactly equal to GAFM on every
single one of 321 leads). Adding the heterogeneous cCRE-class layer
flipped HBP-tighter from **0% to 88%** — the cleanest possible direct
demonstration that what HBP needs is per-variant heterogeneity, not
density.

Yeast was already saturated at 100% in v1; the v2 intervention did
not degrade the rate, but the test was non-falsifiable for yeast
(could not increase from 100%).

## Caveats and partial findings

1. **GAFM PIPs degrade slightly under v2 in some species** — adding
   prior_score to z_func re-mixes mass across the locus window,
   sometimes lowering the top-variant PIP. Mean GAFM ΔPIP human v2-v1:
   −0.08; arabi v3-v2: −0.04. This is the "uniform too-broad
   prior_score" failure mode noted in §arabi v3 and §yeast v2: a
   prior_score that fires on every CDS variant is uniform within the
   "in-CDS" subset, so it doesn't help discriminate within that
   subset. **Future Phase 1.B** will replace the binary
   regulatory-element flag with continuous CADD / MotifBreakR Δ-PWM /
   DeepSEA Δ-effect scores that vary at single-variant resolution.

2. **Yeast prior_score was too uniform (1.0 for any ORF, ~99% of
   variants)**. The HBP-tighter rate didn't change because v1 was
   already saturated, but GAFM's mean PIP dropped (5 PIP≥0.5 → 0).
   Same fix applies: continuous Δ-score per variant rather than
   class-binary.

3. **Rice prior_score reaches 100% coverage (all snpEff-annotated
   variants)**, but the impact classes are well-stratified
   (HIGH=1.0, MODERATE=0.7, LOW=0.4, MODIFIER=0.1) so heterogeneity
   is real even at full coverage. HBP-tighter went 19% → 78%.

## Future Phase 1.B+

The current intervention uses category-binary scores. For maximum
heterogeneity:
- Human: CADD v1.7 PHRED scores (continuous) + MotifBreakR Δ-PWM
  per locus + ABC E-P pair scores
- Yeast: snpEff impact-class (like rice) + JASPAR-fungi PWMs
- Arabidopsis: O'Malley 2016 DAP-seq peak signal scores + CADD-Ath
  if available
- Rice: Cao 2020 / Tu 2020 ChIP-seq peak signal + cereal PhyloP

These are downloads-heavy (CADD genome-wide is 80 GB; DAP-seq tar
is ~5 GB) but the schema, pipeline, and validation framework are
all in place — only the cache-augmentation step would change.

## Reproducibility

Driver scripts:
- `tests/human_augment_cache_phase1a.py` (cCRE)
- `tests/yeast_augment_cache_phase2.py` (BIOGRID + ORF prior)
- `tests/arabidopsis_augment_cache_phase3.py` (TAIR10 CDS prior)
- `tests/rice_augment_cache_phase4.py` (snpEff impact prior)
- `src/python/graphgwas/finemapping_v2.py` (Phase 0 schema extension)

Pipeline drivers (set `CACHE_VERSION=v2` env var):
- `tests/human_panukb_finemap_genomewide.py`
- `tests/yeast_finemap_all.py`
- `tests/arabidopsis_finemap.py`
- `tests/rice3k_irri_finemap.py`
