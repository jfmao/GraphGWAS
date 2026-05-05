# `graphgwas finemap-sumstats` — sumstats-only fine-mapping (v0.1.5)

## Description

Sumstats-only entry path: runs GAFM, HBP, or one of the v0.1.5 mixture-prior
variants (GAFM-MX, HBP-MX, ENS) on a pre-computed locus given a per-variant
z-score and an LD matrix. No Neo4j, BGEN, or genotype panel required.
Useful for:

- Pan-UKB and other public summary-statistics releases.
- Pre-computed in-sample LD matrices stored as `.npz` files.
- Embedding fine-mapping inside larger pipelines that already produce
  z-scores and LD.

## Methods

| Method | Description |
|---|---|
| `gafm` | Graph-Augmented Fine-Mapping (single-Gaussian prior on the LD-deconvolved statistic; α-blended functional prior) |
| `hbp` | Hierarchical Belief Propagation on a variant→gene→pathway factor graph |
| `gafm-mx` | **v0.1.5**: GAFM + λ_GC deflation + SBayesRC-style 4-component Wakefield mixture-BF posterior reweighting on the LD-deconvolved z-scores |
| `hbp-mx` | **v0.1.5**: HBP + same deflation + mixture-prior posterior |
| `ens` | **v0.1.5**: Mean-of-PIPs ensemble of GAFM-MX and HBP-MX |

The mixture prior uses SBayesRC's 4-component defaults: π = (0.005, 0.003,
0.001, 0.001) (renormalised within non-zero components) and γ = (0.001, 0.01,
0.1, 1.0). The mixture BF is applied to LD-deconvolved unique_stats, not raw
marginal z, so GAFM/HBP's LD-aware ranking is preserved during reweighting.

## Usage

```bash
graphgwas finemap-sumstats \
    --variants-tsv variants.tsv --z-tsv z.tsv --ld-npz ld.npz \
    --n-samples N --method gafm-mx [--lambda-gc LAMBDA] \
    --out-tsv credible_set.tsv
```

## Inputs

| Argument | Description |
|---|---|
| `--variants-tsv` | TSV with columns `variant_id, chr, pos, ref, alt, af` — one row per variant in the locus window. |
| `--z-tsv` | TSV with columns `variant_id, z`, where z = BETA / SE. |
| `--ld-npz` | `.npz` file with key `R_sq` = n × n squared-correlation matrix in the same row order as `variants-tsv`. |
| `--n-samples` | GWAS sample size (used to scale the mixture variance V_k = N · γ_k). |
| `--method` | One of `gafm, hbp, gafm-mx, hbp-mx, ens`. Default: `gafm-mx`. |
| `--lambda-gc` | Genomic-control inflation factor for z-deflation (z → z/√λ_GC). Default: no deflation. Use the trait's actual λ_GC for inflated panels (e.g. 1.4–1.8 on the rice 3kRG XI/GJ panel). |
| `--alpha` | GAFM functional/statistical blend; default 0.7 (regime-specific default; use 0.9 for strong-signal panels). |
| `--coverage` | Credible-set coverage; default 0.95. |

## Output

TSV with columns `variant_id, pip, in_credible_set, rank` sorted by descending PIP.

## Calibration caveats

The mixture-prior methods (GAFM-MX, HBP-MX, ENS) are **anti-conservative at
high PIP bins**: empirical TDR ≈ 0.62 in the [0.9, 1.0] PIP bin vs the
diagonal-target ≈ 0.95 (50 reps × 3 h² on chr22 single-causal F1 sims). PIPs
should be interpreted as an operational ranking score (sharper credible sets,
higher rank-1 rate) rather than as posterior probabilities. For applications
needing calibrated PIPs, run the base `gafm` or `hbp` method (or SuSiE) and
report those PIPs alongside.

Null FPR (100 H0 loci, y = ε):
- GAFM, HBP: 0/100 at PIP ≥ 0.5
- SuSiE / SuSiE-inf / FINEMAP-inf: 1/100 at PIP ≥ 0.5
- GAFM-MX, ENS: 6/100 at PIP ≥ 0.5
- HBP-MX: 10/100 at PIP ≥ 0.5
- All methods ≤ 1/100 at PIP ≥ 0.9

## Performance

Median wall-clock per locus on a single CPU thread:

| Method | Time / locus |
|---|---|
| GAFM | 0.02–0.08 s |
| HBP | 0.02–0.08 s |
| GAFM-MX, HBP-MX | 0.03–0.05 s |
| ENS | 0.05–0.08 s |
| SuSiE (R subprocess) | 5–60 s |
| SuSiE-inf, FINEMAP-inf | 1–10 s |
| SBayesRC (per-trait MCMC) | ~10 min genome-wide |

## Example

```bash
# Fine-map the rice 3kRG GS3 locus (Chr3:16,733,441) for grain weight (TGW)
# at λ_GC = 1.41 (trait inflation factor)
graphgwas finemap-sumstats \
    --variants-tsv data/rice_3k/results/grain_finemap/Chr3_GS3_variants.tsv \
    --z-tsv data/rice_3k/results/grain_finemap/Chr3_GS3_z.tsv \
    --ld-npz data/rice_3k/results/grain_finemap/Chr3_GS3_R_sq.npz \
    --n-samples 1847 --method ens --lambda-gc 1.41 \
    --out-tsv results/Chr3_GS3_ens_credible_set.tsv
```

## Reproducibility

- Cross-species v0.1.5 benchmarks: `tests/yeast_finemap_all.py`,
  `tests/arabidopsis_finemap.py`, `tests/rice3k_irri_finemap.py`,
  `tests/rice3k_grain_shape_finemap.py`.
- chr22 head-to-head simulations: `tests/benchmark_v15_chr22.py`.
- Calibration + null FPR: `tests/benchmark_v15_calibration_null.py`.
- Hyperparameter sensitivity: `tests/benchmark_v15_hyperparam_sensitivity.py`.

See [`docs/REPRODUCIBILITY.md`](../../REPRODUCIBILITY.md) §6.3 for full
reproduction commands.

## See also

- [`graphgwas finemap`](finemap.md) — Neo4j / BGEN / Pan-UKB fine-mapping
  with the GAFM, HBP, CLGF, and GLEM method classes
- [`vignettes/fine-mapping-quickstart.md`](../../../vignettes/fine-mapping-quickstart.md)
- [`docs/MATHEMATICAL_PROOFS.md`](../../MATHEMATICAL_PROOFS.md) — theorems for
  HBP, GAFM, and CLGF
