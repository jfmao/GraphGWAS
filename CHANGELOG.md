# Changelog

All notable changes to the GraphGWAS package are documented here.
The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and the project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.5] — 2026-05-04

### Added
- **`graphgwas.finemapping_v2.deflate_z_for_lambda_gc(z, lambda_gc)`**:
  standard genomic-control z-score deflation for inflated panels.
- **`graphgwas.finemapping_v2.apply_mixture_posterior(pips, z, n_samples, ...)`**:
  SBayesRC-style 4-component Wakefield mixture-prior posterior reweighting.
  Defaults: `pi=(0.005, 0.003, 0.001, 0.001)`, `gamma=(0.001, 0.01, 0.1, 1.0)`.
- **`graphgwas.finemapping_v2.credible_set_from_pips(pips, coverage)`**:
  standalone 95% credible-set helper.
- **`graphgwas.finemapping_v2.gafm_mx_from_sumstats(...)`** (paper-facing
  name: GAFM-MX): GAFM + λ_GC deflation + LD-deconvolved mixture posterior.
- **`graphgwas.finemapping_v2.hbp_mx_from_sumstats(...)`** (paper-facing
  name: HBP-MX): HBP + λ_GC deflation + LD-deconvolved mixture posterior.
- **`graphgwas.finemapping_v2.ensemble_from_sumstats(...)`** (paper-facing
  name: ENS): mean-of-PIPs ensemble of GAFM-MX and HBP-MX.
- New CLI subcommand `graphgwas finemap-sumstats --method
  {gafm,hbp,gafm-mx,hbp-mx,ens}` for sumstats-only fine-mapping.
- Unit tests in `tests/test_finemap_v15_mixture.py` (9 tests).
- Cross-species v0.1.5 benchmarks:
  `tests/benchmark_v15_chr22.py` (single-causal F1 across `h² ∈ {0.02, 0.05, 0.10}`),
  `tests/benchmark_v15_calibration_null.py` (PIP calibration + null FPR),
  `tests/benchmark_v15_hyperparam_sensitivity.py` (mixture π/γ robustness).

### Changed
- The mixture-prior reweighting in `gafm_mx_from_sumstats` and
  `hbp_mx_from_sumstats` is applied to the **LD-deconvolved**
  `unique_stats` (output of `_ld_deconvolve`), not to the raw marginal
  z-scores. This preserves GAFM/HBP's LD-aware ranking; otherwise
  LD-correlated noise variants with marginally larger marginal `|z|`
  outcompete the causal during multiplicative reweighting.
- `tests/rice3k_grain_shape_finemap.py` now calls the library
  `gafm_mx_from_sumstats` / `hbp_mx_from_sumstats` /
  `ensemble_from_sumstats` rather than its previous inline mixture
  posterior.
- IRRI, yeast, and Arabidopsis fine-map scripts (`rice3k_irri_finemap.py`,
  `yeast_finemap_all.py`, `arabidopsis_finemap.py`) extended to report
  GAFM-MX, HBP-MX, ENS alongside base GAFM/HBP.
- Pan-UKB benchmark (`benchmark_panukb_finemap.py`) extended to report
  v0.1.5 mixture-prior methods.

### Fixed
- **CRITICAL** alignment bug in `gafm_mx_from_sumstats` and
  `hbp_mx_from_sumstats`: `_build_candidates` returns the candidate list
  sorted by `combined_score`, not in input variants order; the previous
  wrappers extracted PIPs in sorted order then multiplied element-wise by
  z-scores in variants order, misaligning the arrays. Wrappers now
  re-index PIPs by `variant_id` before the mixture step. All v0.1.5
  numbers reported on the rice 3kRG grain pass before this fix should be
  regenerated; see `feedback_v15_misalignment_bug.md` in session memory.
- `hbp_finemap_from_sumstats` now defensively handles `graph_cache=None`
  (treated as empty cache) instead of raising on `.get()`.

### Documentation
- Methods §"v0.1.5 GAFM/HBP enhancements" updated to specify
  LD-deconvolved input to the mixture; rationale documented.
- Supplementary Figure S4 caption extended with v0.1.5 calibration
  panels (e–g), reporting empirical TDR ≈ 0.62 in the [0.9, 1.0] PIP
  bin (anti-conservative tail), null FPR (GAFM-MX/ENS 6%, HBP-MX 10%
  at PIP ≥ 0.5; all ≤ 1% at PIP ≥ 0.9), and π/γ hyperparameter
  insensitivity (Δrank-1 < 0.005 across 5 perturbations).
- New Supplementary Table S9: cross-species CS=1 sharpening (yeast
  0% → 24.5%, Arabidopsis 3.7% → 64.8%).
- New Supplementary Table S10: chr22 cross-h² head-to-head
  (mixture-prior delivers 2–3× confidence sharpening at same rank
  parity).

## [0.1.3] — 2026-05-02

Frozen PyPI release accompanying the original fine-mapping paper
submission.

### Added
- HBP and GAFM core fine-mapping (`finemapping_v2.py`).
- Multi-omics graph cache loaders for yeast, Arabidopsis, rice, human
  Pan-UKB.
- Three-tier ground-truth recovery scorecards.
