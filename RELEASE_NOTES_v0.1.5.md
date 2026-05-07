# GraphGWAS v0.1.5 — Release Notes

**Released:** 2026-05-07
**Tag:** `v0.1.5`
**PyPI:** <https://pypi.org/project/graphgwas/0.1.5/>
**Zenodo DOI:** [10.5281/zenodo.20065705](https://doi.org/10.5281/zenodo.20065705)

This release accompanies the Nature Genetics submission *"Relational biological
structure improves fine-mapping of causal GWAS variants under weak signal"*
(Estaji et al., 2026, paper #1).

## Headline additions

- **Mixture-prior posterior reweighting** — three new methods that add a
  SBayesRC-style 4-component Wakefield Bayes-factor reweighting on the
  LD-deconvolved z-scores, plus standard λ_GC deflation:
  - `gafm_mx_from_sumstats` → paper-facing **GAFM-MX**
  - `hbp_mx_from_sumstats` → paper-facing **HBP-MX**
  - `ensemble_from_sumstats` → paper-facing **ENS**
- On the panel-matched 3000 Rice Genomes grain weight + shape pass, each of
  the three reaches **10/21 (47.6%) top-1-PIP exact-position recovery** of the
  Niu 2021 21-stable-QTN catalogue — the highest of any method tested,
  beating SuSiE (28.6%) and SBayesRC (14.3%) while remaining 200–700×
  faster per locus than SuSiE.
- **SBayesRC integration** — custom 3kRG eigen-decomposed LD reference built
  from scratch with GCTB 2.05beta + R `LDstep1`–`LDstep4` over 23 targeted
  ±250 kb blocks (188 K MAF ≥ 0.05 variants; ~10 min compute), enabling
  head-to-head benchmarking.

## Other changes since v0.1.3

- **CLI grew from ~44 to 53 commands across 15 functional groups**
  (fine-mapping, single-locus, multi-locus, population structure, GNN,
  heritability, multivariate, MR, PRS, infrastructure).
- **Pre-registered control experiments** wired into the benchmark suite:
  placebo-prior (permuted-marginal control), null-overlap permutation,
  hyperparameter sensitivity sweeps over (α, λ_GC, π, γ).
- **3kRG IRRI 18-trait whole-genome scan** + 4-trait grain weight + shape
  pass with Niu 2021 ground truth.
- **Cross-species replication** at v0.1.5 across yeast (1011 Genomes,
  35 traits), Arabidopsis (1001 Genomes, 14 traits), rice (3kRG IRRI 18,
  3kRG grain 4) and Pan-UK Biobank (4 ancestries, 4 traits).
- **PIP calibration + null FPR** redone under v0.1.5 with `α = 0.7` default.
- **Practitioner documentation** — new `docs/INPUT_OUTPUT_GUIDE.md`,
  `docs/manual/index.md` (per-command reference at 53 commands),
  `vignettes/fine-mapping-quickstart.md` (~15-min FTO/BMI Pan-UKB walkthrough).

## Bugfix series leading up to v0.1.5

- **v0.1.4 → v0.1.5**: fixed mixture-prior alignment bug (sorted-PIPs were
  multiplied with variants-order z-scores; LD double-counting on the
  precision factor). Rerun affects rice 3kRG numbers; cross-species and
  Pan-UKB results unaffected.
- Removed v0.1.4 wheels that contained the alignment bug — install v0.1.5.

## Paper-facing → Python module mapping

Paper-facing names map to Python module identifiers as documented in
[`paper/finemapping_v1/methods.tex` Code Availability](paper/finemapping_v1/methods.tex):

| Paper | Python module | CLI |
|---|---|---|
| GAFM | `l1_finemap_from_sumstats` | `graphgwas finemap-sumstats --method gafm` |
| HBP | `hbp_finemap_from_sumstats` | `graphgwas finemap-sumstats --method hbp` |
| GAFM-MX | `gafm_mx_from_sumstats` | `graphgwas finemap-sumstats --method gafm-mx` |
| HBP-MX | `hbp_mx_from_sumstats` | `graphgwas finemap-sumstats --method hbp-mx` |
| ENS | `ensemble_from_sumstats` | `graphgwas finemap-sumstats --method ens` |
| CLGF | `cross_locus_graph_finemap` | (Python only) |
| LPCE | `graphgwas.epistasis_v2.lpce_*` | `graphgwas epistasis scan` |

The `l1_*` and `m1_*` prefixes are historical from the platform's internal
layer-numbering and are retained for back-compatibility; paper-facing names
are canonical.

## Reproducibility

Every figure and table in the paper is regeneratable from a single command
documented in [`docs/REPRODUCIBILITY.md`](docs/REPRODUCIBILITY.md).
Per-benchmark α values used for the headline runs:

- `α = 0.7`: chr22 head-to-head (`benchmark_v15_chr22.py`),
  cross-species replication, 3kRG grain weight + shape pass
  (`rice3k_grain_shape_finemap.py`), v0.1.5 PIP calibration / null FPR
  (`benchmark_v15_calibration_null.py`)
- `α = 0.5`: weak-signal tissue-specific eQTL regime
  (`benchmark_weak_signal.py`; the 27–2 GAFM-vs-SuSiE result)
- `α = 1.0`: sumstats-only Pan-UK Biobank fine-mapping
  (`benchmark_panukb_finemap.py`)

## Install

```bash
pip install graphgwas==0.1.5
```

## Citing

Citation metadata is in `CITATION.cff` and the Zenodo deposit metadata
(`zenodo/.zenodo.json`). For software citation use the Zenodo DOI
[10.5281/zenodo.20065705](https://doi.org/10.5281/zenodo.20065705); for
the accompanying paper, see the manuscript reference once published.
