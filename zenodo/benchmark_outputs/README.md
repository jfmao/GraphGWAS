# Benchmark outputs

JSON and TSV outputs backing every benchmark-derived figure and
table in the paper. Structure: one file per benchmark scenario.
All files are regeneratable from the source code at this software
record's tagged commit; the copies in this deposit are the
specific runs cited in the paper.

## File schema (common)

Every JSON file at minimum contains:

```json
{
  "scenario": "F1_weak_h2_0.02",
  "seed_base": 42,
  "n_replicates": 79,
  "date_run": "2026-04-20",
  "hardware": "Intel i9-13900K / 64 GB DDR5 / Ubuntu 24.04",
  "records": [
    {
      "replicate_id": 0,
      "chr": "chr22", "start": 17500000, "end": 17550000,
      "causal_variant_id": "chr22:17515432:A:G",
      "methods": {
        "l1": {"rank_1": true, "causal_rank": 1, "pip": 0.81, "runtime_s": 0.069, "credible_set_size": 3},
        "susie": {...},
        "finemap": {...},
        "susieinf": {...},
        "finemapinf": {...},
        "hbp": {...}
      }
    },
    ...
  ]
}
```

## Files

| File | Figure / Table in paper | N | Scenario |
|---|---|---:|---|
| `weak_signal_l1_vs_susie_79rep.json` | Main Fig 3 | 79 | F1-eQTL; β=0.15; h²=0.01 |
| `weak_signal_l1_vs_susie_100rep.json` | §2.3 honest replication | 100 | F1; β=0.20; h²=0.02; chr22 17–46 Mb |
| `wu2026_baselines_comparison.json` | Main Fig 4, Supp Table S2 | 30 | F1; h²=0.05; all 6 baselines + HBP + GAFM |
| `polyfun_proxy_comparison.json` | Main Fig 4 (Polyfun-proxy panel) | 20 | F1-eQTL; β=0.15 |
| `pip_calibration_200rep.json` | Main Fig 5 (calibration) | 200 | F1 × 4 h² levels |
| `null_fpr_100rep.json` | Main Fig 5 (null FPR) | 100 | Null (β=0) |
| `power_vs_N_30rep_per_N.json` | Main Fig 6 | 4 × 30 | F1; h²=0.10; N ∈ {500, 1k, 2k, 3k} |
| `cross_ancestry_1kg_30rep_per_ancestry.json` | Main Fig 7a | 3 × 30 | F1; h²=0.05; EUR/AFR/EAS |
| `arabidopsis_ft10_hbp.json` | Main text §2.8 | — | Real phenotype; FT10 at chr4 |
| `yeast_35traits_fine_mapping.json` | Main text §2.8 | 35 | Real phenotypes on 1011 yeast |
| `epistasis_m1_s1_5rep.json` | Supp Fig S1 | 5 | S1 pure interaction; β_I=1.5 |

## Regeneration

Every file above is regenerable from a single command documented
in `docs/REPRODUCIBILITY.md` of the source code. Typical command:

```bash
python tests/benchmark_weak_signal.py --n-rep 79 --out weak_signal_l1_vs_susie_79rep.json
python tests/benchmark_inf_methods.py --n-rep 30 --out wu2026_baselines_comparison.json
python tests/benchmark_pip_calibration.py --n-rep 200 --out pip_calibration_200rep.json
# ...
```
