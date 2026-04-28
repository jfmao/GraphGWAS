# Pan-UK Biobank cross-ancestry fine-mapping outputs

Per-locus JSON outputs from the real-biobank demonstration in
main-text §2.7 (Pan-UKB) and Supplementary Fig. S2.

## Files

| File | Locus | Trait (Pan-UKB phenocode) | Lead (GRCh37) |
|---|---|---|---|
| `fto_bmi_4ancestries.json` | FTO | Body mass index (21001) | chr16:53,820,527 |
| `apoa5_triglycerides_4ancestries.json` | APOA5 | Triglycerides (30870) | chr11:116,662,407 |
| `ldlr_ldl_4ancestries.json` | LDLR | LDL cholesterol (30780) | chr19:11,200,038 |
| `hmga2_height_4ancestries.json` | HMGA2 | Standing height (50) | chr12:66,360,071 |
| `summary_all_16_cells.json` | — | 4 loci × 4 ancestries | — |

## JSON schema (per file)

```json
{
  "locus": "FTO",
  "trait": "BMI",
  "phenocode": "21001",
  "trait_type": "continuous",
  "modifier": "irnt",
  "chr": "16",
  "center": 53820527,
  "window": 100000,
  "t_sumstats_s": 2.5,
  "ancestries": {
    "EUR": {
      "N_sumstats": 420531,
      "n_sumstats_variants": 1574,
      "n_ld_variants": 1621,
      "n_intersect": 742,
      "lead_sumstats_variant": "16:53768582:C:T",
      "top_hbp_variant": "16:53805443:GT:G",
      "top_hbp_pip": 1.000,
      "n_cs_hbp": 1,
      "top_l1_variant": "16:53805443:GT:G",
      "top_l1_pip": 1.000,
      "n_cs_l1": 1,
      "t_ld_s": 2.3,
      "t_hbp_s": 0.004,
      "t_l1_s": 0.003
    },
    "CSA": {...},
    "AFR": {...},
    "EAS": {...}
  }
}
```

## Reproducibility

All 16 cells regenerable in ~60 s from a single command:

```bash
cd /path/to/GraphGWAS
python tests/benchmark_panukb_finemap.py --all --out results/panukb/
```

The script streams Pan-UKB sumstats via `tabix` over HTTPS against
the public `pan-ukb-us-east-1` S3 bucket (no authentication
required), lifts GRCh37 variants to GRCh38 using `pyliftover`, and
uses ancestry-matched 1KG BGEN as the LD reference. The primary
software record at `../code_snapshot/` contains the complete
`src/python/graphgwas/panukb.py` module and
`tests/benchmark_panukb_finemap.py` driver.

## Notes on LD reference

The LD matrices in these results come from ancestry-matched
subsets of the 1000 Genomes high-coverage panel. Pan-UKB's own
in-sample LD BlockMatrices (47.6 TB total across six ancestries)
are publicly available on Amazon S3 but require a Hail
environment with `hadoop-aws` Spark classpath configuration to
read directly — an environment-specific setup outside the scope
of this deposit. Replacing the 1KG LD with in-sample Pan-UKB LD
is expected to tighten credible-set sizes for the smaller-N
non-European ancestries (AFR, EAS) and is flagged as post-
submission polish.
