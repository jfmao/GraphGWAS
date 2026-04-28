# irri_3krg_gwas/

The 18-trait whole-genome GWAS scan on the 3{,}000 Rice Genomes
panel (3kRG; Wang et al. 2018) using PLINK 2.0, plus per-trait
fine-mapping outputs. This is the substrate for Sup Table S5.

## Files (staged at upload time by `zenodo/prepare_upload.sh`)

| File | Rows (incl. header) | Description |
|---|---|---|
| `irri_phenotype_summary.tsv` | 19 | Per-trait `N`, mean, sd, min, max, missingness over the 2,266 IRRI-phenotyped accessions. |
| `irri_lambda_gc.tsv` | 19 | Per-trait $\lambda_{GC}$ on \~500{,}000 random 3kRG variants per trait. |
| `irri_gwas_summary.tsv` | 19 | Per-trait min p-value, count of GW-significant ($p \leq 5\times10^{-8}$) variants, count of independent leads after greedy 250 kb-window clustering. |
| `irri_lead_loci.tsv` | varies | Per-trait independent leads, lead position, lead p-value, locus window. |
| `irri_finemap_v2_summary.tsv` | varies | Per-trait fine-mapping results: HBP CS size, HBP top PIP, GAFM CS size, GAFM top PIP. |
| `gwas_irri.<TRAIT>.glm.linear` (×18, gzipped at upload time) | varies | PLINK 2.0 `--linear hide-covar` output per trait. Whole-genome \~29.6M biallelic SNPs each. Compressed total ≈ 500 MB. |

## Trait codes

Defined in Sup Table S5 of the manuscript. 18 IRRI Standard Evaluation
System codes (e.g. APCO_REV_REPRO = apiculus colour at reproductive
stage; SCCO_REV = spikelet/grain pericarp colour). The full
human-readable definitions are in the manuscript caption.

## PLINK 2.0 invocation

```
plink2 --bfile rice_3krg \
       --maf 0.01 --geno 0.05 \
       --pheno irri_18traits.tsv \
       --linear hide-covar \
       --covar rice_10pcs.tsv \
       --out gwas_irri.<TRAIT>
```

10 PCs computed via `plink2 --pca 10` on the 3kRG panel.

## Sample-ID linkage

The IRRI 18-trait phenotype catalogue uses IRGC accession numbers
(`IRGC_xxxxxxx`). The 3kRG VCF uses 3kRG accession identifiers.
The crosswalk is provided in the source repository under
`data/rice_3k/pheno/irri_irgc_to_3krg.tsv`; this Zenodo deposit
includes a copy of that crosswalk (staged at upload time) so
downstream users can re-link without recomputing.

## Reproducibility

Pre-requisites: 3kRG pseudo-canonical VCF (publicly released), the
phenotype TSV (`data/rice_3k/pheno/irri_18traits.tsv` in the
source repo), 10-PC covariate file, plink2 v2.0.0-a.6.5LM. Re-running
all 18 GWAS scans on the workstation hardware described in Methods
takes \~12 hours.

## Licence

MIT (this aggregation). Underlying 3kRG genotype data is publicly
released by the 3000 Rice Genomes Project; users must cite Wang
et al. 2018 (`10.1101/gr.240234.118`).
