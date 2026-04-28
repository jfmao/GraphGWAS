# ld_references/

Pre-computed ancestry-matched LD reference matrices used for the
controlled cross-ancestry simulation (Fig 5 in the paper) and for
Pan-UKB sumstats fine-mapping when the Hail in-sample LD path is not
available. Avoids re-downloading 1KG and recomputing LD for downstream
users.

## Files (staged at upload time by `zenodo/prepare_upload.sh`)

| File | Size (approx) | Population | Content |
|---|---|---|---|
| `1kg_chr22_eur_n503_ld.npz` | \~150 MB | EUR ($n{=}503$) | Chr 22 sparse LD matrix on common ($\text{MAF} > 0.01$) variants. |
| `1kg_chr22_afr_n661_ld.npz` | \~200 MB | AFR ($n{=}661$) | Same. |
| `1kg_chr22_eas_n504_ld.npz` | \~140 MB | EAS ($n{=}504$) | Same. |
| `1kg_nygc_grch38_chr1..22_eur_ld.npz` (×22) | \~5 GB total | EUR ($n{=}503$) | Genome-wide pre-computed LD matrices used by the genome-wide Pan-UKB EUR fine-mapping pipeline. |
| `panukb_ancestry_to_1kg_subset.tsv` | < 1 KB | – | Mapping from Pan-UKB ancestry codes (EUR/CSA/AFR/EAS) to the 1KG subpopulations used as their LD reference. |

## Format

`.npz` (numpy savez) files containing:
- `r2`: float32 matrix (sparse upper-triangular) of squared Pearson correlations between variants
- `variant_ids`: array of `chr:pos:ref:alt` strings
- `positions`: array of int64 positions
- `metadata`: dict with `population`, `n_samples`, `build` ("GRCh38"), `source_release` ("NYGC 30x")

## Reproducibility

LD matrices are computed by `tests/build_ld_reference.py` from the
1000 Genomes Project NYGC 30× release (publicly released). For each
population subset, the script restricts to common variants
($\text{MAF} > 0.01$) and stores the upper-triangular $r^2$ values
above a 0.05 cutoff to keep file size manageable.

## Licence

MIT (this aggregation). The underlying genotype data is governed by
the IGSR data policy (unrestricted use with citation of
Byrska-Bishop et al. 2022 and the 1000 Genomes Project
Consortium 2015).
