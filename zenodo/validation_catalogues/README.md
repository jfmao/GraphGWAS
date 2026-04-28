# validation_catalogues/

Subsets of upstream catalogues actually used for the 250-kb-window
validation of cross-species fine-mapping leads. Provided as
crosswalk pointers (gene → upstream source DOI / URL) so downstream
users can re-link without re-curating, while honouring upstream
licences.

## Files (staged at upload time by `zenodo/prepare_upload.sh`)

| File | Source | Licence | Content |
|---|---|---|---|
| `ren2023_rice_grain_genes.tsv` | Ren, Ding & Qian 2023 (`10.1016/j.molp.2023.04.003`) | TBD — pointer-only if not redistributable | 269 rice grain-quality genes. Columns: `gene_id`, `gene_symbol`, `chromosome`, `start`, `end`, `trait_class`, `source_table`. |
| `aragwas_bonferroni_subset.tsv` | AraGWAS (Togninalli 2018) | open access (NAR) | Bonferroni-significant variants restricted to the 14 *Arabidopsis* priority traits used here. |
| `bloom2015_yeast_qtl.tsv` | Bloom 2015 (`10.1038/ncomms9712`) | CC-BY | Yeast QTL gene table. |
| `peter2018_yeast_growth_genes.tsv` | Peter 2018 (`10.1038/s41586-018-0030-5`) | open access (Nature) | Yeast growth-trait gene table. |
| `human_lipid_canon_subset.tsv` | GIANT 2018 + GLGC 2013 + Klarin 2018 + Yengo 2022 + GWAS-Catalog | open access | 36 loci across BMI, height, LDL, triglycerides — the curated catalogue used for the chr-22 lipid canon validation (Methods §"Cross-species fine-mapping protocol"). |
| `flowering_canon_arabidopsis.tsv` | curated literature pointers (FRI, FLC, CONSTANS, GIGANTEA, VIN3, FT, DOG1, HKT1) | – | Canonical flowering-time genes used for *Arabidopsis* validation when AraGWAS Bonferroni hits do not cover the trait. |

## Schema

```
gene_id          STRING   upstream identifier (Os locus, AGI, SGD, gene symbol)
chromosome       STRING
start            INT
end              INT
trait_class      STRING   short trait/category label
source_doi       STRING   DOI of the originating catalogue or paper
source_table     STRING   table number / supplementary file in the source
notes            STRING   free-text (optional)
```

## Reproducibility

These TSVs are the *exact* gene/variant lists used by
`tests/multispecies_groundtruth_validation.py` to compute the
"validation rate" column of Sup Table S6 and the per-species
validation tallies in the manuscript text.

## Licence

This aggregation is MIT-licenced. **Each row points to an upstream
source whose own licence applies** to redistribution of that row's
content. For Ren 2023 we provide pointer-only entries (gene IDs +
genomic coordinates from public databases) rather than redistributing
the supplementary table; users should consult the original publication
for full annotations.
