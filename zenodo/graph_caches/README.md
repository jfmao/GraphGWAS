# graph_caches/

The per-species multi-omics graph caches that drive every fine-mapping
result in the paper. Each cache is a JSON dictionary keyed by
`variantId` (`chr:pos:ref:alt`) → `{genes, pathways, ppi, eqtl, ...}`.
HBP and GAFM both consume these caches at fine-mapping time;
re-running fine-mapping without them requires re-traversing the
graph database (\~minutes per locus instead of milliseconds).

The intervention extension is the per-variant `prior_score` field
that drives Table 1 (the heterogeneity-not-density intervention).
Caches without `prior_score` reproduce the *baseline* row of Table 1;
caches with `prior_score` reproduce the *+heterogeneous prior* row.

## Files (staged at upload time by `zenodo/prepare_upload.sh`)

### Human (1KG Phase 3 + GTEx v8 + STRING v12 + ENCODE cCRE v4)

| File | Size (approx) | Content |
|---|---|---|
| `gtex_chr22_enhanced_cache.json` | 10 MB | Chr 22 baseline cache: variant → {genes, pathways, ppi, eqtl}. |
| `human_1kg_cache_chr1..22.json` (×22, baseline) | \~700 MB | Per-chromosome baseline caches (full autosome). |
| `human_1kg_cache_chr1..22_intervention.json` (×22) | \~750 MB | Same caches with `prior_score` populated from ENCODE cCRE class flags (PLS, pELS, dELS, CTCF-bound, DNase-H3K4me3 — 9-class weighting documented in Methods §"Heterogeneity-by-intervention experiment"). |

### Rice (3kRG + RAP-DB + Ren 2023 + RicePPINet)

| File | Size (approx) | Content |
|---|---|---|
| `rice_graph_cache_v2_Chr1..12.json` (×12) | \~70 MB | Per-chromosome rice caches (RAP-DB / MSU genes, Ren 2023 grain-quality pathways, RicePPINet PPI ≥ 0.7). |
| `rice_graph_cache_v2_Chr1..12_intervention.json` (×12) | \~75 MB | Same caches with `prior_score` populated from snpEff impact class (HIGH=1.0, MODERATE=0.7, LOW=0.4, MODIFIER=0.1) on the 3kRG pseudo-canonical VCF. |

### *Arabidopsis* (1001 Genomes + TAIR10 + Plant Reactome + STRING-Ath)

| File | Size (approx) | Content |
|---|---|---|
| `arabidopsis_graph_cache_v3_chr1..5.json` (×5) | \~25 MB | Per-chromosome cache: TAIR10 genes, Plant Reactome pathways, STRING-Ath v12 PPI ≥ 700, UniProt `ARATH_3702` resolution. |
| `arabidopsis_graph_cache_v3_chr1..5_intervention.json` (×5) | \~28 MB | Same caches with `prior_score` from TAIR10 feature class (CDS=1.0, exon=0.8, gene+2kb-promoter=0.5). |

### Yeast (1011 Genomes + SGD + GO-slim + BIOGRID)

| File | Size (approx) | Content |
|---|---|---|
| `yeast_graph_cache_v2.json` | 250 MB | Single-file cache (yeast genome small enough): SGD ORFs, GO-slim, BIOGRID 5.0.256 *S. cerevisiae* S288c physical interactions. |
| `yeast_graph_cache_v2_intervention.json` | 270 MB | Same cache with `prior_score` from ORF-feature class (verified=1.0, tRNA/ncRNA=0.7, pseudogene=0.4) plus the BIOGRID PPI layer. |

## Cache schema (per-variant entries)

```jsonc
{
  "chr:pos:ref:alt": {
    "genes":     [{"gene_id": "...", "consequence": "missense_variant"}, ...],
    "pathways":  ["pathway_id_1", "pathway_id_2"],
    "ppi":       [{"partner_gene": "...", "score": 0.85}, ...],
    "eqtl":      [{"tissue": "Liver", "slope": 0.42, "pvalue": 1e-12}, ...],
    "regulatory_class": "pELS-CTCF",        // human only
    "snpeff_impact":    "MODERATE",         // rice only
    "tair10_class":     "CDS",              // Arabidopsis only
    "orf_class":        "verified",         // yeast only
    "prior_score":      0.85                // intervention-only field
  },
  ...
}
```

## Reproducibility

The baseline caches are built by the `tests/build_*_cache.py` scripts
in the source repository, which traverse the species-specific graph
database and write each variant's annotation neighbourhood. The
intervention extension is built by `tests/inject_prior_score_*.py`,
which adds the `prior_score` field from the species-specific
heterogeneity layer.

## Licence

MIT (this aggregation). The annotation layers themselves derive from
upstream consortia (GTEx v8, STRING v12, ENCODE cCRE v4, RAP-DB,
RicePPINet, TAIR10, Plant Reactome, SGD, BIOGRID); each must be
cited per its own licence terms.
