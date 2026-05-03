# Paper #2 §Y.8 (extended) — Human pipeline + binary-vs-weighted A/B

**Question:** can the §Y.8 catalogue extension be exercised on a real
whole-genome human cohort, and does the §Y.7 weighted multi-source
substrate (Billmann qGI + DepMap ED + STRING) actually improve recall
over binary STRING-PPI alone?

**Scripts:**
- `tests/multi_pair_catalogue_benchmark.py --species human` (binary)
- `tests/multi_pair_catalogue_benchmark.py --species human --substrate weighted` (weighted)

**Outputs:**
- `results/paper2_epistasis/multi_pair_catalogue_benchmark_human_binary.json`
- `results/paper2_epistasis/multi_pair_catalogue_benchmark_human_weighted.json`

## Headline grid

| Substrate | Tested | M5★ rank-1 | M5★ NF/0 |
|---|---|---|---|
| Binary STRING-PPI (cache.ppi) | 70/79 | **50/70 (71%)** | 20 |
| Weighted Billmann + DepMap + STRING | 70/79 | **51/70 (73%)** | 19 |

**One pair flips NF/0 → rank 1**: `Condensin_I-PARP1-XRCC1_complex`
(NCAPD3 × PARP1 representative pair) recovers at rank 1 of 21,730 cross-pairs
once the weighted Billmann+DepMap edges augment STRING. *No regressions.*

## Untestable breakdown (9/79)

- 1 HGNC unresolvable: `HLA-B × KIR3DS1` (KIR3DS1 is a polymorphic locus in
  the LRC region; not in the cache's gene-name vocabulary).
- 8 panel-MAF gaps (all CORUM complexes from Billmann tier):
  CCC-Wash, MSL, BRCC, **FA-core (anchor!)**, GPI-GnT, Respiratory Chain I,
  TIM, GARP. These complexes have <3 MAF≥0.05 variants per gene
  in 1KG-3202; some include chrX hemizygous-male MAF artefacts (FANCB on chrX).

The FA-core anchor being untestable is unfortunate — it was clinically
prominent (Fanconi anemia gene family) AND in the AUPRC top-15. The
issue is panel-coverage, not biology: 1KG samples don't reach MAF≥0.05
for many FA-pathway genes because pathogenic FA variants are rare and
heterogenous. A barcoded-segregant or disease-cohort run would recover it.

## Tier breakdown (binary substrate)

| Tier | Tested | Rank-1 | M5★ failed |
|---|---|---|---|
| `A_gwas_pharmacogenetic` | 7 | **7** (100%) | 0 |
| `A_billmann_corum_high` | 9 | 7 (78%) | 2 |
| `B_billmann_corum_supported` | 50 | 34 (68%) | 16 |
| `C_billmann_corum_only` | 4 | 2 (50%) | 2 |

Notable: **all 7 testable GWAS/pharmacogenetic pairs** (HLA-B×ERAP1,
TPMT×NUDT15, VKORC1×CYP2C9, BRCA1×PARP1, HFE×HFE, HLA-DQA1×DQB1,
HLA-DR3×DR4) recover at rank 1. The Billmann tier-A pairs are 7/9
(BAF, E2F-6, PBAF, BRCA1-A, etc.). The 20 NF/0 pairs are concentrated
in the broader Billmann tier-B set, where some CORUM complexes lack
sufficient STRING-PPI density between member-gene HAP1-screen
representatives.

## Why is the binary→weighted gain so small?

The diagnostic `weighted-substrate: loaded 119,318 edges, N added to G`
shows N is typically 0–2 across the catalogue gene-windows we load.
Interpretation: STRING v12 already covers the bulk of biochemically-curated
complex co-membership for the gene-pairs we're testing. The weighted
substrate doesn't help when STRING is already comprehensive.

Where it WOULD help:
- Population-scale disease-modifier scans (genes with no biochemical-PPI
  but with cellular-screen GIs in HAP1).
- Novel pairs from biobank summary statistics (where neither STRING nor
  curated CORUM has an edge yet).
- Cancer-cell-line-specific dependencies (DepMap captures these; STRING
  doesn't).

For the §Y.7 paper claim, this is honest evidence:
- The weighted substrate is a **strict generalisation** of binary STRING (no regressions, +1 rescue across 70 pairs).
- The marginal recall gain on a biochemically-curated catalogue is small
  precisely **because that catalogue overlaps STRING by construction**.
- The expected payoff is larger for queries outside curated complexes.

## What's still missing for the human arm

- The **8 panel-MAF gaps** would need either a deeper cohort
  (UK Biobank, AllOfUs) or relaxed MAF threshold (1% or below).
- The **20 NF/0 binary failures** that don't flip with the weighted
  substrate would benefit from richer human-PPI sources (e.g.,
  BioGRID, HuRI human reference interactome) being folded into
  `data/weighted_edges/human.tsv`.
- The cache's `pathways` field is currently empty for all human
  variants (per `human_graph_cache_v2_summary.tsv`); populating it from
  Reactome would add a `same_pathway` motif source for M2 on human.

## Reproduction commands

```bash
# Step 0: ensure gene index includes human (one-time)
python scripts/build_gene_coordinates_index.py

# Step 1: build weighted edge files (auto-downloads File_S4 + File_S21)
python scripts/build_weighted_edge_files.py

# Step 2a: binary baseline on human
python tests/multi_pair_catalogue_benchmark.py --species human \
    --output-suffix _human_binary

# Step 2b: weighted variant on human
python tests/multi_pair_catalogue_benchmark.py --species human \
    --substrate weighted --output-suffix _human_weighted

# Step 3: consolidated 4-species run (yeast + Arabidopsis + rice + human)
python tests/multi_pair_catalogue_benchmark.py --species all \
    --output-suffix _4species
```

Required data: per-chromosome `human_graph_cache_v2_chr{1..22,X}.json`
and `tests/data/human/1kGP_3202/1kGP_high_coverage_Illumina.chr{1..22}*.vcf.gz`
must be present (already on local disk; total ~22 GB VCF + ~1.4 GB cache).
