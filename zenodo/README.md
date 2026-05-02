# GraphGWAS — Zenodo deposit

This folder contains all artifacts to be uploaded to Zenodo alongside
the v0.1.3 PyPI release of GraphGWAS (https://pypi.org/project/graphgwas/0.1.3/)
and the GitHub release tag v0.1.3, accompanying the Nature Genetics
manuscript submission:

> Estaji, Zhao, Chen, Nie, Mao (2026). *Relational biological structure
> improves fine-mapping of causal GWAS variants under weak signal*
> (submitted).

The target Zenodo record is a **community-of-record companion**:
versioned, citable, and linked to the GitHub repo at
<https://github.com/jfmao/GraphGWAS> via `.zenodo.json` metadata.

---

## Contents (twelve top-level folders)

```
zenodo/
├── README.md                          # this file
├── CITATION.cff                       # author / version / license / pre-print pointer
├── .zenodo.json                       # Zenodo metadata (authors, licence, grants)
├── MANIFEST.md                        # SHA-256 hashes + file sizes of all artifacts
│
├── graph_dumps/                       # pre-built graph-database dumps (Neo4j 5.26)
│   ├── yeast_1011_v0.1.dump          # 0.5 GB — 1011 Yeast Genomes + SGD + GO
│   ├── human_1kg_multiomics_v0.1.dump # 17 GB — 1000 Genomes + GENCODE/GTEx/STRING/ENCODE
│   └── README.md                      # restore instructions, schema version
│
├── graph_caches/                      # NEW: per-species multi-omics caches
│   │                                  #   that drive HBP / GAFM fine-mapping;
│   │                                  #   includes intervention extension with
│   │                                  #   per-variant prior_score (Table 1)
│   ├── gtex_chr22_enhanced_cache.json # human chr22 baseline
│   ├── human_1kg_cache_chr*.json      # human autosomes baseline + intervention
│   ├── rice_graph_cache_v2_Chr*.json  # rice 12 chrs (intervention)
│   ├── arabidopsis_graph_cache_v3_chr*.json
│   ├── yeast_graph_cache_v2.json      # single-file (small genome)
│   └── README.md
│
├── benchmark_outputs/                 # JSON outputs backing every figure
│   ├── hbp_h2h_50rep.json              # 4-method × 3-scenario benchmark
│   ├── inf_methods_h2h.json            # SuSiE-inf / FINEMAP-inf × 3-scenario
│   ├── polyfun_proxy.json              # weak-signal Polyfun-proxy benchmark
│   ├── pip_calibration.json            # PIP calibration (Sup Fig S4)
│   ├── null_calibration.json           # null FPR (Sup Fig S4)
│   ├── power_vs_N.json                 # sample-size scaling (Sup Fig S5)
│   ├── cross_ancestry.json             # 1KG cross-ancestry simulation
│   ├── 100rep_l1_vs_susie_weak.json    # weak-signal headline (Fig 3)
│   └── README.md
│
├── panukb_results/                    # real-biobank fine-mapping outputs
│   ├── fto_bmi_4ancestries.json
│   ├── apoa5_triglycerides_4ancestries.json
│   ├── ldlr_ldl_4ancestries.json
│   ├── hmga2_height_4ancestries.json
│   ├── summary_all_16_cells.json
│   └── README.md
│
├── figure_source_data/                # raw data tables behind each main figure
│   ├── fig1_architecture/  fig2_hbp_schematic/
│   ├── fig3_weak_signal/   fig4_baselines/
│   ├── fig5_cross_ancestry/  fig6_ablation/  fig7_method_selection/
│   ├── tab1_intervention/
│   └── README.md
│
├── simulation_seeds/                  # per-replicate seed + sampled-window logs
│   ├── F1_strong_h2_0.10.tsv  F1_weak_h2_0.02.tsv  F1_eQTL_weak_beta0.15.tsv
│   ├── null_beta0_100rep.tsv  cross_ancestry_1kg.tsv  power_vs_N.tsv
│   └── README.md
│
├── irri_3krg_gwas/                    # NEW: 18-trait whole-genome GWAS scan
│   │                                  #   on the 3kRG panel (substrate for Sup Table S5)
│   ├── gwas_irri.<TRAIT>.glm.linear.gz  # 18 PLINK 2.0 sumstats (compressed)
│   ├── irri_phenotype_summary.tsv
│   ├── irri_lambda_gc.tsv
│   ├── irri_gwas_summary.tsv
│   ├── irri_lead_loci.tsv
│   ├── irri_finemap_v2_summary.tsv
│   ├── irri_irgc_to_3krg.tsv          # IRGC ↔ 3kRG sample-ID crosswalk
│   └── README.md
│
├── multispecies_leads/                # NEW: 692 leads × 4 species with
│   │                                  #   baseline & intervention CS sizes
│   │                                  #   and validation-catalogue overlaps
│   ├── cross_species_summary.tsv      # matches Sup Table S6
│   ├── leads_yeast_1011.tsv           # 245 leads
│   ├── leads_arabidopsis_1001g.tsv    # 54 leads
│   ├── leads_rice_3krg.tsv            # 72 leads
│   ├── leads_human_chr22.tsv          # 321 leads (chr22 subset)
│   ├── leads_human_genome_wide.tsv    # 321 leads (genome-wide)
│   ├── heterogeneity_intervention_report.md  # matches Table 1
│   ├── cross_species_report.md
│   └── README.md
│
├── validation_catalogues/             # NEW: subsets of upstream catalogues
│   │                                  #   used for 250-kb validation
│   ├── ren2023_rice_grain_genes.tsv
│   ├── aragwas_bonferroni_subset.tsv
│   ├── bloom2015_yeast_qtl.tsv
│   ├── peter2018_yeast_growth_genes.tsv
│   ├── human_lipid_canon_subset.tsv
│   ├── flowering_canon_arabidopsis.tsv
│   └── README.md
│
├── ld_references/                     # NEW: pre-computed ancestry-matched LD
│   │                                  #   for 1KG cross-ancestry sim + Pan-UKB
│   ├── 1kg_chr22_{eur,afr,eas}_*_ld.npz
│   ├── 1kg_nygc_grch38_chr*_eur_ld.npz
│   ├── panukb_ancestry_to_1kg_subset.tsv
│   └── README.md
│
├── manuscript/                        # NEW: pre-print PDF
│   ├── manuscript_preprint.pdf
│   └── README.md
│
└── code_snapshot/                     # tagged v0.1.3 release (git archive)
    ├── graphgwas_v0.1.3.tar.gz        # source code at submission
    ├── commit_sha.txt
    └── README.md
```

After running `bash zenodo/prepare_upload.sh`, an additional `upload/`
directory is created with one `.tar.gz` per top-level subfolder plus
standalone PDFs and dumps. **This is what you drag into Zenodo's
"New upload" form.** The hierarchical layout above is the working
tree; the flat `upload/` bundle is what reviewers actually receive.

---

## Upload order and DOI strategy

1. Upload the **`code_snapshot/`** tarball + **`CITATION.cff`** first
   (fastest, ~5 MB); this receives the primary *software* DOI and is
   what `CITATION.cff` references.
2. Upload `graph_dumps/` (17.5 GB total) as separate deposits, each
   with its own DOI. Cross-reference both to the software DOI.
3. Upload `graph_caches/` (~2 GB) as a separate deposit (it is
   updated more frequently than the dumps and benefits from its own
   DOI for citation in follow-up papers).
4. Upload `benchmark_outputs/`, `panukb_results/`, `figure_source_data/`,
   `simulation_seeds/`, `multispecies_leads/`, `validation_catalogues/`,
   and `manuscript/` as a single *data* DOI (all together < 1 GB).
5. Upload `irri_3krg_gwas/` (~500 MB compressed) as its own *data*
   DOI — researchers reusing the rice GWAS sumstats will want to
   cite this independently.
6. Upload `ld_references/` (~5 GB) as its own *data* DOI for the
   same reason.
7. Update the paper's `Data Availability` statement and the
   top-level `CITATION.cff` with all assigned DOIs before
   camera-ready.

---

## Status at deposit time (2026-04-30)

| Item | Status | Size (approx) | Ready? |
|---|---|---:|---|
| Source-code tarball (v0.1.3) | Tagged + on PyPI | 5 MB | ✓ |
| Human 1KG + multi-omics graph dump | Generated; in `backups/` | 17 GB | ✓ (local) |
| Yeast 1011 graph dump | Generated; in `backups/` | 0.5 GB | ✓ (local) |
| Benchmark JSONs (all 10+ files) | Generated from `results/benchmark_v2/` | ~50 MB | ✓ (local) |
| Pan-UKB per-locus JSONs | Generated from `results/panukb/` | ~10 MB | ✓ (local) |
| Multi-omics graph caches (4 species, baseline + intervention) | Generated from `data/{annotations,rice_3k,yeast,arabidopsis}/` | ~2 GB | ✓ (local) |
| 3kRG IRRI sumstats (18 traits) | Generated from `data/rice_3k/results/irri_gwas/` | ~500 MB compressed | ✓ (local) |
| 3kRG IRRI summary tables | Already staged (small) | < 2 MB | ✓ |
| Cross-species lead lists (4 species, 692 leads) | Already staged (small) | < 200 KB | ✓ |
| Heterogeneity intervention report | Already staged (small) | < 10 KB | ✓ |
| Manuscript pre-print PDF | Already staged | 1.1 MB | ✓ |
| `CITATION.cff` | Already staged | < 5 KB | ✓ |
| Validation catalogue subsets | Pending — to curate per upstream licence | < 1 MB | ⏳ |
| LD references (chr22 + genome-wide EUR) | Pending — to extract from 1KG NYGC | ~5 GB | ⏳ |
| Figure source data (TSV per figure) | Pending — to extract on acceptance | ~100 MB | ⏳ |
| Simulation seeds | Pending — to extract from simulate.py logs | ~5 MB | ⏳ |
| Container image (Docker / Apptainer) | Deferred to camera-ready | ~2-4 GB | ⏳ (deferred) |

Items marked ⏳ are finalised by `bash zenodo/prepare_upload.sh`
(graph dumps, caches, sumstats, LD refs) or by the post-acceptance
camera-ready window (validation catalogues, container image).
Items marked ✓ are already present in this folder and ready to upload.

---

## Licence

All artifacts are released under the MIT Licence (same as the source
code). The Neo4j dumps embed data derived from:

- 1000 Genomes Project Phase 3 and high-coverage NYGC release
  (IGSR data policy: unrestricted use with citation)
- 1011 Yeast Genomes Project (Peter *et al.* 2018)
- GENCODE v47 (Wellcome Sanger; open licence)
- GTEx v8 (NIH; open-access sumstats)
- STRING v12 (CC-BY 4.0)
- ENCODE cCRE v4 (CC-BY 4.0)
- Reactome (CC-BY 4.0)

Downstream users must cite the original consortia in addition to this
deposit.

---

## Citation

Citing this Zenodo record:

> Estaji, E., Zhao, S.-W., Chen, Z.-Y., Nie, S. & Mao, J.-F. (2026).
> *GraphGWAS v0.1.3: relational fine-mapping of causal GWAS variants
> on a multi-omics knowledge graph.* Zenodo.
> https://doi.org/10.5281/zenodo.[DOI-ASSIGNED-ON-UPLOAD]

Citing the accompanying paper:

> Estaji, E., Zhao, S.-W., Chen, Z.-Y., Nie, S. & Mao, J.-F. (2026).
> Relational biological structure improves fine-mapping of causal GWAS
> variants under weak signal. *Nature Genetics* (submitted).
