# Multi-species all-phenotype GWAS → fine-mapping → ablation plan
**Date:** 2026-04-25
**Goal:** Replicate the 3kRG IRRI workflow for human, yeast, and Arabidopsis, with proper literature-derived ground truth.
**Hardware:** 32 cores, 62 GB RAM (49 GB available).

## Resource inventory (verified 2026-04-25)

| Species | Phenotypes | Genotypes | Multi-omics cache |
|---------|------------|-----------|-------------------|
| **Yeast** (1011-genome) | 35+ grammar-corrected GWAS sumstats already done in `results/grammar_corrected/gwas_*.tsv` | Neo4j dump available (currently DOWN) | Built from Neo4j (`_cache_chr_graph_structure`); SGD genes + BIOGRID PPI + KEGG pathways implied by yeast_biological_discovery.py |
| **Arabidopsis** (1001G) | 633 phenotype CSVs in `tests/data/arabidopsis/phenotypes/` + 3 named (FT10, FT16, FLC) + AraGWAS bonferroni associations file (= ground truth) | 18.3 GB master VCF on disk; BGEN files at `/tmp/arabi/` | TO BUILD (Araport11 + AraPath + STRING-Ath at score≥700) |
| **Human** | 4 Pan-UKB tbi files (BMI=21001, height=50, biomarkers 30780/30870) — sumstats need download; chr22 simulated GWAS exists | 1KG chr22 phased VCF (~430k samples in Pan-UKB) | GTEx chr22 cache exists; STRING v12 9606 done; need full-genome version |

## Plan (3 species × 4 steps; parallelizable across 32 cores)

### Step 0 — Plan + setup (this file)

### Step 1 — Yeast (lowest-effort first)
- 1.1 Start Neo4j yeast DB (10 min)
- 1.2 Build / cache yeast multi-omics graph from Neo4j (30 min)
  → `data/yeast/yeast_graph_cache.json`
- 1.3 Fine-map top-N leads per trait across all 35 GWAS sumstats files (sumstats-only path, no genotype refetch needed if LD ref preserved)
  → `results/yeast_finemap_summary.tsv`
- 1.4 Ablation on 5 best-resolved loci × 5 fractions × 10 seeds
  → `results/ablation_yeast.{png,pdf,md}`
- 1.5 Ground-truth from Bloom 2015 + Peter 2018 (already partially cited; codify table)

### Step 2 — Arabidopsis
- 2.1 Build Arabidopsis multi-omics cache (1-2 h)
  - Variant→Gene: TAIR10 GFF3 / Araport11 (download)
  - Gene→Pathway: AraCyc / KEGG_ath (download or use Plant Reactome already in `/mnt/data/GraphPop/data/raw/plant_reactome/`)
  - Gene↔Gene PPI: STRING v12 species 3702
  → `data/arabidopsis/arabidopsis_graph_cache_chr{1..5}.json`
- 2.2 Pick ~30 priority phenotypes:
  - 3 named (FT10, FT16, FLC)
  - Top 30 traits by AraGWAS-bonferroni-association count
- 2.3 GWAS via plink2 multi-trait pass (use pgen built from VCF; create if absent)
  → `results/arabidopsis/gwas_arabi.<TRAIT>.glm.linear`
- 2.4 Fine-map top-3 GW + top-1 suggestive per trait with L1 + HBP using built cache
  → `results/arabidopsis_finemap_summary.tsv`
- 2.5 Ablation on 5 best-resolved loci
- 2.6 Ground-truth from AraGWAS (already on disk!) + literature for FT/FLC

### Step 3 — Human (Pan-UKB)
- 3.1 Download or locate Pan-UKB sumstats tsv.bgz files for 4 phenotypes (~2-3 GB each)
- 3.2 For each phenotype:
  - Identify GW-sig leads (p ≤ 5e-8) genome-wide (not just chr22)
  - Greedy 500 kb clustering → top-100 leads per trait
  - Use 1KG EUR / chr22 LD ref (or Pan-UKB in-sample LD when available)
- 3.3 Fine-map top-20 leads per trait per chromosome with L1 (sumstats path) + HBP using genome-wide GTEx + STRING cache (need to extend chr22-only cache to full genome)
- 3.4 Ablation
- 3.5 Ground-truth from GWAS Catalog (already-known SNPs for BMI, height, glucose, urate)

### Step 4 — Cross-species summary
- 4.1 Aggregate single-variant CS counts per species per trait
- 4.2 Compare ablation rank-improvement curves across species
- 4.3 Update paper Results section with 4-species roll-up table

## Realistic time budget (parallelism via 32 cores)

| Step | Wall-clock | Why |
|------|-----------:|-----|
| 1.1–1.2 yeast Neo4j + cache | 30–45 min | one-time DB warmup + Cypher query |
| 1.3 yeast fine-map 35 traits | 1–2 h | parallel across 8 traits at a time |
| 1.4 yeast ablation | 30 min | 5 loci × 5 fractions × 10 seeds |
| 2.1 arabidopsis cache | 1–2 h | I/O dominated, single thread |
| 2.3 arabidopsis multi-trait GWAS | 2–4 h | plink2 single pass × 30 traits |
| 2.4 arabidopsis fine-map | 1 h | parallel across 8 traits |
| 2.5 arabidopsis ablation | 30 min | |
| 3.1 download human sumstats | 1–3 h | network-bound |
| 3.3 human fine-map | 2–4 h | many leads, parallel |
| 3.4 human ablation | 30 min | |
| 4.x summary + paper update | 1 h | |
| **Total** | **10–18 h** | mostly background |

## Ground-truth strategy

Three tiers per species:
1. **Database-validated** (highest confidence): AraGWAS bonferroni hits, GWAS Catalog associations, Bloom 2015 / Peter 2018 yeast QTL maps
2. **Literature-canonical**: e.g. FRI/FLC for FT10 in Arabidopsis, BMI heritability genes (FTO, MC4R, etc.) for human
3. **Pathway-level support**: top variant in known trait-relevant pathway

Validation metric per locus: **fraction of GW-sig fine-mapped CS that contains at least one gene in the ground-truth list within 250 kb of the lead variant**.

## Honest limitations

- Yeast Neo4j DB is required for graph cache; if it fails to start, fall back to a SGD+BIOGRID-built static cache.
- Pan-UKB sumstats may not be locally downloaded — only `.tbi` indexes exist; tsv.bgz files are 1-3 GB each. May need to download or use tabix over HTTPS.
- Arabidopsis pgen may not exist (only VCF). plink2 conversion VCF→pgen on 18 GB takes ~30-60 min.
- 633 Arabidopsis phenotypes is too many to fine-map all of; cap at top 30 by AraGWAS hit count.
- Human chr22-only Pan-UKB demo done previously; extending to genome-wide adds significant compute.

## Execution order (start with shortest path to results)

1. ✅ Plan written (this file)
2. → Yeast Neo4j start + cache build (background, while writing other code)
3. → Arabidopsis cache build (background, while yeast runs)
4. → Yeast fine-mapping (35 traits) once cache exists
5. → Arabidopsis pgen + multi-trait GWAS
6. → Human Pan-UKB download check
7. → Arabidopsis fine-mapping
8. → Human fine-mapping
9. → All ablation tests
10. → Cross-species summary report + paper update
