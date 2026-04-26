# Multi-week functional-genomics expansion of GraphGWAS multi-omics cache

**Goal:** Add Tier-A/B/C functional annotations (TFBS, ATAC, histone marks, 3D contacts, regulatory elements, conservation, composite scores) on top of the existing pathway+PPI cache for **human, yeast, Arabidopsis, rice**. Re-run fine-mapping. Predict + verify the heterogeneity-not-density hypothesis.

**Falsifiable claim**: human GTEx+STRING (uniform per variant) showed **0/321** HBP-tighter loci. Adding heterogeneous per-variant priors (MotifBreakR Δ-PWM, ABC E-P contact strength, cCRE state) should produce a **non-zero** HBP-tighter rate on the same loci.

## Phase 0 — schema extension (1 day)

Extend `graph_cache` from `{genes, pathways, ppi}` to:
```
{variant_id: {
   genes: [...],
   pathways: [...],
   ppi: [...],
   regulatory: ["cCRE:CTCF_K562", "DNase:Liver", ...],     # binary overlaps
   tf_binding: ["FOXA1_disrupt:0.85", "CTCF_create:0.62"], # MotifBreakR ΔPWM with sign
   chromatin: ["chromHMM:Liver_TssA", ...],                # state labels
   contacts: ["ABC:LDL_LDLR_2.4kb_score=0.34", ...],       # E-P pairs
   prior_score: float                                       # composite scalar
}}
```

Also extend `hbp_finemap_from_sumstats` and `l1_finemap_from_sumstats` z_func builders to consume the new fields (treat `regulatory`/`tf_binding`/`contacts` as additional gene-like edge types).

## Phase 1 — Human (1-2 weeks)

### 1.A Local files (fastest, no downloads)
- **cCRE overlap**: `data/annotations/encode/encodeCcreCombined.bed` (already on disk) → binary regulatory-element annotation per variant
- **DNase peaks**: `encode_dnase_peaks.bed.gz` (on disk) → open-chromatin overlap per variant
- **CADD chr22**: `CADD_chr22.tsv.gz` (on disk) → composite Δ-score per variant
- Build augmented chr22 cache as proof-of-concept

### 1.B External downloads
- **CADD v1.6 genome-wide**: ~80 GB, https://cadd.bihealth.org → tabix-queryable scoring per variant
- **ABC model E-P pairs**: Nasser 2021 (~5 GB) → variant→target-gene assignment
- **ENCODE rE2G v3 maps**: ~3 GB → enhancer→gene predictions
- **Roadmap chromHMM 25-state** consolidated → per-variant state per cell type
- **JASPAR 2024 PWMs** + **MotifBreakR**: per-locus Δ-PWM scoring (compute-on-fly)
- **DeepSEA / Sei chromatin Δ-score**: precomputed for known variants if available

### 1.C Re-run + verify
- Re-run genome-wide Pan-UKB pipeline with augmented cache (4 phenotypes × 22 chr)
- Re-run ablation on best-resolved loci with new cache
- Compare HBP-tighter rate v1 vs v2: target ≥ 5% (vs current 0/321)

## Phase 2 — Yeast (3-5 days)

### 2.A Add layers
- **BIOGRID 559292** physical PPI (deferred from v1)
- **MacIsaac 2006 + Harbison 2004** ChIP-chip TF binding (~200 TFs, genome-wide)
- **YEASTRACT v8** TF→target regulatory network
- **JASPAR-fungi + YeTFaSCo PWMs** for per-variant Δ-PWM scoring
- **Pokholok 2005** histone modification ChIP-chip data
- **Brogaard 2012** nucleosome positioning

### 2.B Re-run
- Re-run all 35 trait sumstats with augmented cache
- Ablation on top-5 best-resolved loci with new layers

## Phase 3 — Arabidopsis (1 week)

### 3.A Add layers
- **DAP-seq 529 TF cistromes** (O'Malley 2016 Cell) — the gold standard, ~30M peaks
- **Lu 2019 ATAC-seq atlas** (root/shoot/leaf chromatin accessibility)
- **Marand 2021 scATAC** root cell-type atlas
- **Jin 2017 / Fukushima 2020** histone modification ChIP-seq (H3K4me3, H3K27me3, H3K9ac)
- **Wang 2015 / Liu 2017** Arabidopsis Hi-C 3D contacts
- **1001 Epigenomes** CG/CHG/CHH methylation per accession (Kawakatsu 2016)

### 3.B Re-run
- Re-run 14 priority phenotype GWAS finemap with augmented cache
- Ablation

## Phase 4 — Rice (1 week)

### 4.A Add layers
- **Cao 2020 + Tu 2020** rice ChIP-seq atlas (TFs in seed/leaf tissues)
- **Lu 2019 rice ATAC** (cross-tissue chromatin accessibility)
- **PlantPAN 4.0** motif binding-site scans (rice-specific)
- **Liu 2017 + Dong 2018** rice Hi-C contacts (indica/japonica root)
- **Cereal-conservation** PhyloP across rice/sorghum/maize/Brachypodium

### 4.B Re-run
- Re-run 18 IRRI trait fine-mapping
- Re-run rice ablation (BG1, TGW6, GS5, GW2, GW8) with new layers

## Phase 5 — cross-species consolidation (3 days)

- Aggregate HBP-tighter rates v1 vs v2 for each species
- Plot: heterogeneity score (computed from cache density/diversity per variant) vs HBP-tighter rate
- Update paper: prediction-verification table; revised cross-species summary
- Commit final manuscript + supplementary

## Realistic time budget (with 32-core parallelism)

| Phase | Wall-clock | Notes |
|---|---|---|
| 0 schema extension | 1 day | code + tests |
| 1.A human local layers | 1 day | binary overlaps |
| 1.B human downloads | 2-3 days | CADD genome-wide is biggest |
| 1.C human re-run + verify | 2 days | 321 loci re-fine-map |
| 2 yeast | 3-5 days | data smaller |
| 3 arabidopsis | 1 week | DAP-seq is heavy |
| 4 rice | 1 week | downloads + builds |
| 5 consolidation | 3 days | report + paper |
| **Total** | **3-4 weeks** | mostly background |

## Honest expected outcome

Per the heterogeneity-not-density principle, I predict:
1. **Human**: HBP-tighter rate goes from 0/321 to ≥10% (maybe 30-50% with MotifBreakR + ABC) — single biggest win, since the GTEx-only baseline was a uniform cache
2. **Yeast**: HBP-tighter goes from 100% (already saturated) → still 100% but with tighter HBP CS sizes (mean ~880 → ~600?)
3. **Arabidopsis**: HBP-tighter goes from 54% → 70-80% (DAP-seq is highly heterogeneous)
4. **Rice**: HBP-tighter goes from 19% → 35-50% (real ChIP-seq adds variant-scale information)

If these predictions hold, the paper's heterogeneity-not-density claim moves from inferential to empirically demonstrated by direct intervention.

## Execution principles

- Commit per phase
- Memory-persist phase milestones
- Predict before re-running each phase
- Run all heavy downloads in background; do code work in parallel
- Don't promise more than can be delivered in 4 weeks

## Status (live)

- [ ] Phase 0 schema extension
- [ ] Phase 1.A human local layers
- [ ] Phase 1.B human downloads
- [ ] Phase 1.C human re-run + verify
- [ ] Phase 2 yeast
- [ ] Phase 3 arabidopsis
- [ ] Phase 4 rice
- [ ] Phase 5 consolidation + paper
