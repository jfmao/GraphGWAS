# GraphGWAS — Task List

## Platform Status: Validated on 1011 Yeast Genomes

**Environment:** Ubuntu 64GB + RTX 4090, Neo4j 5.26 (port 7688), Python 3.13
**Database:** 1,919,547 variants (SNPs + indels), 1,011 samples, 34 populations, 16 chromosomes
**Annotations:** 6,613 genes (SGD), 679 GO pathway terms, 1.18M variant→gene edges, 50K gene→pathway edges

---

## COMPLETED: Phases 0-11 — Core Platform (from prior sessions)
- Phases 0-5: Data foundation, association, graph-native methods, GNN, AI agent
- Phases 6-11: Parallel computing, API/MCP, heritability, multivariate, MET, PRS, MR
- 30 Python modules, 60+ CLI commands

## COMPLETED: Environment Migration (2026-03-31)
- [x] WSL → Native Ubuntu migration (all paths updated)
- [x] Two Neo4j instances: GraphMana/GraphPop (7687) + GraphGWAS (7688)
- [x] Python dependencies installed, graphgwas CLI operational
- [x] N_SAMPLES auto-detection from Neo4j (`detect_n_samples()`)

## COMPLETED: Yeast Data Import (2026-03-31)
- [x] Downloaded original gVCF from 1002genomes (5.8 GB, real nucleotides)
- [x] Minimal filtering: split multi-allelic, AC>0, AN≥100 → 1,919,547 biallelic variants
- [x] GraphMana bulk import into Neo4j (5.5 seconds)
- [x] Genotype decode validated: AF matches to 1e-5 precision across 2,383 test variants
- [x] 35 growth phenotypes loaded (971/1011 strains matched)
- [x] SGD gene annotations loaded (6,613 ORFs, 679 GO terms, 1.18M HAS_CONSEQUENCE edges)
- [x] Expression data loaded (60 genes for cis-eQTL)
- [x] Phenotype data inventory: 35 growth + 6,517 gene expression + 630 protein abundances

## COMPLETED: GWAS Validation (2026-03-31)
- [x] Genome-wide GWAS on YPETHANOL: 1.92M variants in 67s (parallel, 8 workers)
- [x] Validated against published: 7/7 hits matched, ADH4/ALD6/SSU1 in top 100
- [x] 35-trait GWAS scan: 51 min total (uncorrected)
- [x] cis-eQTL mapping: 60 genes, all with significant cis-eQTLs, ADH2 eQTL at 259bp

## COMPLETED: Graph-Native Methods on Yeast (2026-03-31)
- [x] Epistasis: 3 modules on chr6 (2,000 nodes, 326K co-occurrence edges)
- [x] LD fine-mapping: credible set via betweenness centrality (3,000 nodes, 71K LD edges)
- [x] Spectral heritability: h²=0.70 on chr4 (97K variants, 486 samples)
- [x] Max-flow: 218 pathways scored on chr6 (transport, mitochondrion, biosynthetic process)
- [x] Conductance heritability: working on single chromosomes

## COMPLETED: Population Structure Correction (2026-03-31/04-01)
- [x] Graph spectral PCs: Laplacian eigenvectors from 1.2M rare variants → stored on Sample nodes
- [x] GRM PCs: Yang et al. GRM from 151K common variants → stored on Sample nodes
- [x] Lambda computation fixed: now uses common variants only (AF>5%)
- [x] **GRAMMAR+ mixed model integrated**: GRM → REML h² → residualize → calibrate
- [x] **35-trait GRAMMAR-corrected GWAS**: ALL traits lambda 0.90-1.02 (down from 0.85-5.04)
  - Results in: results/grammar_corrected/

## COMPLETED: Performance & Infrastructure (2026-03-31/04-01)
- [x] Parallel GWAS: 5× speedup (8 workers, 67s for 1.92M variants)
- [x] Fixed parallel mode: staggered worker connections, Neo4j pool config
- [x] TSV output default (`-o file.tsv`), `--store` optional for Neo4j
- [x] GPU acceleration: PyTorch CUDA, 39× speedup on GRM matmul (RTX 4090)
- [x] Streaming result writes via store_results_fn callback
- [x] GRM + genotype matrix cached to disk (.npy) for reuse

## COMPLETED: Testing (2026-04-01)
- [x] 33/33 statistics tests pass
- [x] 8/8 simulation GWAS tests pass
- [x] 7/7 yeast offline tests pass

## COMPLETED: Heritability Comparison (2026-04-01)
- [x] Spectral h² computed for all 35 traits
- [x] Compared to published SNP-h² and total h²
- [ ] **Formula needs fixing**: mean spectral h²=0.25 vs published SNP-h²=0.63 (underestimates)
- [ ] Correlation with published: r ≈ -0.07 (should be positive)

---

## COMPLETED (2026-04-22): Wu et al. 2026 baseline integration + paper polish

### Track 3 (head-to-head with state-of-art): DONE
- [x] 3.1 Polyfun-proxy (eQTL→prior_weights to SuSiE), 20 reps, 1KG chr22
- [x] 3.2 Power vs N curves (N=500/1k/2k/3k, 30 reps), 1KG chr22
- [x] 3.3 Cross-ancestry (EUR/AFR/EAS, 30 reps each), 1KG chr22
- [x] 3.4 100-rep weak-signal benchmark, 1KG chr22, all 4 R-baselines
- [x] **NEW: SuSiE-inf + FINEMAP-inf** (Cui et al. 2024) integrated
- [x] **NEW: SBayesRC** (Wu et al. 2026) — pipeline verified on canonical UKB sumstats
- [x] **NEW: Cross-species** Arabidopsis FT10 — PIP=0.989 single variant

### Paper polish: DONE
- [x] Section renumbering 2.1-2.12 (clean integers, no 2.3a/2.3b)
- [x] Refreshed abstract (~230 words covering all 6 baselines + 3 species)
- [x] Updated figure inventory (9 figures all regenerable)
- [x] Updated table inventory (9 tables all regenerable)
- [x] Added refs 5-7 (Cui 2024, Zheng 2024, Wu 2026)
- [x] Updated cover letter (13-claim headline table)
- [x] All cross-references audited

## REMAINING for NG submission

- [ ] Methods §4 prose expansion (~1500 words; currently terse stubs)
- [ ] Author list + affiliations (placeholders in cover letter + manuscript)
- [ ] Reviewer suggestions (3-5 names for cover letter)
- [ ] Final manual proofread

## ORTHOGONAL track

- [ ] UK Biobank application (paperwork; ~4-6 weeks approval)
- [ ] Optional: SBayesRC head-to-head on our 1KG chr22 aggregated sumstats
- [ ] Optional: tighter Methods-level pseudocode for HBP/M1

## COMPLETED (2026-04-20/21): UKB-Ready Hybrid Architecture + Paper Assets

### Track 1 (hybrid genotype): DONE
- [x] 22 BGEN files built from 1KG VCFs (~14 GB)
- [x] `src/python/graphgwas/bgen_reader.py` (AF matches Neo4j to 1e-4)
- [x] `load_locus_variants()` dispatcher in finemapping_v2.py (Neo4j or BGEN)
- [x] `src/python/graphgwas/summary_import.py` (PLINK2 + regenie parsers)
- [x] End-to-end: PLINK2 GWAS → sumstats → BGEN HBP → graph query (1569/1681 hits)
- [ ] M1 epistasis BGEN refactor (still Neo4j-only)

### Track 2 (multi-omics annotation graph): DONE
- [x] GENCODE v47: 20,092 protein-coding genes
- [x] 38.4M HAS_CONSEQUENCE proximity edges (67s with 6 workers)
- [x] GTEx v8 49 tissues: **43.2M eQTL edges** (8.4 min with 4 workers)
- [x] STRING PPI score≥700: 230,850 INTERACTS_WITH edges
- [x] ENCODE cCRE: 370K RegulatoryElement nodes (partial, 40%)
- [ ] PhyloP conservation (deferred)

### Paper artifacts: DONE
- [x] 7 figures (fig1..7) in results/benchmark_v2/paper_figures/
- [x] 5 tables (table1..5) in results/benchmark_v2/paper_tables/
- [x] docs/PAPER_DRAFT_MANUSCRIPT.md (~2000 words, abstract + 6 result sections)
- [x] docs/PAPER_COVER_LETTER.md (~400-word NG cover letter)
- [x] docs/REPRODUCIBILITY.md (regen every artifact from scratch)
- [x] docs/MATHEMATICAL_PROOFS.md (5 theorems, supplement S1)
- [x] 60/60 unit tests pass

## IN PROGRESS / TODO

### Science (completed in this session)
- [x] Multi-trait genetic correlation on 35 GRAMMAR-corrected traits → mean |r_G|=0.69
- [x] GRAMMAR-corrected eQTL (48 genes) → 40/48 with significant cis-eQTL
- [x] Cross-method integrated analysis on YPDCUSO410MM → GWAS+epistasis+flow+finemap
- [x] PLINK2 benchmark → r(-log10p)=1.0000 on 82K shared variants
- [x] Human chr1 GWAS (5.76M variants, 3202 samples) → 9/9 simulated causal recovered
- [x] Full 1KG import (70.7M variants, 3202 samples, 22 chromosomes)
- [x] Full-genome human GWAS (70.7M variants, 3.5 hours, 14/21 causal GW-sig)
- [x] Epistasis v2: M1-M4 implemented, benchmarked, M1 rank #1 on all replicates
- [x] Fine-mapping v2: L1/L4 implemented, L1 CS=4.4, L4 multi-signal
- [x] 10 comprehensive benchmarks completed
- [x] PLINK interaction comparison (42,000× search reduction demonstrated)
- [x] v1 vs v2 comparison (L1 beats v1 centrality, M1 220× more precise than v1 epistasis)
- [x] Git commit + push to github.com/jfmao/GraphGWAS

---

## NEXT PHASE: Preparing for Nature Genetics Submission

### Priority 1: Head-to-Head with State-of-Art Tools (2-3 weeks)

**Status: PARTIALLY DONE**

**Fine-mapping — L1 vs SuSiE (COMPLETED, 20 replicates):**
- [x] SuSiE installed (R 4.5.3, susieR 0.14.2 via conda)
- [x] FINEMAP v1.4.2 installed (binary)
- [x] L1 (basic) vs SuSiE benchmark: SuSiE wins (rank 1.2 vs 13.1, CS 3.8 vs 9.7)
- [x] Multi-omics annotations loaded (STRING PPI 518 edges, eQTL 48K variants, conservation 1M variants)
- [x] L1 enhanced with multi-layer functional scoring (eQTL + PPI + conservation + gene/pathway)
- [x] **L1 vs SuSiE vs FINEMAP 3-way comparison (2026-04-08)**
  - FINEMAP bug fixed: was using r² matrix, FINEMAP needs signed correlation
  - Strong signal (β=0.5, 30 reps): L1 63% rank#1, FINEMAP 63%, SuSiE 57%
  - Weak signal (β=0.2, 30 reps): L1 50% rank#1, FINEMAP 57%, SuSiE 50%
  - **L1 is 25-35× faster** (0.07s vs 1.8-2.5s per locus)
  - **L1 wins on weak signal + annotations** (Section 8: 13.5:1 win ratio)
  - All results in results/benchmark_v2/tool_comparison/
- [x] **Annotation alpha tuning**: optimal α=0.9 (90% stat + 10% annotation)
  - α=0.5 too aggressive — annotations dilute statistical signal
- [ ] Polyfun+SuSiE comparison — not started

**Epistasis — install and compare to:**
- [x] FAME assessment: FAME estimates variance components (σ²), not per-variant PIPs.
  Not directly comparable to fine-mapping tools. Would compare to LDSC/S-LDSC instead.
- [ ] BOOST — boolean pairwise epistasis (not yet installed)
- [x] PLINK --glm interaction (done — GraphGWAS M1 discovers pairs, PLINK only confirms)
- [ ] All on same 1KG chr22 S1 simulations (100 replicates)

**Metrics (matching GWFM paper standard):**
- [x] PIP calibration (TDR vs PIP bins) — 200 sims × 4 h² levels, all methods calibrated
- [x] Null simulation FPR — 100 nulls, 0% FPR all methods
- [x] Runtime comparison — HBP 0.08s vs SuSiE 1.8s vs FINEMAP 2.5s
- [ ] Power at same FPR (optional — covered by rank-based metrics)
- [ ] Credible set size at same coverage (optional)

### Priority 2: Mathematical Theory — COMPLETED (2026-04-10)

- [x] 5 theorems with full proofs in `docs/MATHEMATICAL_PROOFS.md`
- [x] Theorem 1: M1 search space reduction ≤ k(τ)² (proves 42,000×)
- [x] Theorem 2: HBP convergence via Banach fixed-point (geometric rate)
- [x] Theorem 3: L1 causal variant ranking under LD decay
- [x] Theorem 4: Null PIP bound O(√(log n)/n) (proves 0% FPR)
- [x] Theorem 5: CLGF EM convergence (monotone, O(log 1/ε))

### Graph-Native Fine-Mapping v3 — COMPLETED (2026-04-09)

- [x] **HBP** (Hierarchical Belief Propagation) — best graph-native method
  - Matches SuSiE accuracy (4:4:22 on strong signal), 20× faster
  - fast_hbp_finemap() with pre-cached graph: 0.08s/locus
- [x] **GRSD** (Graph-Regularized Sparse Deconvolution) — dropped (underperforms)
- [x] **CLGF** (Cross-Locus Graph Fine-Mapping) — inconclusive on sparse chr22 graph
- [x] F6 multi-locus simulation added to simulate.py
- [x] 90-rep HBP benchmark, 100-rep null FPR, 200-rep PIP calibration

### Real-Data Biological Discovery — COMPLETED (2026-04-09)

- [x] HBP fine-mapping on 50 yeast loci across 5 traits
- [x] Known genes recovered: TUP1, SPT7 (ethanol), PPG1, BMH2 (heat), GCR1, HXT (galactose), YRR1 (caffeine)
- [x] Novel candidates: RNQ1 (copper), GET4 (heat), MOT3 (galactose), ARC35 (caffeine)
- [x] Report in results/yeast_discovery/YEAST_DISCOVERY_REPORT.md

### Priority 3: Multi-Omics Integration for L1 (3-4 weeks)

**Load public multi-omics annotations into Neo4j:**
- [ ] GTEx eQTL data (~2M eQTLs across 49 tissues)
- [ ] ENCODE regulatory elements (enhancers, promoters)
- [ ] STRING protein-protein interactions (~600K edges in human)
- [ ] OpenTargets drug targets
- [ ] PhyloP/GERP conservation scores

**Enhance L1 fine-mapping:**
- [ ] Multi-layer functional scoring via graph traversal depth
- [ ] Tissue-specific traversal weighting
- [ ] Benchmark L1-multiomics vs SuSiE vs SBayesRC on same simulations
- [ ] **Show L1 improves specifically when annotations are rich** (key result)

### Priority 4: Rigorous Simulations (2 weeks)

- [ ] Expand to 100 replicates per scenario (from current 3-5)
- [ ] PIP calibration plots (TDR vs PIP, matching GWFM paper Fig. 2)
- [ ] Power vs coverage curves (matching GWFM paper Fig. 3)
- [ ] Null simulation FPR (1000 null replicates)
- [ ] Vary sample sizes: N = 1K, 5K, 10K, 50K, 100K
- [ ] Multiple genetic architectures (sparse, large-effects, LDMS)

### Priority 5: Code Polish + Documentation (1-2 weeks)

- [ ] pip-installable package
- [ ] README with quick-start guide
- [ ] Reproducible benchmark scripts in repo
- [ ] Integration tests for v2 methods

---

## UK Biobank Preparation (start application NOW, ~4-6 weeks approval)

### Hybrid Architecture for UKB Scale

GraphGWAS cannot store 500K-sample genotypes in Neo4j (1.6TB gt_packed).
**Solution: hybrid file + graph architecture.**

**Neo4j stores (small, fast):**
- GWAS summary statistics (imported from PLINK2/regenie)
- Multi-omics annotations (eQTL, chromatin, Hi-C, PPI, drug targets)
- Gene / Pathway / GO / RegulatoryElement nodes and relationships
- LD structure for candidate loci (precomputed)

**External files store (large, streamed on demand):**
- Genotypes in .bgen/.bed format
- Loaded per-locus for M1-M4 epistasis and L1/L4 fine-mapping
- Never stored in Neo4j

**Implementation needed:**
- [ ] `graphgwas/bgen_reader.py` — load genotypes from .bgen files per-locus
- [ ] Add `genotype_source` parameter to M1, M3, M4, L1, L4 (default "neo4j", option "bgen")
- [ ] `graphgwas/summary_import.py` — import PLINK2/regenie summary stats into Neo4j
- [ ] Test hybrid pipeline on 1KG .bgen files before UKB access

**Hardware requirements:**
- Current machine (64GB): works for candidate-region analyses
- Cloud/HPC (256GB): recommended for full-scale UKB analysis
- Storage: 2TB NVMe SSD for UKB .bgen files

### What To Do Before UKB Access
1. [ ] Load GTEx eQTL data into Neo4j (~2M eQTL across 49 tissues)
2. [ ] Load ENCODE regulatory elements (enhancers/promoters)
3. [ ] Load STRING protein-protein interactions (~600K edges)
4. [ ] Enhance L1 multi-layer functional scoring via graph traversal
5. [ ] Benchmark L1-multiomics vs SuSiE vs SBayesRC on 1KG simulations
6. [ ] Show that L1 improves specifically when annotations are rich
7. [ ] Implement bgen_reader.py module for hybrid architecture
8. [ ] Test hybrid pipeline on 1KG .bgen files

### When UKB Access Arrives
- [ ] Run PLINK2/regenie GWAS on 48 traits (external, fast)
- [ ] Import summary stats into Neo4j
- [ ] Run M1 epistasis on significant loci → discover novel interactions
- [ ] Run L1-multiomics fine-mapping → find causal variants missed by SuSiE
- [ ] Replicate findings in independent UKB subsample
- [ ] Write paper

### Pending from Earlier Phases
- [ ] GNN training on yeast data (requires PyTorch Geometric + CUDA)
- [ ] Agent interactive testing (requires ANTHROPIC_API_KEY)
- [ ] API server integration test
- [ ] Implement remaining epistasis methods: M5 (random walk), M6 (GNN), M7 (topological)
- [ ] Implement remaining fine-mapping methods: L2 (spectral), L3 (haplotype)
- [ ] MCP server integration test

---

## CRITICAL: Validation Gaps & Benchmark Plan

### Current Validation Status

| What | Validated? | Method |
|------|-----------|--------|
| Single-locus GWAS correctness | **YES** | r=1.000 vs PLINK2, simulated causal recovery |
| GRAMMAR+ calibration | **YES** | Lambda 0.90-1.02 across 35 yeast traits |
| Epistasis detection | **UNVALIDATED — likely broken** | Modules are huge (687 variants each), LD not pruned, no ground truth test |
| Max-flow pathway scoring | **Functional only** | Biologically plausible, NOT compared to MAGMA |
| LD fine-mapping | **Functional only** | Credible set produced, NOT compared to SuSiE |
| Spectral heritability | **Broken** | Underestimates (h²=0.04 vs published 0.63) |
| Rare variant power | **UNVALIDATED — never tested** | All simulated "causal" variants had AF > 2%; single-locus is underpowered for truly rare variants |
| Burden/aggregation tests | **NOT run** | SKAT, MPAT not benchmarked on rare variants |
| Power advantage | **NOT tested** | No comparison of detection power vs PLINK |
| Runtime | **Measured** | 3× faster with cursor pagination; still slower than PLINK for single-locus |

### Benchmark Plan (Priority Order)

#### Priority 1 — Proves the core thesis (graph-native methods outperform)

**1a. Epistasis benchmark with known ground truth (CRITICAL)**
- Simulate phenotype with known epistatic pairs (two variants with no marginal effect but strong joint effect)
- **Must include LD pruning step** (r² < 0.2) — current modules are likely capturing LD blocks not interactions
- Test: does GraphGWAS community detection recover the interacting pair?
- Compute interaction p-value via permutation of case labels
- Compare to: PLINK --epistasis (pairwise), BOOST, MDR
- Metrics: sensitivity, false positive rate, runtime
- **KNOWN GAPS in current implementation:**
  - No LD pruning → modules include LD blocks
  - No interaction term test (enrichment ≠ epistasis)
  - No statistical significance for modules
  - Module sizes are 100-700 variants (should be 2-10 for real interactions)

**1b. Rare variant power comparison (CRITICAL — never tested)**
- Simulate phenotypes with causal variants at AF ∈ {0.001, 0.005, 0.01, 0.05}
- Effect sizes inversely proportional to AF (rare → larger β)
- Compare detection methods:
  (a) Single-locus linear regression (baseline — expected to fail at AF<1%)
  (b) Burden test per gene (aggregate rare variants)
  (c) Our 1/AF-weighted max-flow (graph-native)
  (d) SKAT / SKAT-O (gold standard)
- Metrics: power at each AF bin, false positive rate, minimum sample size for 80% power
- **KNOWN GAPS:**
  - All our simulated "causal" variants have been AF > 2% (not truly rare)
  - Never tested burden tests on human data
  - Never benchmarked 1/AF max-flow weighting
  - Never compared to SKAT

**1c. Pathway benchmark with known architecture**
- Simulate phenotype where all causal variants are in one pathway
- Test: does max-flow rank the correct pathway #1?
- Compare to: MAGMA gene-set analysis, PASCAL
- Metrics: correct pathway rank, AUC for pathway recovery

#### Priority 2 — Correctness of novel methods

**2a. Fine-mapping benchmark**
- Simulate causal variant in an LD block on real 1KG data
- Test: does centrality place the causal variant in the credible set?
- Compare to: SuSiE, FINEMAP
- Metrics: credible set size, causal variant rank, PIP comparison

**2b. Heritability benchmark**
- Fix spectral h² formula
- Compare GRM-REML h² to GCTA-GREML on same data
- Also compare spectral h² to LDSC
- Metrics: h² correlation with published, bias, MSE

#### Priority 3 — Practical value

**3a. Integrated analysis value**
- Run GWAS+epistasis+flow+finemap on same trait with known architecture
- Show that the integrated analysis recovers the true architecture
- Compare to: sequential PLINK → MAGMA → SuSiE pipeline

---

## The 7 Novel Claims — Full Validation Scorecard

Must be validated AFTER current full-genome scan completes. **IMPORTANT: PLAN each benchmark before running — design the simulation, specify metrics, define comparison tools, estimate runtime. Do NOT start runs before planning is reviewed.**

| # | Novel Claim | Status | Blocking Benchmark |
|---|-------------|--------|-------------------|
| 1 | **Indels/rare variants as first-class** | Partial (indels ✓, rare variants ✗) | → Task A |
| 2 | **Graph-queryable results** | Demonstrated ✓ | Complete |
| 3 | **Unified GRAMMAR+ calibration** | Demonstrated ✓ (but not graph-unique) | Complete |
| 4 | **Epistasis via community detection** | Runs but unvalidated | → Task B |
| 5 | **Max-flow disease architecture** | Runs but unvalidated | → Task C |
| 6 | **LD fine-mapping via centrality** | Runs but unvalidated | → Task D |
| 7 | **Spectral heritability from rare-variant graph** | BROKEN (h²=0.04 vs 0.63 pub) | → Task E |

### Task A: Rare Variant Power Benchmark (Claim 1)
**Plan before running.** Design:
- Simulate phenotype with causal variants at AF ∈ {0.001, 0.005, 0.01, 0.05}
- Effect sizes inverse to AF (rare → larger β)
- Methods to compare: single-locus regression, burden test per gene, 1/AF-weighted max-flow, SKAT
- Sample sizes: 1KG (N=3,202) and downsampled (N=500, 1000, 2000)
- Metrics: power at each AF bin, minimum N for 80% power, FPR at p<5e-8
- Runtime estimate: ~4-8 hours

### Task B: Epistasis Ground Truth Benchmark (Claim 4) — HIGHEST PRIORITY
**Plan before running.** Design:
- Simulate phenotype Y = β × (G1 × G2) + noise (pure interaction, no marginals)
- Place causal pair in independent genomic regions (no LD)
- **Add LD pruning step to epistasis pipeline** (r² < 0.2) — current code doesn't have this
- **Add interaction term statistical test** — current code tests enrichment, not interaction
- Test: does community detection place G1 and G2 in the same module?
- Compute p-value via case-label permutation
- Compare to: PLINK --epistasis pairwise testing, BOOST, MDR
- Metrics: sensitivity, FPR, runtime, module size distribution (should be 2-10, not 100-700)
- **Implementation work needed:** LD pruning + interaction term testing BEFORE benchmark
- Runtime estimate: ~1-2 days including code updates

### Task C: Pathway Architecture Benchmark (Claim 5)
**Plan before running.** Design:
- Simulate phenotype where ALL causal variants are in one known pathway
- Example: KEGG glycolysis pathway on yeast
- Test: does max-flow rank the causal pathway #1?
- Compare to: MAGMA gene-set analysis, PASCAL
- Metrics: rank of causal pathway, AUC for pathway recovery, top-10 pathway overlap
- Runtime estimate: ~2-4 hours

### Task D: Fine-Mapping Benchmark (Claim 6)
**Plan before running.** Design:
- Simulate causal variant in real 1KG LD block (use chr22 which we already have)
- Multiple replicates with different causal positions
- Compare centrality credible set vs SuSiE posterior inclusion probability
- Metrics: credible set size, causal variant rank in credible set, PIP comparison, runtime
- Runtime estimate: ~2-3 hours

### Task E: Spectral Heritability Fix (Claim 7)
**Plan before implementing.** Known issues to investigate:
- Formula gives h² ≈ 0.04 vs published 0.63
- GRM normalization fixed but REML still underestimates
- Rare-variant-only approach may be fundamentally limited with SNP-only data
- Must compare to GCTA-GREML and LDSC with same variant set and model
- Need to decide: is the spectral approach conceptually sound, or should we recommend GRAMMAR h² instead?
- Implementation: likely needs finer REML grid (0.001 step), different variant weighting, or different kernel
- Runtime estimate: ~4-8 hours investigation + implementation

### Meta-Benchmark: Integrated Pipeline Demonstration
**Plan before running.** Only after Tasks B-D pass individually:
- Single trait simulation with architecture: marginal effects + epistasis + pathway structure
- Run GWAS → epistasis → max-flow → fine-mapping on same data
- Show integrated pipeline recovers full architecture
- Compare to sequential PLINK → MAGMA → SuSiE pipeline
- This is the **paper-ready showcase** demonstrating GraphGWAS's unique value

---

## CRITICAL HONEST ASSESSMENT

Two central claims of GraphGWAS remain UNPROVEN after extensive validation:

### Claim 1: "Graph-native epistasis detection finds interactions matrix tools miss"
**Evidence so far:** Functional demonstration only. Modules detected, but:
- Module sizes (100-700 variants) are inconsistent with tight epistatic pairs
- No LD pruning — detected modules likely reflect LD blocks, not interactions
- Never tested with known epistatic ground truth
- Never computed interaction term vs marginal effects

### Claim 2: "Rare variants and indels are first-class citizens with stronger signal"
**Evidence so far:** Database correctly includes rare variants. But:
- Single-locus regression on rare variants has inherently low power (same as PLINK)
- All simulated "causal" variants have been AF > 2% (not truly rare)
- Burden tests / SKAT / MPAT never benchmarked on rare variants
- Our "rare variant hits" (AF 2-5%) are low-frequency common, not rare

**Until Priority 1a and 1b benchmarks are completed, both claims should be marked as aspirational, not demonstrated.**

The claims that ARE demonstrated:
- Single-locus GWAS produces identical results to PLINK (r=1.000)
- GRAMMAR+ mixed model achieves perfect calibration on all 35 yeast traits
- Max-flow + fine-mapping + epistasis run on real data and produce biologically plausible output
- Full 1KG genome (70M variants) can be processed through the pipeline

#### Priority 4 — Runtime and scalability

**4a. Runtime comparison table (honest)**
| Operation | GraphGWAS | PLINK2/GCTA | Ratio | Justification |
|-----------|-----------|-------------|-------|---------------|
| Single-locus GWAS | 57 min (chr1) | ~2 sec | 1,700× slower | Graph overhead; enables novel methods |
| Epistasis (100kb) | 5 sec | hours (pairwise) | 100-1000× faster | Graph sparsity |
| Fine-mapping (40kb) | 10 sec | ~1 hour (SuSiE) | 360× faster | Centrality vs Bayesian |
| Pathway scoring | 145 sec | ~30 sec (MAGMA) | 5× slower | But unified with GWAS |

**4b. Scalability test**
- Run on 1KG full genome (70M variants) when import completes
- Measure: total wall time, peak memory, per-variant throughput

---

## Key Findings

### GRAMMAR+ Calibration (all 35 traits)
| Metric | Uncorrected | GRAMMAR+ Corrected |
|--------|------------|-------------------|
| Lambda range | 0.85-5.04 | 0.90-1.02 |
| Traits with λ>1.1 | 28/35 (80%) | 0/35 (0%) |
| Mean lambda | 1.62 | 0.98 |

### Top Heritability Estimates (GRAMMAR h²)
| Trait | h² | Biology |
|-------|-----|---------|
| YPDCUSO410MM | 0.28 | Copper sulfate resistance |
| YPD42 | 0.20 | Heat tolerance (42°C) |
| YPACETATE | 0.19 | Acetate utilization |
| YPGALACTOSE | 0.19 | Galactose utilization |
| YPD40 | 0.19 | Heat tolerance (40°C) |

---

## Module Inventory: 30 modules

| Module | Purpose |
|--------|---------|
| config.py | Constants, N_SAMPLES auto-detection |
| db.py | Neo4j connection (pool config, connection timeout) |
| schema.py | Indexes, schema audit |
| phenotype.py | Phenotype import/activation |
| genotype.py | gt_packed decode, PCA import, covariates |
| qc.py | HWE, MAC, missingness, het_rate |
| assoc.py | Chi2, Fisher, logistic, Firth, linear (quantitative trait fix) |
| results.py | AssociationResult nodes, link_study() |
| mpat.py | Gene-level MPAT |
| plots.py | Manhattan, QQ, TSV export |
| epistasis.py | Co-occurrence network, Louvain, permutation |
| flow.py | Max-flow/min-cut disease architecture |
| finemapping.py | LD graph centrality fine-mapping |
| spectral.py | Rare-variant similarity, Laplacian (fixed: uses all samples) |
| gnn.py | HeteroGNN export/train/explain/import |
| agent.py | LangGraph ReAct agent |
| parallel.py | ProcessPoolExecutor, staggered connections |
| api.py | FastAPI REST server |
| mcp_server.py | MCP tool server |
| heritability.py | 6 graph-native heritability estimators |
| multivariate.py | Cross-trait: r_G, G-matrix, coherence, pleiotropy |
| met.py | MET: G×E, env similarity |
| prs.py | PRS: classical, graph-pruned, pathway |
| mr.py | MR: IVW, Egger, weighted median |
| popstruct.py | GRM, PCA, GRAMMAR+, graph-spectral PCs (GPU-accelerated) |
| sumstats.py | Summary statistics |
| sv.py | Structural variant analysis, fast_gwas |
| cli.py | Click CLI (60+ commands, popstruct group added) |
