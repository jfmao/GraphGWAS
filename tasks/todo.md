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

## IN PROGRESS / TODO

### Science (completed in this session)
- [x] Multi-trait genetic correlation on 35 GRAMMAR-corrected traits → mean |r_G|=0.69
- [x] GRAMMAR-corrected eQTL (48 genes) → 40/48 with significant cis-eQTL
- [x] Cross-method integrated analysis on YPDCUSO410MM → GWAS+epistasis+flow+finemap
- [x] PLINK2 benchmark → r(-log10p)=1.0000 on 82K shared variants
- [x] Human chr1 GWAS (5.76M variants, 3202 samples) → 9/9 simulated causal recovered
- [x] Full 1KG import (70.7M variants, 3202 samples, 22 chromosomes)
- [x] Batch genotype loading (cursor pagination) → 3× speedup (4200→5965 variants/sec)
- [ ] Full-genome human GWAS (in progress, ~26M of 70M done, estimated ~1.5 hr remaining)

### Validation Gaps — MUST PLAN BEFORE RUNNING (see bottom of file)
- [ ] **Task B: Epistasis ground truth benchmark** (HIGHEST PRIORITY)
- [ ] **Task A: Rare variant power benchmark**
- [ ] **Task C: Pathway architecture benchmark**
- [ ] **Task D: Fine-mapping benchmark vs SuSiE**
- [ ] **Task E: Fix spectral heritability formula**
- [ ] **Meta: Integrated pipeline showcase** (after B, C, D pass)

### Engineering
- [ ] Git commit (still zero commits — CRITICAL)
- [ ] VEP/SnpEff functional annotations for yeast and human (improve max-flow specificity)
- [ ] GPU-accelerated batch linear regression (test many variants simultaneously)
- [ ] LD pruning implementation (needed for epistasis benchmark)

### Pending from Earlier Phases
- [ ] GNN training on yeast data (requires PyTorch Geometric + CUDA)
- [ ] Agent interactive testing (requires ANTHROPIC_API_KEY)
- [ ] API server integration test
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
