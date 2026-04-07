# GraphGWAS — Nature Genetics Gap Analysis

## What Nature Genetics Demands (from reviewing recent papers)

### Paper 1: Guan et al. 2025 — "Family-based GWAS designs" (Technical Report)
**What they did:** New statistical estimator for family-based GWAS
**What made it NG-worthy:**
- Mathematical proof of increased effective sample size (46.9% to 106.5%)
- Applied to UK Biobank (408K individuals, 19 phenotypes)
- Superior out-of-sample polygenic prediction vs alternatives
- Released as software package (snipar)

### Paper 2: Wu et al. 2026 — "Genome-wide fine-mapping" (Article)
**What they did:** Genome-wide fine-mapping using SBayesRC
**What made it NG-worthy:**
- Applied to 599 complex traits × 13 million SNPs in UK Biobank
- Identified 19,863 credible sets across 48 traits
- 30% of credible sets were OUTSIDE genome-wide significant loci (novel finding)
- Power prediction framework for future studies
- Found new causal variants missed by SuSiE/FINEMAP (FTO locus, ACTRIB for SCZ)
- Replication in independent samples

### Paper 3: FAME 2025 — "Marginal epistasis" (Nature Genetics)
**What they did:** Biobank-scale epistasis test
**What made it NG-worthy:**
- Computationally efficient (300K individuals)
- First evidence of interaction effects between individual variants and polygenic background
- 16 significant epistasis signals across 12 traits
- Epistatic variance 12× larger than GWAS effects

### Common patterns across these NG papers:

| Requirement | Paper 1 | Paper 2 | Paper 3 |
|-------------|---------|---------|---------|
| **Biobank-scale data** | UK Biobank (408K) | UK Biobank (13M SNPs, 599 traits) | UK Biobank (300K) |
| **Mathematical/statistical novelty** | New estimators with proofs | New α-LCS construction | New MoM estimator |
| **Extensive simulations** | Multiple population structure scenarios | 3 genetic architectures × 100 replicates | Detailed power analysis |
| **Real-data biological discovery** | 19 phenotype results | New causal variants (FTO, ACTRIB) | 16 epistasis signals |
| **Comparison to existing methods** | vs sib-difference, Young et al. | vs SuSiE, FINEMAP, Polyfun+SuSiE | vs previous epistasis tests |
| **Independent replication** | Millennium Cohort Study | Independent replication sample | Cross-trait validation |
| **Software release** | snipar package | SBayesRC (public) | FAME (public) |
| **Practical impact** | Increases effective sample size | Finds variants outside GWAS loci | First systematic epistasis signals |

---

## Honest Assessment: Where GraphGWAS Stands

### What We HAVE

| Requirement | GraphGWAS Status | NG Standard |
|-------------|-----------------|-------------|
| Novel methods | 6 graph-native methods (M1-M4, L1, L4) | ✅ Multiple novel approaches |
| Mathematical novelty | Co-occurrence + community detection, LD deconvolution + functional traversal | ⚠️ Methods are heuristic, not mathematically proven |
| Simulations | 10 benchmarks, 5 replicates, S1-S5 + F1-F2 | ⚠️ Need more replicates, more scenarios |
| Comparison to baselines | vs PLINK (r=1.000), vs v1 (220× more precise) | ⚠️ Missing: vs SuSiE, vs FINEMAP, vs FAME |
| Software release | GitHub repo, 30 modules, CLI | ✅ Comprehensive |
| Real-data biological discovery | Yeast ethanol (ADH4/ALD6/SSU1), copper (transport/membrane) | ❌ **No novel human biological finding** |
| Biobank-scale data | 1KG 70.7M variants (3,202 samples) | ❌ **Too small — NG expects 100K-500K samples** |
| Independent replication | None | ❌ **Critical gap** |
| Cross-ancestry validation | None | ❌ **Expected for NG** |

### What We're MISSING (Critical Gaps for NG)

#### 1. BIOBANK-SCALE APPLICATION (the biggest gap)
All recent NG methods papers are applied to **UK Biobank (100K-500K samples)** or equivalent.
Our largest dataset is 1KG with 3,202 samples — a population genetics reference panel, not a GWAS cohort.

**To fix:** Apply GraphGWAS to UK Biobank (requires data access application, ~4 weeks).
Alternatively: apply to a large publicly-accessible GWAS dataset with summary statistics.

#### 2. NOVEL BIOLOGICAL DISCOVERY IN HUMAN DATA
NG requires discovering something biologically new — not just showing the method works on simulations.
The GWFM paper found new causal variants at FTO and ACTRIB. FAME found 16 epistasis signals.

**To fix:** Apply epistasis M1/M2 to a well-powered human dataset and find novel interacting loci that explain missing heritability or elucidate disease mechanisms.

#### 3. INDEPENDENT REPLICATION
Every NG paper validates findings in an independent cohort.

**To fix:** Find an interaction in one dataset, replicate in another (e.g., discovery in UKB White British, replication in UKB non-British or another cohort).

#### 4. HEAD-TO-HEAD WITH STATE-OF-ART
GWFM paper compared SBayesRC to SuSiE, FINEMAP, Polyfun+SuSiE on the same simulations.
We compared M1 to PLINK, but never to FAME, BOOST, or SuSiE.

**To fix:** Install SuSiE/FINEMAP, run on same simulations. Install FAME, compare epistasis detection.

#### 5. MATHEMATICAL RIGOR
NG methods papers include theoretical proofs (convergence, bias, consistency).
Our methods are algorithmic/heuristic — no formal statistical theory.

**To fix:** Develop theoretical bounds for:
- M1: expected false positive rate under null
- L1: conditions under which functional annotations improve credible set precision
- Power analysis: minimum sample size for epistasis detection at given effect size

---

## What Would Make GraphGWAS NG-Worthy

### Option A: Methods paper (Technical Report format)

**Title:** "Graph-native epistasis detection and fine-mapping with 42,000-fold search reduction"

**Story:**
1. Graph representation enables novel epistasis detection impossible in matrix tools
2. LD-aware co-occurrence + interaction regression finds true pairs as #1 hit
3. Dual-graph fine-mapping (LD + functional) achieves CS of 4.4 variants
4. Dark matter method detects synthetic lethality invisible to all existing tools
5. Applied to UK Biobank → N novel epistatic interactions discovered

**What's needed:**
- [ ] UK Biobank application (or large public dataset)
- [ ] Head-to-head with FAME on same data
- [ ] Head-to-head with SuSiE/FINEMAP on same loci
- [ ] Independent replication of top findings
- [ ] Theoretical bounds on search space reduction

### Option B: Resource/Application paper (Article format)

**Title:** "GraphGWAS: a graph-native platform for multi-scale genetic architecture analysis"

**Story:**
1. Complete platform: GWAS + epistasis + fine-mapping + pathway + heritability in one graph
2. Validated identical to PLINK for single-locus (r=1.000)
3. Enables analyses impossible in matrix tools (co-occurrence, max-flow, spectral)
4. Applied to 1011 yeast genomes: 35 traits, 60 eQTLs, 35 genetic correlations
5. Discovered novel epistatic modules in copper resistance
6. Applied to 1KG human genome: 70.7M variants, all causal recovered

**What's needed:**
- [ ] Deeper yeast biological story (validate epistatic modules experimentally or computationally)
- [ ] KEGG/Reactome pathways for max-flow (not just GO slim)
- [ ] Cross-species comparison (same method on yeast + human)
- [ ] More polished analysis report with publication-quality figures

### Option C: Biological discovery paper (Article format)

**Title:** "Graph-native analysis reveals epistatic architecture of complex traits in yeast"

**Story:** Focus on BIOLOGY, not methods. Use GraphGWAS as the tool, not the subject.
1. Applied novel graph-native methods to 1011 yeast genomes
2. Found epistatic modules for copper resistance, ethanol tolerance
3. Max-flow reveals pathway architecture (transport, membrane)
4. Spectral analysis shows rare-variant structure captures more heritability
5. Multi-trait analysis reveals correlated genetic architecture across drug resistances

**What's needed:**
- [ ] Validate findings experimentally (or with independent yeast data)
- [ ] Connect findings to known yeast biology literature
- [ ] Go deeper on 2-3 traits rather than shallow on 35
- [ ] Publication-quality figures

---

## My Recommendation

**Option A (Methods + UK Biobank) is the strongest path to Nature Genetics**, but requires UK Biobank access.

**Option C (Biological discovery in yeast) is the most achievable with current data**, but may target Nature Communications or Genome Biology rather than NG.

**Option B (Platform paper) is comprehensive but risks being "too broad, not deep enough"** for NG. Better suited for Genome Research or Bioinformatics.

### Realistic timeline to NG-quality:

| Phase | Duration | Deliverable |
|-------|----------|-------------|
| UK Biobank application | 4-6 weeks | Data access |
| Implement on UKB data | 4 weeks | Epistasis + fine-mapping on 48 traits |
| Novel biological discoveries | 4 weeks | 5-10 validated epistatic loci |
| Head-to-head benchmarks | 2 weeks | vs FAME, SuSiE, FINEMAP |
| Independent replication | 2 weeks | Replication in second cohort |
| Mathematical theory | 4 weeks | Power bounds, FPR proofs |
| Writing + figures | 4 weeks | Manuscript |
| **Total** | **~6 months** | NG submission |

### Minimum viable path (without UKB):

Apply to the largest publicly accessible individual-level GWAS dataset:
- **FinnGen** (summary statistics only — insufficient)
- **OpenGWAS/MRC IEU** (summary statistics)
- **1KG + Pan-UKBB summary stats** (combine individual genotypes + external summary stats)
- **Large yeast dataset (8,391 traits)** → deep yeast biological story

The yeast path (Option C) could work for **Genome Biology** or **Nature Communications** — still high-impact, and achievable with current data.
