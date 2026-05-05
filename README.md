# GraphGWAS

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Python 3.11+](https://img.shields.io/badge/python-3.11+-blue.svg)](https://www.python.org/)
[![Citation](https://img.shields.io/badge/cite-CITATION.cff-green.svg)](CITATION.cff)

**Relational fine-mapping of causal GWAS variants on a multi-omics knowledge graph**

GraphGWAS is a graph-native fine-mapping platform built on Neo4j. It carries
multi-omics biological structure — genes, tissue-specific eQTLs, pathways,
protein–protein interactions — *through* the fine-mapping inference as a typed
factor graph, rather than collapsing it to flat per-variant annotation priors
as existing Bayesian fine-mappers do. This relational prior matches the accuracy
of SuSiE / FINEMAP / SuSiE-inf / FINEMAP-inf / SBayesRC at 6–60× the speed under
strong signal, and wins **27–2 head-to-head against SuSiE at weak signal with
tissue-specific eQTL priors**.

## Key features

- **Two new fine-mapping algorithms** with theoretical guarantees
  - **HBP** — hierarchical belief propagation on a variant→gene→pathway factor
    graph with PPI coupling; proved Banach contraction (Theorem 2); 0.02–0.08 s
    per locus
  - **GAFM** (Graph-Augmented Fine-Mapping) — LD-deconvolved evidence combined
    with a graph functional score via adaptive α; proved causal-variant ranking
    under mild LD-decay assumptions (Theorem 3)
- **v0.1.5 mixture-prior variants** — **GAFM-MX**, **HBP-MX** and their
  ensemble **ENS** add a SBayesRC-style 4-component Wakefield mixture-BF
  posterior reweighting on the LD-deconvolved z-scores, plus standard
  λ_GC deflation. On the 3kRG grain weight + shape panel they each reach
  10/21 (47.6%) top-1-PIP exact-position recovery against the Niu 2021
  21-QTN catalogue — the highest of any method tested, beating SuSiE
  (28.6%) and SBayesRC (14.3%) — while remaining 200–700× faster than
  SuSiE per locus
- **Six head-to-head baselines** integrated into a common interface — SuSiE,
  FINEMAP, SuSiE-inf, FINEMAP-inf, PolyFun-proxy, SBayesRC
- **Calibrated PIPs** for base GAFM/HBP with 0% null false-positive rate
  across 100 simulations (mixture-prior variants are operational ranking
  scores; null FPR ≤1% at PIP ≥ 0.9, ≤10% at PIP ≥ 0.5)
- **Multi-omics graph** — 70.7 M variants, 20,092 GENCODE genes, 43.2 M
  GTEx v8 tissue eQTLs, 230,850 STRING interactions (combined score ≥ 700),
  370,000 ENCODE cCREs
- **Biobank-scale** — sumstats-only entry path consumes Pan-UK Biobank summary
  statistics directly via tabix over HTTPS; demonstrated on 4 ancestries
  (EUR N = 420,531; CSA, AFR, EAS)
- **Cross-species** — same codebase applies to yeast, human, *Arabidopsis*
- **Unified package** with 52-command CLI, 37-endpoint FastAPI server, and
  16-tool MCP server for AI-agent access

## Quick start

```bash
# Install
git clone https://github.com/jfmao/GraphGWAS.git
cd GraphGWAS/src/python && pip install -e '.[all]'

# Run fine-mapping from Pan-UKB summary statistics (no Neo4j required)
python -c "
from graphgwas.panukb import fetch_sumstats_locus
from graphgwas.finemapping_v2 import hbp_finemap_from_sumstats
# Fetch BMI sumstats near FTO (GRCh37)
sumstats = fetch_sumstats_locus(
    phenocode='21001', chr='16',
    start=53720000, end=53920000,
    trait_type='continuous', modifier='irnt',
    ancestries=['EUR', 'CSA', 'AFR', 'EAS'],
)
print({anc: len(s.variants) for anc, s in sumstats.items()})
"

# Full pipeline with Neo4j + multi-omics graph:
# (1) Start Neo4j with the pre-built human dump (17 GB, from Zenodo)
# (2) Run GAFM fine-mapping on a lead variant
graphgwas finemap --chr 16 --pos 53820527 --window 100000 \
    --phenotype BMI --method l1 -o credible_set.tsv
```

## The graph schema

```
 Variant ──HAS_CONSEQUENCE──> Gene ──IN_PATHWAY──> Pathway
    │                           │
    ├── (af, qual, gt_packed)   ├── INTERACTS_WITH (STRING PPI ≥ 700)
    ├── eQTL ─────────────> Gene (tissue-specific, GTEx v8)
    ├── IN_REGULATORY ─────> RegulatoryElement  (ENCODE cCRE)
    └── FOR_VARIANT <─── AssociationResult ──IN_STUDY──> GWASStudy
```

The credible-set output is itself a graph object: each reported variant is
co-queryable with its gene, tissue and pathway neighbours in a single Cypher
traversal, eliminating the post-hoc enrichment step that flat-prior pipelines
require.

## Three interfaces

| Interface | Use case | Entry point |
|---|---|---|
| **CLI** (52 commands, 15 groups) | interactive analysis, scripted pipelines | `graphgwas ...` |
| **REST API** (FastAPI, 37 endpoints) | web integration, programmatic access | `graphgwas api serve` |
| **MCP server** (FastMCP, 16 tools) | AI-agent access via any MCP-compatible client | `graphgwas mcp` |

Full documentation in [`docs/manual/`](docs/manual/index.md); end-to-end
walkthrough in [`vignettes/fine-mapping-quickstart.md`](vignettes/fine-mapping-quickstart.md).

## Fine-mapping methods at a glance

| Method | Complexity | Typical runtime / locus | Wins vs SuSiE at |
|---|---|---|---|
| **HBP** (three-layer factor graph + Banach contraction) | O(E × T) | 0.02–0.08 s | accuracy parity; 6–60× faster |
| **GAFM** (LD-deconvolved + adaptive α + graph prior) | O(n²) | 0.07 s | 27–2 at weak signal + tissue-specific eQTL priors |
| **GAFM-MX / HBP-MX** (v0.1.5: + λ_GC deflation + LD-deconvolved 4-component mixture BF) | O(n²) | 0.03 s | 10/21 (47.6%) top-1-PIP exact on rice 21-QTN panel |
| **ENS** (v0.1.5: mean-of-PIPs of GAFM-MX and HBP-MX) | O(n²) | 0.05 s | matches GAFM-MX/HBP-MX |
| **CLGF** (cross-locus EM) | O(L × T) | locus-dependent | multi-locus shared-pathway evidence |
| **GLEM** (graph-latent-embedding fine-mapping) | O(n² + n d) | 0.1 s | multi-signal detection |

## Documentation

- [`docs/INSTALL.md`](docs/INSTALL.md) — detailed installation guide
  (Neo4j, Python env, Hail for Pan-UKB LD, optional GNN deps)
- [`docs/INPUT_OUTPUT_GUIDE.md`](docs/INPUT_OUTPUT_GUIDE.md) — practical
  how-to for end users: prepare inputs, run fine-mapping, check outputs,
  interpret PIPs and credible sets, and diagnose problems
- [`docs/manual/index.md`](docs/manual/index.md) — full CLI reference
  (52 commands across 15 groups)
- [`vignettes/fine-mapping-quickstart.md`](vignettes/fine-mapping-quickstart.md) — 15-min Pan-UKB sumstats → credible set
- [`vignettes/full-1kg-pipeline.md`](vignettes/full-1kg-pipeline.md) — 4–6 h end-to-end: raw 1000 Genomes VCF → GWAS → fine-mapping → graph-queryable credible set
- [`docs/MATHEMATICAL_PROOFS.md`](docs/MATHEMATICAL_PROOFS.md) — theorems 1–5
- [`docs/REPRODUCIBILITY.md`](docs/REPRODUCIBILITY.md) — regenerate every
  paper figure and table from a single command

## Platform scope beyond fine-mapping

GraphGWAS is a platform of which fine-mapping is the first method class
rigorously benchmarked (see the accompanying Nature Genetics paper). The
codebase additionally implements:

- Epistasis (M1 LD-pruned, M2 motif-filtered, M3 differential-subgraph,
  M4 dark-matter pairs) — companion manuscript in preparation
- Heritability (6 estimators including spectral, GRM-REML, conductance)
- Multivariate cross-trait analysis (r_G, G-matrix, coherence, pleiotropy)
- Polygenic risk scores (classical + pathway-weighted)
- Mendelian randomisation (IVW, Egger, weighted median)
- Gene–environment interactions (multi-environment trials)
- Heterogeneous GNN (PyTorch Geometric) and LangGraph AI-agent interface

Honest benchmark-status table in Supplementary Note S3 of the manuscript.

## Data

Pre-built Neo4j graph databases on Zenodo (DOIs assigned on acceptance):

| Dataset | Size | Contents |
|---|---|---|
| Human 1KG + multi-omics | 17 GB | 70.7 M variants, 3,202 samples, 20,092 genes, 43.2 M GTEx eQTLs, 230 K STRING PPIs, 370 K ENCODE cCREs |
| Yeast 1011 Genomes | 0.5 GB | 1.92 M variants, 1,011 strains, SGD gene annotations, 35 growth-trait phenotypes |

Pan-UKB summary statistics are streamed on demand via tabix over HTTPS from
the public Amazon S3 bucket `pan-ukb-us-east-1`; no authentication or bulk
download required.

## Citation

If you use GraphGWAS, please cite the accompanying Nature Genetics
manuscript (*Relational biological structure improves fine-mapping of causal
GWAS variants under weak signal*, submitted 2026) and the Zenodo-versioned
software release. See [`CITATION.cff`](CITATION.cff).

## License

MIT — see [`LICENSE`](LICENSE).
