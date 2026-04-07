# GraphGWAS Statistical Engine Design

## Overview

GraphGWAS implements a layered statistical engine spanning 5 phases,
from classical GWAS (PLINK-equivalent) to genuinely novel graph-native methods.

---

## Layer 2: Single-Locus Association (Phase 1-2)

| Method | Module | Function | When Used |
|--------|--------|----------|-----------|
| Chi-squared allelic | assoc.py | `chi2_allelic_test()` | MAC ≥ 20, all cells ≥ 5 |
| Fisher's exact | assoc.py | `fisher_exact_test()` | MAC < 5 |
| Firth penalized logistic | assoc.py | `firth_logistic_regression()` | 5 ≤ MAC < 20 |
| Standard logistic | assoc.py | `logistic_regression()` | With covariates, common variants |
| Linear regression | assoc.py | `linear_regression()` | Quantitative traits |
| Score test (fallback) | assoc.py | `_score_test_logistic()` | When MLE fails to converge |

### Auto Method Selection
```
MAC ≥ 20 and all cells ≥ 5 → chi-squared (fast screening)
5 ≤ MAC < 20               → Firth penalized logistic
MAC < 5                     → Fisher's exact
```

---

## Layer 2.5: Gene-Level Association (Phase 2)

| Method | Module | Function | Test Statistic |
|--------|--------|----------|----------------|
| MPAT directed | mpat.py | `mpat_gene_test(method='directed')` | T = w^T z / √(w^T Σ w) ~ N(0,1) |
| MPAT undirected | mpat.py | `mpat_gene_test(method='undirected')` | Q = z^T z ~ weighted χ² (Liu et al.) |

---

## Layer 3: Graph-Native Methods (Phase 3)

### Epistasis via Co-Occurrence Networks
| Step | Function | Algorithm |
|------|----------|-----------|
| Build graph | `build_cooccurrence_graph()` | Pairwise carrier-set AND, enrichment scoring |
| Detect modules | `detect_epistatic_modules()` | Louvain community detection (networkx) |
| Test modules | `test_epistatic_module()` | Permutation test (1000 label shuffles) |
| Genome scan | `epistasis_scan()` | Sliding window across chromosomes |

### Max-Flow / Min-Cut Disease Architecture
| Step | Function | Algorithm |
|------|----------|-----------|
| Build network | `build_flow_network()` | Directed: Source→Sample→Variant→Gene→Pathway |
| Compute flows | `compute_pathway_flows()` | Ford-Fulkerson max-flow (networkx) |
| Test pathways | `test_pathway_flow()` | Permutation of case labels |
| Extract architecture | `extract_disease_architecture()` | Union of significant min-cuts |

### LD Fine-Mapping via Centrality
| Step | Function | Algorithm |
|------|----------|-----------|
| Build LD graph | `build_ld_graph()` | Pairwise r² > threshold |
| Fine-map | `centrality_finemapping()` | Betweenness / PageRank / degree centrality |

### Graph Spectral Phenotype Decomposition
| Step | Function | Algorithm |
|------|----------|-----------|
| Build similarity | `build_sample_similarity_graph()` | Rare variant co-carriage (1/AF weighted) |
| Eigendecompose | `spectral_decompose()` | Graph Laplacian eigsh (scipy sparse) |
| Filter phenotype | `filter_phenotype_signal()` | Low-frequency component retention |

---

## Layer 4: GNN Association (Phase 4)

| Component | Function | Architecture |
|-----------|----------|-------------|
| Graph export | `export_to_pyg()` | Neo4j → PyG HeteroData |
| Model | `GraphGWASModel` | HeteroConv (SAGEConv per edge type), 3 layers |
| Training | `train_gnn()` | BCEWithLogitsLoss, class-weighted, Adam |
| Attribution | `explain_variants()` | Gradient-based feature importance |
| Import | `import_embeddings()` | GNN scores → Sample.gnn_phenotype_score |

---

## Layer 5: AI Agent (Phase 5)

| Component | Function | Technology |
|-----------|----------|-----------|
| Agent creation | `create_agent()` | LangGraph ReAct + Claude |
| Query interface | `run_query()` | Natural language → tool calls → response |
| Offline interpretation | `interpret_results()` | Structured analysis without LLM |
| Tools | 7 tools | scan, epistasis, flow, finemap, mpat, cypher, status |

---

## Complete Module Inventory

```
graphgwas/
├── __init__.py          Phase 0  (version)
├── config.py            Phase 0  (constants, thresholds)
├── db.py                Phase 0  (Neo4j connection)
├── schema.py            Phase 1  (indexes, schema audit)
├── phenotype.py         Phase 1  (phenotype import/activation)
├── genotype.py          Phase 1-2 (gt_packed decode, 7 output formats, PCA import, covariates)
├── qc.py                Phase 1  (HWE, MAC, missingness, het_rate)
├── assoc.py             Phase 1-2 (chi2, fisher, logistic, firth, linear, z-scores)
├── results.py           Phase 1  (AssociationResult nodes, queries)
├── mpat.py              Phase 2  (gene-level MPAT, LD matrix, Liu et al.)
├── plots.py             Phase 2  (Manhattan, QQ, TSV export)
├── epistasis.py         Phase 3  (co-occurrence network, Louvain, permutation)
├── flow.py              Phase 3  (max-flow/min-cut, disease architecture)
├── finemapping.py       Phase 3  (LD graph, centrality credible set)
├── spectral.py          Phase 3  (rare-variant similarity, Laplacian, filtering)
├── gnn.py               Phase 4  (HeteroGNN, export/train/explain/import)
├── agent.py             Phase 5  (LangGraph ReAct agent, 7 tools)
└── cli.py               All      (Click CLI, 25+ commands)
```
