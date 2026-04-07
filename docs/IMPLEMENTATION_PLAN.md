# GraphGWAS — Implementation Plan & Roadmap

## Document Purpose
This is the master implementation plan for GraphGWAS — a graph-native GWAS platform
built on Neo4j. It is structured for systematic execution and covers all phases
from schema audit to AI agent integration.

**READ FIRST:** the project instructions, `docs/GWAS_SKILL.md`, `docs/GRAPH_SCHEMA_SKILL.md`

---

## Scientific Context & Motivation

### Why GraphGWAS Exists

Classical GWAS tools (PLINK, SAIGE, Hail) perform massively parallel *univariate*
regression. They answer: "does variant X associate with phenotype Y?" — 10 million
times, independently. This is structurally blind to:

1. **Epistasis** — variants whose effect depends on co-occurrence with other variants
2. **Multi-locus architecture** — diseases caused by combinations of rare variants
3. **Functional context** — whether a variant sits in a gene, pathway, regulatory element
4. **Imbalanced cohorts** — rare diseases with very few cases
5. **Relational structure** — kinship, haplotype sharing, population history

GraphGWAS addresses all five by treating association analysis as **graph topology
computation** rather than matrix regression. Complex trait architecture is a graph
property. GraphGWAS makes it measurable as such.

### Relationship to GraphPop
GraphGWAS runs on the same Neo4j instance as GraphPop, sharing the Variant/Sample/
Gene/Pathway schema. GraphPop-computed statistics (Tajima's D, F_ST, iHS, embeddings)
become covariates in GraphGWAS association tests — a novel integration unavailable
anywhere else.

---

## PHASE 0: Schema Audit & Environment Setup
**Duration: 1–3 sessions | No code yet — discovery only**
**Blocker: Everything else depends on this phase.**

### 0.1 Connect to Neo4j and Run Full Schema Audit

```bash
# Verify Neo4j is running
neo4j status  # or: docker ps | grep neo4j

# Run schema audit procedure (see GRAPH_SCHEMA_SKILL.md section 1)
```

**Deliverable:** `docs/CARRIES_STATUS.md` containing:
- [ ] List of all node labels present
- [ ] List of all relationship types present
- [ ] Properties on Variant nodes (with sample values)
- [ ] Properties on Sample nodes (with sample values)
- [ ] CARRIES edge count (or explicit confirmation of absence)
- [ ] Gene/Pathway node counts
- [ ] LD edge count
- [ ] Phenotype presence on Sample nodes (count)
- [ ] Decision on genotype representation strategy (Options 1–4)

### 0.2 Clarify the CARRIES Decision

Based on the audit, document in `docs/CARRIES_STATUS.md`:

**If CARRIES edges exist (even partially):**
→ Proceed to Phase 1 with full CARRIES-based implementation
→ Document scope: all samples? phenotyped only? all variants? candidate only?

**If CARRIES edges are absent:**
→ Document the reason if discoverable from git history, comments, or task logs
→ Select the recovery strategy (see CLAUDE.md "CARRIES Recovery Strategy")
→ The choice between Options 1–4 depends on:
   - Number of phenotyped samples currently in the database
   - Available disk space for new edges
   - Whether import pipeline can be re-run
→ **Recommend Option 1 (Scoped CARRIES) if ≤ 50K phenotyped samples**
→ **Recommend Option 2 (Carrier Arrays) if > 50K samples or import pipeline unavailable**

### 0.3 Environment Setup

- [ ] Verify Java 21+ and Maven available
- [ ] Verify Python 3.10+ with pip available
- [ ] Verify Neo4j version (`CALL dbms.components()`)
- [ ] Verify Neo4j GDS plugin installed (`CALL gds.version()`)
- [ ] Verify APOC plugin if available
- [ ] Set up project directory structure (see CLAUDE.md)
- [ ] Create `tasks/todo.md` and `tasks/lessons.md`

### 0.4 Test Data Preparation

- [ ] Obtain or simulate a small test VCF (500–5000 samples, chr22 or single gene region)
  with known phenotype–genotype relationships
- [ ] Create synthetic phenotype file: `tests/data/phenotypes_test.csv`
  with columns: `sample_id, case_control, BMI, age, sex, PC1, PC2, PC3`
- [ ] Suggested tool: `msprime` simulation with `slim` selection or
  manual spiking of known associations
- [ ] Verify test data loads cleanly into the existing schema

**Phase 0 Exit Criteria:**
- `docs/CARRIES_STATUS.md` is complete and reviewed
- Genotype representation strategy is decided and documented
- Test data is available and loadable
- All environment dependencies verified

---

## PHASE 1: Data Foundation — Genotype & Phenotype Layer
**Duration: 2–4 sessions**
**Depends on: Phase 0 complete, CARRIES_STATUS.md finalized**

### 1.1 Phenotype Import Pipeline

**File:** `src/python/graphgwas/importer.py`

Implement `graphgwas.phenotype.load()`:
- Accepts CSV with sample_id column + one or more phenotype columns
- Matches Sample nodes by ID
- Writes phenotype MAP to `Sample.phenotypes`
- Sets `Sample.is_case` and `Sample.is_control` derived from primary trait column
- Applies QC filters: missingness, sex check, outlier exclusion
- Reports: N samples matched, N cases, N controls, N unmatched

```python
# CLI interface
python -m graphgwas.importer \
  --neo4j-uri bolt://localhost:7688 \
  --phenotype-file tests/data/phenotypes_test.csv \
  --sample-id-col sample_id \
  --primary-trait case_control \
  --case-value case \
  --control-value control \
  --additional-traits BMI,age,sex
```

**Tests:**
- [ ] Load synthetic phenotype file, verify all samples updated
- [ ] Verify `is_case` / `is_control` correctly set
- [ ] Test with mismatched sample IDs (partial match) — should warn, not fail
- [ ] Verify MAP structure on Sample.phenotypes

### 1.2 Genotype Layer Implementation

**Depends on:** Phase 0 genotype strategy decision

#### Option 1: Scoped CARRIES Import (Preferred)

Extend the GraphPop VCF import pipeline to optionally emit CARRIES edges
only for phenotyped samples:

```python
# src/python/graphgwas/carries_import.py
def import_carries_scoped(vcf_path, neo4j_uri, 
                          phenotyped_sample_ids: set,
                          af_threshold: float = 0.05,
                          consequence_filter: list = ['HIGH', 'MODERATE']):
    """
    Import CARRIES edges only for:
    - Samples in phenotyped_sample_ids
    - Variants with af_total < af_threshold OR consequence in consequence_filter
    """
```

Emit bulk import CSV: `CARRIES_scoped.csv` with columns:
`:START_ID(Sample), :END_ID(Variant), gt:int, phase:int`

Run via `neo4j-admin database import --mode=incremental` or `LOAD CSV`.

**Scale estimate:** 5K samples × 500K candidate variants = 2.5B edges max.
Actual CARRIES per sample ~= 2 × (n_heterozygous_sites). At 5% AF threshold,
roughly 1–5% of variants per sample = 5K–50K CARRIES per sample.
For 5K samples: 25M–250M edges — feasible on Tier 1 hardware.

#### Option 2: Carrier Arrays on Variant Nodes (Fallback)

```cypher
// For each Variant, compute carrier arrays from GWAS cohort
MATCH (s:Sample {is_case: true})-[:CARRIES]->(v:Variant)
WITH v, collect(s.id) AS case_ids
SET v.case_carrier_ids = case_ids,
    v.case_carrier_count = size(case_ids)

MATCH (s:Sample {is_control: true})-[:CARRIES]->(v:Variant)
WITH v, collect(s.id) AS ctrl_ids
SET v.control_carrier_ids = ctrl_ids,
    v.control_carrier_count = size(ctrl_ids)
```

**Note:** If CARRIES don't exist at all, this must be populated from external VCF.
Implement `graphgwas.genotype.populate_carrier_arrays(vcf_path, sample_ids)`.

### 1.3 Functional Annotation Verification

If Gene/Pathway nodes are absent (check from Phase 0):

```python
# src/python/graphgwas/annotation.py
def import_vep_annotations(vep_output_file, neo4j_uri):
    """
    Parse VEP output and create:
    - Gene nodes with symbol, ensembl_id
    - HAS_CONSEQUENCE edges: Variant → Gene
    - IN_PATHWAY edges (from Reactome/KEGG): Gene → Pathway
    - HAS_GO_TERM edges: Gene → GOTerm
    """
```

Minimum viable annotation: HAS_CONSEQUENCE with impact and consequence type.
Pathway and GO annotation can be Phase 2 if absent.

### 1.4 Index Creation

Run all indexes from `docs/GRAPH_SCHEMA_SKILL.md` section 5.
Verify with `SHOW INDEXES` and time a test query before/after.

### 1.5 QC Procedure

Implement `graphgwas.qc.run(phenotype_key)`:
- Hardy-Weinberg filter (exclude HWE p < 1e-10 in controls)
- Call rate filter (exclude variants with call_rate < 0.95)
- MAC filter (flag variants with MAC < 5 in cases)
- Sample missingness (flag samples with genotype_missingness > 0.05)
- Relatedness check: use existing KINSHIP edges, flag pairs > 0.125 for exclusion

**Phase 1 Exit Criteria:**
- [ ] Phenotype data loaded, `is_case`/`is_control` set on all samples
- [ ] Genotype representation confirmed and documented
- [ ] Functional annotation (at minimum HAS_CONSEQUENCE) verified present
- [ ] All GWAS indexes created
- [ ] QC procedure runs without error
- [ ] Test data passes QC with expected outcome

---

## PHASE 2: Single-Locus Association Engine
**Duration: 2–3 sessions**
**Depends on: Phase 1 complete**

### 2.1 Core Association Procedure

**File:** `src/java/graphgwas/procedures/AssociationProcedures.java`

Implement `graphgwas.assoc.single_locus(chr, start, end, phenotype_key, method, covariates[])`:

**Algorithm:**
1. For each Variant in range, retrieve carrier counts (from CARRIES or carrier arrays)
2. Retrieve covariate values from Sample nodes (PCA coords, sex, age, etc.)
3. Build 2×2 table or design matrix
4. Apply chosen method:
   - `method = "firth"`: Firth penalized logistic regression (default for case-control)
   - `method = "linear"`: linear regression (for quantitative traits)
   - `method = "chi2"`: basic chi-squared (fast screening, not for rare variants)
5. Store result as `(:AssociationResult)` node connected to Variant
6. Return results ordered by p_value

**Firth Regression Implementation:**
Use Apache Commons Math or import a Java implementation of Firth's penalized
likelihood correction. The key formula:
```
L*_F(β) = L(β) + 0.5 × log|I(β)|
```
Where `I(β)` is the Fisher information matrix.

Reference implementation: `logistf` R package (port algorithm to Java).

### 2.2 Population Stratification Covariates

Before running GWAS, verify PCA coordinates exist on Sample nodes:
```cypher
MATCH (s:Sample) WHERE s.pca_coords IS NOT NULL RETURN count(s)
```

If absent: run `graphpop.structure.pca()` (from GraphPop layer) to compute
and store PC coordinates. Use first 10 PCs as default covariates.

### 2.3 Genome-Wide Scan Procedure

Implement `graphgwas.assoc.genome_scan(chr, phenotype_key, window_size, method)`:
- Calls `single_locus` in windows across chromosome
- Materializes results as AssociationResult nodes
- Returns summary statistics + significant hit list
- Stores `(:GWASStudy)` node with run metadata

### 2.4 Results Query Interface

```cypher
// Manhattan plot data query
MATCH (ar:AssociationResult {run_id: $run_id})-[:FOR_VARIANT]->(v:Variant)
RETURN v.chr, v.pos, ar.p_value_log10, ar.beta, v.gene_symbol
ORDER BY v.chr, v.pos

// QQ-plot data
MATCH (ar:AssociationResult {run_id: $run_id})
RETURN ar.p_value ORDER BY ar.p_value

// Significant hits with functional context
MATCH (ar:AssociationResult {run_id: $run_id})-[:FOR_VARIANT]->(v:Variant)
WHERE ar.p_value < 5e-8
OPTIONAL MATCH (v)-[:HAS_CONSEQUENCE]->(g:Gene)-[:IN_PATHWAY]->(p:Pathway)
RETURN ar.variant_id, ar.p_value, ar.beta, v.consequence, g.symbol, p.name
ORDER BY ar.p_value
```

**Validation:**
- [ ] Run single-locus GWAS on test data with known associations
- [ ] Compare p-values to PLINK2 output on same data (expect Pearson r > 0.99)
- [ ] Verify Firth correction reduces inflation for rare variants (compare lambda GC)
- [ ] Verify AssociationResult nodes created and queryable

**Phase 2 Exit Criteria:**
- [ ] Single-locus GWAS runs on test data
- [ ] Results validated against PLINK2
- [ ] AssociationResult nodes stored and queryable
- [ ] GWASStudy audit node created per run
- [ ] Manhattan and QQ-plot data extractable in single Cypher query

---

## PHASE 3: Multi-Locus & Rare Disease Engine
**Duration: 3–5 sessions**
**Depends on: Phase 2 complete, CARRIES or carrier arrays confirmed available**
**This phase delivers the primary scientific novelty of GraphGWAS.**

### 3.1 Variant Co-Occurrence Graph

Implement `graphgwas.epistasis.cooccurrence_graph(chr, start, end, phenotype_key, min_cocarriers)`:

1. For each pair of variants in region carried by ≥ min_cocarriers case samples:
   - Compute case co-carriage count
   - Compute control co-carriage count
   - Compute enrichment = case_count / (control_count_normalized)
   - Write `CO_OCCURS` edge if enrichment > threshold
2. Returns: edge count written, top enriched pairs

**Performance note:** This is O(V²) per sample in the worst case. Bound it:
- Limit to small genomic windows (≤ 1 Mb)
- Minimum co-carrier threshold of 3 eliminates most pairs
- Use set intersection on carrier_ids arrays for fast filtering

### 3.2 Epistasis Enrichment Test

Implement `graphgwas.epistasis.enriched_pairs(phenotype_key, fold_enrichment, p_threshold)`:

1. Query all CO_OCCURS edges above enrichment threshold
2. For each pair, compute statistical significance:
   - Fisher's exact test: 2×2 table of case/control co-carriage vs independent carriage
   - Or: permutation test (shuffle case/control labels N=1000 times, empirical null)
3. Apply multiple testing correction (Bonferroni or Benjamini-Hochberg)
4. Store significant pairs as `(:AssociationResult)` nodes with `method = "cooccurrence"`
   and `partner_variant_id` set

### 3.3 Gene Burden Test

Implement `graphgwas.rare.burden_test(gene_id, phenotype_key, af_threshold, consequence_types[])`:

1. For each sample: count rare, qualifying variants in gene
2. Compare distribution between cases and controls
3. Methods: sum test (BURDEN), SKAT (optional, complex), SKAT-O (combination)
4. For BURDEN: Mann-Whitney U or logistic regression on burden count
5. Store result as AssociationResult with `method = "burden_gene"`

**Validation:** Compare burden test p-values to SKAT-O on same data.

### 3.4 Disease Subgraph Extraction (Rare Disease)

Implement `graphgwas.rare.disease_subgraph(phenotype_key, af_threshold, min_cases, min_enrichment)`:

**Algorithm:**
1. Score each Variant: enrichment = case_carrier_count / (control_carrier_count_normalized)
   Weight by 1/af_total (upweight ultra-rare variants)
2. Propagate variant scores to Gene nodes:
   gene_score = sum(variant_scores for variants with HAS_CONSEQUENCE to gene)
3. Propagate gene scores to Pathway nodes:
   pathway_score = sum(gene_scores) / sqrt(pathway_gene_count)  (size-normalized)
4. Extract connected subgraph: all variants, genes, pathways above minimum score
5. Statistical validation: permute case/control labels 1000 times, compute empirical p-value
6. Store as `(:DiseaseSubgraph)` node with IN_DISEASE_SUBGRAPH edges

**This is the key rare disease contribution** — a statistically validated network of
variants, genes, and pathways constituting the genetic architecture of the disease.

### 3.5 Phenotype Label Propagation

Implement `graphgwas.diffusion.phenotype_propagate(phenotype_key, n_iterations, restart_prob)`:

Use Neo4j GDS label propagation or implement custom diffusion:

1. Initialize: case Sample nodes score = 1.0, control Sample nodes score = 0.0
2. Propagate through CARRIES edges → Variant nodes (weight by 1/af to upweight rare)
3. Propagate through HAS_CONSEQUENCE → Gene nodes
4. Propagate through IN_PATHWAY → Pathway nodes
5. With restart probability α: at each step, fraction α returns to initial Sample scores
6. Iterate until convergence (delta < 0.001) or n_iterations
7. Store diffusion_score on each node
8. Permutation test: shuffle Sample labels 1000×, compute null distribution of scores
9. Z-score normalize against permuted null

**Pathway-level diffusion score** = key output: which pathways are enriched for
phenotype-associated rare variants even when no individual variant is significant.

### 3.6 Kinship-Aware Cosegregation (Rare Disease)

Implement `graphgwas.rare.cosegregation_scan(phenotype_key, min_kinship, af_threshold)`:

1. Find pairs of related cases (KINSHIP coefficient > min_kinship)
2. For each pair: find variants co-carried by both
3. Count: how many related case pairs share each variant?
4. Statistical test: compare to expected sharing given kinship coefficient
5. Sort by cosegregation score

**Phase 3 Exit Criteria:**
- [ ] Co-occurrence graph builds in < 10 min on test data
- [ ] Epistasis pairs detected in simulated epistatic data (validation required)
- [ ] Burden test results match expected direction on test data
- [ ] Disease subgraph extracted and statistically validated for rare disease test case
- [ ] Phenotype diffusion runs and pathway scores computed
- [ ] Permutation-based p-values computed for all methods
- [ ] All results stored as queryable AssociationResult or DiseaseSubgraph nodes

---

## PHASE 4: GNN Association Engine
**Duration: 3–4 sessions**
**Depends on: Phase 3 complete, PyTorch Geometric available**

### 4.1 Graph Export to PyG Format

Implement `graphgwas.gnn.export_graph(phenotype_key, chr, start, end, output_path)`:

Export heterogeneous graph from Neo4j to PyTorch Geometric HeteroData:

```python
# Node types and features
hetero_data['Variant'].x = [af_total, cadd_score, consequence_encoded, ld_degree, ...]
hetero_data['Sample'].x = [pca_coords, ancestry_proportions, sex_encoded, ...]
hetero_data['Gene'].x = [gene_length, n_variants, pathway_count, ...]

# Edge types
hetero_data['Sample', 'CARRIES', 'Variant'].edge_index = ...  # from CARRIES edges
hetero_data['Variant', 'LD', 'Variant'].edge_index = ...      # from LD edges
hetero_data['Variant', 'IN_GENE', 'Gene'].edge_index = ...    # from HAS_CONSEQUENCE

# Labels (for training)
hetero_data['Sample'].y = is_case (binary) or phenotype_value (continuous)
```

### 4.2 Heterogeneous GNN Architecture

**File:** `src/python/graphgwas/gnn.py`

Architecture: **HGT (Heterogeneous Graph Transformer)** or **R-GCN**

```python
class GraphGWASModel(nn.Module):
    def __init__(self, metadata, hidden_channels=64, num_layers=3):
        super().__init__()
        # Encoder: HGT with 3 message-passing layers
        self.encoder = HGT(hidden_channels, num_heads=4, num_layers=num_layers, metadata=metadata)
        # Phenotype predictor: MLP on Sample node embeddings
        self.predictor = MLP([hidden_channels, 32, 1])
    
    def forward(self, x_dict, edge_index_dict):
        h = self.encoder(x_dict, edge_index_dict)
        return self.predictor(h['Sample'])
```

**Training:**
- Loss: binary cross-entropy (case-control) or MSE (quantitative)
- Optimizer: Adam, lr=0.001
- Regularization: dropout 0.3, weight decay 1e-4
- Validation: 5-fold cross-validation on samples
- Imbalance handling: class weighting (n_control/n_case) in loss function

### 4.3 Feature Attribution (GNN Explainability)

After training, use GNN explainability to identify important variants:
- **GNNExplainer** (PyG built-in): per-node subgraph importance
- **SHAP values** on GNN output: variant-level attribution
- Import top-attributed variants back to Neo4j as AssociationResult nodes

```python
from torch_geometric.explain import GNNExplainer
explainer = GNNExplainer(model, epochs=200)
node_feat_mask, edge_mask = explainer.explain_node(node_idx, x_dict, edge_index_dict)
```

### 4.4 Import Embeddings Back to Graph

Store GNN-learned Sample embeddings back to Neo4j:
```cypher
MATCH (s:Sample {id: $sample_id})
SET s.gnn_embedding = $embedding,
    s.gnn_phenotype_score = $score
```

These embeddings capture multi-locus genetic effects on phenotype — a richer
representation than PCA, usable as population structure covariate in Phase 2 analyses.

**Phase 4 Exit Criteria:**
- [ ] Graph exports to PyG HeteroData without errors
- [ ] GNN trains without NaN loss
- [ ] AUROC > 0.6 on held-out test samples (better than chance)
- [ ] Feature attribution identifies known causal variants in simulated data
- [ ] GNN embeddings imported back to Neo4j and queryable

---

## PHASE 5: AI Agent Layer
**Duration: 2–3 sessions**
**Depends on: Phases 2–4 complete**

### 5.1 GraphRAG over Association Results

Set up vector index on Variant and AssociationResult embeddings:
```cypher
CREATE VECTOR INDEX assoc_embedding IF NOT EXISTS
FOR (ar:AssociationResult) ON (ar.embedding)
OPTIONS {indexConfig: {`vector.dimensions`: 256, `vector.similarity_function`: 'cosine'}}
```

Implement `graphgwas.graphrag.query(question)` combining:
- Vector similarity search on AssociationResult/Variant embeddings
- Cypher-based structured graph retrieval for relational context
- LLM synthesis of retrieved results

### 5.2 LangGraph GWAS Agent

**File:** `src/python/graphgwas/agent.py`

Implement a ReAct-style agent with tools:
- `run_single_locus_gwas(phenotype, region)` → calls Phase 2 procedure
- `find_epistasis(phenotype, region)` → calls Phase 3 procedure
- `extract_disease_subgraph(phenotype)` → calls Phase 3 procedure
- `query_graph(cypher)` → direct Cypher execution
- `interpret_result(result_dict)` → LLM interpretation of graph query result

**Example agent workflow:**
```
User: "What genetic variants are associated with this rare kidney disease?
       My cohort has 45 cases and 2000 controls."

Agent:
1. [run_single_locus_gwas] → p-values for all variants (likely nothing significant at 5e-8)
2. [extract_disease_subgraph] → pathway-level enrichment despite marginal insignificance
3. [query_graph] → "find variants in top pathways shared by related cases"
4. [interpret_result] → "Three variants in the COL4A3/COL4A4 pathway show consistent
   enrichment in cases (OR=8.2, permutation p=0.003). Two are shared between
   related cases (3rd-degree relatives), suggesting cosegregation with disease."
```

**Phase 5 Exit Criteria:**
- [ ] GraphRAG returns relevant results for natural language GWAS queries
- [ ] LangGraph agent successfully executes multi-step GWAS workflows
- [ ] Agent correctly handles "no significant single-locus hits" case
  by falling back to rare disease subgraph analysis
- [ ] Agent outputs include effect sizes, p-values, and biological interpretation

---

## Summary: Capability by Phase

| Capability | Phase | Novelty |
|---|---|---|
| Phenotype loading into graph | 1 | Foundational |
| CARRIES / carrier arrays confirmed | 1 | Foundational |
| Firth-corrected single-locus GWAS | 2 | Parity with PLINK2 |
| Results as queryable graph nodes | 2 | Graph-native advantage |
| Functional conditioning in single query | 2 | Graph-native advantage |
| Variant co-occurrence / epistasis screen | 3 | **Novel** |
| Rare disease subgraph extraction | 3 | **Novel** |
| Phenotype label diffusion | 3 | **Novel** |
| Kinship-aware cosegregation | 3 | **Novel** |
| Heterogeneous GNN phenotype prediction | 4 | **Novel** |
| GNN feature attribution for multi-locus effects | 4 | **Novel** |
| Natural language GWAS interface | 5 | **Novel** |
| GraphPop stats as GWAS covariates | 2–5 | **Novel** (requires GraphPop) |

---

## Unresolved Questions — Must Answer in Phase 0

Record answers in `docs/CARRIES_STATUS.md`:

1. **What is the current genotype representation in the GraphPop database?**
   (CARRIES edges? Dosage arrays? External file? Nothing at all?)

2. **Were Sample nodes implemented? Do they have phenotype properties?**

3. **Was the functional annotation layer (Gene/Pathway) implemented?**

4. **Were LD edges computed and stored?**

5. **What was the actual reason CARRIES was dropped?**
   (Scale/performance? Import pipeline issues? Deliberate architectural decision?)

6. **Is there an existing GraphPop database we are extending, or starting fresh?**

7. **What is the target cohort? (Rice? Human? Crop population?)**
   This affects: phenotype type (binary disease vs quantitative trait),
   expected case-control ratio, population structure complexity.

8. **What genome build and annotation version?** (GRCh38 + Ensembl 110?)

9. **Is VEP output available for functional annotation, or does it need to be run?**

10. **What is the hardware tier?** (Tier 1: laptop/workstation, Tier 2: cluster)
    This constrains CARRIES feasibility at scale.

---

## Publication Target

Successful completion of Phases 1–3 supports a methods paper:

**"GraphGWAS: Graph-Native Multi-Locus Association Analysis Enables Detection
of Epistatic Architecture and Rare Disease Subgraphs Invisible to Classical GWAS"**

Target: *Nature Methods*, *Genome Research*, or *Bioinformatics*

**Key results to demonstrate:**
1. Single-locus results match PLINK2 (validates implementation correctness)
2. Epistasis detection: simulated epistatic pairs detected with FDR < 0.05 when
   marginal effects are non-significant (demonstrates novel power)
3. Rare disease subgraph: pathway enrichment p < 0.01 with 50 cases when no
   single variant reaches p < 5e-8 (demonstrates imbalanced-data advantage)
4. Phenotype diffusion: pathway-level signal detectable with 20 cases
5. GNN AUROC improvement over single-locus logistic regression on simulated data
6. Benchmark: comparative analysis vs PLINK2, SAIGE, SKAT-O on identical data
