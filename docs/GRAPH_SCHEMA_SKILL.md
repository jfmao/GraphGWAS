# GRAPH_SCHEMA_SKILL.md — Neo4j Schema Patterns for GraphGWAS

## Purpose
This file documents the target graph schema for GraphGWAS, the Cypher patterns
for all key operations, and the relationship to the existing GraphPop schema.
Read this before writing any Cypher or schema migration code.

---

## 1. Schema Discovery Protocol (Always Run First)

```cypher
// Step 1: What node labels exist?
CALL db.labels() YIELD label RETURN label ORDER BY label

// Step 2: What relationship types exist?
CALL db.relationshipTypes() YIELD relationshipType RETURN relationshipType

// Step 3: Sample Variant node properties
MATCH (v:Variant) RETURN keys(v) LIMIT 1

// Step 4: Sample Sample node properties
MATCH (s:Sample) RETURN keys(s) LIMIT 1

// Step 5: CRITICAL — CARRIES status
MATCH ()-[r:CARRIES]->() RETURN count(r) AS carries_count LIMIT 1

// Step 6: Functional annotation presence
MATCH (g:Gene) RETURN count(g) AS gene_count
MATCH (p:Pathway) RETURN count(p) AS pathway_count

// Step 7: LD edges
MATCH ()-[r:LD]->() RETURN count(r) AS ld_count LIMIT 1

// Step 8: Phenotype presence on Samples
MATCH (s:Sample) WHERE s.phenotypes IS NOT NULL RETURN count(s) AS phenotyped_samples
```

**Always write the results to `docs/CARRIES_STATUS.md` before proceeding.**

---

## 2. Target Node Schema

### 2.1 Variant Node (from GraphPop, extended for GWAS)

```cypher
(:Variant {
    // Identity (GraphPop)
    id: STRING,              // "chr1:12345:A:T"
    chr: STRING,
    pos: LONG,
    ref: STRING,
    alt: STRING,
    variant_type: STRING,    // "SNP", "INDEL", "SV"

    // Population allele counts (GraphPop fast-path)
    pop_ids: [STRING],
    ac: [INT],
    an: [INT],
    af: [FLOAT],
    af_total: FLOAT,
    ac_total: INT,
    call_rate: FLOAT,

    // Functional annotation (GraphPop, from VEP)
    consequence: STRING,
    impact: STRING,          // "HIGH", "MODERATE", "LOW", "MODIFIER"
    gene_symbol: STRING,
    cadd_score: FLOAT,

    // GraphGWAS ADDITIONS:
    // Option 2 fallback (if CARRIES absent):
    case_carrier_ids: [STRING],     // sample IDs that are cases and carry this variant
    control_carrier_ids: [STRING],  // sample IDs that are controls and carry this variant
    case_carrier_count: INT,
    control_carrier_count: INT,

    // Pre-computed GWAS statistics (stored after analysis runs)
    gwas_p_value: FLOAT,           // best single-locus p-value across all phenotypes
    gwas_beta: FLOAT,
    gwas_phenotype: STRING,
    gwas_run_id: STRING,

    // Diffusion scores (stored after propagation)
    diffusion_score: FLOAT,        // phenotype diffusion score
    diffusion_run_id: STRING,
})
```

### 2.2 Sample Node (from GraphPop, extended for GWAS)

```cypher
(:Sample {
    // Identity (GraphPop)
    id: STRING,
    population: STRING,
    sex: STRING,
    ancestry_proportions: [FLOAT],
    pca_coords: [FLOAT],
    embedding: [FLOAT],

    // GraphGWAS ADDITIONS:
    phenotypes: MAP,         // {trait_name: value} — already in GraphPop design
                             // e.g., {T2D: 'case', BMI: 24.5, age: 45, hypertension: 'yes'}
    phenotype_source: STRING, // which file/cohort
    is_case: BOOLEAN,        // derived from phenotypes for primary trait
    is_control: BOOLEAN,

    // Quality flags
    genotype_missingness: FLOAT,
    is_excluded_qc: BOOLEAN,
    exclusion_reason: STRING,
})
```

### 2.3 AssociationResult Node (GraphGWAS new)

```cypher
(:AssociationResult {
    id: STRING,              // "gwas_{run_id}_{variant_id}"
    run_id: STRING,
    variant_id: STRING,
    phenotype_key: STRING,
    method: STRING,          // "single_locus_firth", "burden_skat", "cooccurrence",
                             // "diffusion_pathway", "gnn_message_passing"

    // Statistical results
    beta: FLOAT,
    se: FLOAT,
    p_value: FLOAT,
    p_value_log10: FLOAT,    // -log10(p), for Manhattan plot
    odds_ratio: FLOAT,
    ci_lower: FLOAT,
    ci_upper: FLOAT,

    // Context
    n_cases: INT,
    n_controls: INT,
    maf: FLOAT,
    mac: INT,
    population: STRING,
    covariates: [STRING],    // e.g., ["PC1","PC2","PC3","sex","age"]
    af_threshold: FLOAT,     // for burden/rare variant tests

    // Multi-locus specific
    partner_variant_id: STRING,   // for epistasis pairs
    interaction_p_value: FLOAT,   // for epistasis
    joint_effect: FLOAT,          // joint beta minus sum of marginals

    timestamp: DATETIME,
    software_version: STRING,
})
```

### 2.4 GWASStudy Node (GraphGWAS new)

```cypher
(:GWASStudy {
    id: STRING,              // "study_{run_id}"
    run_id: STRING,
    name: STRING,
    phenotype_key: STRING,
    n_variants_tested: INT,
    n_cases: INT,
    n_controls: INT,
    method: STRING,
    genome_build: STRING,    // "GRCh38", "GRCh37"
    significance_threshold: FLOAT,
    n_significant: INT,
    timestamp: DATETIME,
    notes: STRING,
})
```

### 2.5 DiseaseSubgraph Node (GraphGWAS new — rare disease)

```cypher
(:DiseaseSubgraph {
    id: STRING,
    phenotype_key: STRING,
    n_variants: INT,
    n_genes: INT,
    n_pathways: INT,
    enrichment_score: FLOAT,
    permutation_p_value: FLOAT,
    run_id: STRING,
    timestamp: DATETIME,
})
```

---

## 3. Target Relationship Schema

### 3.1 CARRIES (from GraphPop — the critical individual-level genotype edge)

```cypher
(:Sample)-[:CARRIES {
    gt: INT,         // 1 = heterozygous, 2 = homozygous alt
    phase: INT,      // 0 = unphased, 1 = haplotype 1, 2 = haplotype 2
    gq: INT,         // genotype quality (optional)
    dp: INT,         // read depth (optional)
}]->(:Variant)

// HOM_REF is implicit (absence of CARRIES edge = reference homozygote)
// This is the key design decision: sparse representation
```

**CARRIES enables all individual-level GWAS queries.**
If absent, see CARRIES_STATUS.md for the alternative representation.

### 3.2 FOR_VARIANT (connects AssociationResult to Variant)

```cypher
(:AssociationResult)-[:FOR_VARIANT]->(:Variant)
(:AssociationResult)-[:IN_STUDY]->(:GWASStudy)
(:AssociationResult)-[:EPISTATIC_WITH]->(:AssociationResult)  // for interaction pairs
```

### 3.3 HAS_CONSEQUENCE (from GraphPop)

```cypher
(:Variant)-[:HAS_CONSEQUENCE {
    type: STRING,    // "missense_variant", "stop_gained", "regulatory_region_variant"
    impact: STRING,  // "HIGH", "MODERATE", "LOW", "MODIFIER"
    hgvs: STRING,
}]->(:Gene)
```

### 3.4 IN_PATHWAY / HAS_GO_TERM (from GraphPop)

```cypher
(:Gene)-[:IN_PATHWAY]->(:Pathway)
(:Gene)-[:HAS_GO_TERM]->(:GOTerm)
(:Pathway)-[:PART_OF]->(:Pathway)    // pathway hierarchy
(:GOTerm)-[:IS_A]->(:GOTerm)         // GO hierarchy
```

### 3.5 IN_DISEASE_SUBGRAPH (new)

```cypher
(:Variant)-[:IN_DISEASE_SUBGRAPH {enrichment_score: FLOAT}]->(:DiseaseSubgraph)
(:Gene)-[:IN_DISEASE_SUBGRAPH {aggregated_score: FLOAT}]->(:DiseaseSubgraph)
(:Pathway)-[:IN_DISEASE_SUBGRAPH {pathway_score: FLOAT}]->(:DiseaseSubgraph)
```

### 3.6 CO_OCCURS (new — variant epistasis graph)

```cypher
(:Variant)-[:CO_OCCURS {
    count: INT,              // number of samples co-carrying both variants
    expected_count: FLOAT,   // expected under independence
    enrichment: FLOAT,       // count / expected_count
    case_count: INT,         // co-carriage count in cases only
    control_count: INT,      // co-carriage count in controls only
    case_enrichment: FLOAT,  // case_count / control_count (normalized)
    population: STRING,
}]->(:Variant)
```

---

## 4. Key Cypher Query Patterns

### 4.1 Load Phenotype into Sample Nodes

```cypher
// Bulk update from CSV (run from python importer)
LOAD CSV WITH HEADERS FROM 'file:///phenotypes.csv' AS row
MATCH (s:Sample {id: row.sample_id})
SET s.phenotypes = {
    primary_trait: row.primary_trait,
    BMI: toFloat(row.BMI),
    age: toInteger(row.age),
    sex: row.sex
},
s.is_case = (row.primary_trait = 'case'),
s.is_control = (row.primary_trait = 'control')
```

### 4.2 Single-Locus Association Query (via CARRIES)

```cypher
// Count cases and controls carrying a specific variant
MATCH (v:Variant {id: $variant_id})
OPTIONAL MATCH (s_case:Sample {is_case: true})-[:CARRIES]->(v)
OPTIONAL MATCH (s_ctrl:Sample {is_control: true})-[:CARRIES]->(v)
RETURN v.id,
       count(DISTINCT s_case) AS n_case_carriers,
       count(DISTINCT s_ctrl) AS n_ctrl_carriers,
       v.af_total
// Then pass to Java procedure for Firth logistic regression
```

### 4.3 Gene Burden Query

```cypher
// Count rare, high-impact variant burden per sample
MATCH (s:Sample)-[:CARRIES]->(v:Variant)-[:HAS_CONSEQUENCE {impact: 'HIGH'}]->(g:Gene {symbol: $gene})
WHERE v.af_total < $af_threshold
WITH s, count(v) AS burden
RETURN s.id, s.is_case, burden
ORDER BY burden DESC
```

### 4.4 Pathway Burden Query

```cypher
// Burden across an entire pathway
MATCH (s:Sample)-[:CARRIES]->(v:Variant)-[:HAS_CONSEQUENCE]->(g:Gene)-[:IN_PATHWAY]->(p:Pathway {name: $pathway})
WHERE v.af_total < $af_threshold
  AND v.impact IN ['HIGH', 'MODERATE']
WITH s, p, count(DISTINCT v) AS pathway_burden
RETURN s.id, s.is_case, p.name, pathway_burden
```

### 4.5 Epistasis Co-Occurrence Query

```cypher
// Find variant pairs enriched in cases vs controls
MATCH (s:Sample {is_case: true})-[:CARRIES]->(v1:Variant),
      (s)-[:CARRIES]->(v2:Variant)
WHERE v1.id < v2.id
  AND v1.af_total < 0.05 AND v2.af_total < 0.05
WITH v1, v2, count(DISTINCT s) AS case_co_count

MATCH (s2:Sample {is_control: true})-[:CARRIES]->(v1),
      (s2)-[:CARRIES]->(v2)
WITH v1, v2, case_co_count, count(DISTINCT s2) AS ctrl_co_count

WHERE case_co_count >= 3
  AND (case_co_count * 1.0) / (ctrl_co_count + 0.5) > $enrichment_threshold
RETURN v1.id, v2.id, case_co_count, ctrl_co_count,
       (case_co_count * 1.0) / (ctrl_co_count + 0.5) AS enrichment
ORDER BY enrichment DESC
LIMIT 1000
```

### 4.6 Kinship-Aware Cosegregation (Rare Disease)

```cypher
// Find rare variants shared by related cases
MATCH (s1:Sample {is_case: true})-[k:KINSHIP]-(s2:Sample {is_case: true})
WHERE k.coefficient > 0.125
MATCH (s1)-[:CARRIES]->(v:Variant)<-[:CARRIES]-(s2)
WHERE v.af_total < 0.01
  AND v.impact IN ['HIGH', 'MODERATE']
WITH v, collect(DISTINCT s1.id) + collect(DISTINCT s2.id) AS cosegregating_cases
RETURN v.id, v.gene_symbol, v.consequence, size(cosegregating_cases) AS n_cosegregating
ORDER BY n_cosegregating DESC
```

### 4.7 Query Association Results as Graph Objects

```cypher
// Find significant hits in drug-targetable pathways
MATCH (ar:AssociationResult)-[:FOR_VARIANT]->(v:Variant)-[:HAS_CONSEQUENCE]->(g:Gene)-[:IN_PATHWAY]->(p:Pathway)
WHERE ar.p_value < 5e-8
  AND ar.phenotype_key = $phenotype
  AND p.has_drug_target = true
RETURN ar.variant_id, ar.p_value, ar.beta, g.symbol, p.name
ORDER BY ar.p_value ASC
```

### 4.8 Combined Population Genetics + GWAS Query (GraphPop integration)

```cypher
// Find GWAS hits that overlap regions under selection (from GraphPop)
MATCH (ar:AssociationResult {phenotype_key: $phenotype})-[:FOR_VARIANT]->(v:Variant)
WHERE ar.p_value < 5e-8
MATCH (v)-[:IN_WINDOW]->(w:GenomicWindow)
WHERE w.ihs_score > 2.0 OR w.tajima_d < -2.0
RETURN v.id, ar.p_value, ar.beta, w.ihs_score, w.tajima_d, v.gene_symbol
ORDER BY ar.p_value ASC
```

---

## 5. Index Strategy for GWAS

Add these indexes before running any GWAS analysis:

```cypher
// Phenotype index (critical for case/control separation)
CREATE INDEX sample_is_case IF NOT EXISTS FOR (s:Sample) ON (s.is_case)
CREATE INDEX sample_is_control IF NOT EXISTS FOR (s:Sample) ON (s.is_control)
CREATE INDEX sample_population IF NOT EXISTS FOR (s:Sample) ON (s.population)

// Variant frequency index (for rare variant filtering)
CREATE INDEX variant_af_total IF NOT EXISTS FOR (v:Variant) ON (v.af_total)
CREATE INDEX variant_impact IF NOT EXISTS FOR (v:Variant) ON (v.impact)
CREATE INDEX variant_gene IF NOT EXISTS FOR (v:Variant) ON (v.gene_symbol)

// Association result indexes
CREATE INDEX assoc_pvalue IF NOT EXISTS FOR (ar:AssociationResult) ON (ar.p_value)
CREATE INDEX assoc_phenotype IF NOT EXISTS FOR (ar:AssociationResult) ON (ar.phenotype_key)
CREATE INDEX assoc_runid IF NOT EXISTS FOR (ar:AssociationResult) ON (ar.run_id)

// Study index
CREATE CONSTRAINT gwas_study_id IF NOT EXISTS FOR (gs:GWASStudy) REQUIRE gs.id IS UNIQUE
```

---

## 6. Schema Migration: Adding GraphGWAS to Existing GraphPop Instance

If GraphPop is already running, add GraphGWAS schema without disrupting existing data:

```cypher
// Step 1: Add GWAS-specific properties to existing Sample nodes
MATCH (s:Sample) WHERE s.phenotypes IS NULL
SET s.phenotypes = {}
SET s.is_case = false
SET s.is_control = false
SET s.is_excluded_qc = false

// Step 2: Add carrier count properties to Variant nodes (Option 2 fallback)
MATCH (v:Variant) WHERE v.case_carrier_count IS NULL
SET v.case_carrier_count = 0
SET v.control_carrier_count = 0
SET v.case_carrier_ids = []
SET v.control_carrier_ids = []

// Step 3: Create GWASStudy node for audit trail
CREATE (:GWASStudy {
    id: 'schema_migration_' + toString(datetime()),
    run_id: 'migration_v1',
    name: 'Schema migration from GraphPop to GraphGWAS',
    timestamp: datetime()
})
```

---

## 7. CARRIES vs Alternative Representations — Decision Matrix

| Scenario | Recommended Representation | Rationale |
|---|---|---|
| ≤ 10K samples, full GWAS cohort | Full CARRIES edges | Manageable scale, maximum capability |
| 10K–100K samples, phenotyped subset | Scoped CARRIES (phenotyped samples only) | Feasible scale, preserves traversal semantics |
| > 100K samples, common variants only | CARRIES for rare (AF < 0.05) + DenseWindow cache for common | Optimal performance |
| Any scale, feasibility concerns | Carrier arrays on Variant nodes | Fallback: loses traversal but enables set operations |
| External cohort data, read-only | External Parquet + graph index | Last resort: no graph-native GWAS semantics |

**RULE:** Always prefer the highest-capability option that is feasible at your data scale.
Do not pre-emptively drop to a lower option without benchmarking the higher option first.

---

## 8. Common Schema Mistakes to Avoid

1. **Do not store genotype dosage vectors as Variant properties.**
   `v.dosages = [0, 1, 2, 0, 1, ...]` does not scale beyond ~10K samples.
   Use CARRIES edges or carrier arrays instead.

2. **Do not use STRING for boolean phenotype values.**
   Use `is_case: BOOLEAN` not `phenotype: 'case'` — the boolean is indexable and
   enables fast GWAS filtering.

3. **Do not hardcode phenotype key names in procedures.**
   Always pass `phenotype_key` as a parameter. Cohorts have different trait names.

4. **Do not store raw p-values only.**
   Always also store `p_value_log10 = -log10(p_value)` for Manhattan plot queries.
   Comparing very small floats (1e-300) is numerically unstable.

5. **Do not forget `run_id` on all AssociationResult nodes.**
   Without run_id, you cannot distinguish results from different analysis runs or
   parameter settings. Multiple GWAS runs will accumulate in the graph over time.

6. **Do not create CO_OCCURS edges for all variant pairs.**
   This is O(V²) and will explode the graph. Only materialize CO_OCCURS edges
   above a minimum co-carrier count (e.g., ≥ 3 cases) and in bounded genomic windows.
