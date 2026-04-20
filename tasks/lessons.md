# GraphGWAS — Lessons Learned

## How to Use This File
After ANY correction from the user, add an entry here describing:
- What went wrong
- Why it went wrong
- The rule to prevent it next time

Review this file at the start of every session.

---

## Session 1 Lessons (Pre-populated from Design Phase)

### LESSON 001: Always audit schema before writing graph code
**Pattern:** Attempted to write Cypher assuming CARRIES edges exist.
**Cause:** Design documents specify CARRIES but implementation dropped it for unknown reasons.
**Rule:** ALWAYS run the schema discovery protocol from GRAPH_SCHEMA_SKILL.md section 1
before writing any Cypher or Java stored procedure. Never assume a relationship type
or property exists without verifying it in the actual database.

### LESSON 002: CARRIES is the architecturally critical dependency
**Pattern:** Proceeded with multi-locus GWAS implementation without confirming genotype representation.
**Cause:** Assumed the design-specified schema was implemented as designed.
**Rule:** The CARRIES edge status determines which GWAS capabilities are possible.
Single-locus GWAS via fast-path (allele counts only) works without CARRIES.
All individual-level operations (epistasis, cosegregation, disease subgraph) require CARRIES
or a confirmed equivalent. Do not implement Phase 3 features before Phase 1 confirms
the genotype layer is in place.

### LESSON 003: Rare disease ≠ just "add Firth correction"
**Pattern:** Treating rare disease as "same GWAS with small N, use Firth."
**Cause:** Classical GWAS framing applied to graph context.
**Rule:** Rare disease in GraphGWAS requires a fundamentally different approach:
kinship-aware cosegregation (exploit family structure, don't correct it away),
pathway diffusion (aggregate signal across many rare variants), and disease subgraph
extraction (structural rather than statistical signal). Implement all three before
declaring rare disease support complete.

### LESSON 004: AssociationResult nodes must have run_id
**Pattern:** Stored GWAS results without run_id, making it impossible to distinguish runs.
**Cause:** Did not include versioning metadata.
**Rule:** Every AssociationResult node MUST have run_id, timestamp, method, n_cases,
n_controls, and covariates. Without these, graph accumulates results from different
runs that cannot be compared or cleaned up.

### LESSON 005: CO_OCCURS edges must be bounded
**Pattern:** Attempted to build co-occurrence graph for entire chromosome.
**Cause:** Did not appreciate O(V²) complexity.
**Rule:** Always bound CO_OCCURS edge computation to:
(1) minimum co-carrier count ≥ 3,
(2) maximum genomic window ≤ 1 Mb,
(3) only variants above minimum case carrier count.
Full-genome epistasis screen requires pre-filtering to candidate loci.

---

## Session 2 Lessons (Phase 1 Implementation)

### LESSON 006: Neo4j cannot store MAP properties
**Pattern:** Tried to SET s.phenotypes = {key: value, ...} as a MAP property.
**Cause:** Neo4j property values can only be primitives or arrays, NOT nested maps.
**Rule:** Store each phenotype trait as a separate property: `s.pheno_case_control`,
`s.pheno_age`, etc. Use a `phenotypes_loaded: BOOLEAN` flag instead of checking
for MAP existence. GraphMana schema.md says `phenotypes: MAP` but this is aspirational,
not implementable as a native Neo4j property.

### LESSON 007: gt_packed decode is the foundation — validate it first
**Pattern:** Built entire association pipeline before validating that genotype decode matches stored allele counts.
**Cause:** Assumed decode was correct without verification.
**Rule:** Before any statistical computation, verify that unpack_genotypes() on gt_packed
produces AC/AN that exactly matches the stored ac_total/an_total on the same Variant node.
This single check validates the entire genotype access layer. Do it first.

---

## Session 3 Lessons (Phase 2 Implementation)

### LESSON 008: Firth beta matches logistic for common variants
**Pattern:** Firth and standard logistic give nearly identical beta estimates for common
variants (diff < 0.001), confirming the implementation is correct.
**Rule:** Always validate a new regression method by comparing to the existing method on
common variants where both should agree. Divergence on common variants = bug.

### LESSON 009: MPAT directed test misses mixed-direction effects within genes
**Pattern:** APOE gene shows p=0.07 with directed MPAT despite strong single-variant signal.
**Cause:** Variants within APOE have mixed effect directions; the directed sum cancels.
**Rule:** Always run BOTH directed (preserves direction, better for pathways) and undirected
(captures mixed effects, better for genes) MPAT. Report the more significant of the two
for gene-level tests.

---

## Session 4 Lessons (Phase 3-5 Implementation)

### LESSON 010: Neo4j may crash under heavy query load
**Pattern:** Neo4j went down during flow network construction on chr19.
**Cause:** Large result sets from queries returning gt_packed for many variants.
**Rule:** Always use SKIP/LIMIT batching for variant queries. Monitor Neo4j process
during heavy scans. Implement connection retry logic in production code.

### LESSON 011: Graph-native methods require careful scoping
**Pattern:** Co-occurrence graph on 100kb window produced 205K edges — manageable.
**Cause:** All-pairs carrier-set AND is O(V²) but bounded by window + min_cocarrier filter.
**Rule:** Always enforce: (1) window ≤ 1Mb, (2) min_cocarriers ≥ 3, (3) max_variants cap.
Without bounds, a chromosome-wide scan would OOM.

---

## Session 5 Lessons (Yeast 1011 Validation)

### LESSON 012: PLINK .bed/.bim/.fam loses nucleotide identity
**Pattern:** PLINK GWAS matrix encoded alleles as "1"/"2" instead of A/C/G/T.
**Cause:** PLINK binary format only stores allele indices, not actual bases.
**Rule:** Always use the original VCF/gVCF when nucleotide identity matters.
Never rely on PLINK→VCF conversion for data that will be loaded into Neo4j.
bcftools norm can set REF from FASTA but cannot recover ALT from numeric codes.

### LESSON 013: Linear regression for quantitative traits needs different sample selection
**Pattern:** Parallel GWAS hung because quantitative traits set is_case=null,
and get_phenotype_indices() returned empty arrays.
**Cause:** single_locus_scan always called get_phenotype_indices() even for linear method.
**Rule:** For method="linear", use get_all_indices() + get_phenotype_values() instead
of the case/control pathway. Fixed in assoc.py.

### LESSON 014: Population structure inflation is the #1 priority for GWAS calibration
**Pattern:** Many yeast traits showed lambda > 2.0 (caffeine=2.45, benomyl=2.43).
**Cause:** 1011 yeast strains span 34 clades with strong population structure. No PCA covariates.
**Rule:** Always check genomic inflation (lambda) before interpreting GWAS results.
Traits with lambda > 1.1 need population structure correction. GraphGWAS should use
graph Laplacian eigenvectors (not imported PCA) for graph-native correction.

### LESSON 015: Spectral h² captures more variance than SNP-h²
**Pattern:** Spectral h²=0.70 vs published SNP-h²=0.18-0.29 for YPETHANOL.
**Cause:** Laplacian eigendecomposition on rare-variant similarity graph captures
structure that common-variant GRM cannot represent.
**Rule:** This gap is GraphGWAS's strongest scientific argument. Validate across
multiple traits before claiming as a methodological advance.

### LESSON 016: Neo4j result storage is too slow for genome-wide writes
**Pattern:** Storing 1.9M AssociationResult nodes took > 30 min and was killed.
**Cause:** Each batch creates nodes + MATCH for variant links + MATCH for study links.
**Rule:** Default to TSV output (like conventional GWAS tools). Neo4j storage is
optional (--store flag). For graph-native queries on results, load TSV post-hoc.

### LESSON 017: FINEMAP needs signed correlation, not r²
**Pattern:** FINEMAP produced empty .snp files (0.001s runtime) on 200+ variant loci.
**Cause:** LD matrix was r² (squared correlation) but FINEMAP expects signed Pearson
correlation. With r² (all positive), the Bayesian model becomes degenerate.
**Rule:** Always check what LD matrix format a tool expects. FINEMAP: signed correlation.
SuSiE: dosage matrix (computes LD internally). Different tools, different formats.

### LESSON 018: Annotation prior weight must be gentle (α≈0.9)
**Pattern:** L1 with α=0.5 (50% stat + 50% annotation) HURT ranking performance.
Annotations boosted many non-causal variants equally, diluting the statistical signal.
**Cause:** eQTL annotations cover ~26% of variants in a locus. Non-causal eQTL variants
get the same annotation boost as the true causal variant.
**Rule:** Set α=0.85-0.95. Annotations should be a gentle nudge, not equal weight.
Only reduce α (more annotation weight) when annotations strongly discriminate (high
variance across variants in the locus).

### LESSON 020: LASSO/ADMM over-shrinks in fine-mapping; use ridge + Bayes factors
**Pattern:** GRSD with L1+graph-TV penalty via ADMM produced rank 151-517 for causal variants.
After switching to ridge + graph Laplacian, still rank 4-53 (better but not competitive).
**Cause:** LASSO aggressively shrinks most betas to zero, including the causal variant among
hundreds of correlated proxies. Ridge preserves more signal but over-smooths.
**Rule:** For fine-mapping, proper Bayesian variable selection (SuSiE/FINEMAP) beats
penalized regression. Graph structure should augment priors, not replace the likelihood.
HBP works because it uses graph structure for the prior only, keeping statistical evidence
as the primary signal.

### LESSON 021: Pre-cache graph structure to avoid Neo4j bottleneck
**Pattern:** HBP took 223s/locus because each call queried Neo4j for variant→gene→pathway.
After pre-caching chr22 graph structure (9.4s once), HBP runs at 0.08s/locus.
**Rule:** For benchmarks or batch operations, cache the full chromosome's graph structure
in a Python dict once, then pass it to per-locus functions. The cache for chr22 has 514K
entries and loads in ~9s.

### LESSON 019: FAME is not a fine-mapping tool
**Pattern:** Tried to compare FAME to SuSiE/FINEMAP/L1 for per-variant ranking.
**Cause:** FAME estimates epistatic variance components (σ²) partitioned by annotation
categories. It produces no per-variant PIPs.
**Rule:** Match tools by output type. Fine-mapping → PIPs/ranks (SuSiE, FINEMAP, L1).
Heritability partitioning → variance components (FAME, LDSC, S-LDSC).
