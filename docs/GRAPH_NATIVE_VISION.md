# Graph-Native GWAS Vision — Phase 3+ Game-Changing Methods

## The Thesis

Most of what classical GWAS does (logistic regression, chi-squared, burden tests) treats
the genome as a matrix. GraphGWAS's unique value comes from methods that are **impossible
in any matrix tool** — methods where the graph topology IS the computation.

---

## Method 1: Variant Co-Occurrence Network → Epistasis via Community Detection

**The problem:** Testing all M² variant pairs for epistasis is O(10¹⁴) — intractable.
The field has essentially given up on systematic epistasis detection.

**The graph solution:**
1. Build variant co-occurrence graph: edge (v1, v2) weighted by co-carrier count in cases
2. This graph is SPARSE — most pairs don't co-occur at high frequency
3. Run community detection (Leiden/Louvain) on the case-enriched co-occurrence graph
4. Each community = a candidate epistatic module (variants whose joint effect exceeds marginals)
5. Statistical test: compare community enrichment to permuted null (shuffle case labels)

**Complexity:** O(E) where E = sparse edges, not O(M²). Reduces search from 10¹⁴ to ~10⁶-10⁸.

**Why no matrix tool can do this:** The co-occurrence structure is a graph property. In a matrix,
you'd need to compute all pairwise correlations — which IS the O(M²) step you're trying to avoid.

---

## Method 2: Max-Flow / Min-Cut from Phenotype to Biology

**The problem:** GWAS, gene mapping, and pathway analysis are done sequentially with separate
tools, compounding multiple testing and losing power at each step.

**The graph solution:**
1. Model as a flow network: Sample → (CARRIES) → Variant → (HAS_CONSEQUENCE) → Gene → (IN_PATHWAY) → Pathway
2. Source: phenotype signal on case Sample nodes (capacity = 1)
3. Edge capacities: CARRIES weighted by 1/AF (rare variants carry more signal)
4. Compute max-flow from cases to each pathway
5. Min-cut = the most parsimonious variant-gene set connecting phenotype to biology

**Statistical test:** ONE test — permute phenotype labels, recompute flow, empirical null.
Not 10M variant tests + gene enrichment + pathway enrichment.

**Why this is game-changing:** Unifies three analysis layers into a single graph computation.
The min-cut IS the disease architecture.

---

## Method 3: Graph Spectral Phenotype Decomposition

**The problem:** PCA captures global ancestry from common variants but misses rare-variant-based
structure, local ancestry, and pathway-based similarity.

**The graph solution:**
1. Build sample-sample similarity graph from shared rare variant carriage
2. Compute graph Laplacian eigendecomposition
3. Decompose phenotype signal into graph frequencies
4. Low-frequency = genetic component; high-frequency = noise
5. Graph-filtered phenotype = denoised genetic signal

**Why this generalizes PCA:** PCA = eigenvectors of common-variant GRM (a special case).
Graph spectral filtering works with any topology — rare variants, pathways, protein interactions.

---

## Method 4: LD Fine-Mapping via Graph Centrality

**The problem:** Fine-mapping with SuSiE/FINEMAP takes hours per locus and requires
LD reference panels.

**The graph solution:**
1. Build LD graph within each significant locus (edges = r² > threshold)
2. The causal variant is the HUB — highest betweenness centrality
3. All proxy associations "flow through" the causal variant in the LD graph
4. Rank by centrality → credible set in milliseconds

**Why this works:** The LD graph structure IS the prior that Bayesian fine-mapping
methods spend hours estimating.

---

## Method 5: Multi-Phenotype Diffusion for Pleiotropy

**The problem:** Diseases share genetic architecture, but GWAS analyzes one phenotype at a time.

**The graph solution:**
1. Run phenotype diffusion for multiple diseases simultaneously
2. Pathways accumulating signal from multiple phenotypes = pleiotropy hubs
3. Pleiotropy hubs = shared mechanism → drug repurposing targets

**Connection to GNExT:** GNExT does this post-hoc via Drugst.One network medicine.
GraphGWAS does it natively as a multi-signal diffusion problem.

---

## Implementation Priority

| Method | Phase | Novel Power | Feasibility |
|--------|-------|-------------|-------------|
| Co-occurrence + community detection | 3 | Epistasis at genome-wide scale | High (gt_packed + Leiden) |
| Max-flow / min-cut | 3 | Unified variant→gene→pathway test | Medium (flow algorithms well-studied) |
| Graph spectral filtering | 4 | Missing heritability recovery | Medium (Laplacian computation) |
| LD centrality fine-mapping | 3 | Millisecond fine-mapping | High (simple after LD computation) |
| Multi-phenotype diffusion | 4 | Drug repurposing | High (extension of single-phenotype) |
