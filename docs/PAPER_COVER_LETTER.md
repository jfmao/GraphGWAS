# Cover Letter — Nature Genetics Submission

**Manuscript title:** Graph-native fine-mapping and epistasis detection for
genome-wide association studies

**Corresponding author:** [TBD]

**Submission type:** Technical Report (Methods)

---

Dear Editor,

We submit for consideration at *Nature Genetics* a Technical Report presenting
**GraphGWAS**, an end-to-end graph-native platform for fine-mapping and
epistasis detection in genome-wide association studies.

Standard fine-mappers (SuSiE, FINEMAP, PolyFun) treat every locus as an
isolated correlation matrix. Biology — genes, tissue-specific eQTL networks,
pathway membership, protein-protein interactions — is absent from the
inference, appearing only in post-hoc enrichment analyses. This decoupling
limits the precision of credible sets in the regime that matters most for
novel discovery: sub-genome-wide-significant loci in dense LD regions where
statistical signal alone cannot resolve the causal variant.

We introduce two methods that perform fine-mapping *within* the biology graph
rather than on top of it. **Hierarchical Belief Propagation (HBP)** runs
message passing on a three-layer factor graph of variants, genes, and
pathways, with PPI edges coupling the gene layer. Under default
hyperparameters, HBP is a contraction mapping on the probability simplex
(Banach fixed-point theorem, Supplement S1, Theorem 2), converging in five
rounds to a unique posterior. **L1 Bayesian fine-mapping** treats annotations
as a structural prior in a single-effect regression with empirically learned
weights.

Our principal empirical finding is that L1 **decisively outperforms SuSiE in
the weak-signal regime with informative annotations**:

- 79 independent simulated loci (chr22, β = 0.15, h² = 0.01, tissue-specific
  eQTL causal variants)
- Head-to-head: **L1 wins 27, ties 50, loses 2** (13.5 : 1 ratio, sign-test
  p < 10⁻⁶)
- Rank-#1 rate: **77% vs 57%** for SuSiE
- Mean causal rank: **1.65 vs 2.84**

Separately, HBP *matches* SuSiE on strong and weak signal (70% and 60%
rank-#1 respectively; 22/30 ties vs SuSiE head-to-head under strong signal)
while running at **0.08 s per locus — 20–30× faster** than SuSiE (1.8 s)
and FINEMAP (2.2 s). All four methods maintain **zero false positives**
across 100 null simulations at PIP > 0.5, and are calibrated across the
PIP scale (Supplement Table 3).

Methodologically, the platform includes:

- **LD-pruned co-occurrence epistasis (M1)** — Theorem 1 proves search space
  reduction by 1/k(τ)²; on 1KG chromosome 22 this is 42,000-fold (10.5 × 10⁹
  → 250 × 10³ pairs) without loss of sensitivity on ground-truth
  interactions.
- **A hybrid BGEN + graph architecture** — genotypes stream per-locus from
  BGEN files while annotations live in Neo4j. Validated end-to-end on 1000
  Genomes chromosome 22 (BGEN AF matches graph-stored AF to < 1 × 10⁻⁴).
  The same architecture scales to UK Biobank's 500 K × 93 M matrix without
  code change; application awaits data access.
- **A preloaded multi-omics graph** covering 70.7 M variants, 20,092 genes,
  49 GTEx tissues (43.2 M eQTL edges), and 230 K high-confidence STRING
  interactions — released publicly as a 17 GB Neo4j dump.

We believe this work is appropriate for *Nature Genetics* because it addresses
a standing gap in the fine-mapping literature — the absence of methods that
integrate biological context as a structural prior rather than a flat
annotation vector — and presents a reusable open-source platform whose
performance has been systematically benchmarked against the field's existing
gold standards. The decisive 27-to-2 advantage in the weak-signal regime
maps onto the specific discovery regime where novel causal variants are
most often missed, including rare-variant effects and genomic regions of
dense LD in non-European ancestries.

Software, all benchmark scripts, and graph dumps are available at
github.com/jfmao/GraphGWAS. The benchmark figures and tables in the
manuscript are each reproducible from a single command.

We suggest the following reviewers with deep expertise in fine-mapping and
biobank-scale statistical genetics: [names redacted]. We declare no
competing interests.

Thank you for considering this manuscript.

Sincerely,
[TBD]

---

## Appendix — headline numbers at a glance

| Claim | Value | Source |
|---|---|---|
| Match SuSiE on strong signal | 22/30 ties, 4 wins, 4 losses | Table 2, Fig 3d |
| 20–30× speedup over SuSiE | HBP 0.08 s vs SuSiE 1.8 s | Table 2, Fig 3c |
| Beat SuSiE at weak signal + annotations | 27-2 wins (13.5:1) | Table 5, Fig 6 |
| Zero FPR across 100 nulls | FPR @ 0.5 = 0% all methods | Table 4, Fig 4c |
| PIP calibrated | Within-bin TDR matches expected | Table 3, Fig 4a |
| 42,000× epistasis search reduction | 10.5 × 10⁹ → 250 × 10³ | Supplement, Fig 5a |
| Hybrid BGEN pipeline | AF diff < 1e-4 vs graph | Manuscript §4.6 |
| Pre-loaded multi-omics graph | 43.2 M eQTL, 230 K PPI | Table 1 |
