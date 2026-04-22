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

Standard fine-mappers — SuSiE, FINEMAP, their Cui et al. 2024 infinitesimal
extensions (SuSiE-inf, FINEMAP-inf), Polyfun+SuSiE, and the Wu et al. 2026
genome-wide method SBayesRC — are blind to biological context. Each operates
on a per-locus correlation matrix or a genome-wide MCMC fit, with functional
annotations entering only as flat per-variant priors (Polyfun) or
post-hoc enrichment. None integrates the variant→gene→pathway→PPI graph
that links statistical evidence to mechanistic interpretation.

We introduce two methods that perform fine-mapping *within* the biology graph.
**Hierarchical Belief Propagation (HBP)** runs message passing on a
three-layer factor graph (variants, genes, pathways), with PPI coupling at
the gene layer. Under default hyperparameters HBP is a contraction mapping
on the probability simplex (Banach fixed-point theorem, Theorem 2,
Supplement S1) — converging in five rounds to a unique posterior in
~80 ms per locus. **L1 Bayesian fine-mapping** treats annotations as a
structural prior in a single-effect regression with empirically learned
weights.

We benchmarked GraphGWAS comprehensively against **all six methods in
the canonical Wu et al. 2026 NG fine-mapping reference (Fig. 4b)**:

- **HBP matches state-of-art accuracy at 6–60× the speed.** On 30 1KG
  chr22 F1 simulations (h² = 0.05): HBP rank-#1 20/30 vs SuSiE-inf
  22/30, SuSiE 21/30, FINEMAP-inf 21/30 — within 1–2 of the field
  leaders, at 0.016 s vs 0.16–2.5 s.
- **L1 wins decisively in the weak-signal + informative-annotation
  regime.** On 79 yeast loci (β = 0.15, h² = 0.01, tissue-specific
  eQTL causal): L1 wins 27, ties 50, loses 2 vs SuSiE — a 13.5:1
  ratio (sign-test p < 10⁻⁶), 77% rank-#1 vs 57%, mean rank 1.65
  vs 2.84.
- **L1 ties Polyfun-style annotation-informed SuSiE on rank** (11/20
  vs 11/20) and **wins 26× on speed** (0.055 s vs 1.45 s).
- **All methods FPR-controlled** (0% FPR at PIP > 0.5 across 100
  null replicates) and **PIP-calibrated** across four h² levels.
- **SBayesRC integration verified end-to-end** on canonical UKB
  example sumstats (1.15 M HM3 SNPs, 473 PIP > 0.9 in 51 s);
  positioned honestly as complementary (genome-wide GWFM) rather
  than competitive with our region-specific methods.

The platform also delivers:

- **42,000-fold epistasis search-space reduction** (Theorem 1):
  10.5 × 10⁹ exhaustive pairs → 250 × 10³ LD-pruned pairs on 1KG
  chr22, with rank-#1 ground-truth detection in 10/10 simulated
  interactions.
- **Three-kingdom species generalisation**: identical codebase
  applied to fungi (yeast 1011 Genomes, 35 traits), animals
  (1000 Genomes), and plants (*Arabidopsis* 1001 Genomes). On
  Arabidopsis FT10, HBP narrows the chr4 flowering-time peak to
  a single variant at PIP = 0.989 in 0.1 s.
- **Cross-ancestry generalisation**: HBP on EUR/AFR/EAS
  superpopulations of 1KG chr22 — best on AFR (53% rank-#1, mean
  PIP 0.346), consistent with AFR's less-correlated LD; method
  does not break on non-European LD.
- **Hybrid BGEN + graph architecture**: 2-bit packed genotypes
  stream per locus from BGEN files (1.6 TB at UKB scale) while the
  Neo4j graph stores ~82 M biological edges. End-to-end pipeline
  validated on 1KG chr22 (BGEN AF matches graph AF to < 1 × 10⁻⁴).
- **A pre-loaded multi-omics graph** released as a 17 GB Neo4j
  dump: 70.7 M variants, 20,092 genes, 49 GTEx tissues (43.2 M
  eQTL edges), 230 K STRING PPI edges, 370 K ENCODE regulatory
  elements.

We believe this work fits *Nature Genetics* because it addresses a standing
gap in the fine-mapping literature — methods that integrate biology as a
structural prior rather than a flat annotation — and answers it with a
reusable open-source platform comprehensively benchmarked against the
field's existing gold standards. The decisive 27-to-2 advantage in the
weak-signal regime maps onto the discovery regime where novel causal
variants are most often missed, including rare-variant effects and dense
LD regions in non-European ancestries.

Software, all benchmark scripts, all ten paper figures, all nine paper
tables, and the multi-omics graph dumps are available at
github.com/jfmao/GraphGWAS. Every figure and table is reproducible from
a single command (see `docs/REPRODUCIBILITY.md`).

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
| Match SuSiE/SuSiE-inf on accuracy | 20–22/30 rank-#1 (within 1–2) | §2.4, Fig 3 |
| 6–60× faster than all baselines | HBP 0.016 s vs 0.16–2.5 s | §2.4, Fig 3c |
| Beat SuSiE at weak signal + tissue-specific eQTL | 27-2 wins (13.5:1) | §2.3, Table 5, Fig 6 |
| Tie Polyfun-style SuSiE on rank, 26× faster | 11-11 wins, 0.055 s vs 1.45 s | §2.8 |
| Zero FPR across 100 nulls | FPR @ PIP > 0.5 = 0% all methods | §2.6, Table 4, Fig 4c |
| PIP calibrated | Bin TDR matches expected on 800 sims | §2.6, Table 3, Fig 4a |
| 42,000× epistasis search reduction | 10.5 × 10⁹ → 250 × 10³ | §2.7, Fig 5a, Theorem 1 |
| Power scales with N | Mean PIP 0.54 → 0.77 at N = 500 → 3000 | §2.9, Fig 8 |
| Cross-ancestry: works best on AFR | 53% rank-#1, mean PIP 0.346 | §2.10, Fig 9 |
| Cross-species: Arabidopsis FT10 | PIP = 0.989 on chr4 peak in 0.1 s | §2.11 |
| Hybrid BGEN pipeline | AF diff < 1 × 10⁻⁴ vs Neo4j | §4.6 |
| SBayesRC integration verified | 1.15 M SNPs, 473 PIP > 0.9 in 51 s | §2.4 |
| Multi-omics graph released | 43.2 M eQTL + 230 K PPI + 370 K cCRE | §2.1 |
