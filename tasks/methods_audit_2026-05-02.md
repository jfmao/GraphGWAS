# Method-completeness audit — 2026-05-02

Status of every detection method defined in `epistasis_v2.py` and
`epistasis_higher_order.py`, as of commit `e92ab2c`.

| Method | Code | Unit tests | §Y.4 yeast | §Y.5 retry | §Y.6 cross-species | Paper text |
|---|---|---|---|---|---|---|
| **M1** LD-pruned exhaustive  | `ld_pruned_cooccurrence_*` | yes | A/B/C: NF | — | Arabidopsis (partial) | §Y.4 ✓ |
| **M2 P1** same-gene motif    | `motif_filtered_epistasis_*` | yes | A/B/C ✓ | — | Arabidopsis (partial) | §Y.4 ✓ |
| **M2 P2** same-pathway motif | `motif_filtered_epistasis_*` | yes | A/B/C ✓ | — | Arabidopsis (partial) | §Y.4 ✓ |
| **M2 P3** protein-interaction motif | `motif_filtered_epistasis_*` (line 405) | yes | **NEVER** (script only enables P1+P2) | — | NEVER | **GAP** |
| **M3** differential subgraph | `differential_subgraph_*` | yes | A/B/C ✓ | — | Arabidopsis (partial) | §Y.4 ✓ |
| **M4** dark-matter           | `dark_matter_epistasis_*` | yes | A/B/C ✓ | — | Arabidopsis (partial) | §Y.4 ✓ |
| **M5 V-priority** mutual RWR | `mutual_rwr_pair_scores` | yes | A/B/C: NF | — | — | §Y.4 ✓ |
| **M5 G-priority** gene RWR (NEW)  | `gene_pair_rwr_scores` (added today) | unit-test gap | — | A: rank 1522 q=1.08e-5 | — | §Y.5 ✓ |
| **M5 k=3 triplet** RWR       | `mutual_rwr_triplet_scores` | trivial-case test only | **NEVER** on real data | — | — | **GAP** |
| **M6** GNN attention         | `gnn_attention_pairs` (NotImpl) | NotImpl assert | — | — | — | abstract+methods only |
| **M7** persistent homology   | `topological_epistasis` (NotImpl) | NotImpl assert | — | — | — | abstract+methods only |

## Two material gaps for paper #2

### Gap 1 — M2 P3 (protein-interaction motif) never benchmarked
**Where**: `tests/validate_m2_yeast.py:847` and `tests/validate_epistasis_cross_species.py:434`
both pass `motifs=["same_gene", "same_pathway"]`. The P3 motif (variants in
genes that are PPI partners) was implemented at `epistasis_v2.py:405` but
never enabled in any §Y.4 / §Y.6 run.

**Why this matters**: P3 is arguably the most biologically-targeted motif —
direct binding partners, not co-membership in a broad pathway.  For
BCY1×TPK1 specifically (canonical PKA regulatory:catalytic binding pair),
P3 is the natural-fit motif — but only succeeds if the cache's PPI field is
correctly populated. Section §Y.5 found the cache's PPI field is broken
(BCY1 lists itself as its only partner), so P3 will likely also fail on
yeast cache — with or without canonical PPI injection.

**Action**: extend §Y.4 yeast benchmark to include P3, run with both the raw
cache and the BioGRID-injected cache (matching §Y.5's two-cell design).
Document P3's contribution honestly. If P3 fails on raw cache because of the
PPI-curation bug, this is consistent with the §Y.5 finding and reinforces
the paper-#2 message: *graph-priority methods are rate-limited by gene-gene
graph quality*.

### Gap 2 — k=3 (triplet) RWR never benchmarked on real data
**Where**: `mutual_rwr_triplet_scores` in `epistasis_higher_order.py:331`.
Unit test at `tests/test_random_walk_epistasis.py:90` only validates the
trivial "no triangles in 2-pair graph" case.

**Why this matters**: k≥3 detection is the **structurally-unique paper-#2
claim** — pairwise BOOST/MAPIT/MDR cannot formulate it. Currently we have
ZERO real-data benchmark of triplet detection.

**Action**: pick one yeast triplet with strong functional validation
(e.g., the cAMP-PKA core trio BCY1-TPK1-TPK2 or BCY1-TPK1-TPK3, or
SUB1A-SUB1B-SUB1C from rice submergence-tolerance literature). Run
`mutual_rwr_triplet_scores` with appropriate seed pool, measure recovery
rank. Compare runtime + recovery to brute-force enumeration of triplets
(O(n³)). The paper-#2 §Y.5 §Y.5 narrative needs this.

## Out-of-scope for paper #2 (mark as "future work")

- M6 GNN-attention: requires PyTorch Geometric training infrastructure.
  Defer to paper #3 (rule-mining) or future revision.
- M7 persistent homology: requires Gudhi or Giotto-TDA. Defer.

## Recommended order

1. (Now) Extend `validate_m2_yeast.py` to test P3 motif (with and without canonical PPI).
2. (Next) Implement triplet benchmark on a yeast cAMP-PKA k=3 ground-truth
   (BCY1+TPK1+TPK2/3) using `mutual_rwr_triplet_scores`.
3. (After REGENIE) §Y.2 conservativeness rescue benchmark.
