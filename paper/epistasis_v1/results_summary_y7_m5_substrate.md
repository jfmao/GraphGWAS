# Paper #2 §Y.7 — M5 RWR substrate rescue: PPI-bridge + focused seeds

**Question:** the §Y.6 cross-species table reports M5 (variant-priority RWR) as
NF for all three canonical pairs. Is the failure parametric (top-K cap too tight,
hub effects, etc.) or structural (the substrate cannot reach the cross-pair)?

**Script:** `tests/m5_variants_benchmark.py`
**Output:** `results/paper2_epistasis/m5_variants_benchmark.json`
**Data:** identical to §Y.6 (yeast 1011 panel chr9/10/12; Arabidopsis 1001G chr1+chr5+chr3 windows; rice 3K RG chr6 narrow windows).

## Headline grid

| Mode | yeast (BCY1×TPK1) | Arabidopsis (FT×FLC) | rice (Hd1×Hd3a) |
|------|---|---|---|
| 0 baseline (top-2K, mutual) | NF / 1,975 | NF / 1,994 | NF / 2,000 |
| A high-cap (top-50K) | NF / 20,036 | NF / 49,994 | NF / 49,782 |
| B hub-deflated `(P_a+P_b)/(deg_a·deg_b)` | NF / 1,977 | NF / 1,994 | NF / 2,000 |
| C focused seeds (no PPI bridge) | **NF / 0** | **NF / 0** | **NF / 0** |
| D PPI bridge + unfocused seeds | NF / 1,996 | NF / 1,994 | NF / 2,000 |
| **E PPI bridge + focused seeds** | **1 / 10** (q=1.8e-9) | **1 / 14,924** (q→0) | **1 / 25,608** (q→0) |

## Diagnosis (the smoking gun)

Mode C produces **zero cross-pairs** for every species. Cross-pairs are
constructed by iterating `for si in g1_seeds: for sj in g2_seeds:` and
filtering on positive RWR mass. Zero pairs survive this filter, meaning
**the RWR substrate gives the canonical cross-pairs identically zero
mutual mass**.

The reason is that `_bipartite_from_dict` builds the variant↔gene
adjacency only from the cache's `genes` field, ignoring the `ppi` field
entirely. So the variant-variant transition `M_vv = P_vg @ P_gv` is the
2-step variant→gene→variant walk: two variants can be reached from each
other only if they share a common annotated gene (or sit in the same
multi-gene pathway entity). For:

- **BCY1 ↔ TPK1** (yeast): no shared annotation in the cache; only the
  PPI edge connects them.
- **FT ↔ FLC** (Arabidopsis): same.
- **Hd1 ↔ Hd3a** (rice): same.

In all three species the canonical PPI edge is added by
`graphgwas.canonical_ppi.inject_canonical_ppi_edges()` to the cache's
`ppi` field — but the M5 substrate then ignores it. The walk is
structurally trapped within the g1 component, never reaches g2 variants,
and the cross-pair mass is identically zero.

## Fix

**Mode E = PPI-augmented substrate + focused seeds.** Two coupled changes:

1. **Substrate (Mode D component):** extend the variant-variant transition to
   `M_vv = P_vg @ P_gg @ P_gv`, where `P_gg` is the row-normalised
   `I + 1[(g_a, g_b) ∈ PPI]`. This adds a single gene-gene PPI hop into the
   walk: `variant → gene → gene' → variant`. The bridge is symmetric and
   includes self-loops (so the walk can stay at the same gene).

2. **Seeding (Mode C component):** seed RWR ONLY from the variants annotated
   to `g_1 ∪ g_2` (no other annotated variants in the seed set), then score
   every `g_1 × g_2` cross-pair. This concentrates RWR mass in the small
   subspace where the PPI bridge matters; without focusing, the bridge
   contribution gets diluted across thousands of competing seeds (Mode D shows
   this — the substrate alone fails when the seed pool is large).

Both changes are necessary: Mode C alone fails (zero cross-pair mass without
the bridge), and Mode D alone fails (bridge mass diluted by competing seeds).

## Why this matters for paper #2

§Y.5 already noted that the variant-priority M5 fails on yeast and a
gene-priority RWR rescues it via an explicit gene-gene heterograph that
includes PPI. The §Y.7 result generalises that observation: the failure
mode is the same across species, the cause is structural (substrate omits
PPI), and the fix preserves variant-priority granularity (Mode E still
returns variant-pair recovery, not just gene-pair recovery — useful for
downstream variant-level fine-mapping).

The M5★ workflow composes naturally as a two-stage refinement: a coarse
selector (M2 motif filtering, or §Y.5 gene-priority RWR, or M3 case-unique
edges) hands a candidate gene pair `(g_1, g_2)` to M5★, which then ranks
the best variant-pair representative.

## Bug fixed in passing

`_bipartite_from_dict` in `tests/validate_m2_yeast.py` returned variants in
cache iteration order rather than in the input `variant_ids_to_keep` order.
For Arabidopsis/rice the orders coincide (cache built in chromosomal order
matches VCF order), but for yeast they differ, which initially caused
the index lookups in `m5_variants_benchmark.py` to misalign. Fixed to
preserve input order; re-tested on the cross-species runner — same
canonical-pair ranks as before, so Paper #2 §Y.4/§Y.6 numbers are unchanged.

## Reproduction commands

```bash
# Run all 3 species, all 6 modes, with default top-K-high=50000
python tests/m5_variants_benchmark.py --species all --beta 3.0 --seed 2026

# Single species
python tests/m5_variants_benchmark.py --species rice --beta 3.0 --seed 2026
```

Required data layout: identical to `validate_epistasis_cross_species.py`
(see `paper/epistasis_v1/results_summary_y6.md`).

The Mode E kernel is `compute_rwr_matrix(..., use_ppi_bridge=True)` +
`run_focused_seeds(..., use_ppi_bridge=True)`. Both are pure functions
in `tests/m5_variants_benchmark.py`; promoting them to
`src/python/graphgwas/epistasis_v2.py` (e.g. as
`mutual_rwr_focused_with_ppi_bridge`) is a paper-#2 revision target.
