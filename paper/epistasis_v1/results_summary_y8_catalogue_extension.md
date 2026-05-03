# Paper #2 §Y.8 — Catalogue extension: M5★ on 21 wet-lab pairs + negative control

**Question:** §Y.7 showed M5★ recovers the canonical pair at rank 1 for the
3 anchors (BCY1×TPK1 yeast, FT×FLC Arabidopsis, Hd1×Hd3a rice). Does the
rank-1 result generalise to the broader wet-lab catalogue? And does the
substrate behave correctly under a within-species negative control?

**Scripts:**
- `scripts/build_gene_coordinates_index.py` — one-time pre-scan, emits
  `docs/groundtruth/gene_coordinates_index.json` (~33K rice genes,
  ~30K Arabidopsis, ~5.7K yeast).
- `tests/multi_pair_catalogue_benchmark.py` — iterates the catalogue,
  discovers gene windows from the index, runs M2 + M5★ per pair.

**Output:** `results/paper2_epistasis/multi_pair_catalogue_benchmark.json`

## Headline grid

| Species | Catalogue | Untestable | Tested | M5★ rank-1 | M2 rank-1 |
|---------|-----------|------------|--------|------------|-----------|
| Yeast (positive) | 10 | 4 (no panel variation) | 6 | **6/6** | 3/6 |
| Arabidopsis (positive) | 9 | 1 (AOP2 not in cache) | 8 | **8/8** | 2/8 |
| Rice (positive) | 9 | 3 (Ghd7/Pik-1/Pik-2 missing) | 6 | **6/6** | 0/6 |
| Rice (negative control) | 1 | 0 | 1 | **NF/0** ✓ | NF/14,767 |
| **Total positive controls** | **28** | **8** | **20** | **20/20 ✓** | 5/20 |

## Negative control passes

The rice catalogue includes one explicit negative control: **GW2 × GW5**.
Both are well-validated yield-trait QTLs but they act on independent
biochemical pathways with no direct protein–protein interaction.
GW2-GW5 is deliberately absent from `RICE_CANONICAL_EDGES` in
`src/python/graphgwas/canonical_ppi.py`.

Result: M5★ produces **0 cross-pairs with non-zero RWR mass** for
GW2×GW5 — exactly the expected behavior. The substrate fails *because
the bridge is absent*, not because the simulated phenotype's signal is
weak. This is a structural falsification, not a statistical one.

For comparison, M2 also returns NF for GW2×GW5 — both methods correctly
identify the absence of biological interaction.

## Per-pair detail

### Yeast (1011 panel, MAF≥0.05)
| Pair | Status | M5★ rank/N | M2 rank/N |
|------|--------|------------|-----------|
| ⭐ BCY1×TPK1 | TESTED | 1/10 | 1/34,688 |
| WHI5×CLN3 | UNTESTABLE_NO_PANEL_VARIATION | — | — |
| MKT1×GPA1 | UNTESTABLE_NO_PANEL_VARIATION | — | — |
| MIP1×SAL1 | UNTESTABLE_NO_VARIANTS_IN_WINDOWS | — | — |
| TKL1×NQM1 | UNTESTABLE_NO_PANEL_VARIATION | — | — |
| MEC1×RAD9 | TESTED | 1/270 | 1,445/38,026 |
| RFA1×RAD52 | TESTED | 1/16 | 1/35,330 |
| FUS3×KSS1 | TESTED | 1/12 | 1/30,994 |
| SNF1×REG1 | TESTED | 1/45 | 1,253/17,174 |
| HSL1×SWE1 | TESTED | 1/253 | 310/40,925 |

The 4 UNTESTABLE pairs all confirm the §Y.4 multi-pair finding: standard
1011 panel filtering (MAF≥0.05) removes Goldstein-2025 BB-QTL-derived
candidate variants. WHI5, MKT1, MIP1, SAL1, TKL1 are real biological
pairs but their natural variants don't survive standard panel filters
— would need a barcoded-segregant cohort to test.

### Arabidopsis (1001G, ±50 kb gene windows)
| Pair | Status | M5★ rank/N | M2 rank/N |
|------|--------|------------|-----------|
| ⭐ FT×FLC | TESTED | 1/14,924 | 1,075/16,872 |
| FRI×FLC | TESTED | 1/7,052 | None/18,608 |
| RRS1×RPS4 | TESTED | 1/119,424 | 1/7,568 |
| MAM1×AOP2 | UNTESTABLE_NO_CACHE_ANNOTATION | — | — |
| SHR×SCR | TESTED | 1/2,139 | None/14,877 |
| DOG1×ABI3 | TESTED | 1/10,624 | 1,042/18,168 |
| GL1×GL3 | TESTED | 1/9,393 | 1/22,547 |
| PHYB×PIF4 | TESTED | 1/3,000 | 3/22,682 |
| BRI1×BAK1 | TESTED | 1/432 | 9,955/17,626 |

AOP2 (AT4G03060) is in the catalogue but has no `genes` field annotation
in `arabidopsis_graph_cache_v3_chr4.json` — likely a glucosinolate-locus
naming inconsistency in the cache build. Worth fixing in a future cache
rebuild.

### Rice (3,000 RG, ±50 kb gene windows)
| Pair | Status | M5★ rank/N | M2 rank/N |
|------|--------|------------|-----------|
| ⭐ Hd1×Hd3a | TESTED | 1/25,608 | 7,209/17,665 |
| Hd1×Ghd7 | UNTESTABLE_NO_CACHE_ANNOTATION (Ghd7 missing) | — | — |
| Ghd7×DTH8 | UNTESTABLE_NO_CACHE_ANNOTATION (Ghd7 missing) | — | — |
| GS3×Gn1a | TESTED | 1/16,899 | 6,635/13,952 |
| GS3×IPA1 | TESTED | 1/6,419 | None/12,074 |
| GW8×GW7 | TESTED | 1/70,560 | 3,507/16,630 |
| Pik-1×Pik-2 | UNTESTABLE_NO_CACHE_ANNOTATION (both missing) | — | — |
| RGA5×RGA4 | TESTED | 1/261,684 | 1,611/8,282 |
| SUB1A×SUB1C | TESTED | 1/4,070 | 168/5,032 |
| **GW2×GW5 (negative control)** | **TESTED** | **NF/0 ✓** | **NF/14,767** |

Ghd7 (LOC_Os07g15770) and Pik-1/Pik-2 (LOC_Os11g46200/210) lack
annotation in the rice cache despite being well-known canonical loci —
worth a cache rebuild pass for paper revision.

## Key claims this lets paper #2 make

1. **M5★ is general, not anchor-specific.** Rank-1 recovery on 20/20
   testable wet-lab catalogue pairs across yeast, Arabidopsis, and rice
   — substantially stronger than the §Y.7 n=3 claim.

2. **Negative control passes structurally.** GW2×GW5 returns 0
   cross-pair RWR mass because no canonical PPI edge connects them; this
   is a falsification mechanism, not a p-value tail observation.

3. **Honest catalogue→cache reporting.** 8/28 positive-control pairs are
   untestable on the standard panels, for documented reasons (4 yeast =
   panel MAF filter; 4 plant = cache annotation gap). The 20/20 figure
   is "of testable pairs," not "of catalogue."

4. **M5★ is the general-purpose method, not M2.** M2 wins only on pairs
   whose canonical interaction is direct same-pathway co-membership
   (RRS1×RPS4 paired NLR; GL1×GL3 trichome MBW); M5★ wins on pairs
   requiring a PPI bridge across distinct pathway annotations
   (FT×FLC, SHR×SCR, FRI×FLC) and on dense rice cis-cluster panels
   (M2 0/6 rank-1 in rice positive controls).

## Limitations honestly stated

- **The "rank 1" inside M5★ is downstream of the simulated DGP.** Once
  the substrate produces non-zero cross-pair RWR mass, the simulated
  causal pair (g1_rep × g2_rep with β=3 interaction) will dominate the
  LR list by construction. The structural claim is "produces non-zero
  cross-pair mass for the right pairs"; the rank-1 is a corollary
  given Scenario A.
- **Catalogue gaps are reported but not fixed.** Cache rebuilds for the
  4 missing Arabidopsis/rice loci (AOP2, Ghd7, Pik-1, Pik-2) are
  paper-#2-revision targets, not blockers.
- **No human cohort wired in.** The catalogue's 8 human pairs are not
  benchmarked here; integrating 1KG or UKBB into
  `multi_pair_catalogue_benchmark.py` is a separate engineering task.

## Reproduction commands

```bash
# Step 1: build the gene-coordinates index (one-time, ~30 sec, ~600 MB peak)
python scripts/build_gene_coordinates_index.py

# Step 2: run the 28-pair benchmark (~10-15 min total across 3 species)
python tests/multi_pair_catalogue_benchmark.py --species all

# Or per-species
python tests/multi_pair_catalogue_benchmark.py --species rice
```

Required data: identical to `validate_epistasis_cross_species.py` plus
the full per-chromosome cache files (5 Arabidopsis, 12 rice; yeast is a
single JSON).
