# Paper #2 §Y.6 — Cross-species replication: 5-method recovery on FT×FLC and Hd1×Hd3a

**Question:** does the §Y.4 yeast result (M2 motif-filtered dominates) generalise
to other species' canonical epistatic pairs? Specifically: does any single
graph-native method recover all three of (BCY1×TPK1 yeast, FT×FLC Arabidopsis,
Hd1×Hd3a rice) under a single uniform pipeline?

**Script:** `tests/validate_epistasis_cross_species.py`
**Output:** `results/paper2_epistasis/cross_species_validation.json`
**Phenotype DGP:** Scenario A — `y = beta · G_X · G_Y + eps`, beta=3.0, seed=2026.
**Motifs (M2):** `same_gene + same_pathway + protein_interaction` with canonical-PPI injection.

## Headline grid

| Species              | Pair          | n      | Variants (MAF≥0.05) | M1 rank | M2 rank | M3 rank | M5 rank |
|----------------------|---------------|--------|---------------------|---------|---------|---------|---------|
| *S. cerevisiae*      | BCY1×TPK1     | 1,011  | 6,625               | NF      | **1,148** (q=1.7e-5) | NF      | NF      |
| *Arabidopsis*        | FT×FLC        | 1,135  | 13,800              | NF      | **3** (q=2.3e-2, PPI) | NF      | NF      |
| *Oryza sativa*       | Hd1×Hd3a      | 3,024  | 9,779               | **32** (q=6.5e-287, co-occ) | 19,057 (PPI) | 2,405 (q=3.6e-6, case_unique) | NF |

M4 skipped on rice via `--skip-heavy` (3K samples × 500 vars OOMs without it).

## Per-species outcome

### Arabidopsis (FT × FLC, chr1+chr5, ±50 kb each + chr3 control)
- 13,800 MAF≥0.05 SNPs; 3,351 in cache.
- 91 FT-annotated + 164 FLC-annotated representative variants in the loaded panel.
- Canonical-PPI injection: 8/11 edges applied (3 dropped: SHR-SCR chr4, PHYB-PIF4 chr2, BRI1-BAK1 chr4 — outside loaded windows).
- **M2 wins**: rank 3 / 13,980 pairs; q=2.31e-02; detected via `protein_interaction`. Without canonical-PPI injection M2 misses (cache lacks FT↔FLC edge — same curation gap as yeast BCY1↔TPK1).
- M1 NF over 143K tested pairs; M3/M4/M5 all NF.

### Rice (Hd1 × Hd3a, both on chr6)
- Loader narrowed to chr6:9.25-9.45 Mb (Hd1) + chr6:2.75-2.95 Mb (Hd3a) + chr6:5.5-5.7 Mb control. **Critical:** the original chr6:1-12 Mb load yielded 183 K MAF-filtered variants which OOM-killed M3 (M1's per-variant `dosage_list_v` copy = 4.4 GB; 1.7M `_test_interaction` calls accumulated through numpy's memory pool).
- 9,779 MAF≥0.05 SNPs after the narrow load; 7,894 in cache.
- 194 Hd1 + 132 Hd3a representative variants.
- Canonical-PPI injection: 5/13 edges applied (8 dropped: genes outside chr6 windows).
- **M1 wins**: rank 32 / 55,617; q=6.5e-287; beta = -1.86 (LR interaction). The flip vs yeast/Arabidopsis is driven by (i) 3× sample size, (ii) both reps common (MAF 0.48, 0.50) so LD pruning keeps them, (iii) M2 burns its motif-pair budget on dense same_gene cis-clusters in the Hd1 and Hd3a regions.
- M2 ranks the canonical pair at 19,057 / 43,502 (PPI motif fires structurally but is buried by within-Hd3a cis-pairs; best M2 hit is Chr6:2921530 × 2923797, both inside Hd3a, p=5.3e-34).
- M3 contributes a complementary signal: rank 2,405 / 15,840; q=3.6e-6; case_unique edge type. Beta = +0.52.
- M5 still NF (top-2,000 cap on 6.0M non-zero RWR pairs is too tight on the dense rice graph).

## Implication for paper #2

No single graph-native method is uniformly best across the three species:

- **S. cerevisiae**: M2 wins.
- **Arabidopsis**: M2 wins (with mandatory canonical-PPI injection).
- **Rice**: M1 wins; M3 contributes; M2 underperforms.

Paper-#2 framework recommendation: run M1+M2+M3 in parallel and report each
canonical pair at the best-of-three rank with method attribution. Total runtime
on rice = 22+10+10 ≈ 42 s on a single CPU; M5 adds 47 s; M4 adds ~2 min.
The five-method runtime is comparable in cost to a single REGENIE step-2 pass,
naturally parallelisable across pairs.

## Reproduction commands

```bash
# Both species in one JSON (recommended)
python tests/validate_epistasis_cross_species.py --species all --beta 3.0 --seed 2026 --skip-heavy

# Per-species (overwrites the merged JSON; do NOT mix with --species all)
python tests/validate_epistasis_cross_species.py --species arabidopsis --beta 3.0 --seed 2026
python tests/validate_epistasis_cross_species.py --species rice        --beta 3.0 --seed 2026 --skip-heavy
```

Required data layout:
- `tests/data/arabidopsis/1001genomes_snp-short-indel_only_ACGTN.vcf.gz` (18 GB)
- `data/arabidopsis/arabidopsis_graph_cache_v3_chr{1,5}.json`
- `data/rice_3k/rice_3k.{pgen,pvar,psam}` (PLINK 2)
- `data/rice_3k/annotations/rice_graph_cache_v2_Chr6.json`

Per-species canonical PPI edge lists live in `src/python/graphgwas/canonical_ppi.py`
(`ARABIDOPSIS_CANONICAL_EDGES`, `RICE_CANONICAL_EDGES`).
