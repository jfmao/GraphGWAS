# GraphGWAS — Fine-mapping Input/Output Guide

A practical how-to for end users: prepare inputs, run fine-mapping, check
outputs, and diagnose problems. Complements the
[command reference](manual/commands/finemap.md) and the
[Pan-UKB quickstart vignette](../vignettes/fine-mapping-quickstart.md).

This guide covers the eight fine-mapping methods exposed through the GraphGWAS
CLI and Python API:

- **HBP** (`hbp`) — hierarchical belief propagation on the variant–gene–pathway
  factor graph, calibrated PIPs.
- **GAFM** (`l1`) — graph-augmented fine-mapping with adaptive α blend,
  calibrated PIPs.
- **GAFM-MX**, **HBP-MX**, **ENS** — mixture-prior posterior reweightings
  of GAFM/HBP and their ensemble; **operational ranking scores, not
  calibrated posteriors** (see [§5](#5-output-interpretation)).
- **CLGF** (`clgf`), **GLEM** (`glem`) — cross-locus EM and
  graph-latent-embedding fine-mapping; implemented and theoretically
  supported but **not benchmarked at the depth of HBP/GAFM** in the
  current paper. Use with caution and report alongside HBP/GAFM.

---

## 1. Pre-fine-mapping checklist

Before you call any fine-mapper, verify three things on your inputs.

### 1.1 z-score distribution and λ_GC

```python
import numpy as np
import matplotlib.pyplot as plt

z = sumstats["BETA"] / sumstats["SE"]
chi2 = z ** 2
lambda_gc = np.median(chi2) / 0.4549  # 0.4549 = chi-square median at 1 df

print(f"  λ_GC = {lambda_gc:.3f}")
print(f"  |z| > 5 (suggestive): {(np.abs(z) > 5).sum()}")
print(f"  |z| > 5.45 (genome-wide, p<5e-8): {(np.abs(z) > 5.45).sum()}")
```

| λ_GC range | Interpretation | Action |
|---|---|---|
| 0.90–1.05 | Well-calibrated | Proceed |
| 1.05–1.20 | Mild residual stratification | Consider deflating: `--lambda-gc <value>` |
| 1.20–2.00 | Strong stratification (XI/GJ-style) | Always deflate; consider re-running GWAS with more PCs / GRM |
| > 2.00 | Suspect; sumstats may be misformatted | Check pipeline; do not trust fine-mapping output |

If `--lambda-gc` is supplied, GraphGWAS deflates `z → z/√λ_GC` before
fine-mapping. This is standard genomic-control correction; it always preserves
rank but shrinks PIPs toward the centre of the distribution.

### 1.2 LD matrix sanity

```python
import numpy as np

R_sq = np.load("ld.npz")["R_sq"]

assert R_sq.shape[0] == R_sq.shape[1] == n_variants
assert (np.abs(np.diag(R_sq) - 1.0) < 1e-6).all(), "diagonal must be 1"
assert (R_sq >= 0).all() and (R_sq <= 1).all(), "r² ∈ [0, 1]"

# Symmetry
assert np.allclose(R_sq, R_sq.T, atol=1e-6), "R_sq must be symmetric"

# Local LD structure (nearby variants should be more correlated than distant ones)
n = R_sq.shape[0]
near = np.median(np.diag(R_sq, k=1))   # immediate neighbours
far = np.median(R_sq[0, n//2:])        # half-locus apart
print(f"  median r² between neighbours: {near:.3f}")
print(f"  median r² half-locus apart:   {far:.3f}")
# A reasonable locus has near > far. If they're similar, your LD matrix is suspect.
```

Common pitfalls that produce a bad LD matrix:

- **Mismatched samples** between sumstats and LD reference (use the *same*
  cohort for both, or an ancestry-matched reference panel).
- **Mismatched alleles**: if your sumstats encode minor allele but your LD
  matrix encodes ALT, signs of off-diagonal r² will be wrong. GraphGWAS
  squares correlations internally so signs cancel, but this hides the
  underlying mismatch — check upstream.
- **Singular LD** (rank < n) from over-imputed regions: drop monomorphic /
  near-monomorphic variants (MAC ≥ 5 typical).

### 1.3 Variant alignment

Both the variants TSV and the LD matrix must use the **same row order**. The
sumstats-only entry path indexes variants by `variant_id` (format
`chr:pos:ref:alt`), but row order in `R_sq` must match `variants.tsv` exactly.
Mis-ordering is silent and produces nonsense PIPs (this is exactly the bug
that affected v0.1.4 mixture-prior reweighting; see `tasks/lessons.md`
LESSON 022).

---

## 2. Input formats

### 2.1 Sumstats-only entry (`graphgwas finemap-sumstats`)

#### `variants.tsv` (required)

Tab-separated, one row per variant in the locus window. Columns:

| Column | Required | Description |
|---|---|---|
| `variant_id` | yes | `chr:pos:ref:alt` (case-sensitive) |
| `chr` | yes | chromosome (e.g. `16`, `chr16`, `Chr3`) |
| `pos` | yes | base-pair position (1-based) |
| `ref` | recommended | reference allele |
| `alt` | recommended | alternate allele |
| `af` | recommended | alternate allele frequency in the GWAS cohort |

```
variant_id	chr	pos	ref	alt	af
16:53800954:A:T	16	53800954	A	T	0.412
16:53801025:G:C	16	53801025	G	C	0.087
...
```

#### `z.tsv` (required)

Per-variant z-scores. **z = BETA / SE**, signed. Sign convention must match
the reference allele convention in `variants.tsv`.

| Column | Required | Description |
|---|---|---|
| `variant_id` | yes | matching `variant_id` in `variants.tsv` |
| `z` | yes | BETA / SE, signed |

#### `ld.npz` (required)

NumPy `.npz` archive with key `R_sq` = signed r² matrix
(GraphGWAS internally clips and squares to handle sign conventions). Row
order must match `variants.tsv`.

```python
import numpy as np
np.savez_compressed("ld.npz", R_sq=R_sq)  # n × n float64
```

### 2.2 Neo4j entry (`graphgwas finemap --source neo4j`)

The graph database must already contain `Variant` nodes with `gt_packed`
genotypes (2-bit encoded over all samples) and `Sample` nodes with phenotype
values. See [`docs/REPRODUCIBILITY.md`](REPRODUCIBILITY.md) §1–§3 for the
loading pipeline. No TSV / NPZ files needed.

### 2.3 BGEN entry (`graphgwas finemap --source bgen`)

Point `--bgen-dir` at a directory containing per-chromosome `.bgen` and
`.bgi` files. Phenotypes loaded from the graph database. Use this for
biobank-scale or yeast 1011 panels.

### 2.4 Pan-UKB entry (`graphgwas finemap --source panukb`)

No local files needed; per-locus sumstats and in-sample LD are streamed from
the public Amazon S3 bucket via tabix over HTTPS. Specify
`--phenotype <phenocode>` (e.g. `21001` for BMI) and
`--ancestries EUR,CSA,AFR,EAS`.

---

## 3. Choosing the method

| Goal | Method | Why |
|---|---|---|
| Calibrated PIPs (publish posterior probabilities) | **HBP** or **GAFM** | Diagonal calibration on F1 simulations (Sup Fig S4 a) |
| Sharpest single-variant credible set | **GAFM-MX**, **HBP-MX** or **ENS** | 2–3× higher confidence at the same rank parity (Sup Table S11) |
| Pan-UKB / biobank-scale ancestry comparison | **HBP** | Has been validated across four ancestries (Fig 5) |
| Tissue-specific eQTL prior, weak signal | **GAFM** with α=0.5 | Wins 27–2 vs SuSiE on this regime (Fig 3) |
| Two opinions to cross-check | **ENS** | Mean-of-PIPs of GAFM-MX and HBP-MX; most stable |
| Don't know yet | **HBP** as the default | 0.02–0.08 s/locus, calibrated, runs anywhere |

A safe pattern is to **report both HBP and ENS** — HBP gives the calibrated
posterior; ENS gives the operational top-1 ranking. If they agree on the lead
variant, confidence is high; if they disagree, look at the credible sets and
re-check inputs.

---

## 4. Choosing hyperparameters

| Option | Default | When to override |
|---|---|---|
| `--alpha` | 0.7 | 0.9 for strong-signal panels (single dominant locus); 0.5 for weak-signal eQTL-restricted simulations (the 27–2 GAFM-vs-SuSiE result) |
| `--coverage` | 0.95 | 0.99 if you need wider safety margin; 0.50 for "median variant" reporting |
| `--lambda-gc` | none | Set to the trait's actual λ_GC if > 1.05 (see §1.1) |
| `--n-rounds` (HBP) | 5 | 3 for fast preview; 10 if convergence not reached at PIP precision |
| `--damping` (HBP) | 0.5 | Lower (0.3) if oscillating; higher (0.7) for very dense graphs |
| `--r2-smooth` | 0.3 | 0.2 for tighter LD blocks; 0.5 for sparse LD |
| `--n-samples` | required | GWAS sample size; used to scale the mixture variance V_k = N · γ_k |

GraphGWAS-MX is **insensitive** to π and γ within reasonable perturbations
(Sup Fig S4 g): rank-1 and mean PIP are unchanged (Δ < 0.005) under
π-uniform, π-concentrated, and γ shifted ±10× from the SBayesRC defaults.
You generally do not need to tune these.

---

## 5. Output interpretation

### 5.1 Output schema

All fine-mapping commands write a TSV with one row per variant in the input,
sorted by descending `pip`:

| Column | Description |
|---|---|
| `variant_id` | matches input `variant_id` |
| `pip` | posterior inclusion probability ∈ [0, 1] |
| `in_credible_set` | boolean: variant is in the 95% credible set (= cumulative descending PIP ≤ 0.95 + the boundary variant) |
| `rank` | 1 = highest PIP, 2 = second-highest, ..., n |

Verbose outputs (Neo4j and full-pipeline modes) additionally include
`z_stat, z_functional, unique_stat, combined_score, n_ld_neighbors,
annotations`. The latter four are useful for understanding *why* GAFM
ranked a variant where it did.

### 5.2 What does a PIP value mean?

For **calibrated** methods (HBP, GAFM, SuSiE family):

| Reported PIP | Interpretation (≈ empirical TDR) |
|---|---|
| 0.95–1.00 | ~0.95 — the variant is the causal in 95% of similar simulations |
| 0.50–0.95 | the variant is one of a few highly plausible causals |
| 0.05–0.50 | low individual confidence; check whether the signal is split across LD-correlated proxies |
| ≤ 0.05 | not a credible candidate |

For **mixture-prior** methods (GAFM-MX, HBP-MX, ENS) — the same numbers are
**operational ranking scores**, NOT calibrated posteriors. The empirical
TDR in the [0.9, 1.0] PIP bin is ≈ 0.62, not 0.95 (Sup Fig S4 e). When
publishing PIPs as posterior probabilities, prefer base GAFM/HBP. When using
PIPs to *rank* candidates for follow-up (the most common case), the
mixture-prior methods deliver tighter credible sets at the same rank parity.

### 5.3 Credible-set size

A "good" credible-set size depends on signal strength, LD complexity, and
sample size:

| CS size | Interpretation |
|---|---|
| 1 (CS = 1) | Single-variant resolution. Strong signal, tight LD, large N — the ideal outcome. |
| 2–10 | Tight cluster. Typical for biobank-scale (N > 100k) at strongly associated loci. |
| 10–100 | Moderate resolution. Typical for crop/livestock panels (N < 5k) with limited LD-breaking. |
| 100–1000 | Broad LD block; you've localised to a region but not a variant. Try increasing α or narrowing the window. |
| > 1000 | Effectively unresolved — typical at suggestive loci or in panels with strong subpopulation structure. |

### 5.4 Null FPR thresholds

Empirical false-positive rates from 100 H₀ (no causal) chr22 simulations
(Sup Fig S4 f):

| Method | FPR at PIP ≥ 0.5 | FPR at PIP ≥ 0.9 |
|---|---|---|
| GAFM, HBP | 0/100 (0%) | 0/100 (0%) |
| SuSiE, SuSiE-inf, FINEMAP-inf | 1/100 (1%) | 0/100 (0%) |
| GAFM-MX, ENS | 6/100 (6%) | 1/100 (1%) |
| HBP-MX | 10/100 (10%) | 1/100 (1%) |

If you need a PIP threshold above which a positive call is "safe" under any
conceivable null:

- For **base methods** (GAFM, HBP): PIP ≥ 0.5 is already conservative.
- For **mixture-prior methods**: use PIP ≥ 0.9 to keep null FPR ≤ 1%.

---

## 6. Post-fine-mapping diagnostics

### 6.1 Top PIP < 0.5 — should you worry?

Causes, in roughly decreasing prevalence:

1. **Weak signal** (the GWAS lead p-value is at suggestive threshold). Fine-mapping cannot create resolution; expect CS to span the LD block.
2. **LD-correlated proxies** at similar effect size. Look at `n_ld_neighbors` in the verbose output; if the top variant has 50 neighbours at r² > 0.3, the signal is genuinely shared.
3. **Polygenic locus**: multiple causal variants with comparable effects. Try `--n-effects 2` or higher; if SuSiE finds multiple credible sets, GAFM/HBP are correctly assigning low PIPs.
4. **Bad LD matrix** (§1.2). If the matrix doesn't reflect your GWAS cohort's LD, top PIP will be low and unstable across methods.
5. **λ_GC inflation** (§1.1). Re-run with `--lambda-gc <value>`.

### 6.2 Disagreement between methods

If HBP and GAFM produce different rank-1 variants, that is informative:

| Pattern | Interpretation |
|---|---|
| HBP rank-1 = GAFM rank-1, both PIP > 0.5 | Strong consensus. Trust the variant. |
| HBP and GAFM rank-1 differ but are LD-correlated (r² > 0.5) | Same effective signal; report both. |
| HBP and GAFM rank-1 differ and are LD-independent | Two signals at the locus; consider a multi-effect fine-mapper (SuSiE) or `--n-effects 2`. |
| One method gives PIP > 0.9, the other < 0.1 at the same variant | One of them is reading bad inputs. Run the §1 checks. |

### 6.3 Unexpectedly large credible set

If your CS spans 500+ variants:

1. Check λ_GC (§1.1). A λ_GC of 1.4 on a 5000-variant locus is enough to inflate the CS by an order of magnitude.
2. Check LD matrix sanity (§1.2). If neighbour r² is similar to far-locus r², your matrix is likely from the wrong reference panel.
3. Reduce the window: `--window 25000` instead of 100000 strips far-LD noise.
4. Increase α (for GAFM): `--alpha 0.9` weights the statistical evidence higher; useful for strong signals at single causal loci.

### 6.4 Cross-ancestry consistency (Pan-UKB mode)

GraphGWAS reports per-ancestry PIPs side-by-side. Expect:

- Lead variant **conserved** across ancestries: strong cross-ancestry replication, the variant is plausibly causal at the population level.
- Lead variant **different** across ancestries with high PIP in each: ancestry-specific causal variants (real biology).
- Lead variant **different** with low PIP in non-EUR: smaller-N ancestries simply have less power; not necessarily different biology.
- All ancestries report **CS spanning > 500 variants**: LD reference is wrong for at least one of the non-EUR ancestries; check that the BGEN sample list matches the sumstats ancestry filter.

---

## 7. FAQ

**Q: Which mixture-prior method should I use?**
A: ENS as the default — it averages GAFM-MX and HBP-MX, so it inherits the
strengths of both and is the most stable. If you must pick one base, use
HBP-MX (slightly better mean PIP at causal in our cross-species tests).

**Q: My credible set is the entire locus (>1000 variants). What do I do?**
A: In order: re-check λ_GC, re-check the LD matrix, narrow the window,
increase α. If all four are fine and the CS is still huge, your locus
genuinely has no fine-mapping resolution at this sample size; report the
lead variant + CS and move on.

**Q: Should I always use `--lambda-gc`?**
A: No. If λ_GC ∈ [0.90, 1.05], do nothing. If λ_GC > 1.05, deflate. If
λ_GC > 1.4, deflation only partially fixes things — also re-run the GWAS
with more PCs or with GRAMMAR+ to reduce inflation at source.

**Q: Can I trust a PIP of 0.9?**
A: From base GAFM, base HBP, or any SuSiE-family method: yes (calibrated
TDR ≈ 0.95). From GAFM-MX, HBP-MX, or ENS: cautious — the empirical TDR
in the [0.9, 1.0] bin is ≈ 0.62. Use mixture-prior PIPs as ranks, not as
posterior probabilities.

**Q: What's the expected runtime per locus?**
A: HBP, GAFM and their mixture variants: 0.02–0.08 s. ENS: 0.05–0.10 s.
SuSiE (R subprocess): 5–60 s. SuSiE-inf, FINEMAP-inf: 1–10 s. SBayesRC
runs genome-wide and takes ≈ 10 min per trait.

**Q: How do I integrate fine-mapping output back into the multi-omics graph?**
A: Each credible-set variant becomes an `AssociationResult` node connected
to its `Variant` node by `FOR_VARIANT` and to a `GWASStudy` node by
`IN_STUDY`. From there you can Cypher-query: which genes does this variant
regulate? what tissues? which pathways are enriched in the credible set?
See [`vignettes/full-1kg-pipeline.md`](../vignettes/full-1kg-pipeline.md) §10
for an end-to-end example.

---

## See also

- [Command reference: `finemap`](manual/commands/finemap.md) — Neo4j / BGEN / Pan-UKB
- [Command reference: `finemap-sumstats`](manual/commands/finemap-sumstats.md) — sumstats-only
- [Pan-UKB quickstart vignette](../vignettes/fine-mapping-quickstart.md)
- [Full 1KG pipeline vignette](../vignettes/full-1kg-pipeline.md)
- [Reproducibility guide](REPRODUCIBILITY.md) — regenerate every paper figure
- [Mathematical proofs](MATHEMATICAL_PROOFS.md) — Theorems 1–5 (HBP contraction, GAFM ranking)
- [`tasks/lessons.md`](../tasks/lessons.md) — accumulated do/don't lessons
