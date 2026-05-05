# `graphgwas finemap` — fine-map causal variants at a GWAS locus

## Description

Runs GraphGWAS fine-mapping on a specified genomic window, producing a
95% credible set of candidate causal variants ranked by posterior
inclusion probability (PIP). Four methods are available:

- **`hbp`** — Hierarchical Belief Propagation on a variant → gene →
  pathway factor graph with protein-interaction coupling. Matches
  state-of-the-art Bayesian accuracy at 6–60× the speed. Banach
  contraction guarantees geometric convergence (Theorem 2).
- **`l1`** *(default)* — LD-deconvolved Bayesian fine-mapping with a
  multi-layer graph functional score, combined via adaptive α chosen by
  empirical Bayes. Wins 27–2 head-to-head against SuSiE at weak signal
  with tissue-specific eQTL priors (Theorem 3).
- **`clgf`** — Cross-Locus Graph Fine-mapping; EM over multiple loci
  that share pathway evidence.
- **`glem`** — Graph-Latent-Embedding Fine-Mapping; multi-signal
  fine-mapping via low-dimensional embedding of the variant–gene
  factor graph.

Three input sources are supported, auto-detected from the `--source`
argument:

| Source | Requires | Use case |
|---|---|---|
| `neo4j` | Neo4j graph with gt_packed | Individual-level genotypes available |
| `bgen` | BGEN directory | Biobank-scale; local 1KG; yeast |
| `panukb` | Internet + tabix | Pan-UKB summary statistics |

## Usage

```bash
graphgwas finemap --chr CHR --pos POS --window WINDOW \
                  --phenotype PHENO --method METHOD \
                  [OPTIONS] [-o OUTPUT]
```

## Arguments

| Name | Type | Required | Default | Description |
|------|------|----------|---------|-------------|
| `--chr` | string | yes | — | Chromosome (e.g. `chr16`, `16`) |
| `--pos` | int | yes | — | Lead variant position (bp) |
| `--window` | int | no | 100000 | Half-window around `--pos` (bp) |
| `--phenotype` | string | yes | — | Phenotype key (Neo4j mode) or Pan-UKB phenocode |
| `--method` | choice | no | `l1` | `l1`, `hbp`, `clgf`, `glem` |
| `--source` | choice | no | `neo4j` | `neo4j`, `bgen`, `panukb` |

## Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--alpha` | float | 0.5 | Balance: 0 = pure functional, 1 = pure statistical |
| `--r2-smooth` | float | 0.3 | LD-deconvolution threshold |
| `--credible-set` | float | 0.95 | Target coverage for credible set |
| `--n-rounds` | int | 5 | HBP message-passing rounds |
| `--damping` | float | 0.5 | HBP damping parameter λ |
| `--build` | choice | `b38` | `b37` or `b38` (Pan-UKB mode) |
| `--ancestries` | string | `EUR,CSA,AFR,EAS` | Pan-UKB ancestries |
| `-o`, `--output` | path | stdout | Output TSV |

## Output

TSV with one row per variant in the window, sorted by descending PIP:

```
variant_id         chr    pos         ref  alt  af      z_stat   z_functional  unique_stat  combined_score  pip    in_credible_set  annotations
chr16:53805443:GT:G chr16  53805443    GT   G    0.42    8.82     0.85          8.41         6.48            1.000  True             gene:FTO;pathway:energy_metabolism
chr16:53820527:A:T  chr16  53820527    A    T    0.39    8.70     0.40          8.20         6.10            0.000  False            gene:FTO
...
```

## Examples

### Example 1 — Fine-map BMI at FTO using Pan-UKB summary statistics

```bash
graphgwas finemap --source panukb \
    --chr 16 --pos 53820527 --window 100000 \
    --phenotype 21001 \
    --method hbp \
    --ancestries EUR,CSA,AFR,EAS \
    -o results/fto_bmi_hbp.tsv
```

Expected result (from the manuscript, Table/Fig 7): EUR resolves to a
single-variant 95% credible set at PIP = 1.000; CSA PIP = 0.028 with
CS size 628; AFR PIP = 0.012; EAS PIP = 0.011 (smaller-N ancestries
expand credible sets as statistically expected).

### Example 2 — GAFM weak-signal fine-mapping on Neo4j-backed yeast data

```bash
graphgwas finemap --source neo4j \
    --chr chr4 --pos 900000 --window 25000 \
    --phenotype YPETHANOL \
    --method l1 \
    --alpha 0.5 \
    -o results/yeast_ethanol_l1.tsv
```

### Example 3 — HBP on local 1KG BGEN

```bash
graphgwas finemap --source bgen \
    --chr 22 --pos 17200000 --window 50000 \
    --phenotype simulated_causal_1 \
    --method hbp \
    --n-rounds 5 --damping 0.5 \
    -o results/1kg_chr22_hbp.tsv
```

## Theoretical background

HBP convergence is guaranteed by Banach's fixed-point theorem: the update
operator is a strict ℓ₁ contraction on the probability simplex with rate
`L = λ + (1−λ)(1−α)·ρ(M)`, where `ρ(M)` is the spectral radius of the
graph propagation operator. See Theorem 2 in
[`docs/MATHEMATICAL_PROOFS.md`](../../MATHEMATICAL_PROOFS.md).

GAFM's causal-variant ranking is guaranteed under mild LD-decay assumptions:
the LD-deconvolved statistic satisfies `u_c > u_i` for all non-causal
variants `i`. See Theorem 3.

## See also

- [`docs/INPUT_OUTPUT_GUIDE.md`](../../INPUT_OUTPUT_GUIDE.md) — preparing
  inputs, choosing the method, interpreting PIPs and credible sets, post-
  fine-mapping diagnostics, FAQ
- Vignette [`vignettes/fine-mapping-quickstart.md`](../../../vignettes/fine-mapping-quickstart.md)
- Command [`graphgwas plot manhattan`](plot-manhattan.md) to visualise PIPs
- Command [`graphgwas results export`](results-export.md) for further downstream export
- Command [`graphgwas mcp`](mcp.md) to expose fine-mapping to an AI agent
