# Fine-mapping quickstart — from Pan-UKB summary statistics to a credible set

This vignette walks through an end-to-end fine-mapping analysis using
GraphGWAS on **real Pan-UK Biobank summary statistics**, with no Neo4j
installation and no local data required. At the end, you'll have a
cross-ancestry 95% credible set for the BMI/FTO locus replicating the
main result of the accompanying Nature Genetics paper.

**Est. time**: 15 min · **Data**: streamed from public Pan-UKB S3 bucket ·
**Prereqs**: Python 3.11+, `tabix` (htslib ≥ 1.19), internet access.

## 1. Install

```bash
git clone https://github.com/jfmao/GraphGWAS.git
cd GraphGWAS/src/python
pip install -e .                # fine-mapping only, no optional extras
```

Verify:
```bash
python -c "from graphgwas.panukb import fetch_sumstats_locus; print('OK')"
```

## 2. Fetch Pan-UKB sumstats for BMI at FTO

FTO (chr16:53,820,527 in GRCh37) is the canonical BMI locus. Pan-UKB
ships GWAS summary statistics for 7,224 phenotypes across six
ancestries; we pull the BMI sumstats (phenocode `21001`, continuous,
IRNT-transformed) for a ±100 kb window around the FTO lead:

```python
from graphgwas.panukb import fetch_sumstats_locus, RECOMMENDED_ANCESTRIES

sumstats = fetch_sumstats_locus(
    phenocode="21001",          # BMI
    chr="16", start=53720000, end=53920000,
    trait_type="continuous",
    modifier="irnt",
    ancestries=RECOMMENDED_ANCESTRIES,   # EUR, CSA, AFR, EAS
)

for anc, loc in sumstats.items():
    print(f"{anc}: N={loc.n_samples:>8,}  variants={len(loc.variants):>5}")
```

Expected output:
```
EUR: N= 420,531  variants= 1,574
CSA: N=   8,876  variants= 1,199
AFR: N=   6,636  variants= 1,483
EAS: N=   2,709  variants=   812
```

Three points to note:
1. Fetch is via `tabix` over HTTPS against `pan-ukb-us-east-1` — no bulk download, no authentication, no Neo4j.
2. Pan-UKB variants are reported in GRCh37.
3. The four ancestries are those with Pan-UKB `N > 2,500`, where HBP
   fine-mapping is most reliable (MID and AMR are available but too
   small-N for tight credible sets).

## 3. Align to GRCh38 + compute ancestry-matched LD

GraphGWAS's reference LD is a subset of 1000 Genomes high-coverage BGEN
(GRCh38). Two steps: liftover the Pan-UKB positions to GRCh38, and
subset 1KG to ancestry-matched samples for the LD matrix.

The benchmark script `tests/benchmark_panukb_finemap.py` encapsulates
both steps. For the quickstart, just run it for the FTO locus:

```bash
python tests/benchmark_panukb_finemap.py --locus FTO \
    --out results/quickstart
```

The script writes `results/quickstart/FTO.json` with per-ancestry HBP
credible sets. Expected completion: ~10 s cold, ~5 s warm.

## 4. Inspect the credible sets

```python
import json
with open("results/quickstart/FTO.json") as f:
    rec = json.load(f)

print(f"{'Ancestry':9} {'N':>9}  {'intersect':>9}  {'top variant':30} {'PIP':>6} {'CS size':>8}")
for anc, r in rec["ancestries"].items():
    print(f"{anc:9} {r['N_sumstats']:>9,}  {r['n_intersect']:>9}  "
          f"{r['top_hbp_variant']:30} {r['top_hbp_pip']:6.3f} {r['n_cs_hbp']:>8}")
```

Expected output:
```
Ancestry         N   intersect  top variant                       PIP  CS size
EUR        420,531        742  16:53805443:GT:G                1.000        1
CSA          8,876        762  16:53758552:T:C                 0.028      628
AFR          6,636       1081  16:53857466:C:A                 0.012     1012
EAS          2,709        688  16:53861892:A:G                 0.011      633
```

Interpretation:

- **EUR** resolves the FTO signal to a single variant at PIP = 1.000 —
  this is the paper's power-at-$N$ result replicated on real biobank
  data.
- **CSA, AFR, EAS** expand their credible sets to hundreds of variants
  — the expected behaviour at small $N$ and without informative
  tissue-specific eQTL priors for those ancestries.

For contrast, run the same quickstart on **Height at HMGA2**, which
resolves to a single-variant credible set in EUR *and* CSA *and* EAS
(see the `HMGA2` entry in the paper's Figure 7):

```bash
python tests/benchmark_panukb_finemap.py --locus HMGA2 \
    --out results/quickstart
```

Expected: EUR PIP = 1.000, CSA PIP = 1.000, EAS PIP = 0.974 (all
CS-size 1).

## 5. Use HBP programmatically

For custom analyses, skip the benchmark script and call the
sumstats-only entry path directly:

```python
import numpy as np
from graphgwas.panukb import fetch_sumstats_locus
from graphgwas.finemapping_v2 import hbp_finemap_from_sumstats
from tests.benchmark_panukb_finemap import (
    _compute_1kg_ld_matched, _intersect_sumstats_ld, _lift_sumstats_to_b38,
)

# Step A: fetch Pan-UKB sumstats
ss = fetch_sumstats_locus(
    phenocode="21001", chr="16",
    start=53720000, end=53920000,
    trait_type="continuous", modifier="irnt",
    ancestries=["EUR"],
)
ss_b38 = {anc: _lift_sumstats_to_b38(s) for anc, s in ss.items()}

# Step B: ancestry-matched 1KG LD
ld_vars, R_sq = _compute_1kg_ld_matched(
    chr="16", start=53690000, end=53890000,   # GRCh38-ish window
    ancestry="EUR",
    bgen_dir="tests/data/human/1kGP_bgen",
    popmap_path="tests/data/human/1kg_popmap.tsv",
)

# Step C: intersect + fine-map
variants, z, R = _intersect_sumstats_ld(ss_b38["EUR"], ld_vars, R_sq)
candidates = hbp_finemap_from_sumstats(
    variants, z, R, graph_cache={},
    chr_name="chr16",
    n_rounds=5, alpha=0.6, damping=0.5,
    credible_set_coverage=0.95,
)

# Step D: print credible set
print(f"{'variant_id':32} {'PIP':>6}  in_CS")
for c in candidates[:10]:
    print(f"{c.variant_id:32} {c.pip:6.3f}  {c.in_credible_set}")
```

The `graph_cache` parameter above is passed empty (`{}`) because we are
running from Pan-UKB sumstats without the multi-omics annotation graph
loaded. To add relational priors, load the pre-built Neo4j dump
(Zenodo DOIs in `CITATION.cff`) and populate `graph_cache` from
`graphgwas.finemapping_v2._build_graph_cache()`.

## 6. Next steps

- **Try other loci**: `--locus FTO | APOA5 | LDLR | HMGA2` (see
  `LOCI` in `tests/benchmark_panukb_finemap.py`)
- **Try GAFM instead of HBP** with `l1_finemap_from_sumstats` (Python
  prefix `l1_` is the historical name; paper-facing name is GAFM) —
  especially interesting when your causal variant is a tissue-specific
  eQTL, where GAFM beats SuSiE 27–2 (see paper §2.3)
- **Load the full Neo4j multi-omics graph** (17 GB Zenodo dump) and
  use the Cypher query layer for interpretability ("which genes in
  this credible set are drug targets?")
- **Expose the fine-mapping methods to an AI agent** via
  `graphgwas mcp` and query them from any MCP-compatible client
  in natural language
- **Reproduce every figure in the paper** via
  `python tests/generate_paper_figures.py`; see
  [`docs/REPRODUCIBILITY.md`](../docs/REPRODUCIBILITY.md)

## Further reading

- Paper: *Relational biological structure improves fine-mapping of
  causal GWAS variants under weak signal* (submitted, 2026)
- Manual: [`docs/manual/index.md`](../docs/manual/index.md)
- Install: [`docs/INSTALL.md`](../docs/INSTALL.md)
- Math: [`docs/MATHEMATICAL_PROOFS.md`](../docs/MATHEMATICAL_PROOFS.md)
- Platform scope (benchmark-status table): Supplementary Note S3
