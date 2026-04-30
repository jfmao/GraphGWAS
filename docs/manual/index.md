# GraphGWAS CLI — User Manual

**Version 0.1.0** · 40+ commands across 15 functional groups

The `graphgwas` command-line interface is the primary entry point for
interactive analysis and scripted pipelines. All commands share a common
Neo4j connection (or a BGEN / Pan-UKB sumstats fallback where noted) and
write their outputs either to TSV/CSV/JSON on disk or as graph nodes
queryable via Cypher.

> Programmatic access: the same procedures are exposed as a 37-endpoint
> FastAPI REST server (`graphgwas serve`) and a 16-tool MCP server
> (`graphgwas mcp`) for AI-agent use. See [`docs/INSTALL.md`](../INSTALL.md)
> for setup.

## Global options

All commands accept the following connection options, also readable from
environment variables:

| Option | Env var | Default |
|---|---|---|
| `--uri` | `GRAPHGWAS_URI` | `bolt://localhost:7688` |
| `--user` | `GRAPHGWAS_USER` | `neo4j` |
| `--password` | `GRAPHGWAS_PASSWORD` | `graphgwas` |
| `--database` | `GRAPHGWAS_DATABASE` | `neo4j` |

## Table of contents

### Core fine-mapping (covered in the Nature Genetics manuscript)
- [`graphgwas finemap`](commands/finemap.md) — run GAFM, HBP, CLGF or L4 fine-mapping on a locus
- [`graphgwas assoc mpat`](commands/assoc-mpat.md) — gene-level MPAT test
- [`graphgwas assoc top-hits`](commands/assoc-top-hits.md) — top signals from a GWAS run

### Phenotype & data
- [`graphgwas phenotype activate`](commands/phenotype-activate.md) — activate a trait for GWAS

### Population structure
- [`graphgwas popstruct spectral-pcs`](commands/popstruct-spectral-pcs.md) — graph-spectral PCs from rare-variant similarity
- [`graphgwas spectral correction`](commands/spectral-correction.md) — spectral phenotype correction

### Graph-native multi-locus analysis (preview in this paper; paper #2 in preparation)
- [`graphgwas epistasis scan`](commands/epistasis-scan.md) — LD-pruned co-occurrence epistasis (M1–M4)
- [`graphgwas flow architecture`](commands/flow-architecture.md) — max-flow disease architecture

### Heritability (platform methods; not benchmarked in this paper)
- [`graphgwas heritability spectral`](commands/heritability-spectral.md) — Laplacian-eigenspectrum $h^2$
- [`graphgwas heritability conductance`](commands/heritability-conductance.md) — conductance-based $h^2$
- [`graphgwas heritability flow`](commands/heritability-flow.md) — pathway-decomposed $h^2$
- [`graphgwas heritability gnn`](commands/heritability-gnn.md) — GNN-predicted $h^2$
- [`graphgwas heritability multi`](commands/heritability-multi.md) — multi-resolution $h^2$
- [`graphgwas heritability report`](commands/heritability-report.md) — unified $h^2$ report

### Multivariate / cross-trait (platform methods)
- [`graphgwas multivariate correlation`](commands/multivariate-correlation.md) — genetic correlation $r_G$
- [`graphgwas multivariate g-matrix`](commands/multivariate-g-matrix.md) — $T \times T$ G-matrix
- [`graphgwas multivariate multi-resolution`](commands/multivariate-multi-resolution.md) — variant / gene / pathway scale correlation
- [`graphgwas multivariate coherence`](commands/multivariate-coherence.md) — frequency-resolved coherence
- [`graphgwas multivariate pleiotropy`](commands/multivariate-pleiotropy.md) — pleiotropic gene discovery
- [`graphgwas multivariate psi`](commands/multivariate-psi.md) — Gene Pleiotropy Strength Index
- [`graphgwas multivariate report`](commands/multivariate-report.md) — full multivariate report

### Gene × Environment / Multi-Environment Trials (platform methods)
- [`graphgwas met status`](commands/met-status.md) — trial summary
- [`graphgwas met env-similarity`](commands/met-env-similarity.md) — genetic correlation between environments
- [`graphgwas met mega-env`](commands/met-mega-env.md) — mega-environment discovery
- [`graphgwas met gxe`](commands/met-gxe.md) — G×E variance decomposition
- [`graphgwas met diagnostics`](commands/met-diagnostics.md) — experimental-design graph diagnostics
- [`graphgwas met reaction-norm`](commands/met-reaction-norm.md) — variant-level reaction norms
- [`graphgwas met impute`](commands/met-impute.md) — graph-neighbour missing-cell imputation
- [`graphgwas met report`](commands/met-report.md) — comprehensive MET report

### Polygenic risk scores (platform methods)
- [`graphgwas prs classical`](commands/prs-classical.md) — standard weighted-sum PRS
- [`graphgwas prs pathway`](commands/prs-pathway.md) — pathway-partitioned PRS
- [`graphgwas prs report`](commands/prs-report.md) — compare all PRS methods

### Mendelian randomisation (platform methods)
- [`graphgwas mr report`](commands/mr-report.md) — IVW + Egger + weighted median + pleiotropy

### GNN (platform methods, partial)
- [`graphgwas gnn export`](commands/gnn-export.md) — export subgraph to PyTorch Geometric

### Visualisation
- [`graphgwas plot manhattan`](commands/plot-manhattan.md) — Manhattan plot
- [`graphgwas plot qq`](commands/plot-qq.md) — QQ plot with $\lambda_{GC}$

### Results management
- [`graphgwas results export`](commands/results-export.md) — TSV export
- [`graphgwas results delete`](commands/results-delete.md) — delete run results

### Interactive / programmatic interfaces
- [`graphgwas agent`](commands/agent.md) — LangGraph natural-language interface
- [`graphgwas interpret`](commands/interpret.md) — deterministic post-hoc interpretation
- [`graphgwas serve`](commands/serve.md) — start FastAPI REST server
- [`graphgwas mcp`](commands/mcp.md) — start MCP server (for Claude Desktop / Claude Code)

## Command reference template

Every command reference page follows this structure (modelled on GraphPop's
per-command docs):

```markdown
# <short command title>

## Description
One paragraph explaining what the command does and when to use it.

## Usage
graphgwas <group> <command> <ARGS> [OPTIONS]

## Arguments
| Name | Type | Required | Default | Description |
|------|------|----------|---------|-------------|
...

## Options
| Option | Type | Default | Description |
|--------|------|---------|-------------|
...

## Output
What the command returns (schema of TSV / JSON, graph nodes written, etc.)

## Examples
2–3 worked examples with expected output.

## See also
Related commands and vignettes.
```

## Benchmark-status legend

Not every command listed above is rigorously benchmarked in the Nature
Genetics manuscript. See Supplementary Note S3 of the paper
(`paper/manuscript_v1/supplementary.tex`) for the full benchmark-status table.
In short:

- **Fine-mapping** commands (`finemap`, and the underlying HBP / GAFM / CLGF
  methods): **rigorously benchmarked** (90 F1 simulations, 79 weak-signal
  replicates, 100 null replicates, Pan-UKB 4 ancestries).
- **Epistasis** (`epistasis scan`): **previewed**; full benchmark in paper #2.
- **All other platform methods** (heritability, multivariate, PRS, MR, MET,
  GNN, agent): **implemented** but awaiting dedicated benchmarks.
