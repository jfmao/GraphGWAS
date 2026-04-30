# GraphGWAS installation

## Prerequisites

| Component | Version | Notes |
|---|---|---|
| Python | ≥ 3.11 | 3.13 tested |
| Neo4j | 5.26 Community | Optional for sumstats-only workflows |
| Java | 21 | Required by Neo4j |
| htslib / tabix | ≥ 1.19 | For Pan-UKB sumstats streaming |
| R | 4.x | Optional — for SuSiE / FINEMAP baseline comparisons |

External binaries (optional, for baseline comparisons in benchmarks):
`plink2 ≥ 2.0`, `FINEMAP 1.4.2`, `bcftools ≥ 1.19`, `GCTB 2.5`.

## Python package

```bash
git clone https://github.com/jfmao/GraphGWAS.git
cd GraphGWAS/src/python

# Minimal install (fine-mapping + CLI + API + MCP)
pip install -e .

# With optional extras
pip install -e '.[all]'        # GNN + agent + API + MCP + dev
pip install -e '.[gnn]'        # PyTorch Geometric for HeteroGNN
pip install -e '.[agent]'      # LangGraph AI agent
pip install -e '.[mcp]'        # MCP server only
```

After install, the `graphgwas` CLI is on `PATH`:

```bash
graphgwas --help
graphgwas finemap --help
```

## Neo4j setup (optional for sumstats-only workflows)

GraphGWAS supports three data-access modes; **only mode (c) requires Neo4j**.

- **(a) Pan-UKB sumstats** — streamed via tabix from the public S3 bucket
  `pan-ukb-us-east-1`; no Neo4j, no local data, no authentication required
- **(b) Local BGEN** — ships with `bgen` Python package; no Neo4j required
- **(c) Full multi-omics graph** — requires Neo4j 5.26 with pre-built dumps

For mode (c), start Neo4j and restore one of the Zenodo dumps:

```bash
# Download + start Neo4j 5.26
wget https://dist.neo4j.org/neo4j-community-5.26.0-unix.tar.gz
tar -xzf neo4j-community-5.26.0-unix.tar.gz
cd neo4j-community-5.26.0

# Configure heap / page-cache for the 17 GB human dump
echo "server.memory.heap.initial_size=16g"     >> conf/neo4j.conf
echo "server.memory.heap.max_size=16g"          >> conf/neo4j.conf
echo "server.memory.pagecache.size=16g"         >> conf/neo4j.conf
echo "server.default_listen_address=0.0.0.0"    >> conf/neo4j.conf
echo "server.bolt.listen_address=:7688"         >> conf/neo4j.conf  # avoid default port collision

# Restore the Zenodo dump (DOI in CITATION.cff)
bin/neo4j-admin database load --from-path=/path/to/dump.dir neo4j --overwrite-destination

# Start
bin/neo4j start
```

Then point GraphGWAS at it:

```bash
export GRAPHGWAS_URI=bolt://localhost:7688
export GRAPHGWAS_USER=neo4j
export GRAPHGWAS_PASSWORD=graphgwas
graphgwas --help
```

## Hail (optional — for Pan-UKB in-sample LD BlockMatrices)

Pan-UKB's per-ancestry LD matrices (47.6 TB across six ancestries) are shipped
as Hail `BlockMatrix` objects on S3. To read them directly, install Hail in a
dedicated Python 3.11 environment:

```bash
conda create -n hail311 python=3.11 -y
conda activate hail311
pip install hail
```

Note: the Spark `hadoop-aws` classpath must be configured for `s3a://` access
to the Pan-UKB public bucket. GraphGWAS falls back to 1KG-matched LD in the
default sumstats-only pipeline, which is sufficient for demonstration
purposes; the Pan-UKB in-sample LD path is optional polish.

## Verification

```bash
cd GraphGWAS
python -m pytest tests/test_*.py -q
```

Expected: **78 passed**, 2 skipped (offline / optional-data tests).

## Docker (optional)

A containerised setup with Neo4j + GraphGWAS + pre-built yeast dump is
planned for the v0.2 release (post-submission).

## Troubleshooting

**Neo4j runs out of memory on the 17 GB human dump.**
Increase heap to 32 GB (needs ≥ 48 GB RAM machine) or use the 0.5 GB yeast
dump instead:

```bash
export NEO4J_HEAP=32g
```

**`graphgwas finemap` returns no results for a locus.**
(1) Verify the phenotype is loaded: `graphgwas phenotype list`.
(2) Verify the chromosome + position format matches the graph (`chr16` vs
`16`); Pan-UKB sumstats use `16` (GRCh37); human 1KG BGEN in this repo uses
`chr16` (GRCh38). Use `graphgwas --help` to see build-specific flags.

**Pan-UKB tabix query is slow.**
First query per chromosome downloads the tabix index (~2 MB). Warm queries
typical ≤ 3 s. Ensure `tabix` is htslib ≥ 1.19 (`tabix --version`).

**Hail S3 read fails with `ClassNotFoundException: org.apache.hadoop.fs.s3a.S3AFileSystem`.**
This is a Spark-classpath configuration issue for the optional Pan-UKB
in-sample LD path. The default 1KG-matched-LD fallback works without Hail.
