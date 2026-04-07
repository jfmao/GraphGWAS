# GraphGWAS — Mac Mini M1 Development Setup

Migrating from Ubuntu workstation (64GB + RTX 4090) to Mac Mini M1 (16GB).
This document covers everything needed to get a working dev environment.

---

## 1. Prerequisites to Install

### Java 21 (required for Neo4j 5.x)
```bash
brew install openjdk@21
sudo ln -sfn $(brew --prefix openjdk@21)/libexec/openjdk.jdk /Library/Java/JavaVirtualMachines/openjdk-21.jdk
java -version  # should show 21.x
```

### Neo4j 5.x Community Edition
```bash
brew install neo4j
# Or download from https://neo4j.com/download-center/
```

Configure memory for 16GB machine — edit `neo4j.conf`:
```properties
# Location: /opt/homebrew/etc/neo4j/neo4j.conf (Homebrew)
# Or: ~/Library/Application Support/Neo4j/neo4j.conf

server.memory.heap.initial_size=4g
server.memory.heap.max_size=4g
server.memory.pagecache.size=2g
# Leave ~10GB for OS + Python + dev tools
```

Set the password to match config.py default:
```bash
neo4j-admin dbms set-initial-password graphpop
```

Start Neo4j:
```bash
neo4j start
# Verify: neo4j status
# Console: http://localhost:7474
```

### Python 3.11+
```bash
brew install python@3.13
# Or use pyenv:
# brew install pyenv && pyenv install 3.13.12
```

---

## 2. Clone & Install the Project

```bash
# Clone all three sibling projects
cd /path/to/your/workspace
git clone <repo-url> GraphGWAS
git clone <repo-url> GraphMana
git clone <repo-url> GraphPop

# Install GraphGWAS in editable mode
cd GraphGWAS
pip install -e src/python/

# Optional deps (install what you need)
pip install -e "src/python/[dev]"       # pytest, ruff
pip install -e "src/python/[api]"       # FastAPI server
pip install -e "src/python/[mcp]"       # MCP server
pip install -e "src/python/[agent]"     # LangGraph agent (needs ANTHROPIC_API_KEY)

# GNN — see Section 5 for Mac-specific notes
pip install -e "src/python/[gnn]"
```

---

## 3. Critical: N_SAMPLES Constant

**`src/python/graphgwas/config.py` line 29 has `N_SAMPLES = 3202` hardcoded.**

This is the 1000 Genomes (human) sample count. It is used by 15+ modules for
`gt_packed` byte array decoding — the packed genotype format encodes 4 genotypes
per byte, so `N_SAMPLES` determines the expected byte length.

**When working with yeast data, you MUST update this to match your dataset:**
- 1011 Yeast Genomes: `N_SAMPLES = 1011`
- Your actual sample count from Neo4j: `MATCH (s:Sample) RETURN count(s)`

If N_SAMPLES is wrong, every genotype decode will produce garbage.

> **TODO for the project:** This should be auto-detected from the database
> rather than hardcoded. For now, change it manually when switching datasets.

---

## 4. Data Strategy for Mac Mini

### What NOT to transfer (51GB+ of test data on Ubuntu)
| Dataset | Size | Transfer? |
|---------|------|-----------|
| `tests/data/human/` | 33 GB | NO — too large, Ubuntu-only |
| `tests/data/arabidopsis/` | 18 GB | NO — too large, Ubuntu-only |
| `tests/data/yeast/` | 504 MB | YES — primary dev dataset |
| `tests/data/sim/` | 7.3 MB | YES — simulation tests |
| `tests/data/*.csv` | 160 KB | YES — phenotype test files |

### What you need on the Mac
1. **Yeast data** (`tests/data/yeast/`, 504 MB) — copy from Ubuntu or re-download
2. **Simulation data** (`tests/data/sim/`, 7.3 MB) — copy or regenerate with `tests/generate_simulation.py`
3. **Phenotype CSVs** (`tests/data/*.csv`) — copy from Ubuntu

### Neo4j Database
You have two options:

**Option A: Export/import Neo4j dump (recommended for yeast-only)**
On Ubuntu:
```bash
neo4j stop
neo4j-admin database dump neo4j --to-path=/tmp/
# Copy the dump file to Mac
```
On Mac:
```bash
neo4j stop
neo4j-admin database load neo4j --from-path=/path/to/dump/
neo4j start
```

**Option B: Re-import from VCF using GraphMana**
If you want a clean yeast-only database on the Mac, use GraphMana's import
pipeline to load yeast VCF + phenotypes fresh. This is the cleanest approach
but requires GraphMana to be set up first.

> **Note:** The current Ubuntu Neo4j database is 3.4 GB on disk (70.7M human
> variants). A yeast-only DB will be much smaller (~100-200 MB).

---

## 5. GNN on Apple Silicon (Mac-Specific)

### PyTorch + PyG installation
```bash
# PyTorch for Apple Silicon
pip install torch torchvision

# PyTorch Geometric — install from wheels
pip install torch-geometric
# If scatter/sparse deps fail:
pip install torch-scatter torch-sparse -f https://data.pyg.org/whl/torch-2.x.x+cpu.html
```

### MPS (Metal Performance Shaders) status
The `_resolve_device()` in `gnn.py` deliberately defaults to **CPU** even on
Apple Silicon. PyTorch Geometric's scatter/message-passing ops have incomplete
MPS support. The GNN will run on CPU, which is fine for yeast-scale data.

If you want to test MPS manually: `graphgwas gnn train --device mps`
(may crash — use at your own risk).

### Memory limits
With 16GB shared between CPU and GPU:
- Yeast (1K samples × 10K variants): fine on CPU
- If you hit OOM, reduce `--max-variants` or work on a smaller chromosome

---

## 6. Environment Variables

Add to your `~/.zshrc` or `~/.bash_profile`:

```bash
# Neo4j connection (defaults match config.py, only set if non-default)
# export GRAPHGWAS_NEO4J_URI="bolt://localhost:7688"
# export GRAPHGWAS_NEO4J_USER="neo4j"
# export GRAPHGWAS_NEO4J_PASSWORD="graphpop"

# For agent module (Phase 5)
export ANTHROPIC_API_KEY="sk-ant-..."

# Python path for running tests directly
export PYTHONPATH="/path/to/GraphGWAS/src/python:$PYTHONPATH"
```

---

## 7. Verify Installation

```bash
# 1. Check Neo4j connectivity
graphgwas schema audit

# 2. Run offline simulation tests (no Neo4j needed)
PYTHONPATH=src/python python -m pytest tests/test_statistics.py -v
PYTHONPATH=src/python python -m pytest tests/test_simulation_gwas.py -v

# 3. Run yeast offline tests (no Neo4j needed)
PYTHONPATH=src/python python -m pytest tests/test_yeast_offline.py -v

# 4. CLI smoke test
graphgwas --help
```

---

## 8. What Works / What Doesn't on Mac Mini

| Capability | Mac Mini M1 | Notes |
|-----------|-------------|-------|
| All 24 non-GNN modules | Yes | CPU-bound NumPy/SciPy |
| CLI (58+ commands) | Yes | |
| FastAPI server | Yes | `graphgwas serve` |
| MCP server | Yes | `graphgwas mcp` |
| Offline tests | Yes | No Neo4j needed |
| Neo4j + yeast GWAS | Yes | 4GB heap sufficient |
| GNN training (CPU) | Yes | Yeast-scale only |
| GNN training (GPU) | No | MPS unreliable with PyG |
| Full human 1KG (70.7M var) | No | Need 32GB+ Neo4j heap |
| Genome-wide epistasis | Tight | Limit to small windows |

---

## 9. Differences from Ubuntu Workstation

| | Ubuntu (64GB + 4090) | Mac Mini M1 (16GB) |
|---|---|---|
| Neo4j heap | 16 GB | 4 GB |
| N_SAMPLES | 3202 (human) | 1011 (yeast) or per-dataset |
| GNN device | CUDA | CPU |
| Full human data | Yes | No |
| Neo4j DB size | 3.4 GB (70.7M variants) | ~200 MB (yeast) |
| Parallel workers | 16+ cores | 8 cores (4P+4E) |
