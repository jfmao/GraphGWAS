# GraphGWAS — Human 1KG Data Preparation

## Database Status

**Neo4j instance:** `bolt://localhost:7687` (user: `neo4j`, password: `graphmana`)
**Data:** 1000 Genomes Phase 3, imported via GraphMana

| Entity | Count |
|--------|-------|
| Variants | 70,691,875 |
| Samples | 2,734 (of 3,202 total) |
| Populations | 26 (1KG superpopulations + populations) |
| Chromosomes | 22 (autosomes) |
| Gene annotations | 0 (need to load) |
| Pathway annotations | 0 (need to load) |

Chromosome naming: `chr1` through `chr22`
Genotype storage: gt_packed (2-bit packed byte arrays)

## What's Needed to Run GraphGWAS

### 1. N_SAMPLES Configuration
The human database has 2,734 samples. GraphGWAS auto-detects via `detect_n_samples()`.
When connecting to port 7687 (human), specify the URI:
```bash
graphgwas --neo4j-uri bolt://localhost:7687 --neo4j-password graphmana <command>
```

### 2. Gene Annotations
No Gene nodes exist yet. Options:
- **Option A:** Download GENCODE GFF3 for GRCh38, parse gene coordinates, load into Neo4j
- **Option B:** Use VEP-annotated VCF (if available) to create Gene nodes + HAS_CONSEQUENCE edges
- **Option C:** Run VEP/SnpEff on the variant data post-hoc

For the graph-native methods (flow, pathway), Gene + Pathway nodes are required.

### 3. Phenotype Data
No phenotypes are loaded yet. For human GWAS:
- **Simulated phenotypes** (for testing): generate case/control from Alzheimer's simulation
- **Real phenotypes** (UK Biobank, etc.): would need separate data access

Existing test phenotypes in `tests/data/`:
- `phenotypes_alzheimer_sim.csv` — simulated Alzheimer's (3,202 samples, binary)
- `phenotypes_null.csv` — null phenotype (no genetic effect, for calibration)

### 4. Running GraphGWAS on Human Data

```bash
# Connect to human database
export GRAPHGWAS_NEO4J_URI=bolt://localhost:7687
export GRAPHGWAS_NEO4J_PASSWORD=graphmana

# Or use CLI flags
graphgwas --neo4j-uri bolt://localhost:7687 --neo4j-password graphmana status

# Load phenotypes
graphgwas --neo4j-uri bolt://localhost:7687 --neo4j-password graphmana \
  phenotype load --phenotype-file tests/data/phenotypes_alzheimer_sim.csv \
  --sample-id-col sample_id

# Activate and scan
graphgwas ... phenotype activate --trait case_control --case-value case --control-value control
graphgwas ... assoc scan --chr chr19 --method auto -o gwas_alzheimer_chr19.tsv
```

### 5. Scale Considerations

| Operation | Yeast (1.9M variants) | Human (70.7M variants) | Scale factor |
|-----------|----------------------|------------------------|--------------|
| Genome-wide GWAS | 67s (parallel) | ~40 min (estimated) | 37× |
| GRM computation | 800s (151K common) | ~8 hrs (5M common) | ~35× |
| Epistasis (100kb window) | 5s | ~5s (same window) | 1× |
| LD fine-mapping (40kb) | 10s | ~10s | 1× |
| Max-flow (1 chr) | 145s | ~145s | 1× |

Window-based methods (epistasis, fine-mapping, flow) scale with region size, not genome size.
Genome-wide scans scale linearly with variant count.
GRM/GRAMMAR scales quadratically with variant count — will need GPU acceleration.

### 6. Known Issue: packed_index Mismatch

The current human database has a **packed_index mismatch**:
- Sample nodes have packed_indices from the original 3,202-sample VCF (range 0-3201, with gaps)
- gt_packed byte arrays are encoded for the 2,734 imported samples only (684 bytes = 2,736 positions)
- When GraphGWAS tries to decode genotype at index 3201, it reads past the 684-byte array

**Fix options:**
1. **Re-import** with GraphMana ensuring packed_index matches gt_packed position
2. **Remap indices**: query all Sample packed_indices, build a mapping to contiguous 0..N-1
3. **Re-encode gt_packed** for the full 3,202 positions (pad missing samples with 0x03=missing)

This needs to be resolved before running GWAS on human data.

### 7. Memory Budget (64GB RAM)

Both Neo4j instances running simultaneously:
| Component | RAM |
|-----------|-----|
| Neo4j human (port 7687) | 32 GB (heap=16g + pagecache=16g) |
| Neo4j yeast (port 7688) | 8 GB (heap=4g + pagecache=4g) |
| Python + GWAS | 10 GB |
| OS | 4 GB |
| **Total** | **54 GB** (fits in 64 GB) |
