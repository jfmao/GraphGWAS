# Reproducibility Guide

Every figure, table, and benchmark in the GraphGWAS manuscript can be
regenerated from this repository. This document lists the exact commands.

## 1. Environment

**Hardware**
- 64 GB RAM, 32+ CPU cores, 2 TB NVMe SSD (for full 1KG + annotations)
- Optional: NVIDIA GPU with ≥ 16 GB VRAM (for GRM accelerations)

**Software**

```bash
# Neo4j 5.x Community Edition
# Install to ~/neo4j-gwas, set password to 'graphgwas', bolt port 7688
# Heap: 16g, pagecache: 16g (for human 70M-variant graph)

# Python 3.13 with these packages
pip install neo4j==6.1 pandas numpy scipy scikit-learn matplotlib \
            bgen bgen-reader pybgen pyBigWig

# External tools
# plink2 >= v2.0.0-a.6.5LM (for BGEN I/O and GWAS)
# R + susieR (for SuSiE baseline)
# FINEMAP v1.4.2 (for FINEMAP baseline)
```

## 2. Data

**Public inputs (download once):**

| Resource | URL | Size |
|---|---|---|
| 1000 Genomes Phase 3 VCFs | https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/ | 25 GB |
| GENCODE v47 basic GTF | https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.basic.annotation.gtf.gz | 32 MB |
| GTEx v8 eQTL tar | https://gtexportal.org/home/downloads/adult-gtex/qtl | 1.5 GB |
| STRING PPI v12 | https://stringdb-downloads.org/download/protein.links.v12.0/9606.protein.links.v12.0.txt.gz | 79 MB |
| ENCODE cCRE registry | https://hgdownload.soe.ucsc.edu/gbdb/hg38/encode3/ccre/encodeCcreCombined.bb | 47 MB |
| Yeast 1011 genomes VCF | http://1002genomes.u-strasbg.fr/files/1011Matrix.gvcf.gz | 5.8 GB |

**Prebuilt dumps (release artifacts):**

| Dump | Contents | Size |
|---|---|---|
| `yeast_1011_1.92M_variants_with_annotations.dump` | Yeast 1011 × 1.92M variants × 35 traits | 519 MB |
| `human_1kg_70.7M_variants_3202_samples_chr22_annotated.dump` | 1KG + chr22 annotations | 17 GB |
| `human_1kg_70.7M_multiomics_STRING_eQTL_conservation.dump` | Same + STRING + eQTL + conservation | 17 GB |

## 3. One-time: load the annotation graph

```bash
# Restore the human 1KG dump
~/neo4j-gwas/bin/neo4j stop
cp backups/human_1kg_70.7M_multiomics_STRING_eQTL_conservation.dump \
   backups/neo4j.dump
~/neo4j-gwas/bin/neo4j-admin database load neo4j \
   --from-path=backups/ --overwrite-destination=true
~/neo4j-gwas/bin/neo4j start

# Load genome-wide gene annotations (20,092 protein-coding genes)
python -m graphgwas.annotations load-gencode \
   --gtf data/annotations/gencode/gencode.v47.basic.annotation.gtf.gz

# Build 38.4 M HAS_CONSEQUENCE proximity edges (6 parallel workers, ~70 s)
python -m graphgwas.annotations build-proximity-edges --parallelism 6

# Load 49 GTEx tissues → 43.2 M eQTL edges (4 parallel workers, ~8 min)
python -m graphgwas.annotations load-gtex-all \
   --gtex-dir data/annotations/gtex_v8/GTEx_Analysis_v8_eQTL \
   --parallelism 4

# Load STRING PPI @ score ≥ 700 → 230 K INTERACTS_WITH edges
python -m graphgwas.annotations load-string \
   --links-gz data/annotations/9606.protein.links.v12.0.txt.gz \
   --info-gz  data/annotations/9606.protein.info.v12.0.txt.gz \
   --score-threshold 700

# ENCODE cCRE regulatory elements (370 K nodes)
python -m graphgwas.annotations load-encode-ccre \
   --bed data/annotations/encode/encodeCcreCombined.bed
```

After this, the graph has 70.7 M Variant, 20,092 Gene, 370 K
RegulatoryElement, 49 GTEx tissues encoded as edge properties, and
~82 M relationships total.

## 4. One-time: convert 1KG VCFs to BGEN

```bash
# Convert each chromosome VCF to BGEN v1.2 (8-bit)
mkdir -p tests/data/human/1kGP_bgen
for chr in $(seq 1 22); do
  plink2 \
    --vcf tests/data/human/1kGP_3202/1kGP_high_coverage_Illumina.chr${chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz \
    --export bgen-1.2 bits=8 \
    --out tests/data/human/1kGP_bgen/chr${chr} \
    --memory 32000 --threads 8
done
```

Total output: ≈ 14 GB for chr1–22.

## 5. Regenerate each paper artifact

### Figures (all 7 panels)

```bash
python tests/generate_paper_figures.py
# Writes results/benchmark_v2/paper_figures/fig{1..7}.{png,pdf}
```

### Tables (all 5)

```bash
python tests/generate_paper_tables.py
# Writes results/benchmark_v2/paper_tables/table{1..5}.{md,csv}
```

### Benchmarks

| Paper section | Script | Expected runtime |
|---|---|---|
| Figure 3 / Table 2 — HBP vs SuSiE/FINEMAP | `tests/benchmark_hbp_vs_susie_finemap.py` | ≈ 20 min (30 reps × 3 scenarios) |
| Figure 4 / Table 3 — PIP calibration | `tests/benchmark_pip_calibration.py` | ≈ 60 min (200 sims × 4 h² levels) |
| Figure 4 / Table 4 — Null FPR | `tests/benchmark_null_simulations.py` | ≈ 30 min (100 nulls × 4 methods) |
| Figure 5 — M1 epistasis | `tests/benchmark_tool_comparison_fast.py` | ≈ 30 min (5 reps × 2 pairs) |
| Figure 6 / Table 5 — Weak-signal headline | `tests/benchmark_weak_signal.py` | ≈ 5 min for 79 reps |

Each script writes its JSON under `results/benchmark_v2/<name>/`. The figure
and table generators read those JSONs, so the correct sequence is:

```bash
python tests/benchmark_hbp_vs_susie_finemap.py
python tests/benchmark_pip_calibration.py
python tests/benchmark_null_simulations.py
python tests/benchmark_tool_comparison_fast.py
python tests/benchmark_weak_signal.py
python tests/generate_paper_figures.py
python tests/generate_paper_tables.py
```

## 6. Hybrid BGEN pipeline — end-to-end test

Validates that the same methods run from streaming BGEN genotypes, giving
identical results to Neo4j-loaded genotypes. This is the UK-Biobank-readiness
test (§ 4.6 of the manuscript).

```bash
# 1. Prepare a phenotype file (FID IID trait)
awk -F',' 'NR==1 {print "FID IID sim_trait"; next} {print 0, $1, $2}' \
   tests/data/human/simulated_chr22_pheno.csv > /tmp/pheno.txt

# 2. Run PLINK2 GWAS on the BGEN
plink2 --bgen tests/data/human/1kGP_bgen/chr22.bgen ref-first \
       --sample tests/data/human/1kGP_bgen/chr22.sample \
       --pheno /tmp/pheno.txt --pheno-name sim_trait \
       --glm allow-no-covars \
       --out /tmp/chr22 --threads 8 --memory 8000

# 3. Import sumstats into the graph
python -m graphgwas.summary_import plink2 \
       --file /tmp/chr22.sim_trait.glm.linear \
       --run-id ukb_test_chr22 \
       --phenotype sim_trait --p-threshold 1e-4

# 4. BGEN-backed HBP fine-mapping on the top locus
python - <<'PY'
from graphgwas.bgen_reader import BgenReader
from graphgwas.finemapping_v2 import load_locus_variants, fast_hbp_finemap
import numpy as np, pandas as pd
reader = BgenReader("tests/data/human/1kGP_bgen")
variants = load_locus_variants("chr22", 17_004_544, 50_000, source=reader)
pheno_df = pd.read_csv("tests/data/human/simulated_chr22_pheno.csv")
sids = [s.decode() if isinstance(s, bytes) else s for s in reader.samples("22")]
m = dict(zip(pheno_df["sample_id"], pheno_df["sim_trait"]))
y = np.array([m.get(s, np.nan) for s in sids])
valid = ~np.isnan(y)
for v in variants: v["dosage"] = v["dosage"][valid]
cands = fast_hbp_finemap(variants, y[valid], graph_cache={}, n_rounds=3)
cands.sort(key=lambda c: -c.pip)
for c in cands[:5]:
    print(c.variant_id, c.pip, c.in_credible_set)
PY
```

Expected output: top PIP > 0.9 at the causal variant; all 4 steps run in
under 2 min on a laptop-scale machine. Scaling to UK Biobank requires only
swapping the `BgenReader` path; the graph is the same.

## 7. Unit tests

```bash
pytest tests/test_statistics.py \
       tests/test_simulation_gwas.py \
       tests/test_yeast_offline.py \
       tests/test_bgen_reader.py \
       tests/test_hybrid_architecture.py \
       tests/test_summary_import.py
```

Should report 60 passed.

## 8. Provenance

All results in the published manuscript were produced on a single
Ubuntu 22.04 workstation with 64 GB RAM, 32 cores, an NVIDIA RTX 4090,
and Neo4j 5.26.0. Software versions and the exact commit hashes
used for each figure are recorded in `results/benchmark_v2/*/PROVENANCE.json`
where present, and in the git log otherwise.
