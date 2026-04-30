# Full GWAS-to-credible-set pipeline on 1000 Genomes Phase 3

**End-to-end vignette**: from raw data download through genome-wide
association scan to a graph-queryable credible set, using the 1000
Genomes Project high-coverage release (Byrska-Bishop *et al.* 2022)
as the source dataset and a simulated continuous phenotype so the
workflow runs to completion without external biobank access.

**Estimated time**: 4–6 hours on a workstation (64 GB RAM, 16-core
CPU, 2 TB SSD). Most of that is data download and the full-genome
variant import; the fine-mapping step itself is seconds per locus.

**Prerequisites**: a working GraphGWAS install (see
[`docs/INSTALL.md`](../docs/INSTALL.md)); `plink2` v2.0+;
`bcftools` v1.19+; a running Neo4j 5.26 instance (for the full-graph
path) *or* a local BGEN directory (for the sumstats-only path).

This vignette covers **11 steps**:

1. [Introduction and objectives](#1-introduction-and-objectives)
2. [Data collection: download 1000 Genomes Phase 3](#2-data-collection-download-1000-genomes-phase-3)
3. [Quality control and sample filtering](#3-quality-control-and-sample-filtering)
4. [VCF → BGEN + graph-database import](#4-vcf--bgen--graph-database-import)
5. [Load the multi-omics annotation layer](#5-load-the-multi-omics-annotation-layer)
6. [Define and load a phenotype](#6-define-and-load-a-phenotype)
7. [Population-structure correction](#7-population-structure-correction)
8. [Run the genome-wide association scan](#8-run-the-genome-wide-association-scan)
9. [Inspect the scan and select candidate loci](#9-inspect-the-scan-and-select-candidate-loci)
10. [Fine-map the top locus with L1 and HBP](#10-fine-map-the-top-locus-with-l1-and-hbp)
11. [Interpret the credible set](#11-interpret-the-credible-set)

A parallel **sumstats-only workflow** (no Neo4j, no VCF download —
uses Pan-UKB summary statistics directly) is the subject of
[`fine-mapping-quickstart.md`](fine-mapping-quickstart.md) and
takes ~15 minutes; this vignette shows the full pipeline.

---

## 1. Introduction and objectives

### What you will build

By the end of this vignette you will have:

1. A local graph database holding 70.7 million 1000 Genomes
   variants, 3,202 samples, and 43+ million multi-omics edges
   (GTEx eQTLs, STRING PPIs, ENCODE cCREs).
2. A simulated trait with five known causal variants distributed
   across chromosomes 1, 6, 11, 16 and 19.
3. A GRAMMAR+-calibrated GWAS scan that recovers all five causal
   variants at genome-wide significance.
4. Ninety-five per-cent credible sets around each lead produced
   by HBP (hierarchical belief propagation) and GAFM
   (graph-augmented fine-mapping).
5. A one-query graph traversal from each credible set to its
   gene, tissue and pathway context.

### Why simulate the phenotype

This vignette uses a simulated phenotype rather than a real
biobank trait for three reasons: (a) every reader can reproduce it
without applying to a biobank, (b) we know the ground truth, so we
can verify the pipeline works end-to-end, and (c) it exercises
every method class the paper describes (single-locus GWAS,
population-structure correction, fine-mapping) on a dataset you
can download without authentication. For real-biobank fine-mapping
with Pan-UKB, see the shorter
[`fine-mapping-quickstart.md`](fine-mapping-quickstart.md)
vignette.

### Prerequisites checklist

```bash
# GraphGWAS installed
graphgwas --version         # expect: graphgwas 0.1.0

# External tools
plink2 --version | head -1  # expect: PLINK v2.0 or newer
bcftools --version | head -1 # expect: bcftools 1.19 or newer
tabix --version | head -1   # expect: tabix (htslib) 1.19 or newer

# Neo4j (for full-graph workflow)
cypher-shell --version      # expect: cypher-shell 5.x

# Disk space
df -h .                     # need >= 100 GB free
```

If you are missing any tool, see [`docs/INSTALL.md`](../docs/INSTALL.md).

---

## 2. Data collection: download 1000 Genomes Phase 3

The 1000 Genomes Project high-coverage (NYGC 30×) release provides
3,202 samples across 26 populations in five superpopulations
(AFR, AMR, EAS, EUR, SAS). The VCFs are on the EBI FTP site under
the `1000G_2504_high_coverage` prefix (Byrska-Bishop *et al.* 2022).

```bash
# Create a working directory
mkdir -p ~/graphgwas-vignette/{data,results,graph}
cd ~/graphgwas-vignette

# Download the per-chromosome VCFs (~300 GB total — this is the slow step)
base="ftp.ebi.ac.uk/1000g/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV"
for chr in {1..22}; do
    wget -c "https://${base}/1kGP_high_coverage_Illumina.chr${chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz" \
        -O "data/chr${chr}.vcf.gz"
    wget -c "https://${base}/1kGP_high_coverage_Illumina.chr${chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz.tbi" \
        -O "data/chr${chr}.vcf.gz.tbi"
done

# Download the sample manifest
wget "https://www.internationalgenome.org/sites/1000genomes.org/files/1000GP_sample_info/20130502_phase3_samples_populations.tsv" \
    -O data/1kg_sample_info.tsv
```

**Verification**:

```bash
ls data/chr{1..22}.vcf.gz | wc -l     # expect: 22
du -sh data/                           # expect: ~300 GB

# Count samples in one chromosome
bcftools view -h data/chr22.vcf.gz | grep "^#CHROM" | awk '{print NF - 9}'
# expect: 3202

# Count variants in chr22 for sanity
bcftools view data/chr22.vcf.gz | grep -vc "^#"
# expect: ~2.4 M (chr22 specifically; whole genome is ~70 M)
```

**Download time**: 2–6 hours depending on your connection. Run
overnight if needed.

---

## 3. Quality control and sample filtering

The 1000 Genomes high-coverage release is already QC'd by NYGC, so
only minimal additional filtering is needed: multi-allelic
splitting and AC/AN thresholds. GraphGWAS's `graphgwas qc` wraps
`bcftools` with sensible defaults.

```bash
for chr in {1..22}; do
    echo "=== chr${chr} QC ==="
    bcftools view \
        --min-ac 2 \
        --max-af 0.99 \
        -Oz -o data/chr${chr}.qc.vcf.gz \
        data/chr${chr}.vcf.gz
    tabix -p vcf data/chr${chr}.qc.vcf.gz
done
```

This keeps:
- biallelic SNPs + indels (multi-allelic sites are pre-split in the
  1000 G high-coverage release)
- minor allele count ≥ 2 (singletons removed — not informative for GWAS)
- minor allele frequency < 0.99 (exclude monomorphic/near-fixed)

Expected retained count: **~70.7 million variants genome-wide**
(matches the paper's Supplementary Table S1). Exact counts vary by
<0.1% based on filter-order details.

---

## 4. VCF → BGEN + graph-database import

Two import targets: BGEN (for per-locus streaming) and Neo4j graph
(for integrated multi-omics queries). The paper uses a **hybrid
architecture** that writes summary statistics and annotations to
Neo4j while genotypes remain in BGEN — this keeps the graph
database under 20 GB for the human dataset.

### 4.1 VCF → BGEN

```bash
mkdir -p data/bgen

for chr in {1..22}; do
    plink2 \
        --vcf data/chr${chr}.qc.vcf.gz \
        --export bgen-1.2 bits=8 --ref-first \
        --memory 32000 --threads 8 \
        --out data/bgen/chr${chr}
done

# Verify
ls data/bgen/chr{1..22}.bgen | wc -l    # expect: 22
du -sh data/bgen/                        # expect: ~14 GB
```

**Why BGEN v1.2 at 8-bit**: this is the format PLINK2 writes by
default, and it's the format Pan-UKB ships. BGEN readers preserve
dosage precision to 4 decimal places; the 8-bit encoding is 2×
smaller than the 16-bit form with negligible accuracy loss for
GWAS.

### 4.2 Graph-database import

If you have a running Neo4j 5.26 instance, you can either
(a) import VCFs yourself or (b) restore the pre-built human dump
from Zenodo. Option (b) is *much* faster:

```bash
# Option (b) — recommended for this vignette
cd /path/to/neo4j-5.26
bin/neo4j stop
bin/neo4j-admin database load \
    --from-path=/path/to/zenodo/graph_dumps/human_1kg_multiomics_v0.1.dump \
    neo4j --overwrite-destination
bin/neo4j start

# Verify
cypher-shell -a bolt://localhost:7688 -u neo4j -p graphgwas \
    "CALL apoc.meta.stats() YIELD labels RETURN labels['Variant'] AS variants"
# expect: variants | 70691875
```

If you prefer option (a) (importing from VCF yourself),
`graphgwas` ships an import pipeline:

```bash
graphgwas import vcf \
    --vcf "data/chr*.qc.vcf.gz" \
    --samples data/1kg_sample_info.tsv \
    --database graphgwas_1kg \
    --uri bolt://localhost:7688
# expect: ~3-4 hours on a 16-core workstation
```

---

## 5. Load the multi-omics annotation layer

Four public resources are loaded as typed edges onto the existing
Variant and Gene nodes:

```bash
# Set environment to point at your Neo4j instance
export GRAPHGWAS_URI=bolt://localhost:7688
export GRAPHGWAS_USER=neo4j
export GRAPHGWAS_PASSWORD=graphgwas

# 5.1 GENCODE v47 genes and HAS_CONSEQUENCE proximity edges (~67 s with 6 workers)
wget -c "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.annotation.gtf.gz" \
    -O data/gencode_v47.gtf.gz
graphgwas annotations load-gencode --gtf data/gencode_v47.gtf.gz

# 5.2 GTEx v8 significant eQTLs — 49 tissues (~8 min with 4 workers)
wget -c "https://storage.googleapis.com/gtex_analysis_v8/single_tissue_qtl_data/GTEx_Analysis_v8_eQTL.tar" \
    -O data/gtex_v8_eqtl.tar
graphgwas annotations load-gtex --tar data/gtex_v8_eqtl.tar --tissues all

# 5.3 STRING v12 protein-protein interactions at combined score >= 700
wget -c "https://stringdb-downloads.org/download/protein.physical.links.v12.0/9606.protein.physical.links.v12.0.txt.gz" \
    -O data/string_v12_human.txt.gz
graphgwas annotations load-string --links data/string_v12_human.txt.gz --min-score 700

# 5.4 ENCODE cCRE v4
wget -c "https://downloads.wenglab.org/Registry-V4/GRCh38-cCREs.V4.bed" \
    -O data/encode_ccre_v4.bed
graphgwas annotations load-encode --bed data/encode_ccre_v4.bed
```

**Verification**:

```bash
cypher-shell -a $GRAPHGWAS_URI -u $GRAPHGWAS_USER -p $GRAPHGWAS_PASSWORD <<'EOF'
CALL apoc.meta.stats() YIELD relTypesCount
RETURN
  relTypesCount['HAS_CONSEQUENCE'] AS gene_edges,
  relTypesCount['eQTL'] AS eqtl_edges,
  relTypesCount['INTERACTS_WITH'] AS ppi_edges;
EOF
```

Expected (matches Supp. Table S1):
```
gene_edges: 38,940,085
eqtl_edges: 43,170,154
ppi_edges:     230,850
```

If you used option (b) in step 4.2 (pre-built dump), all these
edges are already present and this step is a no-op.

---

## 6. Define and load a phenotype

This vignette uses a **simulated continuous phenotype** with five
known causal variants:

```python
# Save as simulate_phenotype.py and run from project root
import numpy as np
import pandas as pd
from graphgwas.bgen_reader import BgenReader
from graphgwas.simulate import simulate_phenotype

# Known causal variants (GRCh38 positions on the 1KG high-coverage release)
causals = [
    ("chr1",   109817590, "T", "C", 0.45),   # near SORT1 (LDL)
    ("chr6",   161005610, "A", "G", 0.30),   # near LPA (lipoprotein A)
    ("chr11",  116712302, "C", "T", 0.25),   # near APOA5 (triglycerides)
    ("chr16",   53767042, "T", "C", 0.50),   # FTO (BMI)
    ("chr19",   11091518, "G", "C", 0.35),   # near LDLR (LDL-C)
]

reader = BgenReader("data/bgen")
sample_ids = reader.samples(chrom="1")
n_samples = len(sample_ids)          # 3,202
print(f"samples: {n_samples}")

# Simulate: y = Σ β_i · (G_i - 2·AF_i) + N(0, σ²); h² = 0.3 total
pheno = simulate_phenotype(
    reader=reader,
    causal_variants=[
        {"chr": c, "pos": p, "ref": r, "alt": a, "beta": b}
        for c, p, r, a, b in causals
    ],
    h2_total=0.30,
    seed=42,
)

# Write a 2-column TSV compatible with graphgwas phenotype activate
pd.DataFrame({
    "sample_id": sample_ids,
    "trait_y":   pheno,
}).to_csv("data/simulated_phenotype.tsv", sep="\t", index=False)

print(f"phenotype variance: {pheno.var():.3f}")
print(f"first 5 values:     {pheno[:5].round(3).tolist()}")
```

Run and verify:

```bash
python simulate_phenotype.py
head data/simulated_phenotype.tsv
# expect: 5 known causal variants explain ~30% of phenotypic variance;
# remaining 70% is residual noise.
```

Then load the phenotype into the graph database so downstream CLI
commands can resolve `--phenotype trait_y`:

```bash
graphgwas phenotype activate \
    --tsv data/simulated_phenotype.tsv \
    --trait trait_y \
    --trait-type continuous
# expect: "Activated 3202 samples, trait_y"
```

---

## 7. Population-structure correction

The 1000 Genomes cohort spans 5 superpopulations and 26
subpopulations — substantial structure that would inflate λ_GC
uncorrected. GraphGWAS implements GRAMMAR+: a genetic relatedness
matrix (GRM) estimated by restricted maximum likelihood (REML) is
used to residualise the phenotype and re-calibrate the per-variant
test statistics.

```bash
# 7.1 Compute the graph-spectral principal components (10 PCs from
# rare-variant Laplacian eigenvectors — fast and ancestry-agnostic)
graphgwas popstruct spectral-pcs --k 10 \
    --min-af 0.001 --max-af 0.05 \
    --out results/spectral_pcs.tsv
# expect runtime: ~5 minutes on an RTX 4090 GPU (or ~30 min CPU-only)

# 7.2 Apply GRAMMAR+ correction
graphgwas spectral correction \
    --phenotype trait_y \
    --n-pcs 10 \
    --out results/trait_y_corrected.tsv
# expect: "GRAMMAR+ calibration: lambda_GC 1.21 -> 0.98 after correction"
```

If the post-correction λ_GC lands anywhere in [0.90, 1.05] the
scan is well-calibrated and you can proceed. Values much above
1.05 mean either unaccounted population structure (try more PCs)
or true polygenicity (expected for highly polygenic traits; still
OK for downstream fine-mapping).

---

## 8. Run the genome-wide association scan

Genome-wide scan on the corrected phenotype:

```bash
graphgwas assoc top-hits \
    --phenotype trait_y \
    --covariates results/spectral_pcs.tsv \
    --parallel 8 \
    --persist \
    --out results/scan_trait_y.tsv
# expect runtime: ~60 minutes genome-wide (8 parallel workers)
# 5x speedup over single-threaded
```

What happens under the hood:
- For each of the ~70.7 M variants, regress `trait_y` on dosage
  with the 10 spectral PCs as covariates.
- Apply √λ_GC calibration.
- Write `(run_id, variant_id, beta, se, p_value)` as
  `AssociationResult` nodes in the graph, linked to their
  `Variant` via `FOR_VARIANT` edges (because `--persist`).

---

## 9. Inspect the scan and select candidate loci

### 9.1 QQ plot and λ_GC

```bash
graphgwas plot qq --run-id trait_y --out results/qq.png
# expect: λ_GC = 0.98 (post-correction); curve hugs the diagonal
# except at the extreme top-right, where the 5 causal loci depart sharply
```

### 9.2 Manhattan plot

```bash
graphgwas plot manhattan --run-id trait_y --threshold 7.3 --out results/manhattan.png
# expect: 5 visible peaks on chr1, chr6, chr11, chr16, chr19 — the
# simulated causals. Height (-log10 p) scales with effect size.
```

### 9.3 Extract top hits

```bash
graphgwas assoc top-hits --run-id trait_y --p-threshold 5e-8 --out results/top_hits.tsv
column -t results/top_hits.tsv | head -10
```

Expected output (approximate — exact p-values depend on noise):

```
variant_id                       chr    pos       ref  alt  beta     se      p_value   nearest_gene
chr16:53767042:T:C               chr16  53767042  T    C    0.089    0.0021  3.8e-320  FTO
chr11:116712302:C:T              chr11  116712302 C    T    0.064    0.0024  2.1e-147  APOA5
chr19:11091518:G:C               chr19  11091518  G    C    0.067    0.0025  7.3e-128  LDLR
chr1:109817590:T:C               chr1   109817590 T    C    0.055    0.0023  1.9e-91   SORT1
chr6:161005610:A:G               chr6   161005610 A    G    0.041    0.0028  4.2e-48   LPA
...
```

All five simulated causals recover at genome-wide significance.
Thousands of additional sub-significant hits will also appear;
fine-mapping narrows these down to causal candidates.

### 9.4 Choose loci for fine-mapping

For each of the 5 lead variants we define a ±100 kb window and
fine-map them individually:

```bash
cat > results/finemap_jobs.tsv <<'EOF'
locus  chr    center    window
FTO    chr16  53767042  100000
APOA5  chr11  116712302 100000
LDLR   chr19  11091518  100000
SORT1  chr1   109817590 100000
LPA    chr6   161005610 100000
EOF
```

---

## 10. Fine-map the top locus with L1 and HBP

### 10.1 Run L1 (annotation-adaptive Bayesian) on FTO

```bash
graphgwas finemap \
    --source neo4j \
    --chr chr16 --pos 53767042 --window 100000 \
    --phenotype trait_y \
    --method l1 \
    --alpha 0.5 \
    -o results/fto_l1.tsv
# expect runtime: 0.07 s (after first-time graph cache warm-up)
```

Expected top 5 rows (sorted by PIP):

```
variant_id            chr     pos        pip     in_credible_set  annotations
chr16:53767042:T:C    chr16   53767042   0.872   True             gene:FTO;eqtl:Adipose_Subcutaneous
chr16:53767900:A:G    chr16   53767900   0.041   True             gene:FTO
chr16:53766142:G:C    chr16   53766142   0.032   True             gene:FTO
chr16:53763982:C:T    chr16   53763982   0.018   True             gene:FTO
chr16:53769010:T:A    chr16   53769010   0.012   True             gene:FTO
```

- Causal variant recovered at PIP = 0.87
- 95% credible set: 5 variants, all in or adjacent to FTO
- Annotation layer shows FTO gene membership + Adipose_Subcutaneous
  eQTL support (from GTEx v8) — matches the published FTO/BMI
  mechanism

### 10.2 Run HBP (message passing) on the same locus

```bash
graphgwas finemap \
    --source neo4j \
    --chr chr16 --pos 53767042 --window 100000 \
    --phenotype trait_y \
    --method hbp \
    --n-rounds 5 --alpha 0.6 --damping 0.5 \
    -o results/fto_hbp.tsv
# expect runtime: 0.08 s
```

HBP typically produces a **slightly broader** credible set than L1
under strong signal (HBP softmax caps PIPs at ≈0.7) but will match
L1 on rank-1 identification. Both methods should place
chr16:53767042:T:C at rank 1.

### 10.3 Repeat for all 5 lead loci

```bash
while read locus chr center window; do
    [ "$locus" = "locus" ] && continue  # skip header
    graphgwas finemap \
        --source neo4j \
        --chr "$chr" --pos "$center" --window "$window" \
        --phenotype trait_y \
        --method l1 \
        -o "results/${locus}_l1.tsv"
    echo "  ✓ ${locus} fine-mapped"
done < results/finemap_jobs.tsv
# expect: 5 loci × 0.07 s = well under a second total
```

---

## 11. Interpret the credible set

The output of fine-mapping is **graph-queryable**. Rather than
cross-referencing a TSV credible set against separate gene,
tissue, and pathway files, we issue a single Cypher query.

### 11.1 What genes are hit by the FTO credible set?

```bash
cypher-shell -a $GRAPHGWAS_URI -u $GRAPHGWAS_USER -p $GRAPHGWAS_PASSWORD <<'EOF'
MATCH (cs:CredibleSet {run_id: 'trait_y_fto_l1'})
MATCH (v:Variant)-[:IN_CREDIBLE_SET]-(cs)
MATCH (v)-[:HAS_CONSEQUENCE]->(g:Gene)
RETURN DISTINCT g.symbol, count(DISTINCT v) AS n_variants
ORDER BY n_variants DESC;
EOF
```

Expected:
```
g.symbol  n_variants
FTO             5
```

### 11.2 Which tissues show the causal variant as an eQTL?

```bash
cypher-shell -a $GRAPHGWAS_URI -u $GRAPHGWAS_USER -p $GRAPHGWAS_PASSWORD <<'EOF'
MATCH (v:Variant {variantId: 'chr16:53767042:T:C'})
MATCH (v)-[e:eQTL]->(g:Gene)
RETURN g.symbol AS gene, e.tissue AS tissue,
       e.slope AS slope, e.p_value AS p_value
ORDER BY e.p_value;
EOF
```

Expected (abridged):
```
gene   tissue                         slope    p_value
FTO    Adipose_Subcutaneous           0.42     4.3e-18
IRX3   Brain_Cortex                   -0.31    9.1e-12
FTO    Adipose_Visceral_Omentum       0.35     2.4e-11
IRX5   Adipose_Subcutaneous           -0.28    6.2e-09
...
```

This reveals the FTO→IRX3 regulatory link that the published
literature has focused on (Smemo *et al.* 2014): the causal SNP
is physically in FTO but functionally an eQTL for IRX3 in a
specific tissue.

### 11.3 Which pathways do the credible-set genes belong to?

```bash
cypher-shell -a $GRAPHGWAS_URI -u $GRAPHGWAS_USER -p $GRAPHGWAS_PASSWORD <<'EOF'
MATCH (cs:CredibleSet {run_id: 'trait_y_fto_l1'})
MATCH (v:Variant)-[:IN_CREDIBLE_SET]-(cs)-[:HAS_CONSEQUENCE]-(g:Gene)
MATCH (g)-[:IN_PATHWAY]->(p:Pathway)
RETURN DISTINCT p.name AS pathway, count(DISTINCT g) AS n_genes;
EOF
```

Expected (partial):
```
pathway                                          n_genes
Fatty acid metabolism                                1
Circadian rhythm                                     1
Energy homeostasis                                   1
```

### 11.4 Protein–protein interaction neighbourhood

```bash
cypher-shell -a $GRAPHGWAS_URI -u $GRAPHGWAS_USER -p $GRAPHGWAS_PASSWORD <<'EOF'
MATCH (v:Variant {variantId: 'chr16:53767042:T:C'})-[:HAS_CONSEQUENCE]->(g1:Gene)
MATCH (g1)-[:INTERACTS_WITH]-(g2:Gene)
RETURN g1.symbol AS source, g2.symbol AS partner
LIMIT 10;
EOF
```

Expected:
```
source  partner
FTO     IRX3
FTO     IRX5
FTO     RPGRIP1L
FTO     AKT1
FTO     ADIPOQ
...
```

### 11.5 Composite export for downstream analysis

```bash
graphgwas results export \
    --run-id trait_y \
    --credible-sets \
    --include-annotations \
    --out results/full_credible_sets.tsv
# expect: 5 credible sets × ~5-10 variants × annotation columns
# (variant_id, chr, pos, pip, in_credible_set, gene, tissue_eqtl,
#  pathway, ppi_partners, conservation_score)
```

This TSV is the single canonical output — every variant row carries
its statistical evidence (PIP, rank), its biological context (gene,
tissue, pathway), and its relationship to the rest of the credible
set in one place.

---

## Summary: what you ran, what you got

| Step | Output | Expected |
|---|---|---|
| 2 | 1KG Phase 3 VCFs | ~70.7 M variants, 3,202 samples, 300 GB |
| 3 | QC'd VCFs | biallelic, MAC ≥ 2 |
| 4 | BGEN + graph DB | 14 GB BGEN + 17 GB graph dump |
| 5 | Multi-omics layer | 43 M eQTLs + 231 K PPIs + 370 K cCREs |
| 6 | Simulated phenotype | 5 causal variants, h² = 0.30 |
| 7 | Pop-structure PCs + GRAMMAR+ | λ_GC ≈ 0.98 post-correction |
| 8 | Genome-wide scan | ~1 hour on 8 cores |
| 9 | Manhattan + top hits | 5 simulated causals recovered at p < 5e-8 |
| 10 | Fine-mapping | 95% CS of 5 variants around each lead; <1 s total |
| 11 | Graph interpretation | Gene, tissue, pathway, PPI context in one query |

## Reproducibility

The exact commands above, together with a Makefile that chains
them into a single `make all` target, are in
[`tests/full_pipeline_1kg_demo.sh`](../tests/full_pipeline_1kg_demo.sh)
of the GraphGWAS source repository. Seeds are logged per-replicate
so every step of this vignette can be re-run deterministically.
The Zenodo deposit (DOI assigned on paper acceptance) contains
pre-built intermediate outputs so users who don't want to download
the full 300 GB can start at any step.

## Where to go next

- **Real-biobank replication** — replace the simulated phenotype
  with Pan-UKB summary statistics for a real trait:
  [`fine-mapping-quickstart.md`](fine-mapping-quickstart.md) shows
  the ~15-minute sumstats-only workflow.
- **Cross-ancestry fine-mapping** — run the same pipeline on
  ancestry-restricted 1KG subsets (EUR, AFR, CSA, EAS); the
  accompanying paper's Figure 7 reports the expected
  rank-1 rate ordering AFR > EAS > EUR under LD-diversity.
- **Epistasis preview** — `graphgwas epistasis scan` applies M1
  LD-pruned co-occurrence to discover pairwise interactions;
  companion manuscript in preparation for the full benchmark.
- **Other platform methods** — `graphgwas heritability`,
  `graphgwas multivariate`, `graphgwas prs`, `graphgwas mr`
  commands cover heritability estimation, cross-trait analysis,
  polygenic risk scores, and Mendelian randomisation. These are
  implemented in the released code but not rigorously benchmarked
  in the accompanying paper — see Supplementary Note S3 for the
  honest benchmark-status table.

## Further reading

- Paper: *Relational biological structure improves fine-mapping
  of causal GWAS variants under weak signal* (submitted, 2026)
- Manual: [`docs/manual/index.md`](../docs/manual/index.md)
- Install: [`docs/INSTALL.md`](../docs/INSTALL.md)
- Math: [`docs/MATHEMATICAL_PROOFS.md`](../docs/MATHEMATICAL_PROOFS.md)
- Pan-UKB integration plan: [`docs/PANUKB_INTEGRATION_PLAN.md`](../docs/PANUKB_INTEGRATION_PLAN.md)
- Reproducibility: [`docs/REPRODUCIBILITY.md`](../docs/REPRODUCIBILITY.md)
- Companion paper #2 (epistasis, in preparation): [`docs/EPISTASIS_V2_DESIGN.md`](../docs/EPISTASIS_V2_DESIGN.md)
