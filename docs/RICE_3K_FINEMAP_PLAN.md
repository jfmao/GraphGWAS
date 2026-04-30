# Rice 3K GraphGWAS fine-mapping plan

**Objective**: apply GraphGWAS's HBP / GAFM fine-mapping to the 3000 Rice
Genomes Project (3K RG) dataset and validate against textbook
known-causal genes (e.g. *sd1*, *Hd1*, *GS3*, *GW2*, *Waxy*). Extend
the paper's cross-species story beyond *Arabidopsis* (FT10,
single-locus validation) to a systematic crop-scale demonstration.

**Status**: plan only; inputs inventoried; execution queued.

---

## 1. What we have — inventory of `/mnt/data/GraphPop/data/raw/3kRG_data/`

| File | Size | Role |
|---|---|---|
| `NB_final_snp.vcf.gz` | 12.3 GB | **Primary genotypes** — 3,024 accessions × ~29 M SNPs aligned to Nipponbare (IRGSP-1.0); INFO field carries pre-computed snpEff annotations with rice gene IDs (`LOC_OsXgXXXXX`) |
| `NB_final_snp.bed.gz` / `.bim.gz` / `.fam.gz` | 4.7 GB / 174 MB / 16 KB | PLINK binary equivalent (same 3,024 samples) |
| `NB_bialSNP_pseudo_canonical_ALL.vcf.gz` | 614 MB | Consensus-format pseudo-canonical VCF; useful for schema bootstrap but not per-sample analysis |
| `3kRG_PhenotypeData_v20170411.xlsx` | 1.1 MB | IRRI phenotype records — 2,266 accessions × 53 traits (mostly IRRI Standard Evaluation System codes 1–9; morphology-dominated) |
| `41586_2018_63_MOESM3_ESM.xlsx` | 301 KB | Wang *et al.* 2018 Nature supplementary — accession metadata, geographic origin, subpopulation assignments (3,025 rows) |
| `rice_3k_panel.txt` | 70 KB | Sample → population mapping (15 subpopulations; 5 superpopulations: GJ-tmp / GJ-trp / XI-1A / XI-1B / XI-2 / XI-3 etc.) |
| `chr_rename.txt` | 90 B | Numeric-to-Chr prefix mapping (`1` → `Chr1`, …, `12` → `Chr12`) |
| `allgenome.genome.zip` | 101 MB | Nipponbare reference genome + annotations |

### What's missing for a robust fine-mapping demonstration

The included IRRI phenotype file is dominated by ordinal morphological
codes (1–9 scale for colour / pubescence / plant type) that have low
heritability and limited continuous dynamic range. For rigorous
validation we need the **quantitative traits** on the same 3K panel,
which are publicly available from SNP-seek
(<https://snp-seek.irri.org/>) under the same accession IDs:

| Trait | Units | Known causal gene(s) (validation target) |
|---|---|---|
| Plant height (**PH**) | cm | *sd1* (Os01g66100) — Green Revolution dwarfing |
| Days to heading (**DTH**) | days | *Hd1* (Os06g16370); *Hd3a* (Os06g06320); *Ghd7* (Os07g15770) |
| Grain length (**GL**) | mm | *GS3* (Os03g07050) |
| Grain width (**GW**) | mm | *GW2* (Os02g14720); *GW5* / *qSW5* (Os05g09520) |
| Grain weight (**TGW**) | g/1000 | *GS3* + *GW2* jointly |
| Amylose content (**AC**) | % | *Waxy/GBSS1* (Os06g04200) |

These six traits are the **textbook QTL validation set** — any
serious crop fine-mapping method must recover the bolded genes.
SNP-seek provides them as `phenotype_data.csv` indexed by the same
IRIS IDs in `rice_3k_panel.txt`.

### Rice-specific annotation layer (to be loaded)

| Resource | Rice equivalent of human… | Coverage / status |
|---|---|---|
| RAP-DB / MSU v7 gene models | GENCODE | ~55,986 protein-coding transcripts for Nipponbare |
| Gramene (Reactome Plant) pathways | Reactome | ~400 annotated rice pathways |
| STRING v12 rice subset (*Oryza sativa*, TaxID 39947) | STRING human | ~7,500 PPI edges at combined score ≥ 700 |
| **RicePPINet** (Liu *et al.* 2017 *Plant J*) | — (rice-specific) | **708,819 PPI edges across 16,895 proteins**, 71% TP rate (<https://netbio.sjtu.edu.cn/riceppinet/>) |
| **PRIN** (Gu *et al.* 2011 *BMC Bioinformatics*) | — (rice-specific) | 4 tiers: experimental + predicted-HC + predicted-HCov + multi-species (<https://bis.zju.edu.cn/prin/>) |
| PlantGSEA / Oryzabase GO | GO | standard GO BP/MF/CC terms |
| RiceXPro tissue-expression atlas | GTEx | 7 tissues × developmental stages; lacks formal eQTLs |
| PlantPAN cis-element / TF binding | ENCODE cCRE | partial; ~1,000 rice TFs characterised |

**Key gap** (which this analysis will highlight): the closest rice
analogue to GTEx v8 doesn't exist. RiceXPro gives tissue expression
but not cis-eQTLs with effect sizes. This is exactly the "non-model
species annotation gap" that GraphGWAS is designed to accommodate.

---

## 2. The eight-step plan

Estimated total runtime: **8–12 hours** on a 16-core workstation
with 64 GB RAM and the full 12 GB VCF.

### Step 1 — Phenotype + genotype alignment (1 h)

```bash
# 1.1 Pull the SNP-seek quantitative phenotypes (~50 MB)
mkdir -p data/rice_3k/pheno
curl -L "https://snp-seek.irri.org/_download.zul" -o data/rice_3k/pheno/snpseek_pheno.csv
# (or download manually from https://snp-seek.irri.org/_download.zul
#  selecting "3K RG 404k phenotypes" checkbox)

# 1.2 Map SNP-seek IRIS IDs onto the VCF sample IDs
python - <<'EOF'
import pandas as pd
panel = pd.read_csv("/mnt/data/GraphPop/data/raw/3kRG_data/rice_3k_panel.txt",
                    sep="\t")
pheno = pd.read_csv("data/rice_3k/pheno/snpseek_pheno.csv")
merged = panel.merge(pheno, left_on="SampleID",
                     right_on="IRIS_ID", how="inner")
print(f"Matched accessions: {len(merged)}")
# Expect: ~2,500-3,000 with non-null quantitative traits
merged[["SampleID", "Population", "Superpopulation",
        "PH", "DTH", "GL", "GW", "TGW", "AC"]].to_csv(
    "data/rice_3k/pheno/matched.tsv", sep="\t", index=False)
EOF

# 1.3 Apply chromosome rename (numeric → Chr prefix)
bcftools annotate \
    --rename-chrs /mnt/data/GraphPop/data/raw/3kRG_data/chr_rename.txt \
    /mnt/data/GraphPop/data/raw/3kRG_data/NB_final_snp.vcf.gz \
    -Oz -o data/rice_3k/NB_final_snp.renamed.vcf.gz
tabix -p vcf data/rice_3k/NB_final_snp.renamed.vcf.gz
```

### Step 2 — VCF → BGEN for hybrid GraphGWAS access (2 h)

```bash
mkdir -p data/rice_3k/bgen
for chr in Chr{1..12}; do
    plink2 \
        --vcf data/rice_3k/NB_final_snp.renamed.vcf.gz \
        --chr ${chr#Chr} \
        --export bgen-1.2 bits=8 --ref-first \
        --memory 32000 --threads 8 \
        --out data/rice_3k/bgen/${chr}
done
# Expected output: 12 BGEN files, ~2-3 GB total
```

### Step 3 — Build the rice graph database (2 h)

Two options:

**(a) From scratch** — import VCF via graphgwas:

```bash
graphgwas import vcf \
    --vcf data/rice_3k/NB_final_snp.renamed.vcf.gz \
    --samples /mnt/data/GraphPop/data/raw/3kRG_data/rice_3k_panel.txt \
    --database graphgwas_rice3k \
    --reference IRGSP-1.0 \
    --uri bolt://localhost:7689
# Runtime: ~2 hours on 16-core
```

**(b) Reuse the existing GraphPop rice Neo4j dump** (faster if
available) — see `/mnt/data/GraphPop/backups/` for the
rice-3k-with-annotations dump; restore via `neo4j-admin load`.

### Step 4 — Load rice annotation layer (1 h)

```bash
# 4.1 Gene models: RAP-DB GFF
wget https://rapdb.dna.affrc.go.jp/download/archive/RAP-DB_gff.zip \
    -O data/rice_3k/anno/rapdb.gff.zip
graphgwas annotations load-gencode \
    --gtf data/rice_3k/anno/rapdb_IRGSP-1.0.gff \
    --species rice
# Adapts the GENCODE loader to RAP-DB format; outputs HAS_CONSEQUENCE edges

# 4.2 Rice pathways from Gramene / Reactome Plant
wget "https://reactome.org/download/current/ReactomePathwaysRelation.txt" \
    -O data/rice_3k/anno/reactome_plant.tsv
graphgwas annotations load-pathway \
    --tsv data/rice_3k/anno/reactome_plant.tsv \
    --species oryza_sativa

# 4.3 Rice PPI — three orthogonal sources, load all three as
#     separate edge types so we can compare coverage and run
#     ablations later.

# 4.3a RicePPINet (Liu et al. 2017, Plant J) — highest-coverage
#      rice-specific computational network: 708,819 interactions
#      across 16,895 proteins at 71% true-positive rate,
#      ML-predicted from 11 features (GO, co-expression,
#      interolog, Rosetta-Stone, phylogenetic profile, structural
#      similarity). ~94× more edges than the STRING rice subset.
wget "https://netbio.sjtu.edu.cn/riceppinet/download/riceppinet_full_dataset.dat.gz" \
    -O data/rice_3k/anno/riceppinet_full.dat.gz
graphgwas annotations load-ppi \
    --tsv data/rice_3k/anno/riceppinet_full.dat.gz \
    --source RicePPINet \
    --id-mapping RAP-DB \
    --edge-type INTERACTS_WITH
# Expected: ~709,000 INTERACTS_WITH edges at default threshold

# 4.3b PRIN (Gu et al. 2011, BMC Bioinformatics) — Zhejiang
#      University resource with four download tiers:
#        - Experimental (literature-curated interactions)
#        - Predicted High-Confidence (strict confidence filter)
#        - Predicted High-Coverage (lenient confidence filter)
#        - Multi-Species Confidence (interolog support across species)
#      Use Experimental + High-Confidence as a second independent
#      source for confirmatory ablation experiments.
# Downloads at https://bis.zju.edu.cn/prin/download.do (manual UI)
graphgwas annotations load-ppi \
    --tsv data/rice_3k/anno/prin_experimental.xls \
    --source PRIN_exp \
    --edge-type INTERACTS_WITH
graphgwas annotations load-ppi \
    --tsv data/rice_3k/anno/prin_high_confidence.csv \
    --source PRIN_highconf \
    --edge-type INTERACTS_WITH

# 4.3c STRING v12 rice subset — orthogonal coverage
#      (combines text-mining + co-expression + database evidence).
#      Loaded as a third typed-edge source with its own
#      combined-score cutoff.
wget "https://stringdb-downloads.org/download/protein.physical.links.v12.0/39947.protein.physical.links.v12.0.txt.gz" \
    -O data/rice_3k/anno/string_rice.txt.gz
graphgwas annotations load-string \
    --links data/rice_3k/anno/string_rice.txt.gz \
    --min-score 700
# Expected: ~7,500 INTERACTS_WITH edges

# Rationale for three sources: a reviewer will (rightly) ask
# whether the fine-mapping result is robust to the choice of PPI
# layer. Loading all three lets us answer that with a single
# HBP re-run per locus, toggling the PPI source in the graph_cache.

# 4.4 (Optional) RiceXPro tissue expression as typed `eQTL`-like edges
#     Without formal eQTL analysis, we load tissue-specific
#     expression-enrichment scores as a lighter-weight prior.
graphgwas annotations load-ricexpro \
    --ricexpro-tissues data/rice_3k/anno/ricexpro_tissues.tsv
```

Expected graph state after loading (for the rice subgraph):

| Node / edge | Count |
|---|---:|
| Variant | ~29 M |
| Sample | 3,024 |
| Gene (RAP-DB + MSU v7) | ~56,000 |
| Pathway (Reactome Plant) | ~400 |
| HAS_CONSEQUENCE edges | ~12 M |
| IN_PATHWAY edges | ~25,000 |
| STRING INTERACTS_WITH edges (cross-species) | ~7,500 |
| RicePPINet INTERACTS_WITH edges (rice-specific ML) | ~709,000 |
| PRIN INTERACTS_WITH edges (experimental + HC) | ~50,000 |

### Step 5 — Population-structure correction (30 min)

Rice 3K spans 12–15 subpopulations (XI = *indica*, GJ = *japonica*,
aus, aromatic). Population structure is extreme by human standards
(*F*_ST between XI and GJ ≈ 0.5). Graph-spectral PCs + GRAMMAR+ is
more robust than a few linear PCs.

```bash
# 5.1 Compute 10 graph-spectral PCs from rare variants (MAF < 5%)
graphgwas popstruct spectral-pcs --k 10 \
    --min-af 0.001 --max-af 0.05 \
    --out data/rice_3k/results/spectral_pcs.tsv

# 5.2 GRAMMAR+ on each of the 6 target traits
for trait in PH DTH GL GW TGW AC; do
    graphgwas spectral correction \
        --phenotype $trait \
        --n-pcs 10 \
        --out data/rice_3k/results/${trait}_corrected.tsv
done
# Expected: λ_GC drops from 2-5 pre-correction to ~0.95-1.02 post
```

### Step 6 — GWAS scan (2 h for 6 traits)

```bash
for trait in PH DTH GL GW TGW AC; do
    graphgwas assoc top-hits \
        --phenotype $trait \
        --covariates data/rice_3k/results/spectral_pcs.tsv \
        --parallel 8 --persist \
        --out data/rice_3k/results/scan_${trait}.tsv
    graphgwas plot manhattan --run-id $trait \
        --threshold 7.3 \
        --out data/rice_3k/results/manhattan_${trait}.png
    graphgwas plot qq --run-id $trait \
        --out data/rice_3k/results/qq_${trait}.png
done
```

### Step 7 — Fine-mapping (minutes)

Extract the lead variant at each *known* target locus and run
GAFM + HBP. For each trait, expect 3–8 genome-wide-significant loci
(many of which are the textbook ones listed in §1).

```bash
# Example: plant height around sd1 (Os01g66100)
graphgwas finemap \
    --source bgen --bgen-dir data/rice_3k/bgen \
    --chr Chr1 --pos 38382382 --window 100000 \
    --phenotype PH --method l1 \
    -o data/rice_3k/results/sd1_PH_l1.tsv

graphgwas finemap \
    --source bgen --bgen-dir data/rice_3k/bgen \
    --chr Chr1 --pos 38382382 --window 100000 \
    --phenotype PH --method hbp \
    -o data/rice_3k/results/sd1_PH_hbp.tsv

# Similarly for DTH/Hd1 (Chr6), GL/GS3 (Chr3), GW/GW2 (Chr2),
# AC/Waxy (Chr6), etc.
```

### Step 8 — Validation + interpretation

**Primary ground truth**: the comprehensive review by
Ren, Ding & Qian (2023, *Science Bulletin* 68:314–350,
doi:10.1016/j.scib.2023.01.026) catalogues **269 rice causal
genes** with confirmed effects on grain size, grain quality,
starch biosynthesis, nutritional quality, and related traits.
This catalogue is extracted into the TSV
`data/rice_3k/ground_truth/grain_quality_causal_genes.tsv` with
columns: `gene` (primary symbol + aliases), `LOC_Os_id` (MSU v7),
`Os_id` (RAP-DB), `chr`, `category`, `pathway`,
`molecular_class`, `trait_effect`, and the Ren *et al.*
reference tag. The 269 genes distribute across all 12
chromosomes (44 on Chr2, 44 on Chr3, 30 on Chr6, …) and the
seven trait categories (104 grain-size; 68 nutritional quality;
46 starch biosynthesis; 27 grain; 14 sugar transport/loading;
5 grain shape/chalkiness/HRR; 5 grain quality).

For each trait × known-gene pair:

1. **Recovery check**: is the known causal gene in the 95% credible
   set? For textbook rice QTLs this should be **yes ≥ 90% of the time**
   — these are major-effect loci with h² > 0.3.
2. **Credible-set size**: how many variants does GAFM/HBP return?
   The narrower the better. Some rice QTLs span LD blocks of hundreds
   of variants; HBP's graph prior should collapse these using the
   pathway + PPI context.
3. **Novel candidates**: at sub-genome-wide-significant loci
   (p ∈ [1e-6, 5e-8]), does GAFM's weak-signal advantage (§2.3 of the
   paper) recover additional candidates? This is the "novel
   discovery" deliverable for the rice application.
4. **Cross-subpopulation consistency**: restrict to XI-only vs
   GJ-only and re-fine-map. Does the credible set shift between
   subspecies? (*e.g.*, *Hd1* allelic variation differs between
   indica and japonica.)

```bash
# Automated report across all 6 traits × textbook loci
graphgwas interpret --rice \
    --runs PH,DTH,GL,GW,TGW,AC \
    --known-loci data/rice_3k/known_qtl.tsv \
    --out data/rice_3k/results/validation_report.md
```

---

## 3. Expected outcomes and novel deliverables

### Definitional recovery checks (paper-ready table)

| Trait | Known-gene target (from Ren *et al.* 2023) | LOC_Os ID | Chr | Expected GAFM recovery | Expected HBP recovery |
|---|---|---|---|---|---|
| Grain width | *GW2* | LOC_Os02g14720 | Chr2 | ~80% | ~85% |
| Grain width | *GW5* / *GSE5* / *qDEC5* | LOC_Os05g09520 | Chr5 | ~75% (2 candidates in LD) | ~80% |
| Grain length | *GS3* | LOC_Os03g407050 | Chr3 | ~95% | ~95% |
| Grain length | *qGL3* / *GL3.1* / *OsPPKGAFM* | LOC_Os03g44500 | Chr3 | ~90% | ~90% |
| Grain length/width | *GL7* / *GW7* / *SLG7* | LOC_Os07g41200 | Chr7 | ~85% | ~90% |
| Grain size | *GW8* / *OsSPGAFM6* | LOC_Os08g41940 | Chr8 | ~85% | ~90% |
| Grain size | *GS5* | LOC_Os05g06660 | Chr5 | ~90% | ~90% |
| Grain size | *TGW6* | LOC_Os06g41850 | Chr6 | ~85% | ~85% |
| Grain size | *BG1* | LOC_Os03g07920 | Chr3 | ~80% | ~85% |
| Panicle architecture | *DEP1* / *qPE9-1* | LOC_Os09g26999 | Chr9 | ~85% | ~90% |
| Amylose content | *Waxy* / *GBSS1* | LOC_Os06g04200 | Chr6 | ≥ 95% | ≥ 95% |
| Plant height | *sd1* (non-grain; included for breadth) | LOC_Os01g66100 | Chr1 | ≥ 95% (single-variant CS) | ≥ 95% |

The full ground-truth set (269 genes, 7 trait categories) is in
`data/rice_3k/ground_truth/grain_quality_causal_genes.tsv` — we
run the same recovery check on every gene whose trait category
matches the GWAS phenotype, giving a systematic *rice-wide*
recovery scorecard rather than cherry-picked single-locus hits.

If the aggregate recovery rate across the 269-gene catalogue
exceeds ~70% (i.e. the gene appears in the 95% credible set
around a significant GWAS lead within ±250 kb), the paper gains
a **second rigorous cross-species validation** (after
*Arabidopsis* FT10 at PIP=0.989) and can legitimately claim
crop-scale demonstration — this time with a literature-gold-
standard ground truth, not a simulated phenotype.

### Novel-discovery angle

The 3K RG GWAS literature is rich but not exhausted. Candidate
targets for novel discovery using GAFM's weak-signal advantage:

- *DTH* loci below genome-wide significance in the full panel but
  above in within-subspecies analyses
- *GW*-modifier loci adjacent to *GW2* and *GW5* with weaker
  marginal effect
- Amylose-content modifier loci (e.g. *SSIIa* at Chr6:6.7 Mb) that
  are sub-threshold individually but coherent with *Waxy* via
  starch-biosynthesis pathway coupling

### Paper integration

Two options:

- **Option A (minimal)**: add a one-paragraph §2.8b "Rice 3K
  validation" to the current manuscript. Concrete recovery table
  for 6 textbook loci. ~150 words of Results + one supplementary
  figure.
- **Option B (medium)**: expand the cross-species section (§2.8)
  to three systematic demonstrations: yeast 1011 (already in the
  paper), *Arabidopsis* FT10 (already in §2.8), **rice 3K** (new).
  Recovery-table figure becomes main-text Figure 8 or 9.

Option A is safer for the current submission cycle; Option B is
better for a revised submission if reviewers push for more
cross-species breadth.

---

## 4. Connection to the non-model-species annotation problem

The rice application directly exercises the "awkward annotation
vector" problem identified in the paper's Discussion (see the
"Graph-native priors for non-model species" passage). Rice has:

- Well-curated gene models (RAP-DB ≈ human GENCODE) → ✓
- Pathway databases (Reactome Plant / Gramene) → partial
- PPI (STRING rice) → partial, ~7,500 edges vs 230,850 for human
- Tissue-specific eQTLs → **absent** (no rice GTEx equivalent)
- Regulatory element atlas → partial

Flat-prior methods require a *single* annotation weight vector
built by stratified LD-score regression. S-LDSC requires:

1. Curated functional categories (human Baseline 2.2: 96 categories);
2. Large GWAS N (S-LDSC is unreliable below N ≈ 20,000, rice 3K is
   N ≈ 3,000);
3. Matched population LD for the target cohort.

None of these three requirements are met by rice 3K. PolyFun on
rice would need to fabricate a weight vector from ad-hoc heuristics
(gene distance, chromatin-class guesses, synteny transfer from
human) — fragile and unvalidated.

**GraphGWAS's position**: edges that exist get added; edges that
don't, don't. There is no "ground-truth weight vector" the method
demands. The graph prior degrades gracefully with annotation depth:
for species where only gene models exist, HBP reduces to
variant→gene propagation with a flat intra-gene prior. Adding
pathway annotations enables the second layer; adding STRING PPI
enables coupling on the gene layer; adding tissue-specific
expression enables reweighting of eQTL edges. Each addition is
additive — the method doesn't break without any of them.

This is why the rice 3K application matters beyond its own
scientific result: it is a **proof of portability for the
non-model-species regime**, where fine-mapping methods that require
harmonised scalar priors cannot operate without heavy manual
curation.

---

## 5. Timeline + resources

- **Person-time**: ~2 weeks for Steps 1–7 including debugging
- **Compute**: 1 × 16-core workstation, 64 GB RAM, 100 GB scratch
- **Dependencies to install**: `bcftools`, `plink2`, `tabix`,
  `graphgwas` CLI (already available)
- **External downloads**: SNP-seek quantitative phenotypes (~50 MB),
  RAP-DB GFF (~100 MB), Reactome Plant pathways (~10 MB),
  STRING rice (~30 MB), RiceXPro tissues (~50 MB)
- **Risk items**:
  - SNP-seek may require IRRI account registration; if blocked,
    fall back to the rice quantitative phenotypes released with
    Wang *et al.* 2018 Nature (MOESM3_ESM.xlsx, already in
    `/mnt/data/GraphPop/data/raw/3kRG_data/`)
  - The full VCF (12 GB) may need a 2 TB NVMe for decompressed
    intermediate files
  - Rice chromosome naming: numeric (1–12) in VCF vs prefixed
    (Chr1–Chr12) in some annotation sources — use `chr_rename.txt`
    consistently

---

## 6. Deliverables when the plan executes

1. `docs/RICE_3K_FINEMAP_REPORT.md` — validation report, 6 traits,
   recovery table, per-locus credible sets
2. `results/rice_3k/{traits}/fine_mapping/` — JSON outputs per
   locus (schema matching the paper's Pan-UKB output)
3. `vignettes/rice-3k-plant-height.md` — companion vignette modelled
   on `full-1kg-pipeline.md` but targeting rice
4. Main-text paragraph (Option A) or sub-section (Option B) in the
   paper's cross-species Results
5. One new supplementary figure: "Rice 3K validation — known-QTL
   recovery vs credible-set size"

---

## 7. Why this matters for the paper's cross-species story

Currently the paper demonstrates cross-species portability on two
single-locus showcases: yeast 1011 (fungi) and *Arabidopsis* FT10
(plant, chromosome 4 flowering time). Rice 3K would add a
**systematic, multi-trait, crop-relevant** demonstration — the
species of global agricultural importance whose quantitative
traits drive the Green Revolution and ongoing breeding programmes.
It also lets us make a concrete claim about **non-model species
fine-mapping** (Discussion passage), which is the broader scientific
value of the graph-prior design.
