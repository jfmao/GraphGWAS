Yes — and honestly, yeast may be an even more powerful test system for GraphGWAS than Arabidopsis in several important respects. Let me give you the full picture.

---

## The Yeast Resource: Richer Than You Might Expect

### The Core Dataset: 1,011 Genomes (Peter et al. 2018)

This landmark study reports the whole-genome sequencing and phenotyping of 1,011 *Saccharomyces cerevisiae* isolates, providing an accurate evolutionary picture of the genomic variants shaping the species-wide phenotypic landscape.

SGD hosts all datasets from this paper including the actual genome sequences of the 1,011 isolates; core ORFs plus variable ORFs making up the pangenome; copy number variation (CNV) data; **phenotyping data for 36 conditions**; SNPs and indels relative to the S288C genome; and much more.

The data is at: `http://1011mt.epfl.ch` (EPFL mirror) and mirrored at SGD `https://www.yeastgenome.org`.

The 36 phenotyped conditions cover growth in the presence of various stressors — carbon sources, drugs, temperature extremes, osmotic stress, oxidative stress — all measured as quantitative growth fitness values. These are clean, reproducible, experimentally controlled quantitative traits, ideal for testing GraphGWAS linear regression methods.

### The 2025 Expansion: 1,086 Genomes, 8,391 Traits

This is the most exciting recent development. A Nature paper published in October 2025 presents an extensive genomic and phenotypic resource based on near telomere-to-telomere assemblies of 1,086 natural isolates. By incorporating the full spectrum of genetic variation, GWAS were conducted across **8,391 molecular and organismal traits**. The inclusion of structural variants and small insertion-deletion mutations improved heritability estimates by an average of 14.3% compared to analyses based only on SNPs.

8,391 traits is an extraordinary number — nearly an order of magnitude more than Arabidopsis's ~1,000. These include organismal traits (growth, morphology, stress response) and molecular traits (gene expression, protein abundance), making yeast by far the most phenotypically dense population genomics dataset in any model organism.

### The Expanded Population: 3,034 Genomes

A comprehensive 2024 study assembled 3,034 yeast genomes sampled from 94 countries across six continents, detecting a total of 1,918,693 SNPs and 58,947 InDels. The isolates span winemaking, beer fermentation, spirits production, bakery, clinical, and wild environments.

This gives you scalability options: start with the well-phenotyped 1,011 cohort, then scale graph traversal experiments to the 3,034 genomes for population structure analysis.

### Machine Learning Already Applied: Your Benchmark Exists

A 2025 study used machine learning models to explore phenotype predictions for 223 traits measured across 1,011 genome-sequenced *S. cerevisiae* strains. Both the proteome and transcriptome of the 1,011 *S. cerevisiae* collection have also been quantified, providing molecular traits that can aid predictions of higher-level phenotypes from genetic variants.

This means GraphGWAS's GNN layer (Layer 4) has a direct published benchmark to compare against — you can run your heterogeneous GNN on the same data and compare AUROC to the published ML results.

---

## Why Yeast Is Particularly Compelling for GraphGWAS

Beyond simple data availability, yeast has several properties that make it genuinely superior to Arabidopsis for testing GraphGWAS's novel methods specifically:

### 1. Known Epistasis at Unprecedented Scale

Yeast has the most thoroughly characterized epistatic interaction network of any eukaryote. The synthetic genetic array (SGA) experiments by the Boone lab have mapped genetic interactions for essentially the entire genome — ~400,000 pairwise interactions. This means when GraphGWAS's Layer 3 co-occurrence epistasis test identifies a variant pair, you can cross-reference it against the known genetic interaction network. A true positive is not just statistically significant — it connects to a known biological interaction. This is a validation framework no other organism can match.

### 2. CNV as a Major Phenotypic Driver — Graph Structure Advantage

The largest numbers of variants identified by genome-wide association are copy-number changes, which have a greater phenotypic effect than single nucleotide polymorphisms.

This is critical for GraphGWAS. Classical GWAS tools handle SNPs well but struggle with CNVs, which require representing presence/absence and dosage simultaneously. In GraphGWAS's graph schema, CNVs can be represented as distinct Variant nodes with a `variant_type: "CNV"` property and `gt` values representing copy number. The graph naturally encodes the non-additive, dosage-sensitive effects of CNVs through CARRIES edge properties — something a simple SNP matrix cannot represent cleanly. The 14.3% heritability improvement from including CNVs in the 2025 1,086-genome study is a direct argument for the graph representation.

### 3. Ploidy and Aneuploidy Variation

Domesticated isolates exhibit high variation in ploidy, aneuploidy and genome content. Yeast strains range from haploid to tetraploid with frequent aneuploidies. In a genotype matrix this is a nightmare — the dosage model breaks. In GraphGWAS, ploidy is simply a property on the Sample node (`ploidy: INT`), and CARRIES edge `gt` values naturally encode copy number. Queries can condition on ploidy:

```cypher
MATCH (s:Sample {ploidy: 2})-[:CARRIES]->(v:Variant)
```

This is a natural graph query. In PLINK it requires preprocessing and assumptions.

### 4. Domestication as a "Disease" Analog

One of the most compelling traits in the 1,011 dataset is the wild vs. domesticated classification. Wild and domesticated yeast clades cluster separately despite different yeast clades having been domesticated independently at different time points, supporting that domestication is the major determinant of yeast phenotypic variation.

For GraphGWAS testing, domestication status (wild=0, domesticated=1) is a **binary phenotype with known genetic architecture** — exactly the rare-disease analog. The genetic basis of domestication involves a handful of key loci (genes involved in sugar metabolism, ethanol tolerance, flocculation) but with strong population structure confounding (domesticated yeasts cluster by industry of origin). This is a realistic, biologically grounded test case for GraphGWAS's population stratification correction and binary trait association methods, far more controllable than an actual human disease.

---

## Comparison: Arabidopsis vs. Yeast for GraphGWAS Testing

| Feature | Arabidopsis (1001 Genomes) | Yeast (1011 Genomes) |
|---|---|---|
| Sample size | ~2,000 accessions | 1,011 (full phenotype), 3,034 (genotype only) |
| Number of phenotypes | ~1,038 (AraPheno) | 223 (classical), **8,391** (2025 paper) |
| Molecular phenotypes (eQTL/pQTL) | RNA-Seq for hundreds | Full proteome + transcriptome for 1,011 |
| Known epistasis | FLC-FRI system | Genome-wide SGA map (~400,000 interactions) |
| CNV importance | Moderate | **Major** — drives most phenotypic variance |
| Ploidy variation | Minimal (inbred lines) | High (haploid to tetraploid + aneuploidy) |
| Genome size | ~135 Mb (5 chromosomes) | **~12 Mb (16 chromosomes)** — much smaller |
| Population structure | Strong geographic structure | Wild/domesticated + geographic |
| Inbred/clonal | Yes (selfing) | Mostly clonal (mainly asexual) |
| Binary phenotype analog | Disease resistance | Domestication status |
| ML benchmark published | No | **Yes** (Khaiwal et al. 2025) |
| Graph annotation (SGD) | TAIR — good | **SGD — exceptional, manually curated** |

The genome size difference is important practically: at ~12 Mb, the entire yeast genome is about 1/10th the size of Arabidopsis chr4. You can load the **complete genome** of 1,011 yeast strains into GraphGWAS's Tier 1 (single workstation) trivially, whereas Arabidopsis requires subsetting to one chromosome for initial testing.

---

## Recommended Starting Point for Yeast in GraphGWAS

For your first run, the optimal combination is:

**Genotype:** 1,011 genomes SNP/indel VCF from `http://1011mt.epfl.ch/` — whole genome, ~1.9M SNPs, 1,011 strains. Loads entirely on a workstation.

**Phenotype:** The 36 growth conditions from the Peter et al. 2018 supplementary data (hosted at SGD). Start with ethanol tolerance and temperature sensitivity — both have known genetic architecture and are quantitative, clean traits suitable for linear regression validation.

**Validation ground truth:** Published GWAS results from Peter et al. 2018 and the 2025 1,086-genome paper — your Layer 2 results should reproduce their significant associations.

**Epistasis positive control:** The *FLO* gene family (flocculation) and *IRA1/IRA2* (cAMP pathway) have well-characterized functional interactions relevant to growth phenotypes — your Layer 3 co-occurrence test has known biology to recover.

**Functional annotation:** SGD's Gene Ontology annotations are the most comprehensive and manually curated of any model organism, making GraphGWAS's Gene/Pathway nodes exceptionally well-populated for yeast — far better than most organisms.

In short: yes, yeast has outstanding phenotype+genotype data, and for GraphGWAS specifically it may be the single best test organism available — small genome, many phenotypes, published ML benchmarks, exceptional functional annotation, and a genome-wide epistasis map that no other organism can match.