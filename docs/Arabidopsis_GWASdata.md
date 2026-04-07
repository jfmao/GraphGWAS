The answer is a strong yes — and it is actually one of the best-matched genotype+phenotype resources in all of plant biology for exactly your purpose. Let me give you the full picture.

---

## The Arabidopsis 1001 Genomes Ecosystem for GraphGWAS

### The Core Resource: AraPheno + AraGWAS Catalog

AraPheno is a public phenotype database for *Arabidopsis thaliana*, hosted at `https://arapheno.1001genomes.org`, currently containing 28 studies and **1,038 phenotypes total**, with 536 published, last updated February 2025.

This is not a modest resource. *A. thaliana* provides one of the best and most extensive collections of population-scale phenotype and genotype data that exist, making the species a perfect tool for statistical method development in GWAS. These models are nowadays heavily used even outside the *A. thaliana* community, for example in human genetics. That last sentence is the key point — the field literally uses Arabidopsis as a GWAS methods benchmark because the data is so clean and complete.

### What Phenotypes Are Available

The range is biologically rich and spans exactly the kinds of traits useful for benchmarking GraphGWAS's different analysis modes:

**Quantitative traits** (ideal for linear regression testing): Five quantitative traits from 944 samples including two morphological traits (stem branching number and rosette leaf number) and three flowering time traits (days to first flower opening, days to 1cm inflorescence elongation, days from sowing to visible floral buds). Flowering time is the classic Arabidopsis GWAS phenotype — it has known strong associations at *FLC*, *FRI*, and *CRY2* that serve as positive controls.

Beyond those, AraPheno contains phenotypes covering disease resistance (bacterial and fungal pathogens), abiotic stress responses (drought, cold, salt), metabolite content, epigenetic variation, root growth, and climate adaptation traits — all quantitative or binary, all with known genetic architecture from prior GWAS, giving you ground truth to validate against.

**RNA-Seq expression data**: AraPheno has been extended to provide gene expression data from RNA-Seq experiments for hundreds of accessions and thousands of genes. Gene expression data could be used directly for GWAS as well as for transcription-wide association studies (TWAS), where transcript abundances are treated as explanatory variables. This opens a fascinating angle for GraphGWAS — eQTL analysis using the graph, which maps the functional layer from variant → gene expression → phenotype natively through the graph edges.

### The Genotype Data

The AraGWAS Catalog uses fully imputed genotype data for 2,029 *A. thaliana* lines from the 1001 Genomes Consortium, combined with existing SNP chip data, giving a SNP matrix for 2,029 accessions on **10,709,466 segregating markers**.

For GraphGWAS at Tier 1 scale (your local machine): ~2,000 samples × ~10M variants. That is very manageable. With scoped CARRIES for phenotyped samples this is around 2,000 × a few million candidate variants — well within Neo4j Community Edition on a workstation.

The VCF data is available from `https://1001genomes.org/data/GMI-MPI/releases/v3.1/`.

### The AraGWAS Catalog: Your Validation Ground Truth

The AraGWAS Catalog (`https://aragwas.1001genomes.org`) provides a publicly available, manually curated and standardized collection of marker-trait associations. It currently includes 167 phenotypes and more than 222,000 SNP-trait associations with p < 10⁻⁴, of which 3,887 are significantly associated using permutation-based thresholds.

This is enormously valuable for GraphGWAS. Every time you implement a new method — Layer 2 single-locus regression, Layer 3 burden test, Layer 3 diffusion — you can compare your results against the AraGWAS Catalog's pre-computed associations. If GraphGWAS's single-locus p-values correlate strongly with the catalog's values on the same phenotype, your implementation is validated. If your multi-locus methods find additional signals not in the catalog, those become discovery claims.

All data, including the imputed genotype matrix used for GWAS, are easily downloadable via the respective databases. The catalog also provides downloadable R scripts used for the standardized GWAS pipeline at `https://github.com/arthurkorte/GWAS`, meaning you have the exact classical-tool baseline to compare against.

### How to Download

There are three access routes:

**1. Direct download** — full database dump:
```
https://arapheno.1001genomes.org  →  "Download Database" button
```

**2. REST API** — programmatic access to phenotypes and metadata:
The REST API can be used to retrieve phenotype data and meta-information from AraPheno via URLs.
```python
# Example: fetch phenotype 6 (flowering time) for all accessions
import requests
r = requests.get('https://arapheno.1001genomes.org/rest/phenotype/6/values.json')
df = pd.DataFrame(r.json())
```

**3. Genotype VCF** — from the 1001 Genomes FTP:
```
https://1001genomes.org/data/GMI-MPI/releases/v3.1/
# File: 1001genomes_snp-short-indel_only_ACGTN.vcf.gz
# ~10M variants, 1,135 accessions, chr1-5 + organelles
```

---

## Why Arabidopsis Is Especially Well-Suited for GraphGWAS Specifically

Beyond the data availability, there are three biological properties that make Arabidopsis ideal for testing your novel graph-native methods:

**1. Inbred lines = no heterozygosity complexity.** Due to its selfing nature, individuals sampled from nature are generally inbred lines, homozygous throughout their genome. This allows efficient collection of many different phenotypes from genetically identical plants. For GraphGWAS this means CARRIES edges carry only gt=2 (homozygous alt) or absent (homozygous ref) — no heterozygous ambiguity. The CARRIES graph is simpler and faster to traverse, perfect for prototyping.

**2. Strong population structure = a real test of your stratification correction.** Arabidopsis has well-characterized population structure — Swedish, Spanish, Central Asian, North American accessions cluster distinctly. This means population stratification confounding is real and detectable, giving you a genuine test of whether GraphGWAS's graph-derived PCA covariates correctly control inflation. The AraGWAS Catalog computed kinship matrices from the same data, so you can directly compare.

**3. Known epistasis in flowering time.** The *FLC*–*FRI* interaction is one of the best-characterised epistatic systems in plant biology — *FRI* upregulates *FLC* expression, which represses flowering, but loss-of-function in either gene leads to early flowering. This gives you a **real biological positive control for your epistasis detection layer** (GraphGWAS Layer 3). If your co-occurrence enrichment test detects the *FRI*–*FLC* functional interaction, you have validated the core novel capability of the platform on real data with known biology.

---

## Recommended Starting Dataset for GraphGWAS

Given all this, here is the specific data combination I would recommend for your first test run:

| Component | Source | What to download |
|---|---|---|
| Genotype VCF | 1001genomes.org v3.1 | Chromosome 4 only (~2M variants) — smaller, but contains *FLC* and *FRI* |
| Phenotype | AraPheno via REST API | Phenotype IDs 6, 29, 30 (flowering time traits) — well-known associations |
| Ground truth associations | AraGWAS Catalog download | Associations for those same phenotype IDs — your validation target |
| Functional annotation | TAIR (arabidopsis.org) | VEP/SnpEff annotation for TAIR10 reference — Gene/Pathway nodes |

Chromosome 4 contains both *FLC* (the key flowering time repressor, chr4: ~1Mb) and is manageable in size. Running GraphGWAS Layer 2 on flowering time with Chr4 should reproduce the strong *FLC* association within minutes, directly validating your implementation against the AraGWAS Catalog result. Then Layer 3 epistasis screening between *FRI* (chr4) and *FLC* variants gives you the biological epistasis positive control.

This is genuinely the ideal benchmark dataset for your purposes — well-known biology, clean data, published ground truth, manageable scale, and biologically motivated epistasis to test your novel methods against.