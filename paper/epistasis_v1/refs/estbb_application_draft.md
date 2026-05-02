# Estonian Biobank — Data Access Application Draft

**Submitted by:** Jian-Feng Mao, Ph.D., Umeå University (jianfeng.mao@umu.se)
**Project title:** Graph-native epistasis detection rescues anti-conservative
bias in biobank-scale GWAS
**Date:** [TO FILL — submission date]
**Companion preprint:** Yelmen et al. 2026, bioRxiv 10.1101/2025.11.21.689603

---

## A. Preliminary inquiry (Stage 1 — ~5 business-day response)

Submit to: data.access@geenivaramu.ee  *(verify exact contact at https://genomics.ut.ee/en/content/estonian-biobank)*

**Subject:** Data access inquiry — graph-native epistasis bias-rescue method

Dear Estonian Biobank Data Access Office,

Our group at Umeå University (Sweden) is preparing a methodological
manuscript for *Nature Genetics* on **graph-native detection of
omitted-interaction bias in genome-wide association studies**. The work
extends the bias-quantification framework recently published by Yelmen
*et al.* (2026, bioRxiv 10.1101/2025.11.21.689603 — using EstBB data),
adding (i) a graph-typed maximum-correlation diagnostic (ρ_max) computed
from biological motif pairs, and (ii) a covariate-rescue protocol that
collapses the spurious-significance regime in REGENIE LMM analyses.

To validate the rescue claim on the same data substrate as Yelmen *et
al.*, we request access to:

- **Chromosome 21 + 22 GSA SNP-array genotypes** (the same ~8,170
  bi-allelic SNPs Yelmen *et al.* used) for ~100,000–200,000 unrelated
  participants
- **A small set of quantitative phenotypes** (1–3 traits, registry-derived
  is sufficient — height / BMI / lipids would all serve)

Could you confirm:

1. Approximate cost for this data tier (genotyping array + ≤3 registry
   phenotypes for ~100K individuals)?
2. Whether a local Estonian collaborator is required for this scope?
3. Estimated timeline from preliminary inquiry to data delivery?

We are happy to expand into a full SAC application once feasibility is
confirmed.

Sincerely,
Jian-Feng Mao, Ph.D.

---

## B. Scientific Advisory Committee (SAC) application — research proposal

### B.1 Project title
Graph-native epistasis detection rescues anti-conservative bias in
biobank-scale GWAS — extension of Yelmen et al. 2026 to method
development

### B.2 Principal investigator
Jian-Feng Mao, Ph.D., Professor, Department of Plant Physiology, Umeå
Plant Science Centre, Umeå University
- ORCID: [0000-0000-0000-0000]
- Email: jianfeng.mao@umu.se
- Lab page: [URL]

### B.3 Co-investigators
- Ehsan Estaji (Umeå University)
- Shi-Wei Zhao (Umeå University)
- Zhao-Yang Chen (Umeå University)
- Shuai Nie (Guangdong Academy of Agricultural Sciences)

### B.4 Background

Yelmen *et al.* (2026, in press) recently demonstrated that classical
genome-wide association studies — including state-of-the-art linear
mixed-model implementations such as REGENIE — produce
**anti-conservative bias** in test statistics when the true genetic
architecture involves epistasis but the analytical model omits it. The
expected mean shift of the *t*-statistic is

    μ = ρ · √(λ · n) / √(1 − λ ρ²)

with ρ the correlation between target SNP and realised interaction
signal, λ the variance fraction from interactions, and *n* the sample
size. Critically, the bias **grows with sample size** — biobank-scale
GWAS is more vulnerable, not less.

Yelmen *et al.* used the Estonian Biobank's 210,145-sample chr21 + chr22
genotype panel to validate the framework empirically. Their analysis
revealed thousands of spuriously genome-wide-significant SNPs under
strict no-path null DGPs.

### B.5 Research question
**Can the spurious-significance regime characterised by Yelmen *et al.*
be rescued by adding biology-typed interaction covariates detected by
graph-native methods?**

### B.6 Specific aims

**Aim 1.** Reproduce Yelmen *et al.* Fig 1 using the chr21+chr22
GSA-array genotypes; verify our independently-implemented ρ_max
estimator (graphgwas.bias module, MIT-licenced, https://github.com/jfmao/GraphGWAS)
matches their published distributions within tolerance.

**Aim 2.** Compute the **graph-typed ρ_max** (the key paper #2
contribution): substitute Yelmen *et al.*'s random Z with a biology-typed
Z built from same-gene (motif P1), PPI-partner (P2), and same-pathway
(P3) interaction features derived from public annotations (GENCODE,
STRING, Reactome). Test the hypothesis that motif-typed Z does *not*
inflate ρ_max for null target SNPs.

**Aim 3.** Run REGENIE LMM on simulated phenotypes (Yelmen's strict
no-path null DGP) under three conditions:
(a) baseline LMM,
(b) LMM with random-Z interaction covariates,
(c) LMM with motif-Z interaction covariates from GraphGWAS M2.
Report the spurious-significant-hit count for each. Headline rescue
metric: percent reduction from (a) → (c).

### B.7 Data requested

| Data tier | Sample-size we need | Variant scope | Phenotype |
|-----------|--------------------:|---------------|-----------|
| Genotyping array (GSA) | 100,000–200,000 unrelated | chr21 + chr22 only | n/a |
| Registry phenotypes | same individuals | n/a | 1–3 quantitative traits (height, BMI, lipids) for power calibration |

**No** WGS, exome, methylation, mRNA, or microbiome data needed.

### B.8 Data management plan

- Storage: encrypted analysis on Umeå University HPC cluster (HPC2N,
  ISO 27001 certified).
- Processing: pseudonymised IDs only; no re-identification attempts;
  no linkage to external datasets beyond the registry phenotypes
  released by EstBB.
- Sharing: results published as aggregate summary statistics in the
  *Nature Genetics* manuscript and on Zenodo. Individual-level data
  remain on the Umeå HPC and are deleted on project completion.
- Retention: 5 years from project end, per EstBB Tier-1 standard.

### B.9 Ethics

Analysis is methodological — no clinical interpretation, no return of
results to participants, no contact with participants. All proposed
analyses fall within the participant consent forms granted to EstBB.
Local ethics review at Umeå University waived for de-identified
secondary data; Estonian Committee on Bioethics and Human Research
(EBIN) application will be filed alongside the SAC application per
EstBB procedure.

### B.10 Timeline

| Month | Milestone |
|-------|-----------|
| 1 | SAC application submission |
| 2–3 | EBIN ethics review |
| 3–4 | Data delivery |
| 4–6 | Aim 1 (Yelmen reproduction on EstBB) |
| 6–8 | Aim 2 (motif-typed ρ_max) |
| 8–10 | Aim 3 (REGENIE rescue benchmark) |
| 10–12 | Manuscript revision + Zenodo deposit |

### B.11 Dissemination

- Open-source software: graphgwas Python package on PyPI (already live
  at https://pypi.org/project/graphgwas/), MIT licence, source at
  https://github.com/jfmao/GraphGWAS
- Preprint: bioRxiv at submission
- Peer-reviewed publication: *Nature Genetics* methodology paper
- Zenodo deposit: code snapshot, simulated-phenotype seeds, REGENIE
  configurations, and aggregate summary statistics

### B.12 Funding statement

[TO FILL by user — Wallenberg Initiatives in Forest Research (WIFORCE)
covered paper #1; paper #2 funding source needs confirmation.]

---

## C. Notes for the user before submission

1. **Email contact** at the top is from the EstBB website; double-check
   the exact email/portal URL by visiting
   https://genomics.ut.ee/en/content/estonian-biobank (the page text is
   what we drafted from).
2. **Cost estimate** based on EstBB published tariff:
   - €6,020 processing fee
   - €0.30/person × 100K registry phenotypes ≈ €30,000
   - Total ≈ €36,000–€50,000 for the scope above
   Confirm budget availability before formal SAC application.
3. **Local Estonian collaborator** — not strictly required per the EstBB
   page, but informally accelerates SAC review. Possible candidates:
   Tõnu Esko or Lili Milani (both at the Estonian Genome Centre and
   listed as authors on Yelmen *et al.*); courtesy email proposing
   co-authorship may smooth the path.
4. **EBIN ethics** — Estonia-specific; the SAC office can guide on form
   requirements. Plan ~4–6 weeks for ethics review.
5. **ORCID placeholder** — replace with your actual ORCID before
   submission.
6. **Funding statement** (Section B.12) needs filling.
7. **Author affiliations** assume Mao lab @ Umeå and Nie @ GDAAS as on
   paper #1; adjust if paper-#2 author list differs.
