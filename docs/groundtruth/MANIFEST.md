# Functional ground-truth catalogue for paper #2 §Y.4 / §Y.6 benchmarks

**Compiled**: 2026-05-02
**Source**: 4 parallel literature searches (yeast, Arabidopsis, rice, human),
each via web + bioRxiv MCP + PubMed MCP.
**Purpose**: expand benchmark beyond the single anchor pair per species
that §Y.4 currently uses. Each pair was selected for **wet-lab functional
validation** (double mutant, CRISPR DKO, biochemistry, paired-NLR
cosegregation), not pure GWAS-only signal.

Machine-readable schema: `data/groundtruth/epistasis_pairs.json`

## Counts at a glance

| Species | Anchor | Additional pairs | With-natural-variation in panel | Negative controls |
|---|---|---|---|---|
| Yeast       | BCY1×TPK1   | 9  | 5 (Goldstein 2025-confirmed) | 0 |
| Arabidopsis | FT×FLC      | 8  | 6 (allelic in 1001G)         | 0 |
| Rice        | Hd1×Hd3a    | 9  | 9 (3K RG haplotyped)         | 1 (GW2×GW5 — independent pathways) |
| Human       | (none)      | 8  | 7 (varies in 1KG/UKBB)       | 0 |

## Tier-A pairs by species (most reliably testable)

### Yeast — 1011 Yeast Genomes
1. **BCY1 × TPK1** — anchor (cAMP-PKA reg:cat). Already in §Y.4.
2. **WHI5 × CLN3** — Goldstein 2025 confirmed segregating. G1/S commitment.
3. **MKT1 × hubs (HAP1, GPA1, IRA2, KRE33)** — pleiotropic mRNA-stability
   factor with multiple background-effect partners. Goldstein 2025.
4. **MIP1 × SAL1** — mitonuclear epistasis, validated in recombinant cross.
5. **TKL1 × NQM1** — pentose-phosphate paralog redundancy (CellMap-validated).

Lab-strain pairs (NOT 1011-detectable but useful as positive-control architecture):
MEC1×RAD9 (DDR), RFA1×RAD52 (HR), FUS3×KSS1 (MAPK), SNF1×REG1 (carbon switch), HSL1×SWE1 (morphogenesis).

### Arabidopsis — 1001 Genomes
1. **FT × FLC** — anchor.
2. **FRI × FLC** — vernalization, latitudinal cline. ~58% of accessions LoF for FRI.
3. **RRS1 × RPS4** — paired NLR (obligate heterodimer). P/A polymorphic.
4. **MAM1 × AOP2** — chr5×chr4 glucosinolate chemotype QTL. Strong haplotype variation.
5. **DOG1 × ABI3** — seed dormancy. DOG1 is the major dormancy QTL across 1001G.
6. **GL1 × GL3** — trichome MBW complex. GL1 LoF variants segregate.
7. **PHYB × PIF4** — light/temperature integration. PHYB has natural alleles, PIF4 cis-eQTL.
8. **SHR × SCR** — root patterning. Highly conserved (limited natural variation).
9. **BRI1 × BAK1** — BR co-receptor. Limited natural variation.

### Rice — 3K Rice Genomes
1. **Hd1 × Hd3a** — anchor.
2. **Hd1 × Ghd7** — heading date LD/SD. Both vary in 3K RG.
3. **Ghd7 × DTH8** — heading date NLD; jointly co-repress Ehd1. Both vary.
4. **GS3 × Gn1a** — grain length × grain number. CRISPR-validated.
5. **GS3 × OsSPL14/IPA1** — grain × architecture. CRISPR triple in ZH11.
6. **GW8 × GW7** — grain quality × yield. GW7 11-bp deletion is 3K RG haplotype.
7. **Pik-1 × Pik-2** — paired NLR (rice blast). Pikm/Pikp/Pikh haplotypes in 3K RG.
8. **RGA4 × RGA5** — paired NLR (Pia). Both polymorphic.
9. **SUB1A × SUB1C** — submergence tolerance. SUB1A presence/absence haplotype.
10. **GW2 × GW5** — *negative control*: independent pathways, no expected epistasis.

### Human — 1KG / UKBB / Pan-UKB
*Field caveat: human epistasis literature is sparse and contested. Hemani 2014
trans-eQTL pairs largely retracted. Most long-range LD pair hits from
eMERGE/UKBB 2023 are statistically marginal with no functional follow-up.
Below pairs survived independent functional validation.*

1. **HLA-B × ERAP1** — ankylosing spondylitis. Allele-specific peptidome.
2. **TPMT × NUDT15** — thiopurine toxicity. Prospective dose-titration RCT.
3. **VKORC1 × CYP2C9** — warfarin dose. EU-PACT RCT.
4. **HLA-B Bw4 × KIR3DS1** — HIV control. NK functional assays.
5. **BRCA1/2 × PARP1** — synthetic lethal cancer therapy. CRISPR + clinical RCTs.
6. **HFE C282Y × H63D** — hereditary hemochromatosis (compound het, same gene).
7. **HLA-DQA1 × HLA-DQB1** — celiac (DQ2.5 heterodimer).
8. **HLA-DRB1*15:01 × *03:01** — T1D (DR3/DR4 trans-dimer).

For 1KG **chr21+22 simulation** specifically: only CYP2D6 (chr22) ×
CYP3A4 (chr7) is on-chr22. Most canonical human pairs involve HLA (chr6).
Recommend simulating *functional-architecture-matched* epistasis (paralog
pairs, enhancer×promoter pairs from ExP STARR-seq) rather than insisting
on real chr21+22 hits.

## Action items for paper #2 §Y.4 / §Y.6 expansion

1. Extend `tests/validate_m2_yeast.py` ground-truth list from 1 pair to 5
   Tier-A 1011-segregating pairs (BCY1-TPK1, WHI5-CLN3, MKT1+hubs,
   MIP1-SAL1, TKL1-NQM1). Per-pair recovery rank in 5-method × 3-scenario × 5-pair grid.
2. Extend cross-species runner to include the Tier-A pair set per species.
3. Add the rice **GW2×GW5 negative control** — methods should NOT detect
   epistasis there. Report null-recovery rate as FPR proxy.
4. For human chr21+22 simulation: build 5-10 synthetic-architecture-matched
   epistasis ground-truth from same-gene (HFE-style compound het), paralog
   (CYP family), and ExP-style enhancer×promoter pairs that map to chr21/22.

## References

Yeast: [Goldstein 2025](https://academic.oup.com/genetics/article/231/2/iyaf136/8206093) ·
[De Chiara 2025](https://www.nature.com/articles/s41586-025-09637-0) ·
[Costanzo CellMap](https://thecellmap.org/about/)

Arabidopsis: [Caicedo 2004](https://pmc.ncbi.nlm.nih.gov/articles/PMC524852/) ·
[Saucet 2015](https://www.nature.com/articles/ncomms7338) ·
[Long 2024 SHR×SCR](https://www.nature.com/articles/s41586-023-06971-z) ·
[Brachi/Roux 2022](https://elifesciences.org/articles/67784)

Rice: [Sun 2022 Hd1/Ghd7/DTH8](https://doi.org/10.1016/j.jgg.2022.02.018) ·
[Wang 2018 3K RG](https://www.nature.com/articles/s41586-018-0030-5) ·
[Białas 2022 Pik-1/Pik-2](https://doi.org/10.1073/pnas.2116896119) ·
[Singh 2020 SUB1A](https://doi.org/10.1038/s41598-020-65588-8)