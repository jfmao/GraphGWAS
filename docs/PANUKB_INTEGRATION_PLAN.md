# Pan-UKB Integration Plan

**Goal:** Replace the 1KG-simulated cross-ancestry experiment (manuscript §2.10) with a real-data cross-ancestry fine-mapping demonstration on Pan-UKB summary statistics, using Pan-UKB in-sample LD matrices.

**Scope:** limited — one new Results subsection (§2.10 replacement), one supplementary figure, and a sumstats-only entry path for L1 and HBP. *Not* a full UKB application; we consume only the public Pan-UKB release.

---

## Why Pan-UKB

- Preempts the reviewer question "why didn't you use biobank data?"
- Upgrades §2.10 from simulated F1 on 1KG subsamples (N = 500–661 per ancestry) to real UKB sumstats (N ≈ 2,700–420,000).
- No UKB application required; data is fully public.

## Pan-UKB resources (verified)

Public on Amazon S3 bucket `pan-ukb-us-east-1`, no authentication, no requester-pays:

| Resource | Path prefix | Per-ancestry size |
|---|---|---|
| Per-phenotype sumstats | `sumstats_release/*.tsv.bgz` | ~few GB per trait |
| Phenotype manifest | `sumstats_release/phenotype_manifest.tsv.bgz` | — |
| LD BlockMatrix (r², bias-adj) | `ld_release/UKBB.{ANC}.ldadj.bm/` | 2.8–15.5 TB |
| LD variant manifest (GRCh37) | `ld_release/UKBB.{ANC}.ldadj.variant.ht/` | — |
| LD variant manifest (GRCh38) | `ld_release/UKBB.{ANC}.ldadj.variant.b38.ht/` | — |
| LD scores | `ld_release/UKBB.{ANC}.ldscore.ht/` | — |

Ancestry sample sizes (from Pan-UKB manuscript):

| Ancestry | N | Notes |
|---|---:|---|
| EUR | ~420,000 | primary cohort |
| CSA | ~8,900 | |
| AFR | ~6,600 | |
| EAS | ~2,700 | |
| MID | ~1,600 | noisy LD likely |
| AMR | ~980 | too small for reliable fine-mapping |

**Decision:** report fine-mapping on the four ancestries with N > 2,500 (EUR, CSA, AFR, EAS). MID and AMR are noted but deprioritised.

## Loci selection

Pick 8–10 canonical multi-ancestry GWAS loci with well-characterised causal variants:

| Trait | Locus | Lead SNP (approx) | Why |
|---|---|---|---|
| LDL cholesterol | LDLR (19p13.2) | chr19:11,200,038 | classic monogenic, replicates across ancestries |
| Lipoprotein(a) | LPA (6q25.3) | chr6:161M | strong AFR-specific signal |
| Height | HMGA2 (12q14.3) | chr12:66M | European-discovered, replicates |
| BMI | FTO (16q12.2) | chr16:53.8M | landmark GWAS hit |
| T2D | TCF7L2 (10q25.2) | chr10:114M | ancestry-heterogeneous effect |
| HbA1c | G6PD (Xq28) | chrX:153M | X-linked, AFR signal |
| Triglycerides | APOA5 (11q23.3) | chr11:116.7M | strong EAS signal |
| Serum urate | SLC2A9 (4p16.1) | chr4:10M | multi-ancestry eQTL |

Final list to be confirmed against Pan-UKB phenotype manifest (some codes map differently).

## Manuscript impact

Restructure §2.10 as:

> **2.10 Cross-ancestry fine-mapping on Pan-UKB.** For each of 4 canonical loci × 4 ancestries (EUR, CSA, AFR, EAS), we fetched Pan-UKB sumstats, computed credible sets via HBP and L1 using Pan-UKB in-sample LD matrices, and compared to published fine-mapping maps. Results: …

Old §2.10 (1KG simulated cross-ancestry) becomes supplementary — it still demonstrates methodological robustness, but real UKB data is the headline.

## Engineering checkpoints

- [ ] **C1** — `src/python/graphgwas/panukb.py` with:
  - `PanUKBClient` with phenotype manifest reader
  - `fetch_sumstats_locus(phenocode, chr, start, end, pops)` → pandas DataFrame
  - `fetch_ld_slice(ancestry, chr, start, end)` → (R matrix, variant_ids) using Hail BlockMatrix slice
  - lazy Hail import (module usable without Hail installed)
  - per-ancestry sample-size constants
- [ ] **C2** — `l1_finemap_from_sumstats(z, R, variant_ids, graph_cache)` in `finemapping_v2.py`
- [ ] **C3** — `hbp_finemap_from_sumstats(z, R, variant_ids, graph_cache)` in `finemapping_v2.py`
- [ ] **C4** — `tests/test_sumstats_finemap.py` verifying the sumstats path agrees with the dosage path on a synthetic locus (max PIP difference < 1e-6)
- [ ] **C5** — Hail install + Pan-UKB BlockMatrix smoke test
- [ ] **C6** — Download phenotype manifest + sumstats for 8 candidate loci
- [ ] **C7** — Fetch LD slices for 4 ancestries × 8 loci
- [ ] **C8** — Run HBP/L1 fine-mapping on all (locus × ancestry) cells; aggregate PIPs + credible sets
- [ ] **C9** — New figure: cross-ancestry credible-set comparison panel
- [ ] **C10** — Write §2.10 replacement text; move old §2.10 to supplement

**Budget:** ~2–3 weeks elapsed time.
- C1–C4 (engineering): ~1 week
- C5–C7 (data + Hail): ~3–5 days
- C8–C10 (analysis + write-up): ~1 week

## Risks and fallbacks

| Risk | Probability | Fallback |
|---|---|---|
| Hail install fails on Python 3.13 | Medium | Use separate Python 3.11 env; or shell out to Hail subprocess; or use Alkes lab EUR-only LD (`s3://broad-alkesgroup-ukbb-ld`) |
| LD BlockMatrix slicing too slow | Low | Pre-extract ±2 Mb windows to local parquet; query cached |
| Non-EUR LD matrices noisy at modest N | Medium | Report LD quality metrics, restrict to EUR + CSA + AFR if EAS noise dominates |
| Loci disagree between populations | **This is the point** — fine-mapping should expose it | Report as finding, not failure |

## What stays unchanged

- §2.3 (weak-signal L1 vs SuSiE headline) — yeast
- §2.4 (Wu et al. baselines) — 1KG
- §2.6 (PIP calibration + null FPR) — yeast
- §2.7 (epistasis teaser) — unchanged
- §2.9 (power vs N) — 1KG
- §2.11 (Arabidopsis) — unchanged
- §2.12 (method portfolio) — unchanged

Pan-UKB is an extension, not a rewrite.
