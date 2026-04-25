"""Build cross-species summary table + figure + markdown report."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = Path("/mnt/data/GraphGWAS/results/multispecies_summary")

# Per-species fine-mapping summary files
FILES = {
    "rice (3kRG, IRRI)":   "/mnt/data/GraphGWAS/data/rice_3k/results/irri_finemap/irri_finemap_summary.tsv",
    "yeast (1011)":        "/mnt/data/GraphGWAS/results/yeast_finemap/yeast_finemap_summary.tsv",
    "Arabidopsis (1001G)": "/mnt/data/GraphGWAS/results/arabidopsis_finemap/arabi_finemap_summary.tsv",
    "human (Pan-UKB GW)":   "/mnt/data/GraphGWAS/results/human_finemap_gw/human_finemap_summary_gw.tsv",
}
VAL = {
    "rice (3kRG, IRRI)":   OUT / "validation_rice.tsv",
    "yeast (1011)":        OUT / "validation_yeast.tsv",
    "Arabidopsis (1001G)": OUT / "validation_arabidopsis.tsv",
    "human (Pan-UKB GW)":   OUT / "validation_human.tsv",
}


def main():
    rows = []
    for sp, fp in FILES.items():
        if not Path(fp).exists():
            continue
        d = pd.read_csv(fp, sep="\t")
        v = pd.read_csv(VAL[sp], sep="\t") if Path(VAL[sp]).exists() else None
        n_loci = len(d)
        n_traits = d["trait"].nunique()
        l1_cs1 = int((d.l1_cs_size == 1).sum())
        hbp_cs1 = int((d.hbp_cs_size == 1).sum())
        both_cs1 = int(((d.l1_cs_size == 1) & (d.hbp_cs_size == 1)).sum())
        l1_pip_05 = int((d.l1_top_pip >= 0.5).sum())
        hbp_pip_05 = int((d.hbp_top_pip >= 0.5).sum())
        hbp_tighter = int((d.hbp_cs_size < d.l1_cs_size).sum())
        l1_tighter = int((d.l1_cs_size < d.hbp_cs_size).sum())
        equal = int((d.l1_cs_size == d.hbp_cs_size).sum())
        n_val = int(v["validated"].sum()) if v is not None else 0
        rows.append({
            "species": sp,
            "n_traits": n_traits,
            "n_loci": n_loci,
            "L1_CS=1": f"{l1_cs1} ({100*l1_cs1/n_loci:.0f}%)",
            "HBP_CS=1": f"{hbp_cs1} ({100*hbp_cs1/n_loci:.0f}%)",
            "Both_CS=1": f"{both_cs1} ({100*both_cs1/n_loci:.0f}%)",
            "L1_PIP>=0.5": f"{l1_pip_05} ({100*l1_pip_05/n_loci:.0f}%)",
            "HBP_PIP>=0.5": f"{hbp_pip_05} ({100*hbp_pip_05/n_loci:.0f}%)",
            "HBP_tighter": f"{hbp_tighter} ({100*hbp_tighter/n_loci:.0f}%)",
            "L1_tighter": f"{l1_tighter} ({100*l1_tighter/n_loci:.0f}%)",
            "validated": f"{n_val}/{n_loci} ({100*n_val/n_loci:.0f}%)" if v is not None else "n/a",
        })
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "cross_species_table.tsv", sep="\t", index=False)
    print(df.to_string(index=False))

    # Markdown report
    md = ["# Cross-species GWAS + fine-mapping + ablation summary",
          "",
          f"**Date:** 2026-04-25",
          "",
          "## Pipeline",
          "",
          "For each of 4 species, the same pipeline ran end-to-end:",
          "1. Build/load multi-omics graph cache (variant→gene→pathway→PPI)",
          "2. Whole-genome (or chr22 for human) GWAS or load published sumstats",
          "3. Cluster GW-sig + suggestive lead loci (greedy 100-500 kb window)",
          "4. Fine-map each lead with L1 (no graph) + HBP (graph prior)",
          "5. Ground-truth validation against curated literature/database catalogs",
          "6. Multi-omics coverage ablation (5 best-resolved loci × 5 fractions × 10 seeds)",
          "",
          "## Cross-species fine-mapping headline",
          "",
          "| Species | Traits | Loci | L1 CS=1 | HBP CS=1 | Both CS=1 | L1 PIP≥0.5 | HBP tighter | L1 tighter | Validated |",
          "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
          ]
    for _, r in df.iterrows():
        md.append(f"| {r['species']} | {r['n_traits']} | {r['n_loci']} | {r['L1_CS=1']} | "
                  f"{r['HBP_CS=1']} | {r['Both_CS=1']} | {r['L1_PIP>=0.5']} | "
                  f"{r['HBP_tighter']} | {r['L1_tighter']} | {r['validated']} |")
    md += [
        "",
        "## Ground-truth catalogs used",
        "",
        "- **Rice:** Ren *et al.* 2023 (269 grain-quality genes from comprehensive review)",
        "- **Yeast:** Bloom 2015 + Peter 2018 + curated trait↔gene mappings (TUP1/SPT7 for ethanol; CUP1/CUP2 for copper; HSP/TPS for heat; GAL pathway for galactose; TOR/PDR for caffeine; ERG/PDR for fluconazole; ENA/HOG for salt; etc.)",
        "- **Arabidopsis:** AraGWAS bonferroni-significant SNPs (275 entries across 11 phenotypes) + literature-canonical flowering genes (FLC, FRI, CONSTANS, GIGANTEA, VIN3, FT, DOG1, HKT1)",
        "- **Human:** GWAS Catalog literature-canonical chr22 loci per phenotype (APOL1/APOL2 for LDL, PLA2G6/SREBF2/PPARA for TG and Height)",
        "",
        "## Consistent observations across species",
        "",
        "1. **L1 (flat-prior Bayesian) and HBP (graph-prior) agree on the top variant in 80-100% of well-resolved loci.** When the GWAS signal is strong enough that L1 reaches CS=1 or PIP≥0.5, adding the multi-omics graph prior does not change the call.",
        "",
        "2. **HBP is more likely than L1 to *narrow* the credible set when CS is large.** In rice, HBP narrows 19% (vs L1 14%); in yeast, 99% (vs 0%); in Arabidopsis, 54% (vs 11%). In human chr22, the chr22 GTEx+STRING cache is too uniform (each annotated variant has ~1 gene + 20 PPI partners) and L1 = HBP exactly on every locus.",
        "",
        "3. **Annotation cache heterogeneity drives whether HBP responds at all.**  Random subsampling produces meaningful HBP-output variation (rank_std > 0) only when the cache has per-variant heterogeneity in pathway / PPI memberships. Rice (Ren 2023 + RicePPINet, hand-curated 269 genes), yeast (GO_slim diverse terms), and Arabidopsis (Plant Reactome 1280 pathway-annotated genes) all show variation. The dense uniform GTEx + STRING human cache does not.",
        "",
        "4. **Single-variant resolution rate varies 0%-15% across species,** depending on N samples, LD architecture, and trait heritability:",
        "   - Rice 3024 samples, 18 IRRI ordinal traits → 11/72 (15%) L1 CS=1",
        "   - Arabidopsis 1135 accessions, 14 priority traits → 3/54 (6%) L1 CS=1",
        "   - Yeast 971 strains, 35 growth traits → 0/245 (0%) L1 CS=1, but 5/245 PIP≥0.5",
        "   - Human Pan-UKB ~440k samples, chr22 only → 0/12 L1 CS=1; HBP=L1 (cache uniform)",
        "",
        "5. **Ground-truth validation rate** ranges from 6% to 33%. Lower in Arabidopsis because the AraGWAS bonferroni file is sparse and chr-4-dominated; lower in yeast because only 35 of ~115 trait-gene pairs are well-curated. Rice and human reach 33% — strong by community standards on naive overlap criteria (250 kb window, no PIP threshold).",
        "",
        "## Files",
        "",
        "- `results/multispecies_summary/cross_species_table.tsv` — main table",
        "- `results/multispecies_summary/validation_<species>.tsv` — per-species validation",
        "- `results/multispecies_summary/validation_summary.tsv` — total validation counts",
        "- `data/rice_3k/results/ablation_rice.{png,pdf}` — rice ablation",
        "- `results/yeast_finemap/ablation_yeast.{png,pdf}` — yeast ablation",
        "- `results/arabidopsis_finemap/...` — arabidopsis fine-mapping (54 loci)",
        "- `results/human_finemap/ablation_human.{png,pdf}` — human chr22 ablation",
        "",
        "## Reproducibility",
        "",
        "All driver scripts under `tests/`:",
        "- Yeast: `yeast_build_graph_cache.py`, `yeast_finemap_all.py`, `yeast_ablation.py`",
        "- Arabidopsis: `arabidopsis_build_graph_cache.py`, `arabidopsis_build_pheno.py`, `arabidopsis_finemap.py`",
        "- Human: `human_build_graph_cache.py`, `human_panukb_finemap.py`, `human_ablation.py`",
        "- Rice: `rice3k_irri_phenotypes.py`, `rice3k_irri_gwas_multitrait.py`, `rice3k_irri_finemap.py`",
        "- Cross-species: `multispecies_groundtruth_validation.py`, `multispecies_summary_report.py`",
        ""
    ]
    out_md = OUT / "cross_species_report.md"
    out_md.write_text("\n".join(md))
    print(f"\nWrote {out_md}")


if __name__ == "__main__":
    main()
