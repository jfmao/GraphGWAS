"""Three-tier ground-truth recovery scoring for the 5-method grain fine-mapping run.

Tiers:
  1. Niu 2021 21 stable QTNs — variant-level recovery: is the QTN's chr:pos
     contained in the 95% credible set (or top-K PIP) for any locus of the
     matching trait? Strictest tier.
  2. Niu 2021 7 candidate genes (NEW) — locus-level recovery: does any
     fine-mapped lead for the matching trait fall within ±100 kb of the
     reported QTN coordinates that flagged the new gene? (Equivalent to
     "did the method localise to the right genomic neighbourhood".)
  3. Ren 2023 grain catalogue — gene-level: any catalogue gene of the matching
     biological category within ±250 kb of any fine-mapped lead variant
     (top PIP) for the same trait?

Outputs:
  data/rice_3k/results/grain_finemap/grain_recovery_scorecard.tsv
  data/rice_3k/results/grain_finemap/grain_recovery_summary.md
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

DATA = Path("/mnt/data/GraphGWAS/data/rice_3k")
RES = DATA / "results"
FM_OUT = RES / "grain_finemap"
GT = DATA / "ground_truth"

QTNS = pd.read_csv(GT / "niu2021_qtns.tsv", sep="\t")
REN = pd.read_csv(GT / "grain_quality_causal_genes.tsv", sep="\t")
GENE_IDX = pd.read_csv(GT / "gene_position_index.tsv", sep="\t")

METHODS = ["GAFM", "HBP", "SuSiE", "SuSiE-inf", "FINEMAP-inf", "SBayesRC"]
TRAITS = ["TGW", "GL", "GW", "RLW"]
TIER1_NEAR_WINDOW = 10_000   # ±bp around Niu QTN to count as "in CS (near)"
TIER2_WINDOW = 100_000   # ±bp around Niu 2021 QTN to call its candidate gene "recovered"
TIER3_WINDOW = 250_000   # ±bp around fine-map top PIP to call any Ren gene "recovered"


def _parse_pos(variant_id: str) -> tuple[str, int]:
    chr_, pos, _, _ = variant_id.split(":")
    return chr_, int(pos)


def load_locus_results():
    """Load every per-locus JSON written by rice3k_grain_shape_finemap.py."""
    rows = []
    for jf in sorted(FM_OUT.glob("*_Chr*_*.json")):
        d = json.load(open(jf))
        rows.append(d)
    return rows


def cs_contains_pos(cs_variants, target_chr, target_pos):
    for v in cs_variants:
        c, p, _, _ = v.split(":")
        if c == target_chr and int(p) == target_pos:
            return True
    return False


def cs_contains_window(cs_variants, target_chr, target_pos, window):
    for v in cs_variants:
        c, p, _, _ = v.split(":")
        if c == target_chr and abs(int(p) - target_pos) <= window:
            return True
    return False


def main():
    locus_results = load_locus_results()
    print(f"Loaded {len(locus_results)} per-locus result files")

    # Per-method runtime + CS-size summary
    timing_rows = []
    for r in locus_results:
        for method, mr in r["methods"].items():
            if mr.get("error"):
                continue
            timing_rows.append({
                "method": method,
                "trait": r["trait"],
                "lead_chr": r["lead_chr"],
                "lead_pos": r["lead_pos"],
                "n_variants": r["n_variants"],
                "cs_size": mr.get("cs_size"),
                "top_pip": mr.get("top_pip"),
                "time_s": mr.get("time_s"),
            })
    timing_df = pd.DataFrame(timing_rows)
    if not timing_df.empty:
        agg = timing_df.groupby("method").agg(
            n_loci=("trait", "count"),
            median_cs=("cs_size", "median"),
            mean_cs=("cs_size", "mean"),
            cs_eq_1=("cs_size", lambda s: int((s == 1).sum())),
            median_time_s=("time_s", "median"),
            max_time_s=("time_s", "max"),
        ).reindex(METHODS)
        print("\nPer-method runtime + CS-size summary:")
        print(agg.to_string())
        agg.to_csv(FM_OUT / "grain_method_runtime_cs.tsv", sep="\t")

    # ==== Tier 1: Niu 2021 21 QTNs ====
    tier1_rows = []
    for _, qtn in QTNS.iterrows():
        target_chr = qtn["chr"]
        target_pos = int(qtn["pos"])
        # Find all loci of the matching trait whose lead falls within TIER2_WINDOW
        # of the QTN — these are candidates that *might* contain it.
        loci = [
            r for r in locus_results
            if r["trait"] == qtn["trait"]
            and r["lead_chr"] == target_chr
            and abs(r["lead_pos"] - target_pos) <= TIER2_WINDOW
        ]
        for method in METHODS:
            in_cs_exact = False
            in_cs_near = False
            in_top1_exact = False
            in_top1_near = False
            top_pip_at_qtn = None
            min_distance = None
            for r in loci:
                mr = r["methods"].get(method, {})
                if mr.get("error"):
                    continue
                if mr.get("top_pip") is not None and mr.get("top_pip") < 0.01:
                    continue
                cs_v = mr.get("cs_variants") or []
                top_v = mr.get("top_variant")
                if cs_contains_pos(cs_v, target_chr, target_pos):
                    in_cs_exact = True
                if cs_contains_window(cs_v, target_chr, target_pos, TIER1_NEAR_WINDOW):
                    in_cs_near = True
                if top_v:
                    tc, tp = _parse_pos(top_v)
                    if tc == target_chr and tp == target_pos:
                        in_top1_exact = True
                    if tc == target_chr and abs(tp - target_pos) <= TIER1_NEAR_WINDOW:
                        in_top1_near = True
                for v in cs_v:
                    vc, vp = _parse_pos(v)
                    if vc == target_chr:
                        dist = abs(vp - target_pos)
                        if min_distance is None or dist < min_distance:
                            min_distance = dist
                            if dist == 0:
                                top_pip_at_qtn = mr.get("top_pip")
            tier1_rows.append({
                "tier": 1,
                "tier_name": "Niu 2021 QTN (variant)",
                "trait": qtn["trait"],
                "qtn": qtn["qtn"],
                "qtn_chr": target_chr,
                "qtn_pos": target_pos,
                "novel": bool(qtn["novel"]),
                "method": method,
                "in_credible_set": in_cs_exact,
                "in_credible_set_near": in_cs_near,
                "in_top1_pip": in_top1_exact,
                "in_top1_pip_near": in_top1_near,
                "min_dist_in_cs_bp": min_distance,
                "n_loci_in_window": len(loci),
            })

    tier1 = pd.DataFrame(tier1_rows)

    # Per-method tier 1 summary
    print("\n=== Tier 1 — Niu 2021 21 QTNs (variant-level CS membership) ===")
    pivot1 = tier1.groupby("method").agg(
        n_qtns=("qtn", "count"),
        in_cs_exact=("in_credible_set", "sum"),
        in_cs_near=("in_credible_set_near", "sum"),
        in_top1_exact=("in_top1_pip", "sum"),
        in_top1_near=("in_top1_pip_near", "sum"),
    ).reindex(METHODS)
    pivot1["recovery_cs_pct"] = (100 * pivot1["in_cs_near"] / pivot1["n_qtns"]).round(1)
    pivot1["recovery_top1_pct"] = (100 * pivot1["in_top1_near"] / pivot1["n_qtns"]).round(1)
    print(pivot1.to_string())

    # ==== Tier 2: Niu 2021 7 NEW candidate genes (locus-level recovery) ====
    tier2_rows = []
    novel_qtns = QTNS[QTNS["novel"]].copy()
    for _, qtn in novel_qtns.iterrows():
        target_chr = qtn["chr"]
        target_pos = int(qtn["pos"])
        candidate_genes = qtn["candidate_gene_new"]
        for method in METHODS:
            # locus-level criterion: any fine-mapped locus for this trait whose lead is
            # within TIER2_WINDOW of the QTN AND whose CS spans the QTN window
            recovered = False
            for r in locus_results:
                if r["trait"] != qtn["trait"]:
                    continue
                if r["lead_chr"] != target_chr:
                    continue
                if abs(r["lead_pos"] - target_pos) > TIER2_WINDOW:
                    continue
                mr = r["methods"].get(method, {})
                if mr.get("error"):
                    continue
                if mr.get("top_pip") is not None and mr.get("top_pip") < 0.01:
                    continue
                cs_v = mr.get("cs_variants") or []
                if cs_contains_window(cs_v, target_chr, target_pos, TIER2_WINDOW):
                    recovered = True
                    break
            tier2_rows.append({
                "tier": 2,
                "tier_name": "Niu 2021 NEW candidate gene (locus)",
                "trait": qtn["trait"],
                "qtn": qtn["qtn"],
                "qtn_chr": target_chr,
                "qtn_pos": target_pos,
                "candidate_genes": candidate_genes,
                "method": method,
                "recovered": recovered,
            })
    tier2 = pd.DataFrame(tier2_rows)

    print("\n=== Tier 2 — Niu 2021 7 NEW candidate genes (locus-level recovery) ===")
    pivot2 = tier2.groupby("method").agg(
        n=("qtn", "count"), recovered=("recovered", "sum")
    ).reindex(METHODS)
    pivot2["recovery_pct"] = (100 * pivot2["recovered"] / pivot2["n"]).round(1)
    print(pivot2.to_string())

    # ==== Tier 3: Ren 2023 catalogue (gene-level near top PIP) ====
    # Map Ren genes to coordinates via gene_position_index.
    GENE_IDX["LOC_Os_id"] = GENE_IDX["LOC_Os_id"].astype(str)
    GENE_IDX["chr"] = GENE_IDX["chr"].astype(str)
    GENE_IDX["pos"] = GENE_IDX["pos"].astype(int)
    REN["LOC_Os_id"] = REN["LOC_Os_id"].astype(str)
    # Ren has its own 'chr' column (semantic / chromosome notation may not match
    # gene_position_index format) — drop it before merging so we use the index's
    # canonical Chr1..Chr12 chromosome notation.
    ren_no_chr = REN.drop(columns=[c for c in ["chr", "pos"] if c in REN.columns])
    ren_with_pos = ren_no_chr.merge(GENE_IDX, on="LOC_Os_id", how="left")
    print(f"\nRen 2023 catalogue: {len(REN)} genes; {ren_with_pos['pos'].notna().sum()} mapped to coordinates")

    # Restrict to grain-shape relevant categories (everything in this catalog is grain)
    # but we filter to those whose 'trait_effect' string mentions length/width/weight/size.
    def matches_trait(row, trait):
        t = str(row.get("trait_effect", "")).lower()
        if trait == "TGW":
            return any(k in t for k in ["weight", "size", "yield"])
        if trait == "GL":
            return any(k in t for k in ["length", "size"])
        if trait == "GW":
            return any(k in t for k in ["width", "size"])
        if trait == "RLW":
            return any(k in t for k in ["length", "width", "shape", "size"])
        return False

    tier3_rows = []
    for trait in TRAITS:
        gene_subset = ren_with_pos[
            ren_with_pos.apply(matches_trait, axis=1, args=(trait,))
            & ren_with_pos["pos"].notna()
        ].copy()
        n_genes = len(gene_subset)
        for method in METHODS:
            # for each fine-mapped locus of the trait, count Ren genes within TIER3_WINDOW
            # of the top PIP variant
            recovered_genes = set()
            for r in locus_results:
                if r["trait"] != trait:
                    continue
                mr = r["methods"].get(method, {})
                if mr.get("error") or mr.get("top_variant") is None:
                    continue
                tc, tp = _parse_pos(mr["top_variant"])
                near = gene_subset[
                    (gene_subset["chr"] == tc)
                    & ((gene_subset["pos"] - tp).abs() <= TIER3_WINDOW)
                ]
                for _, g in near.iterrows():
                    recovered_genes.add(g["LOC_Os_id"])
            tier3_rows.append({
                "tier": 3,
                "tier_name": "Ren 2023 catalogue gene (gene)",
                "trait": trait,
                "n_eligible_genes": n_genes,
                "method": method,
                "n_recovered_genes": len(recovered_genes),
                "recovered_genes": ";".join(sorted(recovered_genes)),
            })
    tier3 = pd.DataFrame(tier3_rows)

    print("\n=== Tier 3 — Ren 2023 catalogue (gene within ±250 kb of any top PIP) ===")
    pivot3 = tier3.groupby("method").agg(
        recovered_total=("n_recovered_genes", "sum")
    ).reindex(METHODS)
    print(pivot3.to_string())

    # ==== Save the consolidated scorecard ====
    out_t1 = FM_OUT / "grain_recovery_tier1.tsv"
    out_t2 = FM_OUT / "grain_recovery_tier2.tsv"
    out_t3 = FM_OUT / "grain_recovery_tier3.tsv"
    tier1.to_csv(out_t1, sep="\t", index=False)
    tier2.to_csv(out_t2, sep="\t", index=False)
    tier3.to_csv(out_t3, sep="\t", index=False)

    # Master scorecard
    scorecard_rows = []
    for method in METHODS:
        t1 = pivot1.loc[method] if method in pivot1.index else None
        t2 = pivot2.loc[method] if method in pivot2.index else None
        t3 = pivot3.loc[method] if method in pivot3.index else None
        scorecard_rows.append({
            "method": method,
            "tier1_qtns_in_cs": int(t1["in_cs_near"]) if t1 is not None else None,
            "tier1_qtns_in_cs_exact": int(t1["in_cs_exact"]) if t1 is not None else None,
            "tier1_qtns_total": int(t1["n_qtns"]) if t1 is not None else None,
            "tier1_recovery_pct": float(t1["recovery_cs_pct"]) if t1 is not None else None,
            "tier1_top1_count": int(t1["in_top1_near"]) if t1 is not None else None,
            "tier2_novel_recovered": int(t2["recovered"]) if t2 is not None else None,
            "tier2_novel_total": int(t2["n"]) if t2 is not None else None,
            "tier2_recovery_pct": float(t2["recovery_pct"]) if t2 is not None else None,
            "tier3_ren_genes_recovered": int(t3["recovered_total"]) if t3 is not None else None,
        })
    scorecard = pd.DataFrame(scorecard_rows)
    scorecard.to_csv(FM_OUT / "grain_recovery_scorecard.tsv", sep="\t", index=False)
    print(f"\nWrote {FM_OUT}/grain_recovery_scorecard.tsv")
    print(scorecard.to_string(index=False))

    # Markdown summary
    md = []
    md.append("# Rice 3kRG grain weight + shape — fine-mapping recovery scorecard\n")
    md.append("Ground truth: Niu et al. 2021 *BMC Genomics* (21 stable QTNs, 7 NEW candidate genes) + Ren et al. 2023 *Sci Bull* (104 grain-size genes).\n")
    md.append("## Tier 1 — Niu 2021 21 stable QTNs (variant-level recovery in 95 % CS)\n")
    md.append("```\n" + pivot1.to_string() + "\n```\n")
    md.append("\n## Tier 2 — Niu 2021 7 NEW candidate genes (locus-level recovery)\n")
    md.append("```\n" + pivot2.to_string() + "\n```\n")
    md.append("\n## Tier 3 — Ren 2023 catalogue (gene within ±250 kb of fine-map top PIP)\n")
    md.append("```\n" + pivot3.to_string() + "\n```\n")
    md.append("\n## Method summary\n")
    md.append("```\n" + scorecard.to_string(index=False) + "\n```\n")
    md.append("\n")
    (FM_OUT / "grain_recovery_summary.md").write_text("\n".join(md))
    print(f"Wrote {FM_OUT}/grain_recovery_summary.md")


if __name__ == "__main__":
    main()
