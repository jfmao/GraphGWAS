"""Postprocess plink2 GWAS for the 4 grain traits.

Inputs:
  data/rice_3k/results/grain_gwas/gwas_grain.<TRAIT>.glm.linear

Outputs:
  data/rice_3k/results/grain_gwas/grain_gwas_summary.tsv
  data/rice_3k/results/grain_gwas/grain_lead_loci.tsv      (clustered at 250 kb)
  data/rice_3k/results/grain_gwas/grain_lambda_gc.tsv      (genomic inflation per trait)
  data/rice_3k/results/grain_gwas/grain_niu2021_match.tsv  (direct chr:pos hits to 21 QTNs)
"""
from __future__ import annotations

import glob
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

DATA = Path("/mnt/data/GraphGWAS/data/rice_3k")
GWAS_OUT = DATA / "results" / "grain_gwas"
GT_QTNS = DATA / "ground_truth" / "niu2021_qtns.tsv"

GW_THRESHOLD = 5e-8
SUG_THRESHOLD = 1e-5
WINDOW_BP = 250_000
NIU_MATCH_WINDOW = 100_000  # call a Niu 2021 QTN matched if any lead within this distance

TRAIT_ORDER = ["TGW", "GL", "GW", "RLW"]


def cluster_leads(df: pd.DataFrame, p_thresh: float, window: int) -> list[dict]:
    d = df[df["P"] <= p_thresh].copy()
    if d.empty:
        return []
    d["CHROM"] = d["CHROM"].astype(str)
    d = d.sort_values(["CHROM", "P"])
    leads = []
    for chrom, sub in d.groupby("CHROM", sort=False):
        taken = []
        for _, r in sub.iterrows():
            if any(abs(r["POS"] - t) < window for t in taken):
                continue
            taken.append(r["POS"])
            leads.append({
                "chrom": chrom,
                "pos": int(r["POS"]),
                "p": float(r["P"]),
                "beta": float(r["BETA"]),
                "ref": r["REF"],
                "alt": r["ALT"],
                "variant_id": f"{chrom}:{int(r['POS'])}:{r['REF']}:{r['ALT']}",
            })
    return leads


def main():
    files = sorted(glob.glob(str(GWAS_OUT / "gwas_grain.*.glm.linear")))
    if not files:
        print(f"No glm.linear files at {GWAS_OUT}")
        return
    print(f"Found {len(files)} per-trait glm.linear files")

    summary_rows, all_leads, lambda_rows = [], [], []

    for fp in files:
        trait = Path(fp).stem.replace("gwas_grain.", "").replace(".glm", "")
        if trait.endswith(".glm"):
            trait = trait[:-4]
        print(f"\n=== {trait} ===", flush=True)
        df = pd.read_csv(fp, sep="\t", low_memory=False)
        df.columns = [c.lstrip("#") for c in df.columns]
        df = df[(df["TEST"] == "ADD") & (df["ERRCODE"] == ".")].copy()
        for c in ("BETA", "SE", "P"):
            df[c] = pd.to_numeric(df[c], errors="coerce")
        df = df.dropna(subset=["P", "BETA", "SE"])
        df["CHROM"] = df["CHROM"].astype(str).apply(
            lambda c: c if c.startswith("Chr") else f"Chr{c}"
        )
        n_var = len(df)
        if n_var == 0:
            print("  no valid test rows, skip")
            continue

        # Lambda_GC from chi^2_1 quantile
        chi2 = stats.norm.isf(df["P"].values / 2) ** 2
        lam = float(np.median(chi2) / stats.chi2.ppf(0.5, df=1))

        n_sug = int((df["P"] <= SUG_THRESHOLD).sum())
        n_gw = int((df["P"] <= GW_THRESHOLD).sum())
        leads_gw = cluster_leads(df, GW_THRESHOLD, WINDOW_BP)
        leads_sug = cluster_leads(df, SUG_THRESHOLD, WINDOW_BP)

        for l in leads_gw:
            l["trait"] = trait
            l["sig_level"] = "GW"
        for l in leads_sug:
            already_in_gw = any(
                l2["trait"] == trait
                and l2["chrom"] == l["chrom"]
                and abs(l2["pos"] - l["pos"]) < WINDOW_BP
                for l2 in leads_gw
            )
            if already_in_gw:
                continue
            l["trait"] = trait
            l["sig_level"] = "suggestive"
        all_leads.extend(leads_gw + [l for l in leads_sug if l.get("sig_level") == "suggestive"])

        top = df.loc[df["P"].idxmin()]
        summary_rows.append({
            "trait": trait,
            "N_variants": n_var,
            "lambda_gc": round(lam, 4),
            "n_suggestive_1e-5": n_sug,
            "n_gw_sig_5e-8": n_gw,
            "n_indep_gw_leads": len(leads_gw),
            "n_indep_sug_leads": len(leads_sug),
            "min_p": float(df["P"].min()),
            "top_variant": f"{top['CHROM']}:{int(top['POS'])}:{top['REF']}:{top['ALT']}",
        })
        lambda_rows.append({"trait": trait, "lambda_gc": round(lam, 4), "N_variants": n_var})
        print(
            f"  {n_var:,} variants; lambda_GC={lam:.3f}; min_p={df['P'].min():.2e}; "
            f"{n_gw} GW-sig ({len(leads_gw)} indep); "
            f"{n_sug} suggestive ({len(leads_sug)} indep)",
            flush=True,
        )

    # Order traits canonically
    summary = pd.DataFrame(summary_rows)
    if not summary.empty:
        summary["trait"] = pd.Categorical(summary["trait"], TRAIT_ORDER, ordered=True)
        summary = summary.sort_values("trait")
    summary.to_csv(GWAS_OUT / "grain_gwas_summary.tsv", sep="\t", index=False)
    print(f"\nWrote {GWAS_OUT}/grain_gwas_summary.tsv")
    print(summary.to_string(index=False))

    lam_df = pd.DataFrame(lambda_rows)
    lam_df.to_csv(GWAS_OUT / "grain_lambda_gc.tsv", sep="\t", index=False)

    if all_leads:
        leads_df = pd.DataFrame(all_leads).sort_values(["trait", "p"])
        leads_df.to_csv(GWAS_OUT / "grain_lead_loci.tsv", sep="\t", index=False)
        print(f"Wrote {GWAS_OUT}/grain_lead_loci.tsv ({len(leads_df)} leads total)")

        # ----- Direct match against Niu 2021 21 QTNs -----
        if GT_QTNS.exists():
            qtns = pd.read_csv(GT_QTNS, sep="\t")
            match_rows = []
            for _, q in qtns.iterrows():
                # restrict to leads of the same trait, on the same chromosome, within window
                mask = (
                    (leads_df["trait"] == q["trait"])
                    & (leads_df["chrom"] == q["chr"])
                    & (leads_df["pos"].between(q["pos"] - NIU_MATCH_WINDOW, q["pos"] + NIU_MATCH_WINDOW))
                )
                hits = leads_df[mask]
                if hits.empty:
                    match_rows.append({
                        "trait": q["trait"],
                        "qtn": q["qtn"],
                        "qtn_chr": q["chr"],
                        "qtn_pos": int(q["pos"]),
                        "matched": False,
                        "lead_chr": "",
                        "lead_pos": "",
                        "distance_bp": "",
                        "lead_p": "",
                        "lead_sig_level": "",
                        "known_gene": q["known_gene"] if isinstance(q["known_gene"], str) else "",
                        "candidate_gene_new": q["candidate_gene_new"] if isinstance(q["candidate_gene_new"], str) else "",
                    })
                    continue
                hits = hits.assign(distance=(hits["pos"] - q["pos"]).abs()).sort_values("distance")
                top = hits.iloc[0]
                match_rows.append({
                    "trait": q["trait"],
                    "qtn": q["qtn"],
                    "qtn_chr": q["chr"],
                    "qtn_pos": int(q["pos"]),
                    "matched": True,
                    "lead_chr": top["chrom"],
                    "lead_pos": int(top["pos"]),
                    "distance_bp": int(top["distance"]),
                    "lead_p": top["p"],
                    "lead_sig_level": top["sig_level"],
                    "known_gene": q["known_gene"] if isinstance(q["known_gene"], str) else "",
                    "candidate_gene_new": q["candidate_gene_new"] if isinstance(q["candidate_gene_new"], str) else "",
                })
            mdf = pd.DataFrame(match_rows)
            out = GWAS_OUT / "grain_niu2021_match.tsv"
            mdf.to_csv(out, sep="\t", index=False)
            n_match = int(mdf["matched"].sum())
            print(
                f"\nNiu 2021 QTN match (within ±{NIU_MATCH_WINDOW//1000} kb of any GWAS lead): "
                f"{n_match}/{len(mdf)} recovered"
            )
            print(f"Wrote {out}")
            print(mdf.to_string(index=False))


if __name__ == "__main__":
    main()
