"""Post-process multi-trait plink2 GWAS output to build per-trait summary
and lead-loci table.

Input: data/rice_3k/results/irri_gwas/gwas_irri.<TRAIT>.glm.linear  (one per trait)
Output:
  data/rice_3k/results/irri_gwas/irri_gwas_summary.tsv
  data/rice_3k/results/irri_gwas/irri_lead_loci.tsv
"""
from __future__ import annotations

import glob
from pathlib import Path

import pandas as pd

DATA = Path("/mnt/data/GraphGWAS/data/rice_3k")
GWAS_OUT = DATA / "results" / "irri_gwas"

GW_THRESHOLD = 5e-8
SUG_THRESHOLD = 1e-5
WINDOW_BP = 250_000


def cluster_leads(df: pd.DataFrame, p_thresh: float, window: int) -> list[dict]:
    d = df[df["P"] <= p_thresh].copy()
    if d.empty: return []
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
                "chrom": chrom, "pos": int(r["POS"]),
                "p": float(r["P"]), "beta": float(r["BETA"]),
                "ref": r["REF"], "alt": r["ALT"],
                "variant_id": f"{chrom}:{int(r['POS'])}:{r['REF']}:{r['ALT']}",
            })
    return leads


def main():
    files = sorted(glob.glob(str(GWAS_OUT / "gwas_irri.*.glm.linear")))
    print(f"Found {len(files)} per-trait glm.linear files")

    summary_rows = []
    all_leads = []

    for fp in files:
        trait = Path(fp).stem.replace("gwas_irri.", "").replace(".glm", "")
        # handle case where trait name has dots (unlikely)
        if trait.endswith(".glm"): trait = trait[:-4]
        print(f"\n=== {trait} ===", flush=True)
        df = pd.read_csv(fp, sep="\t", low_memory=False)
        df.columns = [c.lstrip("#") for c in df.columns]
        df = df[(df["TEST"]=="ADD") & (df["ERRCODE"]==".")].copy()
        for c in ("BETA","SE","P"):
            df[c] = pd.to_numeric(df[c], errors="coerce")
        df = df.dropna(subset=["P","BETA","SE"])
        df["CHROM"] = df["CHROM"].astype(str).apply(
            lambda c: c if c.startswith("Chr") else f"Chr{c}")
        n_var = len(df)
        if n_var == 0:
            print(f"  no valid test rows, skip")
            continue
        n_sug = int((df["P"] <= SUG_THRESHOLD).sum())
        n_gw = int((df["P"] <= GW_THRESHOLD).sum())
        leads_gw = cluster_leads(df, GW_THRESHOLD, WINDOW_BP)
        leads_sug = cluster_leads(df, SUG_THRESHOLD, WINDOW_BP)
        for l in leads_gw: l["trait"] = trait; l["sig_level"] = "GW"
        for l in leads_sug:
            if any(l2["trait"]==trait and l2["chrom"]==l["chrom"]
                   and abs(l2["pos"]-l["pos"])<WINDOW_BP for l2 in leads_gw):
                continue
            l["trait"] = trait; l["sig_level"] = "suggestive"
        all_leads.extend(leads_gw + [l for l in leads_sug
                                     if l.get("sig_level")=="suggestive"])
        top = df.loc[df["P"].idxmin()]
        summary_rows.append({
            "trait": trait, "n_variants": n_var,
            "n_suggestive_1e-5": n_sug, "n_gw_sig_5e-8": n_gw,
            "n_indep_gw_leads": len(leads_gw),
            "n_indep_sug_leads": len(leads_sug),
            "min_p": float(df["P"].min()),
            "top_variant": f"{top['CHROM']}:{int(top['POS'])}:{top['REF']}:{top['ALT']}",
        })
        print(f"  {n_var:,} variants; min_p={df['P'].min():.2e}; "
              f"{n_gw} GW-sig ({len(leads_gw)} indep); "
              f"{n_sug} suggestive ({len(leads_sug)} indep)", flush=True)

    summary = pd.DataFrame(summary_rows).sort_values("min_p")
    summary.to_csv(GWAS_OUT / "irri_gwas_summary.tsv", sep="\t", index=False)
    print(f"\nWrote {GWAS_OUT}/irri_gwas_summary.tsv")
    print(summary.to_string(index=False))

    if all_leads:
        leads_df = pd.DataFrame(all_leads)
        leads_df = leads_df.sort_values(["trait", "p"])
        leads_df.to_csv(GWAS_OUT / "irri_lead_loci.tsv", sep="\t", index=False)
        print(f"\nWrote {GWAS_OUT}/irri_lead_loci.tsv  ({len(leads_df)} leads total)")


if __name__ == "__main__":
    main()
