"""Run plink2 --glm linear for each IRRI phenotype on the 3kRG panel.

For each pheno file in data/rice_3k/pheno/irri_*.pheno, runs plink2 with
10 PCs as covariates and the 3kRG pgen. Writes per-trait summary and a
combined genome-wide-significant lead-loci table.
"""
from __future__ import annotations

import glob
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

DATA = Path("/mnt/data/GraphGWAS/data/rice_3k")
RES = DATA / "results"
PHENO_DIR = DATA / "pheno"
PGEN = DATA / "rice_3k"
PCA = DATA / "tmp" / "pca.eigenvec"
GWAS_OUT = RES / "irri_gwas"
GWAS_OUT.mkdir(parents=True, exist_ok=True)

GW_THRESHOLD = 5e-8
SUG_THRESHOLD = 1e-5
WINDOW_BP = 250_000


def run_gwas(pheno_file: Path, trait: str) -> Path | None:
    """Run plink2 --glm linear for this trait. Returns path to output file."""
    out_prefix = GWAS_OUT / f"gwas_{trait}"
    out_file = Path(f"{out_prefix}.{trait}.glm.linear")
    if out_file.exists():
        print(f"  [cache] {out_file.name}", flush=True)
        return out_file
    cmd = [
        "/home/jfmao/bin/plink2",
        "--pfile", str(PGEN),
        "--pheno", str(pheno_file),
        "--covar", str(PCA),
        "--covar-name", "PC1-PC10",
        "--glm", "hide-covar",
        "--out", str(out_prefix),
        "--threads", "8",
    ]
    print(f"  running plink2 on {trait} ...", flush=True)
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
    if r.returncode != 0:
        print(f"  [err] {trait}: plink2 exit {r.returncode}")
        print(r.stderr[-500:])
        return None
    if not out_file.exists():
        # plink sometimes uses slightly different case
        alt = list(GWAS_OUT.glob(f"gwas_{trait}.*.glm.linear"))
        if alt:
            out_file = alt[0]
    return out_file if out_file.exists() else None


def cluster_leads(df: pd.DataFrame, p_thresh: float, window: int) -> list[dict]:
    """Greedy lead-locus clustering within `window` bp on each chromosome."""
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
                "trait": None, "chrom": chrom, "pos": int(r["POS"]),
                "p": float(r["P"]), "beta": float(r["BETA"]),
                "ref": r["REF"], "alt": r["ALT"],
                "variant_id": f"{chrom}:{int(r['POS'])}:{r['REF']}:{r['ALT']}",
            })
    return leads


def main():
    pfiles = sorted(glob.glob(str(PHENO_DIR / "irri_*.pheno")))
    print(f"Found {len(pfiles)} IRRI pheno files\n")

    summary_rows = []
    all_leads = []

    for pf in pfiles:
        trait = Path(pf).stem.replace("irri_", "")
        print(f"\n=== {trait} ===", flush=True)
        gwas = run_gwas(Path(pf), trait)
        if gwas is None:
            print("  SKIP — GWAS failed"); continue
        df = pd.read_csv(gwas, sep="\t", low_memory=False)
        df.columns = [c.lstrip("#") for c in df.columns]
        df = df[(df["TEST"]=="ADD") & (df["ERRCODE"]==".")].copy()
        for c in ("BETA","SE","P"):
            df[c] = pd.to_numeric(df[c], errors="coerce")
        df = df.dropna(subset=["P","BETA","SE"])
        n_var = len(df)
        n_sug = int((df["P"] <= SUG_THRESHOLD).sum())
        n_gw = int((df["P"] <= GW_THRESHOLD).sum())
        # Cluster
        leads_gw = cluster_leads(df, GW_THRESHOLD, WINDOW_BP)
        leads_sug = cluster_leads(df, SUG_THRESHOLD, WINDOW_BP)
        for l in leads_gw: l["trait"] = trait; l["sig_level"] = "GW"
        for l in leads_sug:
            # Dedupe: skip sug leads already captured as GW
            if any(l2["trait"]==trait and l2["chrom"]==l["chrom"]
                   and abs(l2["pos"]-l["pos"])<WINDOW_BP for l2 in leads_gw):
                continue
            l["trait"] = trait; l["sig_level"] = "suggestive"
        all_leads.extend(leads_gw + [l for l in leads_sug
                                     if l.get("sig_level")=="suggestive"])
        min_p = float(df["P"].min())
        top = df.loc[df["P"].idxmin()]
        summary_rows.append({
            "trait": trait, "n_variants": n_var,
            "n_suggestive_5e-5": n_sug, "n_gw_sig_5e-8": n_gw,
            "n_indep_gw_leads": len(leads_gw),
            "n_indep_sug_leads": len(leads_sug),
            "min_p": min_p,
            "top_variant": f"{top['CHROM']}:{int(top['POS'])}:{top['REF']}:{top['ALT']}",
        })
        print(f"  {n_var:,} variants; min_p={min_p:.2e}; "
              f"{n_gw} GW-sig ({len(leads_gw)} indep); "
              f"{n_sug} suggestive ({len(leads_sug)} indep)", flush=True)

    # Write outputs
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
