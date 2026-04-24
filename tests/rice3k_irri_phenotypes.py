"""Build IRRI-trait phenotype files for 3kRG GWAS.

Maps the 3kRG_PhenotypeData_v20170411.xlsx (IRRI morphological codes,
2266 accessions x 54 traits) onto plink2 psam sample IDs via Wang 2018
metadata (IRGC accession cross-reference).

Selects the most-informative quantitative-like traits (those with >=5
distinct values and sufficient N), writes one pheno file per trait,
ready for plink2 --glm.

Outputs:
  data/rice_3k/pheno/irri_{TRAIT}.pheno  (plink2 #IID\tVALUE format)
  data/rice_3k/pheno/irri_phenotype_summary.tsv
"""
from __future__ import annotations

import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

DATA = Path("/mnt/data/GraphGWAS/data/rice_3k")
PHENO_DIR = DATA / "pheno"
PHENO_DIR.mkdir(parents=True, exist_ok=True)

PSAM = DATA / "rice_3k.psam"
PHENO_XLSX = Path("/mnt/data/GraphPop/data/raw/3kRG_data/3kRG_PhenotypeData_v20170411.xlsx")
META_XLSX = Path("/mnt/data/GraphPop/data/raw/3kRG_data/41586_2018_63_MOESM3_ESM.xlsx")

# Heuristic: pick traits with >=5 distinct integer values and <50% missing.
MIN_DISTINCT = 5
MAX_MISSING_FRAC = 0.5
MIN_N = 800


def parse_irgc_n(x):
    """Parse integer IRGC accession number from any form."""
    if pd.isna(x): return pd.NA
    m = re.search(r"(\d+)", str(x))
    return int(m.group(1)) if m else pd.NA


def main():
    # Load psam (3kRG sample IDs — B001, IRIS_313-XXXX, CX prefixes)
    psam = pd.read_csv(PSAM, sep="\t")
    psam.columns = [c.lstrip("#") for c in psam.columns]
    psam["IID"] = psam["IID"].astype(str)
    print(f"psam: {len(psam)} samples ({psam['IID'].iloc[0]} ... {psam['IID'].iloc[-1]})")

    # Load Wang 2018: IRIS_UNIQUE_ID (= psam IID) ↔ IRGC_Accno_source (int)
    meta = pd.read_excel(META_XLSX, sheet_name="Table 1 Metadata", skiprows=1)
    meta["IRGC_Accno_source_n"] = meta["IRGC_Accno_source"].apply(parse_irgc_n)
    meta = meta[["3K_DNA_IRIS_UNIQUE_ID", "IRGC_Accno_source_n"]].dropna()
    meta.columns = ["IID", "IRGC_n"]
    print(f"metadata: {len(meta)} B-code/IRIS ↔ IRGC pairs")

    # Load phenotype xlsx — use Source_Accno (int) as the IRGC number
    ph = pd.read_excel(PHENO_XLSX, sheet_name="Phenotype Data")
    ph["IRGC_n"] = ph["Source_Accno"].astype("Int64")
    ph = ph.dropna(subset=["IRGC_n"])
    print(f"phenotype: {len(ph)} rows with IRGC accessions")

    # Merge via IRGC_n → IID (psam)
    merged = ph.merge(meta, on="IRGC_n", how="inner")
    merged = merged.merge(psam[["IID"]], on="IID", how="inner")
    print(f"merged (in psam): {len(merged)} accessions with phenotypes")

    # Identify candidate trait columns (exclude identifiers)
    id_cols = {"Seqno","STOCK_ID","GS_ACCNO","NAME","Source_Accno","CROPYEAR",
               "IRGC_n","IID"}
    trait_cols = [c for c in merged.columns if c not in id_cols]

    # Filter to numeric traits with enough distinct values and N.
    # IRRI 1-9 ordinal scheme uses 999 as "missing / not assessed" — recode.
    rows = []
    for t in trait_cols:
        x = pd.to_numeric(merged[t], errors="coerce")
        # Recode IRRI missing-sentinel codes
        x = x.where((x != 999) & (x != 99), np.nan)
        if x.isna().all():
            continue
        n_valid = int(x.notna().sum())
        n_distinct = x.dropna().nunique()
        if n_valid < MIN_N or n_distinct < MIN_DISTINCT: continue
        miss = x.isna().mean()
        if miss > MAX_MISSING_FRAC: continue
        # Keep recoded column on merged for later dump
        merged[t] = x
        rows.append({
            "trait": t,
            "n_samples": n_valid,
            "n_distinct": n_distinct,
            "missing_frac": round(float(miss), 3),
            "mean": round(float(x.mean()), 3),
            "std": round(float(x.std()), 3),
            "min": float(x.min()),
            "max": float(x.max()),
        })
    summary = pd.DataFrame(rows).sort_values("n_distinct", ascending=False)
    print("\nQuantitative-like IRRI traits:")
    print(summary.to_string(index=False))

    summary.to_csv(PHENO_DIR / "irri_phenotype_summary.tsv", sep="\t", index=False)

    # Write one .pheno file per selected trait (plink2 #IID\tVALUE)
    for t in summary["trait"]:
        df = merged[["IID", t]].copy()
        df[t] = pd.to_numeric(df[t], errors="coerce")
        df = df.dropna()
        out = PHENO_DIR / f"irri_{t}.pheno"
        out_df = df.rename(columns={"IID": "#IID", t: t.upper()})
        out_df.to_csv(out, sep="\t", index=False)
        print(f"  wrote {out.name}  (N={len(df)})")

    print(f"\nTotal trait phenofiles: {len(summary)}")


if __name__ == "__main__":
    main()
