"""Build grain weight + grain shape phenotype files for the 3kRG fine-mapping rerun.

Reads the two new phenotype sources prepared in docs/:
  - docs/3K.shape.phe.xlsx       (2,453 accessions; columns: 3K_ID, GLWR, GL, GW)
  - docs/thousand_grain_weight.txt  (1,847 accessions; tab-separated: ID, TGW)

Both files use the plink2 psam IID convention directly (B/C/IRIS_313 prefixes),
no IRGC bridge needed (unlike the IRRI ordinal pheno path).

Renames GLWR → RLW for paper alignment with Niu et al. 2021 BMC Genomics, where
the four traits are TGW (g), GL (mm), GW (mm), RLW (= GL/GW ratio).

Outputs:
  data/rice_3k/pheno/grain_TGW.pheno          (one trait per file, plink2 #IID/VALUE)
  data/rice_3k/pheno/grain_GL.pheno
  data/rice_3k/pheno/grain_GW.pheno
  data/rice_3k/pheno/grain_RLW.pheno
  data/rice_3k/pheno/grain_all_traits.pheno   (combined; NA-padded)
  data/rice_3k/pheno/grain_pheno_qc.tsv       (per-trait N, mean, SD, range, vs Niu 2021)
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path("/mnt/data/GraphGWAS")
DATA = REPO / "data" / "rice_3k"
PHENO_DIR = DATA / "pheno"
PHENO_DIR.mkdir(parents=True, exist_ok=True)

PSAM = DATA / "rice_3k.psam"
SHAPE_XLSX = REPO / "docs" / "3K.shape.phe.xlsx"
TGW_TXT = REPO / "docs" / "thousand_grain_weight.txt"

# Niu 2021 Table 1 — Whole-population variance components and h^2.
# We use this for an honest sanity check only; we do not claim to reproduce
# the variance components (their model includes year and replicate).
NIU2021_REFERENCE = {
    "TGW": {"h2": 0.92, "approx_mean_g": 25.0},
    "GL":  {"h2": 0.88, "approx_mean_mm": 8.0},
    "GW":  {"h2": 0.93, "approx_mean_mm": 3.1},
    "RLW": {"h2": 0.92, "approx_mean": 2.7},
}


def main():
    psam = pd.read_csv(PSAM, sep="\t")
    psam.columns = [c.lstrip("#") for c in psam.columns]
    psam_ids = set(psam["IID"].astype(str))
    print(f"psam: {len(psam_ids)} samples", flush=True)

    # Shape file: 3K_ID / GLWR / GL / GW
    shape = pd.read_excel(SHAPE_XLSX)
    assert set(shape.columns) >= {"3K_ID", "GLWR", "GL", "GW"}, shape.columns
    shape = shape.rename(columns={"3K_ID": "IID", "GLWR": "RLW"})
    shape = shape[["IID", "RLW", "GL", "GW"]].copy()
    print(f"shape file: {len(shape)} rows", flush=True)

    # TGW file: tab-separated, no header
    tgw = pd.read_csv(TGW_TXT, sep="\t", header=None, names=["IID", "TGW"])
    tgw["IID"] = tgw["IID"].astype(str)
    print(f"tgw   file: {len(tgw)} rows", flush=True)

    # Filter to psam-resolvable IIDs
    shape_in = shape[shape["IID"].isin(psam_ids)].copy()
    tgw_in = tgw[tgw["IID"].isin(psam_ids)].copy()
    print(f"shape ∩ psam: {len(shape_in)}/{len(shape)}", flush=True)
    print(f"tgw   ∩ psam: {len(tgw_in)}/{len(tgw)}", flush=True)

    # Build each per-trait pheno file
    qc_rows = []
    for trait, src in [("TGW", tgw_in), ("GL", shape_in), ("GW", shape_in), ("RLW", shape_in)]:
        df = src[["IID", trait]].dropna().copy()
        df[trait] = pd.to_numeric(df[trait], errors="coerce")
        df = df.dropna()
        out = PHENO_DIR / f"grain_{trait}.pheno"
        out_df = df.rename(columns={"IID": "#IID"})
        out_df.to_csv(out, sep="\t", index=False)
        qc_rows.append({
            "trait": trait,
            "N": len(df),
            "mean": round(float(df[trait].mean()), 4),
            "sd": round(float(df[trait].std()), 4),
            "min": round(float(df[trait].min()), 4),
            "max": round(float(df[trait].max()), 4),
            "niu2021_h2": NIU2021_REFERENCE[trait]["h2"],
        })
        print(f"  wrote {out.name}  N={len(df)}", flush=True)

    qc = pd.DataFrame(qc_rows)
    qc.to_csv(PHENO_DIR / "grain_pheno_qc.tsv", sep="\t", index=False)
    print("\nPhenotype QC:")
    print(qc.to_string(index=False))

    # Combined wide pheno file (one row per IID; NA where any trait missing).
    # Need outer-join on IID across all four traits.
    parts = []
    for trait in ["TGW", "GL", "GW", "RLW"]:
        d = pd.read_csv(PHENO_DIR / f"grain_{trait}.pheno", sep="\t")
        d.columns = [c.lstrip("#") for c in d.columns]
        parts.append(d.set_index("IID"))
    combined = pd.concat(parts, axis=1, join="outer").reset_index()
    combined = combined.rename(columns={"IID": "#IID"})
    out = PHENO_DIR / "grain_all_traits.pheno"
    combined.to_csv(out, sep="\t", index=False, na_rep="NA")
    print(f"\nCombined: {len(combined)} rows × {combined.shape[1] - 1} traits → {out}")

    # Pairwise N (samples with both traits non-missing)
    print("\nPairwise N (samples with both traits non-missing):")
    cmat = combined.set_index("#IID")
    pairs = pd.DataFrame(
        {a: {b: int((cmat[a].notna() & cmat[b].notna()).sum())
             for b in cmat.columns} for a in cmat.columns}
    )
    print(pairs.to_string())

    # Pairwise correlation (descriptive)
    print("\nPairwise Pearson r:")
    print(cmat.corr(method="pearson").round(3).to_string())


if __name__ == "__main__":
    main()
