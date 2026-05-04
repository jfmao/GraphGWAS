"""Convert plink2 .glm.linear sumstats to SBayesRC .ma format for the grain pass.

.ma format expected by SBayesRC:
    SNP A1 A2 freq b se p N

Restricts to variants present in the rice_3k_sbayes_subset.bim (the 188K
common variants in the 23 targeted blocks) so the .ma file matches the
genotype panel used for LD construction.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

DATA = Path("/mnt/data/GraphGWAS/data/rice_3k")
GWAS_OUT = DATA / "results" / "grain_gwas"
SBAYES_DIR = DATA / "sbayesrc_rice"
SUMSTATS_DIR = SBAYES_DIR / "sumstats"
SUMSTATS_DIR.mkdir(parents=True, exist_ok=True)

BIM = SBAYES_DIR / "rice_3k_sbayes_subset.bim"

bim = pd.read_csv(BIM, sep="\t", header=None,
                  names=["chr", "snp", "cm", "pos", "a1", "a2"])
keep_snps = set(bim["snp"].astype(str))
print(f"BIM: {len(keep_snps):,} variants in subset")

for trait in ["TGW", "GL", "GW", "RLW"]:
    fp = GWAS_OUT / f"gwas_grain.{trait}.glm.linear"
    df = pd.read_csv(fp, sep="\t", low_memory=False,
                     usecols=["#CHROM", "POS", "ID", "REF", "ALT", "A1",
                              "A1_FREQ", "BETA", "SE", "P", "OBS_CT",
                              "TEST", "ERRCODE"])
    df.columns = [c.lstrip("#") for c in df.columns]
    df = df[(df["TEST"] == "ADD") & (df["ERRCODE"] == ".")].copy()
    df["ID"] = df["ID"].astype(str)
    df = df[df["ID"].isin(keep_snps)]

    df["A2"] = df.apply(lambda r: r["REF"] if r["A1"] != r["REF"] else r["ALT"], axis=1)
    out = df[["ID", "A1", "A2", "A1_FREQ", "BETA", "SE", "P", "OBS_CT"]].copy()
    out.columns = ["SNP", "A1", "A2", "freq", "b", "se", "p", "N"]
    out = out.dropna()
    out["N"] = out["N"].astype(int)
    out_path = SUMSTATS_DIR / f"grain_{trait}.ma"
    out.to_csv(out_path, sep="\t", index=False)
    print(f"  wrote {out_path.name}  N_snps={len(out):,}")
