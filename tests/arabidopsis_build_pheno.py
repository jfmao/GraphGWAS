"""Build a multi-trait plink2 pheno file for the 18 priority Arabidopsis phenotypes."""
from __future__ import annotations

import pandas as pd
from pathlib import Path

DATA = Path("/mnt/data/GraphGWAS/data/arabidopsis")
TESTS = Path("/mnt/data/GraphGWAS/tests/data/arabidopsis")

# 18 priority pheno IDs from AraGWAS + named flowering traits
PRIORITY = pd.read_csv(DATA / "priority_phenotypes.tsv", sep="\t")
print(f"Loading {len(PRIORITY)} priority phenotypes")

# psam to ensure intersection
psam_path = DATA / "genotype" / "arabi_1001.psam"
if psam_path.exists():
    psam = pd.read_csv(psam_path, sep="\t")
    psam.columns = [c.lstrip("#") for c in psam.columns]
    sample_ids = set(psam["IID"].astype(str))
else:
    print("psam not yet built; using VCF accession list")
    sample_ids = None

dfs = []
for _, r in PRIORITY.iterrows():
    pid = r["phenotype_id"]
    pname = str(r["name"]).strip()
    fp = TESTS / "phenotypes" / f"pheno_{pid}.csv"
    if not fp.exists():
        print(f"  WARN: {fp} missing"); continue
    d = pd.read_csv(fp)
    d["accession_id"] = d["accession_id"].astype(str)
    d = d[["accession_id", "phenotype_value"]].copy()
    # safe column name
    safe = f"PH{pid}_{pname.replace(' ','').replace('/','').replace('-','')[:12]}"
    d.columns = ["IID", safe]
    if sample_ids is not None:
        d = d[d["IID"].isin(sample_ids)]
    d = d.drop_duplicates("IID")
    dfs.append(d.set_index("IID"))
    print(f"  {safe}: N={len(d)}")

combined = pd.concat(dfs, axis=1, join="outer").reset_index()
combined = combined.rename(columns={"IID":"#IID"})
out = DATA / "arabi_priority.pheno"
combined.to_csv(out, sep="\t", index=False, na_rep="NA")
print(f"\nWrote {out}: {len(combined)} samples × {len(combined.columns)-1} traits")
