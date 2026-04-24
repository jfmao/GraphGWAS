"""Single-run multi-phenotype plink2 --glm scan for all IRRI traits.

Much faster than N separate GWAS runs — plink2 shares genotype I/O
across phenotypes. Produces one <trait>.glm.linear per column in the
combined phenotype file.
"""
from __future__ import annotations

import glob
import subprocess
from pathlib import Path

import pandas as pd

DATA = Path("/mnt/data/GraphGWAS/data/rice_3k")
RES = DATA / "results"
PHENO_DIR = DATA / "pheno"
PGEN = DATA / "rice_3k"
PCA = DATA / "tmp" / "pca.eigenvec"
GWAS_OUT = RES / "irri_gwas"
GWAS_OUT.mkdir(parents=True, exist_ok=True)


def main():
    # Combine all single-trait pheno files into one wide file
    pheno_files = sorted(glob.glob(str(PHENO_DIR / "irri_*.pheno")))
    print(f"Found {len(pheno_files)} single-trait pheno files", flush=True)

    dfs = []
    for pf in pheno_files:
        d = pd.read_csv(pf, sep="\t")
        d.columns = [c.lstrip("#") for c in d.columns]
        dfs.append(d.set_index("IID"))
    combined = pd.concat(dfs, axis=1, join="outer").reset_index()
    combined = combined.rename(columns={"IID": "#IID"})
    out = PHENO_DIR / "irri_all_traits.pheno"
    combined.to_csv(out, sep="\t", index=False, na_rep="NA")
    trait_cols = [c for c in combined.columns if c != "#IID"]
    print(f"  {len(combined)} samples, {len(trait_cols)} traits", flush=True)
    print(f"  wrote {out}", flush=True)

    # Single plink2 call with all phenotypes
    out_prefix = GWAS_OUT / "gwas_irri"
    cmd = [
        "/home/jfmao/bin/plink2",
        "--pfile", str(PGEN),
        "--pheno", str(out),
        "--covar", str(PCA),
        "--covar-name", "PC1-PC10",
        "--glm", "hide-covar",
        "--out", str(out_prefix),
        "--threads", "8",
    ]
    print(f"\nRunning plink2 multi-phenotype GWAS:", flush=True)
    print(" ", " ".join(cmd), flush=True)
    r = subprocess.run(cmd, text=True)
    print(f"\nplink2 exit code: {r.returncode}", flush=True)

    # List outputs
    outs = sorted(GWAS_OUT.glob("gwas_irri.*.glm.linear"))
    print(f"\nProduced {len(outs)} GWAS output files:")
    for o in outs:
        size_mb = o.stat().st_size / 1e6
        print(f"  {o.name}  ({size_mb:.0f} MB)")


if __name__ == "__main__":
    main()
