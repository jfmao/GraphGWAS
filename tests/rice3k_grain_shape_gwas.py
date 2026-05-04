"""Multi-phenotype plink2 --glm scan over the 4 grain weight + shape traits.

Single plink2 call shares genotype I/O across the 4 phenotypes — much faster
than 4 separate runs. Produces one <trait>.glm.linear file per column in the
combined phenotype file.

Reuses the same PCA covariate file as the IRRI scan (PC1-PC10 from
data/rice_3k/tmp/pca.eigenvec, computed in an earlier session).

Outputs:
  data/rice_3k/results/grain_gwas/gwas_grain.<TRAIT>.glm.linear  (4 files)
  data/rice_3k/results/grain_gwas/gwas_grain.log
"""
from __future__ import annotations

import subprocess
from pathlib import Path

DATA = Path("/mnt/data/GraphGWAS/data/rice_3k")
RES = DATA / "results"
PHENO_DIR = DATA / "pheno"
PGEN = DATA / "rice_3k"  # plink2 prefix; .pgen/.pvar/.psam
PCA = DATA / "tmp" / "pca.eigenvec"
PHENO_FILE = PHENO_DIR / "grain_all_traits.pheno"
GWAS_OUT = RES / "grain_gwas"
GWAS_OUT.mkdir(parents=True, exist_ok=True)

PLINK2 = "/home/jfmao/bin/plink2"


def main():
    assert PHENO_FILE.exists(), f"missing {PHENO_FILE}"
    assert PCA.exists(), f"missing {PCA}"
    out_prefix = GWAS_OUT / "gwas_grain"
    cmd = [
        PLINK2,
        "--pfile", str(PGEN),
        "--pheno", str(PHENO_FILE),
        "--covar", str(PCA),
        "--covar-name", "PC1-PC10",
        "--glm", "hide-covar",
        "--out", str(out_prefix),
        "--threads", "8",
    ]
    print("Running:", " ".join(cmd), flush=True)
    r = subprocess.run(cmd, text=True)
    print(f"\nplink2 exit code: {r.returncode}")

    outs = sorted(GWAS_OUT.glob("gwas_grain.*.glm.linear"))
    print(f"\nProduced {len(outs)} GWAS output files:")
    for o in outs:
        size_mb = o.stat().st_size / 1e6
        print(f"  {o.name}  ({size_mb:.0f} MB)")


if __name__ == "__main__":
    main()
