"""Extract SBayesRC PIPs at the 41 fine-mapped loci and merge with the 5-method
recovery scorecard to produce a 6-method final.

Inputs:
  data/rice_3k/sbayesrc_rice/results/sbayesrc_<TRAIT>.txt
      (SNP A1 BETA SE PIP BETAlast — per-variant posteriors)
  data/rice_3k/sbayesrc_rice/rice_3k_sbayes_subset.bim
      (chromosome map for variant ID → chr/pos)
  data/rice_3k/results/grain_finemap/*.json
      (existing 5-method per-locus output)

Outputs:
  data/rice_3k/results/grain_finemap/<TRAIT>_<chr>_<pos>.json   (updated in place
      with a new "SBayesRC" entry under r["methods"])
  data/rice_3k/results/grain_finemap/grain_finemap_summary.tsv  (extended)
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

DATA = Path("/mnt/data/GraphGWAS/data/rice_3k")
SBAYES = DATA / "sbayesrc_rice"
FM_OUT = DATA / "results" / "grain_finemap"
COVERAGE = 0.95
TRAITS = ["TGW", "GL", "GW", "RLW"]
WINDOW_BP = 100_000


def load_chr_pos_map():
    bim = pd.read_csv(SBAYES / "rice_3k_sbayes_subset.bim", sep="\t", header=None,
                      names=["chr", "snp", "cm", "pos", "a1", "a2"])
    bim["snp"] = bim["snp"].astype(str)
    return bim


def main():
    bim = load_chr_pos_map()
    snp2chrpos = dict(zip(bim["snp"].astype(str),
                           zip("Chr" + bim["chr"].astype(str), bim["pos"].astype(int))))

    # Per-trait: load SBayesRC PIPs and join with chr/pos
    sb = {}
    # v0.1.5: prefer the deflated SBayesRC outputs in results_v15/ if they exist
    # (built by tests/rice3k_grain_shape_v15_sbayesrc.R), else fall back to
    # results/ from v0.1.4.
    for trait in TRAITS:
        fp_v15 = SBAYES / "results_v15" / f"sbayesrc_{trait}.txt"
        fp_v14 = SBAYES / "results" / f"sbayesrc_{trait}.txt"
        fp = fp_v15 if fp_v15.exists() else fp_v14
        if not fp.exists():
            print(f"  warning: {fp.name} not found")
            continue
        df = pd.read_csv(fp, sep="\t")
        df["SNP"] = df["SNP"].astype(str)
        df["chr"] = df["SNP"].map(lambda s: snp2chrpos.get(s, (None, None))[0])
        df["pos"] = df["SNP"].map(lambda s: snp2chrpos.get(s, (None, None))[1])
        df = df.dropna(subset=["chr", "pos"])
        df["pos"] = df["pos"].astype(int)
        sb[trait] = df
        print(f"  loaded {trait}: {len(df):,} variants with PIPs")

    # Per-locus update: for each existing JSON, add SBayesRC method block
    summary_rows = []
    for jf in sorted(FM_OUT.glob("*_Chr*_*.json")):
        d = json.load(open(jf))
        trait = d["trait"]
        chrom = d["lead_chr"]
        pos = d["lead_pos"]
        if trait not in sb:
            continue
        df = sb[trait]
        win = df[(df["chr"] == chrom) & (df["pos"].between(pos - WINDOW_BP, pos + WINDOW_BP))].copy()
        if win.empty:
            d["methods"]["SBayesRC"] = {
                "method": "SBayesRC", "error": True,
                "error_msg": "no variants in window",
                "cs_size": None, "top_pip": None, "top_variant": None,
                "cs_variants": [], "top_variants": [], "time_s": None,
                "n_credible_sets": 0, "cs_size_aggregated": None,
            }
        else:
            win = win.sort_values("PIP", ascending=False).reset_index(drop=True)
            cum = np.cumsum(win["PIP"].values)
            total = float(cum[-1])
            if total <= 1e-6:
                k = 1
            else:
                target = COVERAGE * total
                k = int(np.searchsorted(cum, target)) + 1
            cs = win.iloc[:k]
            top = win.iloc[0]
            cs_variants = [f"{r['chr']}:{int(r['pos'])}:?:?" for _, r in cs.iterrows()]
            d["methods"]["SBayesRC"] = {
                "method": "SBayesRC",
                "n_credible_sets": 1,
                "cs_size": int(k),
                "cs_size_aggregated": int(k),
                "top_pip": float(top["PIP"]),
                "top_variant": f"{top['chr']}:{int(top['pos'])}:?:?",
                "cs_variants": cs_variants,
                "top_variants": [
                    (f"{r['chr']}:{int(r['pos'])}:?:?", float(r["PIP"]))
                    for _, r in win.head(5).iterrows()
                ],
                "time_s": None,  # SBayesRC is genome-wide; no per-locus time
                "error": False,
            }
        with open(jf, "w") as f:
            json.dump(d, f, indent=2, default=str)

        # Append summary row for the new method
        mr = d["methods"]["SBayesRC"]
        summary_rows.append({
            "trait": trait,
            "lead_chr": chrom,
            "lead_pos": pos,
            "sig_level": d["sig_level"],
            "lead_p": d["lead_p"],
            "n_variants": d["n_variants"],
            "cache_coverage": d.get("cache_coverage", 0),
            "n_samples": d["n_samples"],
            "method": "SBayesRC",
            "n_credible_sets": mr.get("n_credible_sets"),
            "cs_size": mr.get("cs_size"),
            "cs_size_aggregated": mr.get("cs_size_aggregated"),
            "top_pip": mr.get("top_pip"),
            "top_variant": mr.get("top_variant"),
            "time_s": mr.get("time_s"),
            "error": mr.get("error"),
        })

    # Append SBayesRC rows to grain_finemap_summary.tsv
    sb_df = pd.DataFrame(summary_rows)
    summ_path = FM_OUT / "grain_finemap_summary.tsv"
    if summ_path.exists():
        existing = pd.read_csv(summ_path, sep="\t")
        existing = existing[existing["method"] != "SBayesRC"]
        out = pd.concat([existing, sb_df], ignore_index=True)
        out.to_csv(summ_path, sep="\t", index=False)
        print(f"\nWrote {summ_path} ({len(out)} rows total; +{len(sb_df)} SBayesRC)")
    else:
        sb_df.to_csv(summ_path, sep="\t", index=False)


if __name__ == "__main__":
    main()
