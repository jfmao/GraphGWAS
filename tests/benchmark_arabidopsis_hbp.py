"""Fine-map the top Arabidopsis FT10 locus with BGEN-backed HBP.

Demonstrates GraphGWAS HBP fine-mapping generalizes to plants: same code,
different species. No annotation graph loaded — pure statistical fine-mapping
on ±25 kb around the top GWAS peak.
"""
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

from graphgwas.bgen_reader import BgenReader
from graphgwas.finemapping_v2 import load_locus_variants, fast_hbp_finemap


OUT = Path("/mnt/data/GraphGWAS/results/arabidopsis")


def main():
    t0 = time.time()
    reader = BgenReader("/tmp/arabi", name_template="chr{chr}.bgen")

    # Pick top locus from prior GWAS output
    gwas_summary = json.loads((OUT / "chr4_ft10_summary.json").read_text())
    top_pos = int(gwas_summary["top_pos"])
    print(f"Top GWAS locus: chr4:{top_pos} (p={gwas_summary['top_p']:.2e})")

    # Load ±25 kb window via BGEN reader
    variants = load_locus_variants(
        chr="4", center=top_pos, window=50_000, source=reader
    )
    print(f"[{time.time()-t0:.1f}s] loaded {len(variants)} variants via BgenReader")

    # Align phenotype to BGEN sample order
    pheno_df = pd.read_csv(
        "/mnt/data/GraphGWAS/tests/data/arabidopsis/FT10.csv"
    ).drop_duplicates(subset=["accession_id"])
    pheno_df["accession_id"] = pheno_df["accession_id"].astype(str)
    pmap = dict(zip(pheno_df["accession_id"], pheno_df["phenotype_value"]))
    samples = [str(s) for s in reader.samples("4")]
    y_full = np.array([pmap.get(s, np.nan) for s in samples], dtype=float)
    valid = ~np.isnan(y_full)
    y = y_full[valid]

    # Impute missing genotypes and subset to valid samples
    kept = []
    for v in variants:
        d = v["dosage"][valid]
        # Impute NaN per variant
        col_mean = np.nanmean(d)
        if not np.isfinite(col_mean):
            continue
        d = np.where(np.isnan(d), col_mean, d)
        v["dosage"] = d
        af = float(col_mean / 2.0)
        maf = min(af, 1 - af)
        if maf < 0.02:
            continue
        v["af_total"] = af
        kept.append(v)
    variants = kept
    print(f"[{time.time()-t0:.1f}s] {len(variants)} variants after MAF>0.02 filter")
    if len(variants) < 3:
        print("Too few variants to fine-map")
        return

    # Run HBP fine-mapping (no annotation graph → pure stat-based)
    t1 = time.time()
    cands = fast_hbp_finemap(
        variants, y, graph_cache={}, n_rounds=5, chr_name="chr4"
    )
    print(f"[{time.time()-t0:.1f}s] HBP converged in {time.time()-t1:.3f}s")
    cands.sort(key=lambda c: -c.pip)

    # Write top 20 and summary
    tsv = OUT / "chr4_ft10_hbp_credset.tsv"
    with tsv.open("w") as fh:
        fh.write("rank\tvariant_id\tpos\tz_stat\tpip\tin_CS\n")
        for i, c in enumerate(cands[:30], 1):
            fh.write(f"{i}\t{c.variant_id}\t{c.pos}\t{c.z_stat:.3f}\t"
                     f"{c.pip:.4f}\t{c.in_credible_set}\n")

    cs_size = sum(1 for c in cands if c.in_credible_set)
    summary = {
        "locus_center": top_pos,
        "window_bp": 50_000,
        "n_variants_in_window": len(variants),
        "hbp_runtime_s": time.time() - t1,
        "credible_set_size": cs_size,
        "top_variant": cands[0].variant_id,
        "top_pip": float(cands[0].pip),
        "top_z": float(cands[0].z_stat),
    }
    with (OUT / "chr4_ft10_hbp_summary.json").open("w") as fh:
        json.dump(summary, fh, indent=2)

    print("\n=== HBP fine-mapping summary ===")
    print(json.dumps(summary, indent=2))
    print("\nTop 5 credible variants:")
    for c in cands[:5]:
        print(f"  chr4:{c.pos}  z={c.z_stat:+.2f}  pip={c.pip:.4f}  "
              f"in_CS={c.in_credible_set}")


if __name__ == "__main__":
    main()
