"""Run FT10 GWAS on Arabidopsis chr4 via BgenReader + vectorized regression.

This demonstrates that GraphGWAS's hybrid BGEN architecture applies to plant
genomics without any code change. PLINK2 segfaults on this dataset; we use
Python numpy regression.
"""
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as sp_stats

from graphgwas.bgen_reader import BgenReader


OUT = Path("/mnt/data/GraphGWAS/results/arabidopsis")
OUT.mkdir(exist_ok=True)


def main(chr_name: str = "4", chunk_mb: int = 2):
    t0 = time.time()
    reader = BgenReader("/tmp/arabi", name_template="chr{chr}.bgen")

    # Build phenotype vector aligned to BGEN sample order
    samples = reader.samples(chr_name)
    sample_ids = [(s.decode() if isinstance(s, bytes) else s) for s in samples]
    pheno_df = pd.read_csv(
        "/mnt/data/GraphGWAS/tests/data/arabidopsis/FT10.csv"
    )
    pheno_df = pheno_df.drop_duplicates(subset=["accession_id"])
    pheno_df["accession_id"] = pheno_df["accession_id"].astype(str)
    pmap = dict(zip(pheno_df["accession_id"], pheno_df["phenotype_value"]))
    y_full = np.array([pmap.get(s, np.nan) for s in sample_ids], dtype=float)
    valid = ~np.isnan(y_full)
    y = y_full[valid]
    n = y.shape[0]
    y_centered = y - y.mean()
    y_ss = (y_centered ** 2).sum()
    print(f"[{time.time()-t0:.1f}s] phenotype: n={n} samples with FT10")

    # Chunked GWAS across chr4
    # Use reader's position array to find variants
    positions = reader._pos_array(chr_name)
    min_pos, max_pos = min(positions), max(positions)
    print(f"[{time.time()-t0:.1f}s] chr{chr_name}: {len(positions)} variants, "
          f"range {min_pos}-{max_pos}")

    chunk_bp = chunk_mb * 1_000_000
    n_chunks = (max_pos - min_pos) // chunk_bp + 1
    results = []
    total_tested = 0
    for i in range(n_chunks):
        start = min_pos + i * chunk_bp
        end = start + chunk_bp - 1
        variants, dosage = reader.load_locus(chr_name, start, end, format="dosage")
        if len(variants) == 0:
            continue
        # subset dosage to valid samples
        X = dosage[valid]  # (n, n_variants), may contain NaN (missing GT)
        # Mean-impute missing genotypes per variant, then MAF filter
        col_mean = np.nanmean(X, axis=0)
        nan_mask = np.isnan(X)
        if nan_mask.any():
            X = np.where(nan_mask, col_mean, X)
        af = col_mean / 2.0
        # `af` here = allele frequency of allele 0 (see bgen_reader); MAF is min(af, 1-af)
        maf = np.minimum(af, 1 - af)
        keep = (maf > 0.02) & np.isfinite(col_mean)
        X = X[:, keep]
        vsub = variants.loc[keep].reset_index(drop=True)
        af_keep = af[keep]
        if X.shape[1] == 0:
            continue

        # Vectorized univariate regression:
        # beta = cov(x, y) / var(x);  se = sqrt(rss / (n-2) / ss_x)
        X_c = X - X.mean(axis=0)
        ss_x = (X_c ** 2).sum(axis=0)
        ss_x = np.where(ss_x < 1e-10, np.nan, ss_x)
        beta = (X_c * y_centered[:, None]).sum(axis=0) / ss_x
        fitted = beta[None, :] * X_c
        resid = y_centered[:, None] - fitted
        rss = (resid ** 2).sum(axis=0)
        se = np.sqrt(rss / (n - 2) / ss_x)
        t_stat = beta / se
        pval = 2 * sp_stats.t.sf(np.abs(t_stat), df=n - 2)
        total_tested += X.shape[1]

        # Keep top hits (p < 1e-3) to avoid massive output
        sig = np.where(np.isfinite(pval) & (pval < 1e-3))[0]
        for j in sig:
            results.append(dict(
                chr=vsub.iloc[j]["chr"],
                pos=int(vsub.iloc[j]["pos"]),
                a1=vsub.iloc[j]["a1"],
                a2=vsub.iloc[j]["a2"],
                af=float(af_keep[j]),
                beta=float(beta[j]),
                se=float(se[j]),
                t=float(t_stat[j]),
                p=float(pval[j]),
                log10p=-float(np.log10(pval[j])) if pval[j] > 0 else 999,
            ))

        if i % 3 == 0 or i == n_chunks - 1:
            elapsed = time.time() - t0
            print(f"[{elapsed:.1f}s] chunk {i+1}/{n_chunks} "
                  f"({start//1_000_000}-{end//1_000_000} Mb) "
                  f"tested={total_tested:,} hits<1e-3={len(results)}")

    results.sort(key=lambda r: r["p"])
    # Save top hits
    out_tsv = OUT / f"chr{chr_name}_ft10_top_hits.tsv"
    with out_tsv.open("w") as fh:
        cols = ["chr", "pos", "a1", "a2", "af", "beta", "se", "t", "p", "log10p"]
        fh.write("\t".join(cols) + "\n")
        for r in results:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    summary = {
        "chr": chr_name,
        "n_samples": n,
        "n_variants_tested": total_tested,
        "n_hits_p_lt_1e-3": len(results),
        "n_hits_p_lt_1e-6": sum(1 for r in results if r["p"] < 1e-6),
        "n_hits_p_lt_5e-8": sum(1 for r in results if r["p"] < 5e-8),
        "top_pos": results[0]["pos"] if results else None,
        "top_p": results[0]["p"] if results else None,
        "top_log10p": results[0]["log10p"] if results else None,
        "runtime_s": time.time() - t0,
    }
    with (OUT / f"chr{chr_name}_ft10_summary.json").open("w") as fh:
        json.dump(summary, fh, indent=2)

    print(f"\n=== GWAS complete in {summary['runtime_s']:.1f}s ===")
    print(json.dumps(summary, indent=2))
    if results:
        print("\nTop 5 hits:")
        for r in results[:5]:
            print(f"  chr{r['chr']}:{r['pos']}  beta={r['beta']:.3f}  "
                  f"p={r['p']:.2e}  log10p={r['log10p']:.2f}")


if __name__ == "__main__":
    main()
