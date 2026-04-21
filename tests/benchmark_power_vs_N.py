"""Power vs sample-size curve: HBP and L1 on 1KG chr22 at N ∈ {500, 1000, 2000, 3000}.

Simulates F1-style single-causal-in-LD on random chr22 50 kb windows using the
BGEN reader, subsamples BGEN samples at each N, runs HBP fine-mapping, and
reports rank-1 rate and mean causal rank.

Writes results/benchmark_v2/power_vs_N/power_curve.{json,tsv}.
"""
from __future__ import annotations

import bisect
import json
import time
from pathlib import Path

import numpy as np

from graphgwas.bgen_reader import BgenReader
from graphgwas.finemapping_v2 import fast_hbp_finemap

OUT = Path("/mnt/data/GraphGWAS/results/benchmark_v2/power_vs_N")
OUT.mkdir(parents=True, exist_ok=True)


def simulate_f1(
    reader: BgenReader,
    rng: np.random.Generator,
    window_bp: int = 50_000,
    h2: float = 0.10,
    n_attempts: int = 10,
) -> tuple[list[dict], int, np.ndarray, np.ndarray]:
    """Pick a random chr22 50 kb window with a reasonable causal variant,
    simulate y = β*G_causal + ε with heritability h². Returns
    (variants_list, causal_idx, dosage_matrix_[N×M], y).
    """
    positions = reader._pos_array("22")
    pos_min, pos_max = positions[0] + window_bp, positions[-1] - window_bp
    for _ in range(n_attempts):
        center = rng.integers(pos_min, pos_max)
        lo = bisect.bisect_left(positions, center - window_bp // 2)
        hi = bisect.bisect_right(positions, center + window_bp // 2)
        if hi - lo < 50:
            continue
        vdf, dosage = reader.load_locus("22", center - window_bp // 2, center + window_bp // 2, format="dosage")
        # Filter MAF > 0.05 and impute NaN
        col_mean = np.nanmean(dosage, axis=0)
        dosage = np.where(np.isnan(dosage), col_mean, dosage)
        af = col_mean / 2.0
        maf = np.minimum(af, 1 - af)
        keep = (maf > 0.05) & np.isfinite(col_mean)
        if keep.sum() < 30:
            continue
        dosage = dosage[:, keep]
        vdf = vdf.loc[keep].reset_index(drop=True)
        af_kept = af[keep]
        # Pick a causal variant (prefer AF 10-40%)
        candidates = np.where((af_kept > 0.1) & (af_kept < 0.4))[0]
        if len(candidates) == 0:
            continue
        causal_idx = int(rng.choice(candidates))
        g = dosage[:, causal_idx]
        g_c = g - g.mean()
        var_g = g_c.var()
        # β so that var(β*G) = h² * var(y) and var(ε) = (1-h²) * var(y), with var(y)=1
        beta = np.sqrt(h2 / var_g) if var_g > 0 else 0
        eps = rng.standard_normal(dosage.shape[0]) * np.sqrt(1 - h2)
        y = beta * g_c + eps
        # Pack variants into dict list (chr/pos/af_total are required by fast_hbp_finemap)
        variants = []
        for j in range(dosage.shape[1]):
            row = vdf.iloc[j]
            variants.append({
                "variantId": f"{row['chr']}:{row['pos']}:{row['a1']}:{row['a2']}",
                "chr": row["chr"],
                "pos": int(row["pos"]),
                "ref": row["a1"],
                "alt": row["a2"],
                "af_total": float(af_kept[j]),
                # dosage is set per-sample by the caller after subsampling
            })
        return variants, causal_idx, dosage, y
    raise RuntimeError("Failed to find a suitable locus after attempts")


def run_one_rep(reader: BgenReader, N: int, h2: float, seed: int) -> dict:
    rng = np.random.default_rng(seed)
    variants, causal_idx, dosage, y = simulate_f1(reader, rng, h2=h2)
    # Subsample to N samples
    n_available = dosage.shape[0]
    idx = rng.choice(n_available, size=min(N, n_available), replace=False)
    y_sub = y[idx]
    dosage_sub = dosage[idx, :]
    # Build variants-with-dosage list
    for j, v in enumerate(variants):
        v["dosage"] = dosage_sub[:, j]
    t0 = time.time()
    cands = fast_hbp_finemap(variants, y_sub, graph_cache={}, n_rounds=3, chr_name="chr22")
    dt = time.time() - t0
    # Locate causal variant in cands (HBP returns one per variant)
    cand_order = {id(v): v for v in variants}
    # Match by pos (variants list is the same order as we built)
    causal_pos = variants[causal_idx]["pos"]
    # Sort cands by PIP and find the causal
    cands_sorted = sorted(cands, key=lambda c: -c.pip)
    causal_rank = None
    causal_pip = None
    for rank, c in enumerate(cands_sorted, 1):
        if c.pos == causal_pos:
            causal_rank = rank
            causal_pip = float(c.pip)
            break
    return {
        "N": N,
        "h2": h2,
        "seed": seed,
        "n_variants": len(variants),
        "causal_pos": causal_pos,
        "causal_af": variants[causal_idx]["af_total"],
        "causal_rank": causal_rank,
        "causal_pip": causal_pip,
        "runtime": dt,
    }


def main():
    reader = BgenReader("/mnt/data/GraphGWAS/tests/data/human/1kGP_bgen")
    Ns = [500, 1000, 2000, 3000]
    n_reps = 30
    h2 = 0.10
    all_results = []
    t_start = time.time()
    for N in Ns:
        print(f"\n=== N = {N} ===")
        for rep in range(n_reps):
            r = run_one_rep(reader, N, h2, seed=1000 + N * 100 + rep)
            all_results.append(r)
            if r["causal_rank"] is not None:
                print(f"  rep {rep:2d}: causal_rank={r['causal_rank']:3d} "
                      f"pip={r['causal_pip']:.3f} t={r['runtime']:.2f}s")
            else:
                print(f"  rep {rep:2d}: causal not in output (skipped)")
    # Aggregate
    summary = {}
    for N in Ns:
        sub = [r for r in all_results if r["N"] == N and r["causal_rank"] is not None]
        if not sub:
            continue
        ranks = [r["causal_rank"] for r in sub]
        pips = [r["causal_pip"] for r in sub]
        summary[str(N)] = {
            "n_reps_with_result": len(sub),
            "rank_1_rate_pct": 100 * sum(1 for r in ranks if r == 1) / len(ranks),
            "mean_rank": float(np.mean(ranks)),
            "median_rank": float(np.median(ranks)),
            "mean_pip": float(np.mean(pips)),
            "mean_runtime_s": float(np.mean([r["runtime"] for r in sub])),
        }
    print(f"\n=== Summary (h²={h2}, n_reps={n_reps}) ===")
    print(json.dumps(summary, indent=2))
    print(f"Total wall time: {time.time()-t_start:.1f}s")

    (OUT / "power_curve.json").write_text(
        json.dumps({"summary": summary, "results": all_results}, indent=2))
    # TSV
    with (OUT / "power_curve.tsv").open("w") as fh:
        fh.write("N\trep\tcausal_pos\tcausal_af\tcausal_rank\tcausal_pip\truntime\n")
        for r in all_results:
            fh.write(f"{r['N']}\t{r['seed']}\t{r['causal_pos']}\t"
                     f"{r['causal_af']:.4f}\t{r['causal_rank']}\t"
                     f"{r['causal_pip']}\t{r['runtime']:.3f}\n")
    print(f"\nWrote {OUT}/power_curve.{{json,tsv}}")


if __name__ == "__main__":
    main()
