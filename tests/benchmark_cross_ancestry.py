"""Cross-ancestry F1 fine-mapping on 1KG chr22 (EUR / AFR / EAS superpops).

Subsets the 1KG BGEN samples by superpopulation (503 EUR, 661 AFR, 504 EAS),
runs F1 single-causal-in-LD simulations at h²=0.05 (30 reps per ancestry),
and reports HBP rank-1 rate, mean rank, and mean PIP by ancestry.

Answers the reviewer question: "Does HBP's precision hold across non-European
LD structure?"

Writes results/benchmark_v2/cross_ancestry/cross_ancestry.{json,tsv}.
"""
from __future__ import annotations

import bisect
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

from graphgwas.bgen_reader import BgenReader
from graphgwas.finemapping_v2 import fast_hbp_finemap

OUT = Path("/mnt/data/GraphGWAS/results/benchmark_v2/cross_ancestry")
OUT.mkdir(parents=True, exist_ok=True)
POPMAP = "/mnt/data/GraphGWAS/tests/data/human/1kg_popmap.tsv"


def load_ancestry_masks(reader: BgenReader) -> dict[str, np.ndarray]:
    """Return boolean masks for EUR / AFR / EAS aligned to BGEN sample order."""
    bgen_samples = [str(s) for s in reader.samples("22")]
    df = pd.read_csv(POPMAP, sep="\t")
    sample_to_sp = dict(zip(df["sample"], df["superpopulation"]))
    masks = {}
    for sp in ["EUR", "AFR", "EAS"]:
        masks[sp] = np.array(
            [sample_to_sp.get(s) == sp for s in bgen_samples], dtype=bool
        )
        print(f"  {sp}: {masks[sp].sum()} samples")
    return masks


def run_one_rep(reader: BgenReader, ancestry_mask: np.ndarray, h2: float,
                 seed: int, window_bp: int = 50_000) -> dict | None:
    rng = np.random.default_rng(seed)
    positions = reader._pos_array("22")
    pos_min, pos_max = positions[0] + window_bp, positions[-1] - window_bp
    for _ in range(10):
        center = int(rng.integers(pos_min, pos_max))
        vdf, dosage = reader.load_locus(
            "22", center - window_bp // 2, center + window_bp // 2, format="dosage"
        )
        if len(vdf) == 0:
            continue
        # Restrict to ancestry samples, impute NaN, filter by AF in this ancestry
        dos_anc = dosage[ancestry_mask]
        col_mean = np.nanmean(dos_anc, axis=0)
        dos_anc = np.where(np.isnan(dos_anc), col_mean, dos_anc)
        af_anc = col_mean / 2.0
        maf = np.minimum(af_anc, 1 - af_anc)
        keep = (maf > 0.05) & np.isfinite(col_mean)
        if keep.sum() < 30:
            continue
        dos_keep = dos_anc[:, keep]
        vdf_keep = vdf.loc[keep].reset_index(drop=True)
        af_keep = af_anc[keep]
        # Prefer a causal variant at AF 0.1-0.4 in this ancestry
        cands = np.where((af_keep > 0.1) & (af_keep < 0.4))[0]
        if len(cands) == 0:
            continue
        causal_idx = int(rng.choice(cands))
        break
    else:
        return None
    # Simulate phenotype
    g = dos_keep[:, causal_idx]
    g_c = g - g.mean()
    var_g = g_c.var()
    beta = np.sqrt(h2 / var_g) if var_g > 1e-8 else 0
    eps = rng.standard_normal(dos_keep.shape[0]) * np.sqrt(1 - h2)
    y = beta * g_c + eps
    # Run HBP
    variants = []
    for j in range(dos_keep.shape[1]):
        row = vdf_keep.iloc[j]
        variants.append({
            "variantId": f"chr22:{int(row['pos'])}:{row['a1']}:{row['a2']}",
            "chr": row["chr"],
            "pos": int(row["pos"]),
            "ref": row["a1"],
            "alt": row["a2"],
            "af_total": float(af_keep[j]),
            "dosage": dos_keep[:, j],
        })
    t0 = time.time()
    cands_hbp = fast_hbp_finemap(variants, y, graph_cache={}, n_rounds=3, chr_name="chr22")
    dt = time.time() - t0
    causal_pos = variants[causal_idx]["pos"]
    cands_sorted = sorted(cands_hbp, key=lambda c: -c.pip)
    causal_rank = None
    causal_pip = None
    for r, c in enumerate(cands_sorted, 1):
        if c.pos == causal_pos:
            causal_rank = r
            causal_pip = float(c.pip)
            break
    return {
        "seed": seed,
        "n_samples": int(ancestry_mask.sum()),
        "n_variants": len(variants),
        "causal_pos": causal_pos,
        "causal_af": float(af_keep[causal_idx]),
        "causal_rank": causal_rank,
        "causal_pip": causal_pip,
        "runtime": dt,
    }


def main():
    reader = BgenReader("/mnt/data/GraphGWAS/tests/data/human/1kGP_bgen")
    print("Loading ancestry masks:")
    masks = load_ancestry_masks(reader)

    n_reps = 30
    h2 = 0.05
    all_results = {}
    t_start = time.time()
    for sp, mask in masks.items():
        print(f"\n=== {sp} (n={mask.sum()}) ===")
        sp_results = []
        for rep in range(n_reps):
            r = run_one_rep(reader, mask, h2, seed=3000 + hash(sp) % 1000 + rep)
            if r is None:
                continue
            r["ancestry"] = sp
            sp_results.append(r)
            print(f"  rep {rep:2d}: rank={r['causal_rank']:>3} "
                  f"pip={r['causal_pip']:.3f} af={r['causal_af']:.3f}")
        all_results[sp] = sp_results

    # Summary
    summary = {}
    for sp, results in all_results.items():
        if not results:
            continue
        ranks = [r["causal_rank"] for r in results if r["causal_rank"] is not None]
        pips = [r["causal_pip"] for r in results if r["causal_pip"] is not None]
        summary[sp] = {
            "n_samples": results[0]["n_samples"],
            "n_reps_with_result": len(ranks),
            "rank_1_rate_pct": 100 * sum(1 for r in ranks if r == 1) / len(ranks),
            "mean_rank": float(np.mean(ranks)),
            "median_rank": float(np.median(ranks)),
            "mean_pip": float(np.mean(pips)),
            "mean_runtime_s": float(np.mean([r["runtime"] for r in results])),
        }

    print("\n=== Cross-ancestry summary (h² = {}) ===".format(h2))
    print(json.dumps(summary, indent=2))
    print(f"Total: {time.time()-t_start:.1f}s")

    flat = [r for results in all_results.values() for r in results]
    (OUT / "cross_ancestry.json").write_text(
        json.dumps({"summary": summary, "results": flat,
                    "params": {"n_reps": n_reps, "h2": h2}}, indent=2))
    with (OUT / "cross_ancestry.tsv").open("w") as fh:
        fh.write("ancestry\tseed\tn_samples\tn_variants\tcausal_pos\tcausal_af\t"
                 "causal_rank\tcausal_pip\truntime\n")
        for r in flat:
            fh.write(f"{r['ancestry']}\t{r['seed']}\t{r['n_samples']}\t"
                     f"{r['n_variants']}\t{r['causal_pos']}\t{r['causal_af']:.4f}\t"
                     f"{r['causal_rank']}\t{r['causal_pip']}\t{r['runtime']:.3f}\n")
    print(f"\nWrote {OUT}/cross_ancestry.{{json,tsv}}")


if __name__ == "__main__":
    main()
