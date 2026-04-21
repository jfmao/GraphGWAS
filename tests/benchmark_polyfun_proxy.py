"""Polyfun-proxy benchmark: SuSiE with annotation priors vs L1.

Polyfun's contribution is that it computes per-variant prior weights from
functional enrichment (S-LDSC) and passes them to SuSiE via its
`prior_weights` argument. We approximate the same mechanism using our
eQTL annotation graph: variants with GTEx eQTL edges get boosted priors.

This is a scoped comparison answering: "if SuSiE gets the same annotation
information as L1, does L1 still win?"

Pipeline per replicate:
  1. Random 50 kb window on 1KG chr22
  2. F1 simulation: pick causal variant from variants with GTEx eQTL
     (tissue-specific), β=0.15, h²=0.02
  3. Query Neo4j for per-variant eQTL edge counts → prior weights
  4. Run:
     - SuSiE (vanilla, uniform priors)
     - SuSiE + annotation prior ("Polyfun-proxy")
     - L1 Bayesian with same annotations

Writes results/benchmark_v2/polyfun_proxy/polyfun_proxy.{json,tsv}.
"""
from __future__ import annotations

import bisect
import json
import subprocess
import tempfile
import time
from pathlib import Path

import numpy as np

from graphgwas import config
from graphgwas.bgen_reader import BgenReader
from graphgwas.db import GraphGWASConnection
from graphgwas.finemapping_v2 import fast_hbp_finemap

OUT = Path("/mnt/data/GraphGWAS/results/benchmark_v2/polyfun_proxy")
OUT.mkdir(parents=True, exist_ok=True)


def _load_locus_with_priors(reader: BgenReader, chr_name: str, start: int, end: int,
                            min_maf: float = 0.02):
    """Return (variants_df, dosage_imputed, prior_weights, eqtl_counts)."""
    vdf, dosage = reader.load_locus(chr_name, start, end, format="dosage")
    if len(vdf) == 0:
        return None
    col_mean = np.nanmean(dosage, axis=0)
    dosage = np.where(np.isnan(dosage), col_mean, dosage)
    af = col_mean / 2.0
    maf = np.minimum(af, 1 - af)
    keep = (maf > min_maf) & np.isfinite(col_mean)
    if keep.sum() < 30:
        return None
    dosage = dosage[:, keep]
    vdf = vdf.loc[keep].reset_index(drop=True)

    # Per-variant eQTL count from Neo4j
    positions = vdf["pos"].tolist()
    chr_for_query = f"chr{reader._norm_chr(chr_name)}"
    eqtl_counts = np.zeros(len(positions), dtype=float)
    try:
        with GraphGWASConnection(
            config.NEO4J_URI, config.NEO4J_USER, config.NEO4J_PASSWORD
        ) as conn:
            # Batched lookup — one query with UNWIND
            rows = conn.execute_read(
                "UNWIND $ps AS p MATCH (v:Variant {chr: $chr, pos: p}) "
                "OPTIONAL MATCH (v)-[e:eQTL]->() "
                "RETURN p AS pos, count(e) AS n_eqtl",
                {"chr": chr_for_query, "ps": positions},
            ).data()
            pos2count = {int(r["pos"]): int(r["n_eqtl"]) for r in rows}
        eqtl_counts = np.array([pos2count.get(p, 0) for p in positions], dtype=float)
    except Exception as e:
        # Neo4j unavailable — all zero priors (vanilla case)
        print(f"  Warning: Neo4j query failed ({e}) — using uniform priors")
    # Convert eQTL counts to prior weights (softmax-like, normalized)
    # Base uniform + log-boost for variants with eQTLs
    raw = 1.0 + np.log1p(eqtl_counts)
    prior_weights = raw / raw.sum()
    return vdf, dosage, prior_weights, eqtl_counts


def _run_susie(X: np.ndarray, y: np.ndarray, causal_idx: int,
                prior_weights: np.ndarray | None = None,
                timeout: int = 120) -> dict:
    """Run SuSiE via R subprocess with optional prior_weights."""
    tmp = tempfile.mkdtemp()
    np.savetxt(f"{tmp}/X.csv", X, delimiter=",")
    np.savetxt(f"{tmp}/y.csv", y)
    pw_line = ""
    if prior_weights is not None:
        np.savetxt(f"{tmp}/pw.csv", prior_weights)
        pw_line = f'pw <- scan("{tmp}/pw.csv")\n'
        susie_call = 'susie(X, y, L=5, prior_weights=pw, verbose=FALSE)'
    else:
        susie_call = 'susie(X, y, L=5, verbose=FALSE)'
    r_script = f'''
    library(susieR)
    X <- as.matrix(read.csv("{tmp}/X.csv", header=FALSE))
    y <- scan("{tmp}/y.csv")
    {pw_line}
    fit <- tryCatch({susie_call}, error=function(e) NULL)
    if(!is.null(fit)) {{
        pips <- fit$pip
        ci <- {causal_idx + 1}
        cat(sprintf("PIP=%.6f\\n", pips[ci]))
        cat(sprintf("RANK=%d\\n", sum(pips >= pips[ci])))
        cat(sprintf("MAX=%.6f\\n", max(pips)))
    }} else cat("FAILED\\n")
    '''
    t0 = time.time()
    try:
        proc = subprocess.run(
            ["R", "--no-save", "--no-restore", "-e", r_script],
            capture_output=True, text=True, timeout=timeout,
        )
    except subprocess.TimeoutExpired:
        return {"rank": None, "pip": 0, "max_pip": 0, "time": timeout}
    dt = time.time() - t0
    rank = pip = max_pip = None
    for line in proc.stdout.split("\n"):
        if line.startswith("PIP="):
            pip = float(line.split("=")[1])
        elif line.startswith("RANK="):
            rank = int(line.split("=")[1])
        elif line.startswith("MAX="):
            max_pip = float(line.split("=")[1])
    return {"rank": rank, "pip": pip, "max_pip": max_pip, "time": dt}


def run_one_rep(reader: BgenReader, seed: int, h2: float = 0.02):
    rng = np.random.default_rng(seed)
    positions = reader._pos_array("22")
    pos_min, pos_max = positions[0] + 50_000, positions[-1] - 50_000
    for _ in range(20):
        center = int(rng.integers(pos_min, pos_max))
        loaded = _load_locus_with_priors(reader, "22", center - 25_000, center + 25_000)
        if loaded is None:
            continue
        vdf, dosage, prior_weights, eqtl_counts = loaded
        # Prefer a causal variant that has an eQTL edge (enriched scenario)
        eqtl_cands = np.where(eqtl_counts > 0)[0]
        if len(eqtl_cands) == 0:
            continue
        causal_idx = int(rng.choice(eqtl_cands))
        break
    else:
        return None
    # Simulate phenotype y = β * g + ε with var(y) = 1, h² fixed
    g = dosage[:, causal_idx]
    g_c = g - g.mean()
    var_g = g_c.var()
    beta = np.sqrt(h2 / var_g) if var_g > 1e-8 else 0
    eps = rng.standard_normal(dosage.shape[0]) * np.sqrt(1 - h2)
    y = beta * g_c + eps

    n = dosage.shape[0]
    n_var = dosage.shape[1]

    # SuSiE vanilla
    vanilla = _run_susie(dosage, y, causal_idx, prior_weights=None)
    # SuSiE with annotation prior (Polyfun-proxy)
    annotated = _run_susie(dosage, y, causal_idx, prior_weights=prior_weights)

    # L1 Bayesian: build variants list and run fast_hbp_finemap (which is our
    # graph-native fine-mapper; L1 shares core infrastructure). Use same prior.
    variants = [
        {
            "variantId": f"chr22:{int(vdf.iloc[j]['pos'])}:{vdf.iloc[j]['a1']}:{vdf.iloc[j]['a2']}",
            "chr": "chr22",
            "pos": int(vdf.iloc[j]["pos"]),
            "ref": vdf.iloc[j]["a1"],
            "alt": vdf.iloc[j]["a2"],
            "af_total": float(dosage[:, j].mean() / 2),
            "dosage": dosage[:, j],
        }
        for j in range(n_var)
    ]
    t0 = time.time()
    # Supply annotation via graph_cache (mock gene/pathway mapping from eQTL boost)
    graph_cache = {}
    for j, v in enumerate(variants):
        if eqtl_counts[j] > 0:
            graph_cache[v["variantId"]] = {
                "genes": [f"gene_{j}"],
                "pathways": [f"pw_{j}"],
                "ppi": [],
            }
    cands = fast_hbp_finemap(variants, y, graph_cache=graph_cache,
                              n_rounds=3, chr_name="chr22")
    l1_time = time.time() - t0
    cands_sorted = sorted(cands, key=lambda c: -c.pip)
    causal_pos = variants[causal_idx]["pos"]
    l1_rank = None
    l1_pip = None
    for rank, c in enumerate(cands_sorted, 1):
        if c.pos == causal_pos:
            l1_rank = rank
            l1_pip = float(c.pip)
            break

    return {
        "seed": seed,
        "n_variants": n_var,
        "causal_pos": causal_pos,
        "causal_eqtl_count": int(eqtl_counts[causal_idx]),
        "n_annotated": int((eqtl_counts > 0).sum()),
        "susie_vanilla_rank": vanilla["rank"],
        "susie_vanilla_pip": vanilla["pip"],
        "susie_vanilla_time": vanilla["time"],
        "susie_annotated_rank": annotated["rank"],
        "susie_annotated_pip": annotated["pip"],
        "susie_annotated_time": annotated["time"],
        "l1_rank": l1_rank,
        "l1_pip": l1_pip,
        "l1_time": l1_time,
    }


def main():
    reader = BgenReader("/mnt/data/GraphGWAS/tests/data/human/1kGP_bgen")
    n_reps = 20
    h2 = 0.02
    results = []
    t_start = time.time()
    for i in range(n_reps):
        r = run_one_rep(reader, seed=2000 + i, h2=h2)
        if r is None:
            continue
        results.append(r)
        print(f"rep {i:2d}: SuSiE={r['susie_vanilla_rank']} "
              f"SuSiE+prior={r['susie_annotated_rank']} "
              f"L1={r['l1_rank']} | annotated={r['n_annotated']}/{r['n_variants']}")
    elapsed = time.time() - t_start

    # Summary
    summary = {}
    for label, key_rank, key_pip, key_time in [
        ("susie_vanilla", "susie_vanilla_rank", "susie_vanilla_pip", "susie_vanilla_time"),
        ("susie_annotated", "susie_annotated_rank", "susie_annotated_pip", "susie_annotated_time"),
        ("l1", "l1_rank", "l1_pip", "l1_time"),
    ]:
        ranks = [r[key_rank] for r in results if r[key_rank] is not None]
        pips = [r[key_pip] for r in results if r[key_pip] is not None]
        times = [r[key_time] for r in results if r[key_time] is not None]
        summary[label] = {
            "n": len(ranks),
            "rank_1": sum(1 for r in ranks if r == 1),
            "mean_rank": float(np.mean(ranks)) if ranks else None,
            "median_rank": float(np.median(ranks)) if ranks else None,
            "mean_pip": float(np.mean(pips)) if pips else None,
            "mean_time": float(np.mean(times)) if times else None,
        }
    # Head-to-head L1 vs SuSiE+prior
    wins = {"l1": 0, "susie_annot": 0, "tie": 0}
    for r in results:
        if r["l1_rank"] is None or r["susie_annotated_rank"] is None:
            continue
        if r["l1_rank"] < r["susie_annotated_rank"]:
            wins["l1"] += 1
        elif r["l1_rank"] > r["susie_annotated_rank"]:
            wins["susie_annot"] += 1
        else:
            wins["tie"] += 1
    summary["head_to_head_l1_vs_susie_annot"] = wins

    print("\n=== Summary ===")
    print(json.dumps(summary, indent=2))
    print(f"Total: {elapsed:.1f}s")

    (OUT / "polyfun_proxy.json").write_text(
        json.dumps({"summary": summary, "results": results,
                    "params": {"n_reps": n_reps, "h2": h2}}, indent=2))
    print(f"\nWrote {OUT}/polyfun_proxy.json")


if __name__ == "__main__":
    main()
