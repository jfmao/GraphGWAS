"""Hyperparameter sensitivity analysis for HBP and GAFM.

For each hyperparameter (HBP: alpha, damping/lambda, n_rounds, r2_smooth;
GAFM: alpha, r2_smooth_LD), run 10 simulated loci on the strong-signal
scenario and report the rank-1 rate for each perturbed value.

Output: results/hyperparam_sensitivity/hyperparam_sensitivity.{json,tsv}
        + a Supplementary-Figure-S7 plot.
"""
from __future__ import annotations

import json
import os
import time
from pathlib import Path

import numpy as np

import graphgwas.config as cfg
cfg._N_SAMPLES_DETECTED = False

from graphgwas.db import GraphGWASConnection
from graphgwas.config import detect_n_samples
from graphgwas.simulate import simulate_F1_single_causal_in_ld, load_simulation_as_phenotype
from graphgwas.finemapping_v2 import (
    _load_locus_variants, _compute_association_stats, _compute_ld_matrix,
    _ld_deconvolve, _softmax, _load_eqtl_cache,
    _cache_chr_graph_structure, fast_hbp_finemap,
)
from graphgwas.genotype import get_all_indices, get_phenotype_values

OUT = Path("/mnt/data/GraphGWAS/results/hyperparam_sensitivity")
OUT.mkdir(parents=True, exist_ok=True)
EQTL_PATH = "/mnt/data/GraphGWAS/data/annotations/gtex_chr22_enhanced_cache.json"
WINDOW = 50000
N_REPS = 10
CENTERS = [17000000, 19000000, 21000000, 23000000, 25000000,
           27000000, 29000000, 31000000, 33000000, 35000000]


def _fast_gafm(variants, pheno, eqtl_cache, alpha=0.5, r2_smooth=0.3):
    """GAFM (annotation-adaptive softmax fine-mapper). Local copy with hyperparams exposed."""
    z = _compute_association_stats(variants, pheno)
    R = _compute_ld_matrix(variants)
    u, n = _ld_deconvolve(z, R, r2_smooth)
    zf = np.zeros(len(variants))
    if eqtl_cache:
        for i, v in enumerate(variants):
            eq = eqtl_cache.get(v["variantId"])
            if eq:
                s = eq.get("composite", 0)
                nt = eq.get("n_tissues", 0)
                if 1 <= nt <= 3:
                    s *= 1.5
                elif nt >= 8:
                    s *= 0.3
                zf[i] = np.log1p(s)
        if zf.max() > 0:
            zf = zf / zf.max() * z.max()
    return _softmax(alpha * u + (1 - alpha) * zf)


def _rank_of_causal(pip_array_or_cands, variants, causal_vid):
    """Return the rank of the causal variant (1-based)."""
    if isinstance(pip_array_or_cands, np.ndarray):
        cidx = next((i for i, v in enumerate(variants) if v["variantId"] == causal_vid), -1)
        if cidx < 0:
            return None
        pip_c = pip_array_or_cands[cidx]
        return int(np.sum(pip_array_or_cands >= pip_c))
    else:  # list of FinemapCandidate
        for r, c in enumerate(sorted(pip_array_or_cands, key=lambda c: -c.pip), 1):
            if c.variant_id == causal_vid:
                return r
        return None


def main() -> None:
    print("Loading caches...", flush=True)
    eqtl_cache = _load_eqtl_cache(EQTL_PATH)
    print(f"  eQTL: {len(eqtl_cache)} annotations")

    with GraphGWASConnection("bolt://localhost:7688", "neo4j", "graphgwas") as conn:
        detect_n_samples(conn)
        all_idx = get_all_indices(conn)
        print(f"  N_SAMPLES = {cfg.N_SAMPLES}")

        print("Caching chr22 graph structure...", flush=True)
        graph_cache = _cache_chr_graph_structure(conn, "chr22")
        print(f"  {len(graph_cache)} variants with annotations")

        # Generate the 10 base simulations once
        print(f"\nGenerating {N_REPS} strong-signal simulations (β=0.5, h²=0.10)...", flush=True)
        sims = []
        for rep in range(N_REPS):
            try:
                sim = simulate_F1_single_causal_in_ld(
                    conn, "chr22", CENTERS[rep], beta=0.5, h2_target=0.10,
                    seed=10000 + rep)
            except Exception as e:
                print(f"  rep {rep}: sim failed ({e})")
                continue
            causal_vid = sim.causal_variants[0].variant_id
            causal_pos = sim.causal_variants[0].pos
            load_simulation_as_phenotype(conn, sim, f"hp_sens_{rep}")
            pheno = get_phenotype_values(conn, all_idx)
            variants = _load_locus_variants(conn, "chr22", causal_pos, WINDOW, all_idx)
            cidx = next((i for i, v in enumerate(variants) if v["variantId"] == causal_vid), -1)
            if cidx < 0 or len(variants) < 10:
                continue
            sims.append({
                "rep": rep,
                "causal_vid": causal_vid,
                "causal_pos": causal_pos,
                "variants": variants,
                "pheno": pheno,
            })
            print(f"  rep {rep}: nv={len(variants)}, causal at {causal_pos}", flush=True)

        n_sims = len(sims)
        print(f"\n{n_sims} usable simulations.")

        # ===================================================================
        # Hyperparameter sweeps
        # ===================================================================
        results = {"hbp": {}, "gafm": {}}

        # ----- HBP sweeps -----
        sweeps_hbp = {
            "alpha":     [0.4, 0.5, 0.6, 0.7, 0.8],   # default 0.6
            "damping":   [0.3, 0.4, 0.5, 0.6, 0.7],   # default 0.5
            "n_rounds":  [3, 4, 5, 6, 7],             # default 5
            "r2_smooth": [0.2, 0.25, 0.3, 0.35, 0.4], # default 0.3
        }
        defaults_hbp = dict(alpha=0.6, damping=0.5, n_rounds=5, r2_smooth=0.3)

        for param, values in sweeps_hbp.items():
            print(f"\n=== HBP sweep: {param} ===")
            results["hbp"][param] = []
            for v in values:
                kwargs = dict(defaults_hbp)
                kwargs[param] = v
                ranks = []
                t0 = time.time()
                for s in sims:
                    cands = fast_hbp_finemap(s["variants"], s["pheno"],
                                              graph_cache, eqtl_cache=eqtl_cache,
                                              **kwargs)
                    r = _rank_of_causal(cands, s["variants"], s["causal_vid"])
                    if r is not None:
                        ranks.append(r)
                rank1 = sum(1 for r in ranks if r == 1) / max(len(ranks), 1)
                results["hbp"][param].append({
                    "value": v, "rank_1_rate": round(rank1, 3),
                    "n_eval": len(ranks), "wall_s": round(time.time() - t0, 2),
                })
                print(f"  {param}={v}: rank-1 = {rank1:.2f} ({sum(1 for r in ranks if r == 1)}/{len(ranks)})")

        # ----- GAFM sweeps -----
        sweeps_gafm = {
            "alpha":     [0.3, 0.4, 0.5, 0.6, 0.7],   # default 0.5
            "r2_smooth": [0.2, 0.25, 0.3, 0.35, 0.4], # default 0.3
        }
        defaults_gafm = dict(alpha=0.5, r2_smooth=0.3)

        for param, values in sweeps_gafm.items():
            print(f"\n=== GAFM sweep: {param} ===")
            results["gafm"][param] = []
            for v in values:
                kwargs = dict(defaults_gafm)
                kwargs[param] = v
                ranks = []
                t0 = time.time()
                for s in sims:
                    pip = _fast_gafm(s["variants"], s["pheno"], eqtl_cache, **kwargs)
                    r = _rank_of_causal(pip, s["variants"], s["causal_vid"])
                    if r is not None:
                        ranks.append(r)
                rank1 = sum(1 for r in ranks if r == 1) / max(len(ranks), 1)
                results["gafm"][param].append({
                    "value": v, "rank_1_rate": round(rank1, 3),
                    "n_eval": len(ranks), "wall_s": round(time.time() - t0, 2),
                })
                print(f"  {param}={v}: rank-1 = {rank1:.2f} ({sum(1 for r in ranks if r == 1)}/{len(ranks)})")

        # ===================================================================
        # Save outputs
        # ===================================================================
        out = {
            "n_reps": N_REPS,
            "n_eval": n_sims,
            "scenario": "strong-signal F1 single causal in LD (β=0.5, h²=0.10)",
            "defaults_hbp": defaults_hbp,
            "defaults_gafm": defaults_gafm,
            "results": results,
        }
        with (OUT / "hyperparam_sensitivity.json").open("w") as f:
            json.dump(out, f, indent=2)
        print(f"\nWrote {OUT}/hyperparam_sensitivity.json")

        # TSV
        rows = []
        for method, sweeps in results.items():
            for param, vals in sweeps.items():
                for v in vals:
                    rows.append((method, param, v["value"], v["rank_1_rate"], v["n_eval"]))
        with (OUT / "hyperparam_sensitivity.tsv").open("w") as f:
            f.write("method\tparameter\tvalue\trank_1_rate\tn_eval\n")
            for r in rows:
                f.write("\t".join(str(x) for x in r) + "\n")
        print(f"Wrote {OUT}/hyperparam_sensitivity.tsv")


if __name__ == "__main__":
    main()
