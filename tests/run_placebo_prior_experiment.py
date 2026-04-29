"""Placebo-prior control experiment.

Tests whether the 0% → 88% intervention gain requires *biological* discriminative
information in the per-variant prior_score, or whether any non-uniform prior
distribution would suffice.

Three cache configurations on the same simulated loci:
  (a) baseline    — graph cache without prior_score (uniform prior)
  (b) informative — v2 chr22 cache with cCRE-class prior_score (the original intervention)
  (c) placebo     — same v2 cache, but prior_score values PERMUTED across variants
                    (preserving the marginal distribution of [0.4, 0.5, 0.6, 0.8, 1.0]
                     but randomising which variants receive which class)

Output: results/placebo_prior/placebo_prior.{json,tsv}
"""
from __future__ import annotations

import copy
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

OUT = Path("/mnt/data/GraphGWAS/results/placebo_prior")
OUT.mkdir(parents=True, exist_ok=True)
EQTL_PATH = "/mnt/data/GraphGWAS/data/annotations/gtex_chr22_enhanced_cache.json"
HUMAN_V2_CACHE = "/mnt/data/GraphGWAS/data/annotations/human_graph_cache_v2_chr22.json"
WINDOW = 50000
N_REPS = 30
N_PERMS = 5  # number of placebo permutations to bound the placebo-rate variability
RNG = np.random.default_rng(2026)

CENTERS = [17000000, 19000000, 21000000, 23000000, 25000000,
           27000000, 29000000, 31000000, 33000000, 35000000,
           37000000, 39000000, 41000000, 43000000, 45000000,
           17500000, 20000000, 24000000, 28000000, 36000000,
           18000000, 22000000, 26000000, 30000000, 34000000,
           38000000, 40000000, 42000000, 44000000, 46000000]


def _gafm_with_cache(variants, pheno, eqtl_cache, graph_cache, alpha=0.9, r2_smooth=0.3):
    """GAFM combining LD-deconvolved evidence with eqtl + prior_score from cache."""
    z = _compute_association_stats(variants, pheno)
    R = _compute_ld_matrix(variants)
    u, _ = _ld_deconvolve(z, R, r2_smooth)
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
    # Add prior_score from graph_cache (the intervention layer)
    for i, v in enumerate(variants):
        info = graph_cache.get(v["variantId"], {}) if graph_cache else {}
        ps = info.get("prior_score") if isinstance(info, dict) else None
        if ps is not None:
            zf[i] += float(ps) * 0.5  # 0.5x scaling matches Methods §intervention
    if zf.max() > 0:
        zf = zf / zf.max() * z.max()
    return _softmax(alpha * u + (1 - alpha) * zf)


def _cs_size(pip_array, coverage=0.95):
    """95% credible-set size: minimal sorted-PIP set summing to 0.95."""
    sorted_p = np.sort(pip_array)[::-1]
    cum = np.cumsum(sorted_p)
    return int(np.searchsorted(cum, coverage) + 1)


def _hbp_cs_size(cands, coverage=0.95):
    pips = np.array([c.pip for c in cands])
    if len(pips) == 0:
        return None
    return _cs_size(pips, coverage)


def build_placebo_cache(informative_cache: dict, seed: int) -> dict:
    """Permute prior_score values across variants in the cache."""
    rng = np.random.default_rng(seed)
    placebo = copy.deepcopy(informative_cache)
    # Collect prior_score values
    keys_with_score = [k for k, v in placebo.items() if isinstance(v, dict) and "prior_score" in v]
    scores = [placebo[k]["prior_score"] for k in keys_with_score]
    permuted = rng.permutation(scores)
    for k, s in zip(keys_with_score, permuted):
        placebo[k]["prior_score"] = float(s)
    return placebo


def main() -> None:
    print("Loading caches...", flush=True)
    eqtl_cache = _load_eqtl_cache(EQTL_PATH)
    print(f"  eQTL: {len(eqtl_cache)} annotations")

    print(f"Loading human v2 chr22 cache (the informative cache)...", flush=True)
    informative_cache = json.load(open(HUMAN_V2_CACHE))
    n_with_score = sum(1 for v in informative_cache.values()
                       if isinstance(v, dict) and "prior_score" in v)
    print(f"  total: {len(informative_cache)} variants, {n_with_score} with prior_score")

    # Baseline cache: same structure but with prior_score stripped
    baseline_cache = {k: {kk: vv for kk, vv in v.items() if kk != "prior_score"}
                      for k, v in informative_cache.items()
                      if isinstance(v, dict)}

    with GraphGWASConnection("bolt://localhost:7688", "neo4j", "graphgwas") as conn:
        detect_n_samples(conn)
        all_idx = get_all_indices(conn)
        print(f"  N_SAMPLES = {cfg.N_SAMPLES}")

        # Generate the simulated loci
        print(f"\nGenerating {N_REPS} simulated loci (β=0.5, h²=0.10)...", flush=True)
        sims = []
        for rep in range(N_REPS):
            try:
                sim = simulate_F1_single_causal_in_ld(
                    conn, "chr22", CENTERS[rep], beta=0.5, h2_target=0.10,
                    seed=20000 + rep)
            except Exception as e:
                print(f"  rep {rep}: sim failed ({e})")
                continue
            causal_vid = sim.causal_variants[0].variant_id
            causal_pos = sim.causal_variants[0].pos
            load_simulation_as_phenotype(conn, sim, f"plac_{rep}")
            pheno = get_phenotype_values(conn, all_idx)
            variants = _load_locus_variants(conn, "chr22", causal_pos, WINDOW, all_idx)
            cidx = next((i for i, v in enumerate(variants) if v["variantId"] == causal_vid), -1)
            if cidx < 0 or len(variants) < 10:
                continue
            sims.append({
                "rep": rep,
                "causal_vid": causal_vid,
                "variants": variants,
                "pheno": pheno,
            })
            if rep % 5 == 0:
                print(f"  rep {rep}: nv={len(variants)}", flush=True)

        n_sims = len(sims)
        print(f"\n{n_sims} usable simulations.")

        # ===================================================================
        # Run all three cache configurations
        # ===================================================================
        results_summary = {}

        for cfg_name, cache in [("baseline", baseline_cache),
                                ("informative", informative_cache)]:
            print(f"\n=== {cfg_name} ===", flush=True)
            hbp_cs = []
            gafm_cs = []
            for s in sims:
                hbp = fast_hbp_finemap(s["variants"], s["pheno"], cache,
                                        eqtl_cache=eqtl_cache)
                gafm_pip = _gafm_with_cache(s["variants"], s["pheno"], eqtl_cache, cache)
                hbp_cs.append(_hbp_cs_size(hbp))
                gafm_cs.append(_cs_size(gafm_pip))
            n_hbp_tighter = sum(1 for h, g in zip(hbp_cs, gafm_cs)
                                if h is not None and g is not None and h < g)
            n_eq = sum(1 for h, g in zip(hbp_cs, gafm_cs)
                       if h is not None and g is not None and h == g)
            n_hbp_wider = sum(1 for h, g in zip(hbp_cs, gafm_cs)
                              if h is not None and g is not None and h > g)
            results_summary[cfg_name] = {
                "n": n_sims,
                "hbp_tighter": n_hbp_tighter,
                "hbp_equal": n_eq,
                "hbp_wider": n_hbp_wider,
                "hbp_tighter_rate": round(n_hbp_tighter / n_sims, 3),
                "mean_hbp_cs_size": float(np.mean([c for c in hbp_cs if c is not None])),
                "mean_gafm_cs_size": float(np.mean([c for c in gafm_cs if c is not None])),
            }
            print(f"  HBP-tighter rate: {n_hbp_tighter}/{n_sims} = {n_hbp_tighter/n_sims*100:.0f}%")

        # Placebo: average over N_PERMS permutations
        print(f"\n=== placebo (avg over {N_PERMS} permutations) ===", flush=True)
        placebo_rates = []
        per_perm_details = []
        for perm in range(N_PERMS):
            placebo_cache = build_placebo_cache(informative_cache, seed=30000 + perm)
            hbp_cs = []
            gafm_cs = []
            for s in sims:
                hbp = fast_hbp_finemap(s["variants"], s["pheno"], placebo_cache,
                                        eqtl_cache=eqtl_cache)
                gafm_pip = _gafm_with_cache(s["variants"], s["pheno"], eqtl_cache, placebo_cache)
                hbp_cs.append(_hbp_cs_size(hbp))
                gafm_cs.append(_cs_size(gafm_pip))
            n_hbp_tighter = sum(1 for h, g in zip(hbp_cs, gafm_cs)
                                if h is not None and g is not None and h < g)
            placebo_rates.append(n_hbp_tighter / n_sims)
            per_perm_details.append({
                "perm": perm,
                "hbp_tighter": n_hbp_tighter,
                "rate": round(n_hbp_tighter / n_sims, 3),
            })
            print(f"  perm {perm}: HBP-tighter rate = {n_hbp_tighter}/{n_sims} = {n_hbp_tighter/n_sims*100:.0f}%")

        results_summary["placebo"] = {
            "n_perms": N_PERMS,
            "n_sims": n_sims,
            "placebo_rates": placebo_rates,
            "placebo_rate_mean": round(float(np.mean(placebo_rates)), 3),
            "placebo_rate_std": round(float(np.std(placebo_rates, ddof=1)), 3),
            "placebo_rate_min": round(float(np.min(placebo_rates)), 3),
            "placebo_rate_max": round(float(np.max(placebo_rates)), 3),
            "per_perm": per_perm_details,
        }

        # ===================================================================
        # Save outputs
        # ===================================================================
        out = {
            "scenario": "F1 single causal in LD on chr22, β=0.5, h²=0.10",
            "n_sims_per_cache": n_sims,
            "cache_baseline_size": len(baseline_cache),
            "cache_informative_n_with_prior": n_with_score,
            "cache_human_v2_chr22": HUMAN_V2_CACHE,
            "results": results_summary,
            "interpretation": (
                "Compare HBP-tighter rates across the three cache configurations. "
                "If placebo rate is close to baseline rate, the gain genuinely "
                "requires biological structure in the prior; if placebo rate is "
                "close to informative rate, the gain reflects heterogeneity per se."
            ),
        }
        with (OUT / "placebo_prior.json").open("w") as f:
            json.dump(out, f, indent=2)
        print(f"\nWrote {OUT}/placebo_prior.json")

        # Headline
        print(f"\n=== HEADLINE ===")
        print(f"baseline    HBP-tighter: {results_summary['baseline']['hbp_tighter_rate']*100:.0f}%")
        print(f"informative HBP-tighter: {results_summary['informative']['hbp_tighter_rate']*100:.0f}%")
        print(f"placebo     HBP-tighter: {results_summary['placebo']['placebo_rate_mean']*100:.0f}% "
              f"± {results_summary['placebo']['placebo_rate_std']*100:.0f} (5 perms)")


if __name__ == "__main__":
    main()
