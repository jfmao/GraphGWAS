"""Run SuSiE-inf and FINEMAP-inf on the same 3-scenario campaign as
benchmark_hbp_vs_susie_finemap.py (strong / weak / functional), so that
Fig 4a can be a uniform 6-method x 3-scenario panel.

Reuses the simulation setup from benchmark_hbp_vs_susie_finemap.py exactly
(same CENTERS, same seed_base, same beta and h2 per scenario, same
functional-eQTL window selection) but only invokes the two missing
methods. Output: results/benchmark_v2/inf_methods/inf_methods_h2h.json,
keyed by scenario, mergeable into hbp_h2h_50rep.json at plot time.
"""

import json
import os
import time

import numpy as np

import graphgwas.config as cfg
cfg._N_SAMPLES_DETECTED = False

from graphgwas.db import GraphGWASConnection
from graphgwas.config import detect_n_samples
from graphgwas.simulate import simulate_F1_single_causal_in_ld, load_simulation_as_phenotype
from graphgwas.finemapping_v2 import _load_locus_variants, _load_eqtl_cache
from graphgwas.genotype import get_all_indices, get_phenotype_values

# Reuse inf-method runners
from tests.benchmark_inf_methods import _run_susie_inf, _run_finemap_inf

OUTDIR = "/mnt/data/GraphGWAS/results/benchmark_v2/inf_methods"
os.makedirs(OUTDIR, exist_ok=True)
EQTL_PATH = "/mnt/data/GraphGWAS/data/annotations/gtex_chr22_enhanced_cache.json"
WINDOW = 50000

# Same as the 4-method script
CENTERS = [17000000, 19000000, 21000000, 23000000, 25000000,
           27000000, 29000000, 31000000, 33000000, 35000000,
           37000000, 39000000, 41000000, 43000000, 45000000,
           17500000, 20000000, 24000000, 28000000, 36000000,
           18000000, 22000000, 26000000, 30000000, 34000000,
           38000000, 40000000, 42000000, 44000000, 46000000,
           16500000, 18500000, 20500000, 22500000, 24500000,
           26500000, 28500000, 30500000, 32500000, 34500000,
           36500000, 38500000, 40500000, 42500000, 44500000,
           17200000, 19200000, 21200000, 23200000, 25200000]


def _xy_from_variants(variants, pheno, all_idx, causal_vid):
    """Build (X, y, causal_idx) from the same variants/pheno used by the 4-method script."""
    nv = len(variants)
    X = np.zeros((len(all_idx), nv))
    for i, v in enumerate(variants):
        d = v["dosage"]
        col_mean = np.nanmean(d)
        X[:, i] = np.where(np.isnan(d), col_mean, d)
    y = np.where(np.isnan(pheno), np.nanmean(pheno), pheno)
    cidx = next((i for i, v in enumerate(variants) if v["variantId"] == causal_vid), -1)
    return X, y, cidx


def run_scenario(name, conn, all_idx, eqtl_cache, n_reps, beta, h2,
                 seed_base, functional=False):
    print(f"\n{'='*80}\n{name}\n{'='*80}")
    results = []

    if functional:
        eqtl_pos = sorted([int(k.split(":")[1]) for k in eqtl_cache if k.startswith("chr22:")])
        windows = {}
        for pos in eqtl_pos:
            wk = (pos // WINDOW) * WINDOW
            windows.setdefault(wk, []).append(pos)
        good = [(k, len(v)) for k, v in windows.items() if len(v) >= 10]
        good.sort(key=lambda x: -x[1])
        locus_centers = [w[0] + WINDOW // 2 for w in good[:n_reps]]
    else:
        locus_centers = CENTERS[:n_reps]

    for rep in range(min(n_reps, len(locus_centers))):
        center = locus_centers[rep]
        seed = seed_base + rep

        if functional:
            variants = _load_locus_variants(conn, "chr22", center, WINDOW, all_idx)
            eqtl_vars = [v for v in variants if v["variantId"] in eqtl_cache
                         and 0.02 < v.get("af_total", 0) < 0.98]
            if len(eqtl_vars) < 5 or len(variants) < 20:
                continue
            rng = np.random.default_rng(seed)
            causal_v = rng.choice(eqtl_vars)
            dosage = causal_v["dosage"].copy()
            dosage = np.where(np.isnan(dosage), np.nanmean(dosage), dosage)
            if causal_v.get("af_total", 0) < 0.02:
                continue
            genetic = beta * dosage
            var_g = np.var(genetic)
            if var_g < 1e-10:
                continue
            var_e = var_g * (1 / h2 - 1)
            pheno = genetic + rng.normal(0, np.sqrt(var_e), len(all_idx))
            conn.execute_write("""
                UNWIND $data AS row
                MATCH (s:Sample {packed_index: row.idx})
                SET s.gwas_value = row.val
            """, {"data": [{"idx": int(all_idx[i]), "val": float(pheno[i])}
                            for i in range(len(all_idx))]})
            causal_vid = causal_v["variantId"]
            causal_pos = causal_v["pos"]
            causal_af = causal_v.get("af_total", 0)
        else:
            try:
                sim = simulate_F1_single_causal_in_ld(
                    conn, "chr22", center, beta=beta, h2_target=h2, seed=seed)
            except Exception as e:
                print(f"  rep {rep}: sim failed ({e})")
                continue
            causal_vid = sim.causal_variants[0].variant_id
            causal_pos = sim.causal_variants[0].pos
            causal_af = sim.causal_variants[0].af
            load_simulation_as_phenotype(conn, sim, f"h2h_inf_{seed}")
            pheno = get_phenotype_values(conn, all_idx)
            variants = _load_locus_variants(conn, "chr22", causal_pos, WINDOW, all_idx)

        if len(variants) < 10:
            continue
        X, y, cidx = _xy_from_variants(variants, pheno, all_idx, causal_vid)
        if cidx < 0:
            continue

        susie_inf = _run_susie_inf(X, y, cidx)
        finemap_inf = _run_finemap_inf(X, y, cidx)

        row = {
            "rep": rep,
            "nv": len(variants),
            "af": causal_af,
            "susie_inf_rank": susie_inf["rank"],
            "susie_inf_pip": susie_inf["pip"],
            "susie_inf_time": susie_inf["time"],
            "finemap_inf_rank": finemap_inf["rank"],
            "finemap_inf_pip": finemap_inf["pip"],
            "finemap_inf_time": finemap_inf["time"],
        }
        results.append(row)
        print(f"  rep {rep:2d}: nv={len(variants):4d} | "
              f"SuSiE-inf=#{susie_inf['rank']} ({susie_inf['time']:.2f}s) | "
              f"FINEMAP-inf=#{finemap_inf['rank']} ({finemap_inf['time']:.2f}s)",
              flush=True)
    return results


def main():
    print("Loading eQTL cache...", flush=True)
    eqtl_cache = _load_eqtl_cache(EQTL_PATH)
    print(f"  {len(eqtl_cache)} eQTL annotations")

    with GraphGWASConnection("bolt://localhost:7688", "neo4j", "graphgwas") as conn:
        detect_n_samples(conn)
        all_idx = get_all_indices(conn)
        print(f"  N_SAMPLES={cfg.N_SAMPLES}")

        s_strong = run_scenario(
            "STRONG (random, beta=0.5, h2=0.10)",
            conn, all_idx, eqtl_cache,
            n_reps=30, beta=0.5, h2=0.10, seed_base=1000, functional=False)

        s_weak = run_scenario(
            "WEAK (random, beta=0.2, h2=0.02)",
            conn, all_idx, eqtl_cache,
            n_reps=30, beta=0.2, h2=0.02, seed_base=2000, functional=False)

        s_func = run_scenario(
            "FUNCTIONAL (eQTL causal, beta=0.2, h2=0.02)",
            conn, all_idx, eqtl_cache,
            n_reps=30, beta=0.2, h2=0.02, seed_base=3000, functional=True)

    out = {"strong": s_strong, "weak": s_weak, "functional": s_func}
    fname = f"{OUTDIR}/inf_methods_h2h.json"
    with open(fname, "w") as f:
        json.dump(out, f, indent=2, default=str)
    print(f"\nWrote {fname}")

    # Print per-scenario summary
    for scen, rows in out.items():
        n = len(rows)
        if not n:
            continue
        for k, label in [("susie_inf", "SuSiE-inf"), ("finemap_inf", "FINEMAP-inf")]:
            ranks = [r[f"{k}_rank"] for r in rows if r.get(f"{k}_rank") is not None]
            r1 = sum(1 for r in ranks if r == 1)
            tm = [r[f"{k}_time"] for r in rows if r.get(f"{k}_time") is not None]
            print(f"  {scen:10s} {label:12s}: rank-1 {r1}/{n} ({100*r1/n:.0f}%), "
                  f"mean rank {np.mean(ranks):.2f}, runtime {np.mean(tm):.2f}s")


if __name__ == "__main__":
    main()
