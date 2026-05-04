"""5-method fine-mapping at the GWAS lead loci of the 4 grain traits.

Methods:
  1. GAFM (graph-augmented fine-mapping) — l1_finemap_from_sumstats
  2. HBP  (hierarchical belief propagation) — hbp_finemap_from_sumstats
  3. SuSiE (susieR::susie_rss) via R subprocess
  4. SuSiE-inf (FinucaneLab susieinf) — sumstats with infinitesimal background
  5. FINEMAP-inf (FinucaneLab finemapinf) — SSS sampler with infinitesimal background

For each trait × lead in data/rice_3k/results/grain_gwas/grain_lead_loci.tsv,
loads sumstats + builds in-sample LD from the 3kRG VCF (250 kb window),
runs all 5 methods, and writes:
  data/rice_3k/results/grain_finemap/<trait>_<chr>_<pos>.json     (per-locus full result)
  data/rice_3k/results/grain_finemap/grain_finemap_summary.tsv     (one row per (locus, method))

By default fine-maps the top 5 GW + top 2 suggestive leads per trait
(adjustable via env vars MAX_GW_PER_TRAIT, MAX_SUG_PER_TRAIT).

Uses the existing rice multi-omics graph cache (Chr*.json) for the GAFM/HBP
graph prior. SuSiE/SuSiE-inf/FINEMAP-inf are run with no graph prior.
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
import time
from dataclasses import asdict
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, "/mnt/data/GraphGWAS/src/python")
from graphgwas.finemapping_v2 import (  # noqa: E402
    hbp_finemap_from_sumstats,
    l1_finemap_from_sumstats,
)

import susieinf  # noqa: E402
import finemapinf  # noqa: E402

DATA = Path("/mnt/data/GraphGWAS/data/rice_3k")
RES = DATA / "results"
LEADS = RES / "grain_gwas" / "grain_lead_loci.tsv"
GWAS_OUT = RES / "grain_gwas"
FM_OUT = RES / "grain_finemap"
FM_OUT.mkdir(parents=True, exist_ok=True)
ANNO = DATA / "annotations"

VCF = "/mnt/data/GraphPop/data/raw/3kRG_data/NB_final_snp.vcf.gz"
FM_WINDOW = int(os.environ.get("FM_WINDOW", "100000"))  # ±100 kb fine-map window
MAX_GW_PER_TRAIT = int(os.environ.get("MAX_GW_PER_TRAIT", "5"))
MAX_SUG_PER_TRAIT = int(os.environ.get("MAX_SUG_PER_TRAIT", "2"))
L_EFFECTS = int(os.environ.get("L_EFFECTS", "5"))  # canonical fine-mapping prior
COVERAGE = 0.95


def _parse_gt(g):
    if g in ("0|0", "0/0"):
        return 0.0
    if g in ("1|0", "0|1", "1/0", "0/1"):
        return 1.0
    if g in ("1|1", "1/1"):
        return 2.0
    return np.nan


def load_locus(chrom, lo, hi):
    """bcftools-stream variants in window; impute NaN with column mean; MAF filter."""
    q = subprocess.run(
        ["bcftools", "query", "-r", f"{chrom}:{lo}-{hi}",
         "-f", "%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n", VCF],
        capture_output=True, text=True, timeout=600,
    )
    rows, mats = [], []
    for line in q.stdout.strip().split("\n"):
        if not line:
            continue
        f = line.split("\t")
        if len(f) < 5:
            continue
        d = np.array([_parse_gt(g) for g in f[4:]], dtype=np.float32)
        mu = np.nanmean(d) if np.isfinite(d).any() else 0.0
        d = np.where(np.isnan(d), mu, d)
        af = mu / 2.0 if 0 < mu < 2 else 0.0
        # MAF > 0.05 matches Niu 2021's fine-mapping convention and keeps the
        # 250-kb window tractable for the inf-method baselines (which scale
        # cubically in the number of variants).
        maf = min(af, 1 - af)
        if maf < 0.05:
            continue
        rows.append((f[0], int(f[1]), f[2], f[3], af))
        mats.append(d)
    if not rows:
        return pd.DataFrame(), np.zeros((0, 0))
    return pd.DataFrame(rows, columns=["chr", "pos", "ref", "alt", "af"]), np.vstack(mats)


_TRAIT_SUMSTATS_CACHE: dict = {}


def _load_trait_sumstats(trait: str):
    """Load + cache one trait's full sumstats (CHROM/POS/REF/ALT/BETA/SE/P/OBS_CT/z)."""
    if trait in _TRAIT_SUMSTATS_CACHE:
        return _TRAIT_SUMSTATS_CACHE[trait]
    gw_path = list(GWAS_OUT.glob(f"gwas_grain.{trait}.glm.linear"))
    if not gw_path:
        return None
    print(f"  [warming sumstats cache for {trait}]...", flush=True)
    t0 = time.time()
    # only the columns we need; cuts memory in half
    cols = ["#CHROM", "POS", "REF", "ALT", "TEST", "ERRCODE",
            "BETA", "SE", "P", "OBS_CT"]
    # Read TEST/ERRCODE/CHROM as category to cut memory; keep REF/ALT as str so
    # the variant-key concat works downstream.
    gw = pd.read_csv(gw_path[0], sep="\t", low_memory=False, usecols=cols,
                     dtype={"#CHROM": "category", "TEST": "category",
                            "ERRCODE": "category"})
    gw.columns = [c.lstrip("#") for c in gw.columns]
    gw = gw[(gw["TEST"] == "ADD") & (gw["ERRCODE"] == ".")].copy()
    gw["CHROM"] = gw["CHROM"].astype(str).apply(
        lambda c: c if c.startswith("Chr") else f"Chr{c}"
    )
    gw["REF"] = gw["REF"].astype(str)
    gw["ALT"] = gw["ALT"].astype(str)
    for c in ("BETA", "SE", "P", "OBS_CT"):
        gw[c] = pd.to_numeric(gw[c], errors="coerce")
    gw = gw.dropna(subset=["P", "BETA", "SE"])
    gw["z"] = gw["BETA"] / gw["SE"]
    # sort + index by chromosome for fast slicing
    gw = gw.sort_values(["CHROM", "POS"]).reset_index(drop=True)
    _TRAIT_SUMSTATS_CACHE[trait] = gw
    print(f"  [{trait}] {len(gw):,} variants cached in {time.time()-t0:.1f}s", flush=True)
    return gw


def _build_aligned_inputs(trait, lead_chrom, lead_pos):
    """Returns (variants_df, z_aligned, R_signed, R_squared, n_samples)."""
    lo, hi = lead_pos - FM_WINDOW, lead_pos + FM_WINDOW
    gw = _load_trait_sumstats(trait)
    if gw is None:
        return None
    gw_sub = gw[(gw["CHROM"] == lead_chrom) & (gw["POS"].between(lo, hi))].copy()
    if len(gw_sub) < 20:
        return None
    n_samples = int(gw_sub["OBS_CT"].median()) if "OBS_CT" in gw_sub.columns else 2453

    var_df, D = load_locus(lead_chrom, lo, hi)
    if len(var_df) < 20:
        return None
    mu = D.mean(axis=1, keepdims=True)
    sd = D.std(axis=1, keepdims=True)
    sd = np.where(sd < 1e-10, 1.0, sd)
    Z = (D - mu) / sd
    R = (Z @ Z.T) / Z.shape[1]  # signed correlation, used by SuSiE & inf-methods
    R_sq = np.clip(R ** 2, 0, 1).astype(np.float64)  # used by GAFM/HBP
    R = np.clip(R, -1.0, 1.0)
    np.fill_diagonal(R_sq, 1.0)
    np.fill_diagonal(R, 1.0)

    var_df["variant_id"] = (
        var_df["chr"] + ":" + var_df["pos"].astype(str)
        + ":" + var_df["ref"] + ":" + var_df["alt"]
    )
    gw_sub["key"] = (
        gw_sub["CHROM"].astype(str) + ":" + gw_sub["POS"].astype(str)
        + ":" + gw_sub["REF"] + ":" + gw_sub["ALT"]
    )
    z_map = dict(zip(gw_sub["key"].values, gw_sub["z"].values))
    z_flip = {
        f"{r['CHROM']}:{r['POS']}:{r['ALT']}:{r['REF']}": -r["z"]
        for _, r in gw_sub.iterrows()
    }
    z = np.full(len(var_df), np.nan)
    for i, k in enumerate(var_df["variant_id"].values):
        if k in z_map:
            z[i] = z_map[k]
        elif k in z_flip:
            z[i] = z_flip[k]
    keep = np.isfinite(z)
    if keep.sum() < 20:
        return None
    idx = np.where(keep)[0]
    var_df = var_df.loc[keep].reset_index(drop=True)
    z = z[keep]
    R_sq = R_sq[np.ix_(idx, idx)]
    R = R[np.ix_(idx, idx)]
    return var_df, z, R, R_sq, n_samples


def _run_susie_rss(z, R, n, L=L_EFFECTS, coverage=COVERAGE) -> dict:
    """SuSiE-RSS via R subprocess. Returns per-variant PIP + credible-set membership."""
    tmp = tempfile.mkdtemp(prefix="susie_rss_")
    np.savetxt(f"{tmp}/z.csv", z)
    np.savetxt(f"{tmp}/R.csv", R, delimiter=",")
    r_script = f'''
    library(susieR)
    z <- scan("{tmp}/z.csv")
    R <- as.matrix(read.csv("{tmp}/R.csv", header=FALSE))
    fit <- tryCatch(
        susie_rss(z=z, R=R, n={n}, L={L}, coverage={coverage}, verbose=FALSE),
        error=function(e) NULL
    )
    if(!is.null(fit)) {{
        pips <- fit$pip
        cs <- susie_get_cs(fit, coverage={coverage})
        write.table(data.frame(idx=seq_along(pips), pip=pips),
                    file="{tmp}/pips.tsv", sep="\\t", row.names=FALSE, quote=FALSE)
        if(length(cs$cs) > 0) {{
            cs_lines <- character()
            for(i in seq_along(cs$cs)) {{
                cs_lines <- c(cs_lines,
                    paste(c(i, length(cs$cs[[i]]), cs$cs[[i]]), collapse=","))
            }}
            writeLines(cs_lines, "{tmp}/cs.txt")
        }} else {{
            writeLines(character(0), "{tmp}/cs.txt")
        }}
    }} else {{
        writeLines("FAILED", "{tmp}/status.txt")
    }}
    '''
    t0 = time.time()
    proc = subprocess.run(
        ["R", "--no-save", "--no-restore", "-e", r_script],
        capture_output=True, text=True, timeout=300,
    )
    dt = time.time() - t0
    if not Path(f"{tmp}/pips.tsv").exists():
        return {"pips": None, "cs_indices": [], "time": dt, "error": (proc.stderr or "")[-500:]}
    pips = pd.read_csv(f"{tmp}/pips.tsv", sep="\t")["pip"].values
    cs_indices = []
    if Path(f"{tmp}/cs.txt").exists():
        for line in open(f"{tmp}/cs.txt"):
            line = line.strip()
            if not line:
                continue
            parts = [int(x) for x in line.split(",")]
            cs_indices.append([i - 1 for i in parts[2:]])  # R 1-indexed → Python 0-indexed
    return {"pips": pips, "cs_indices": cs_indices, "time": dt, "error": None}


def _aggregate_pip(raw_pip):
    """Convert a (N,L) per-effect alpha matrix into a (N,) per-variant PIP.
    Handle NaN values by treating them as 0 in the aggregation."""
    if raw_pip.ndim == 2:
        clean = np.nan_to_num(raw_pip, nan=0.0, posinf=1.0, neginf=0.0)
        clean = np.clip(clean, 0.0, 1.0)
        pips = 1 - np.prod(1 - clean, axis=1)
    else:
        pips = np.nan_to_num(raw_pip, nan=0.0)
        pips = np.clip(pips, 0.0, 1.0)
    return pips


def _run_susie_inf(z, R, n, L=L_EFFECTS, coverage=COVERAGE) -> dict:
    t0 = time.time()
    # MLE variance estimation is more robust than moments at extreme z-scores
    # (rice grain QTNs with p<1e-150 generate z>25 that breaks the moments path).
    try:
        fit = susieinf.susie(z=z, meansq=1.0, n=n, L=L, LD=R,
                             method="MLE", verbose=False, maxiter=50)
        raw_pip = fit["PIP"] if isinstance(fit, dict) else fit.PIP
        if raw_pip is None:
            raise ValueError("susieinf returned PIP=None (convergence failure)")
        pips = _aggregate_pip(raw_pip)
        cs_indices = []
        if raw_pip.ndim == 2:
            for li in range(raw_pip.shape[1]):
                col = np.nan_to_num(raw_pip[:, li], nan=0.0)
                col = np.clip(col, 0.0, 1.0)
                order = np.argsort(-col)
                cum = np.cumsum(col[order])
                if cum[-1] < 1e-6:
                    continue
                target = coverage * cum[-1]
                k = int(np.searchsorted(cum, target)) + 1
                cs_indices.append(sorted(order[:k].tolist()))
        if not cs_indices:
            order = np.argsort(-pips)
            cum = np.cumsum(pips[order])
            k = int(np.searchsorted(cum, coverage)) + 1 if cum[-1] > 0 else len(pips)
            cs_indices = [sorted(order[:k].tolist())]
        return {"pips": pips, "cs_indices": cs_indices, "time": time.time() - t0, "error": None}
    except Exception as e:
        return {"pips": None, "cs_indices": [], "time": time.time() - t0, "error": str(e)}


def _run_finemap_inf(z, R, n, L=L_EFFECTS, coverage=COVERAGE) -> dict:
    t0 = time.time()
    try:
        fit = finemapinf.finemap(z=z, meansq=1.0, n=n, L=L, LD=R,
                                 verbose=0, sched_sss=[50, 50, 1000])
        raw_pip = fit["PIP"] if isinstance(fit, dict) else fit.PIP
        if raw_pip is None:
            raise ValueError("finemapinf returned PIP=None")
        pips = _aggregate_pip(raw_pip)
        order = np.argsort(-pips)
        cum = np.cumsum(pips[order])
        if cum[-1] < 1e-6:
            cs_indices = [sorted(order[: min(len(pips), 1)].tolist())]
        else:
            target = coverage * cum[-1]
            k = int(np.searchsorted(cum, target)) + 1
            cs_indices = [sorted(order[:k].tolist())]
        return {"pips": pips, "cs_indices": cs_indices, "time": time.time() - t0, "error": None}
    except Exception as e:
        return {"pips": None, "cs_indices": [], "time": time.time() - t0, "error": str(e)}


def _per_method_rows(r, trait, chrom, pos, sig, lead_p):
    """Helper to flatten a per-locus json into one summary row per method."""
    rows = []
    for m, mr in r["methods"].items():
        rows.append({
            "trait": trait,
            "lead_chr": chrom,
            "lead_pos": pos,
            "sig_level": sig,
            "lead_p": lead_p,
            "n_variants": r["n_variants"],
            "cache_coverage": r["cache_coverage"],
            "n_samples": r["n_samples"],
            "method": m,
            "n_credible_sets": mr.get("n_credible_sets"),
            "cs_size": mr.get("cs_size"),
            "cs_size_aggregated": mr.get("cs_size_aggregated"),
            "top_pip": mr.get("top_pip"),
            "top_variant": mr.get("top_variant"),
            "time_s": mr.get("time_s"),
            "error": mr.get("error"),
        })
    return rows


def _summarize_method(method, variants, pips, cs_indices, dt, top_k=5):
    if pips is None:
        return {
            "method": method, "cs_size": None, "cs_size_aggregated": None,
            "top_pip": None, "top_variant": None,
            "top_variants": [], "time_s": round(dt, 3), "error": True,
        }
    pips = np.asarray(pips, dtype=float)
    order = np.argsort(-pips)
    top = [(variants[i]["variantId"], float(pips[i])) for i in order[:top_k]]
    # Union of all CS indices
    cs_union = sorted(set(j for cs in cs_indices for j in cs))
    return {
        "method": method,
        "n_credible_sets": len(cs_indices),
        "cs_size": len(cs_indices[0]) if cs_indices else None,
        "cs_size_aggregated": len(cs_union),
        "top_pip": float(pips[order[0]]),
        "top_variant": variants[order[0]]["variantId"],
        "top_variants": top,
        "cs_variants": [variants[i]["variantId"] for i in cs_union],
        "time_s": round(dt, 3),
        "error": False,
    }


def run_locus(trait, lead_chrom, lead_pos, lead_p, sig_level, chr_cache):
    t_io = time.time()
    out = _build_aligned_inputs(trait, lead_chrom, lead_pos)
    if out is None:
        return None
    var_df, z, R, R_sq, n_samples = out
    print(f"  loaded sumstats+LD ({len(var_df)} variants, n={n_samples}) in {time.time()-t_io:.1f}s",
          flush=True)

    variants = [
        {
            "variantId": v["variant_id"],
            "chr": v["chr"], "pos": int(v["pos"]),
            "ref": v["ref"], "alt": v["alt"],
            "af_total": float(v["af"]),
        }
        for _, v in var_df.iterrows()
    ]
    window_cache = (
        {k: chr_cache[k] for k in var_df["variant_id"].values if k in chr_cache}
        if chr_cache else {}
    )

    method_results = {}

    def _log(method, mr):
        cs = mr.get("cs_size", "—")
        tp = mr.get("top_pip", None)
        tp_s = f"{tp:.3f}" if isinstance(tp, float) else "—"
        t = mr.get("time_s", "—")
        err = " ERROR" if mr.get("error") else ""
        print(f"    {method:12s} CS={cs} topPIP={tp_s} t={t}s{err}", flush=True)

    # GAFM
    t0 = time.time()
    try:
        gafm = l1_finemap_from_sumstats(
            variants, z, R_sq, z_func=None, alpha=0.7,
            r2_smooth=0.3, credible_set_coverage=COVERAGE, chr_name=lead_chrom,
            graph_cache=window_cache,
        )
        dt = time.time() - t0
        pips = np.array([c.pip for c in gafm])
        cs = [[i for i, c in enumerate(gafm) if c.in_credible_set]]
        method_results["GAFM"] = _summarize_method("GAFM", variants, pips, cs, dt)
    except Exception as e:
        method_results["GAFM"] = {"method": "GAFM", "error": True, "error_msg": str(e),
                                   "time_s": round(time.time() - t0, 3)}
    _log("GAFM", method_results["GAFM"])

    # HBP
    t0 = time.time()
    try:
        hbp = hbp_finemap_from_sumstats(
            variants, z, R_sq, graph_cache=window_cache,
            r2_smooth=0.3, credible_set_coverage=COVERAGE, chr_name=lead_chrom,
        )
        dt = time.time() - t0
        pips = np.array([c.pip for c in hbp])
        cs = [[i for i, c in enumerate(hbp) if c.in_credible_set]]
        method_results["HBP"] = _summarize_method("HBP", variants, pips, cs, dt)
    except Exception as e:
        method_results["HBP"] = {"method": "HBP", "error": True, "error_msg": str(e),
                                  "time_s": round(time.time() - t0, 3)}
    _log("HBP", method_results["HBP"])

    # SuSiE-RSS
    susie = _run_susie_rss(z, R, n_samples, L=L_EFFECTS, coverage=COVERAGE)
    method_results["SuSiE"] = _summarize_method(
        "SuSiE", variants, susie["pips"], susie["cs_indices"], susie["time"]
    )
    if susie["error"]:
        method_results["SuSiE"]["error_msg"] = susie["error"]
    _log("SuSiE", method_results["SuSiE"])

    # SuSiE-inf
    sinf = _run_susie_inf(z, R, n_samples, L=L_EFFECTS, coverage=COVERAGE)
    method_results["SuSiE-inf"] = _summarize_method(
        "SuSiE-inf", variants, sinf["pips"], sinf["cs_indices"], sinf["time"]
    )
    if sinf["error"]:
        method_results["SuSiE-inf"]["error_msg"] = sinf["error"]
    _log("SuSiE-inf", method_results["SuSiE-inf"])

    # FINEMAP-inf
    finf = _run_finemap_inf(z, R, n_samples, L=L_EFFECTS, coverage=COVERAGE)
    method_results["FINEMAP-inf"] = _summarize_method(
        "FINEMAP-inf", variants, finf["pips"], finf["cs_indices"], finf["time"]
    )
    if finf["error"]:
        method_results["FINEMAP-inf"]["error_msg"] = finf["error"]
    _log("FINEMAP-inf", method_results["FINEMAP-inf"])

    return {
        "trait": trait,
        "lead_chr": lead_chrom,
        "lead_pos": lead_pos,
        "lead_p": lead_p,
        "sig_level": sig_level,
        "n_variants": len(variants),
        "cache_coverage": len(window_cache),
        "n_samples": n_samples,
        "methods": method_results,
    }


def main():
    if not LEADS.exists():
        print(f"No lead loci at {LEADS}")
        return
    leads = pd.read_csv(LEADS, sep="\t")
    print(f"Loaded {len(leads)} total leads across {leads['trait'].nunique()} traits")

    def pick(g):
        gw = g[g["sig_level"] == "GW"].head(MAX_GW_PER_TRAIT)
        sug = g[g["sig_level"] == "suggestive"].head(MAX_SUG_PER_TRAIT)
        return pd.concat([gw, sug])

    leads_top = (
        leads.sort_values("p")
        .groupby("trait", group_keys=False)
        .apply(pick)
        .reset_index(drop=True)
    )

    # Augment with leads matched to Niu 2021 QTNs (any lead within ±100 kb of a
    # NIU QTN of the same trait). This ensures we fine-map the literature-relevant
    # loci even if they don't make the top-K marginal-p ranking.
    niu_match_path = GWAS_OUT / "grain_niu2021_match.tsv"
    if niu_match_path.exists():
        nm = pd.read_csv(niu_match_path, sep="\t")
        nm = nm[nm["matched"].astype(str).str.lower().isin(["true", "1"])]
        # Pick the lead variant from leads.tsv: same trait + same chr + within 100 kb of QTN
        niu_leads_rows = []
        for _, q in nm.iterrows():
            trait = q["trait"]; qc = q["qtn_chr"]; qp = int(q["qtn_pos"])
            cand = leads[(leads["trait"] == trait) & (leads["chrom"] == qc)
                         & (leads["pos"].between(qp - 100_000, qp + 100_000))]
            if cand.empty:
                continue
            best = cand.sort_values("p").iloc[0]
            niu_leads_rows.append(best.to_dict())
        niu_leads = pd.DataFrame(niu_leads_rows)
        # Merge: keep top-K leads + Niu-matched, dedupe by (trait, chrom, pos)
        combined = pd.concat([leads_top, niu_leads], ignore_index=True)
        combined = combined.drop_duplicates(subset=["trait", "chrom", "pos"], keep="first")
        leads_final = combined.sort_values(["trait", "p"]).reset_index(drop=True)
        n_added = len(leads_final) - len(leads_top)
        print(f"After cap ({MAX_GW_PER_TRAIT} GW + {MAX_SUG_PER_TRAIT} sug + Niu-matched): "
              f"{len(leads_final)} ({len(leads_top)} top-K + {n_added} Niu-augment)")
    else:
        leads_final = leads_top
        print(f"After cap (max {MAX_GW_PER_TRAIT} GW + {MAX_SUG_PER_TRAIT} sug per trait): "
              f"{len(leads_final)}")
    leads = leads_final

    chr_caches = {}

    def get_chr_cache(chrom):
        if chrom not in chr_caches:
            p = ANNO / f"rice_graph_cache_{chrom}.json"
            if p.exists():
                t0 = time.time()
                print(f"  loading {p.name} ({p.stat().st_size/1e6:.0f} MB)...", flush=True)
                chr_caches[chrom] = json.load(open(p))
                dt = time.time() - t0
                print(f"  cache loaded ({len(chr_caches[chrom]):,} entries) in {dt:.0f}s", flush=True)
            else:
                chr_caches[chrom] = {}
        return chr_caches[chrom]

    summary_rows = []
    all_results = []

    for _, lead in leads.iterrows():
        trait = lead["trait"]
        chrom = lead["chrom"]
        pos = int(lead["pos"])
        sig = lead["sig_level"]
        lead_p = float(lead["p"])
        out = FM_OUT / f"{trait}_{chrom}_{pos}.json"
        if out.exists() and out.stat().st_size > 200:
            # already fine-mapped — load and continue
            try:
                r = json.load(open(out))
                all_results.append(r)
                print(f"\n=== {trait} {chrom}:{pos:,} [{sig}, p={lead_p:.2e}] (cached) ===", flush=True)
            except Exception:
                out.unlink()
            else:
                summary_rows.extend(_per_method_rows(r, trait, chrom, pos, sig, lead_p))
                continue
        print(f"\n=== {trait} {chrom}:{pos:,} [{sig}, p={lead_p:.2e}] ===", flush=True)
        try:
            r = run_locus(trait, chrom, pos, lead_p, sig, get_chr_cache(chrom))
        except Exception as e:
            print(f"  ERROR: {e}", flush=True)
            continue
        if r is None:
            print("  skipped (insufficient data)")
            continue
        all_results.append(r)
        # write per-locus json
        with open(out, "w") as f:
            json.dump(r, f, indent=2, default=str)
        # console summary line per method
        for m, mr in r["methods"].items():
            cs = mr.get("cs_size", "—")
            tp = mr.get("top_pip", None)
            tp_s = f"{tp:.3f}" if isinstance(tp, float) else "—"
            t = mr.get("time_s", "—")
            print(f"  {m:12s}  CS={cs}  topPIP={tp_s}  t={t}s", flush=True)
            summary_rows.append({
                "trait": trait,
                "lead_chr": chrom,
                "lead_pos": pos,
                "sig_level": sig,
                "lead_p": lead_p,
                "n_variants": r["n_variants"],
                "cache_coverage": r["cache_coverage"],
                "n_samples": r["n_samples"],
                "method": m,
                "n_credible_sets": mr.get("n_credible_sets"),
                "cs_size": mr.get("cs_size"),
                "cs_size_aggregated": mr.get("cs_size_aggregated"),
                "top_pip": mr.get("top_pip"),
                "top_variant": mr.get("top_variant"),
                "time_s": mr.get("time_s"),
                "error": mr.get("error"),
            })

    df = pd.DataFrame(summary_rows)
    out = FM_OUT / "grain_finemap_summary.tsv"
    df.to_csv(out, sep="\t", index=False)
    print(f"\nWrote {out}  ({len(df)} (locus, method) rows; {len(all_results)} loci)")


if __name__ == "__main__":
    main()
