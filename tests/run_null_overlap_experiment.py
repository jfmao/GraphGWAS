"""Null-overlap permutation experiment for the cross-species 250 kb validation criterion.

For each species: sample 1000 random "leads" from the same chromosome distribution
as the observed lead set, and compute the 250 kb-window overlap rate against the
species' validation catalogue. Compare observed validation rate vs null distribution.

Output: results/null_overlap/null_overlap.{json,tsv}
"""
from __future__ import annotations

import json
import os
import random
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

OUT = Path("/mnt/data/GraphGWAS/results/null_overlap")
OUT.mkdir(parents=True, exist_ok=True)

VALIDATION_DIR = Path("/mnt/data/GraphGWAS/results/multispecies_summary")
RICE_GT = Path("/mnt/data/GraphGWAS/data/rice_3k/ground_truth/grain_quality_causal_genes.tsv")
ARABI_GT = Path("/mnt/data/GraphGWAS/tests/data/arabidopsis/aragwas_bonf_associations.csv")

WINDOW_KB = 250
N_PERMS = 1000
RNG = np.random.default_rng(2026)

# Approximate chromosome lengths (bp) from each species reference build
CHR_LENGTHS = {
    "yeast": {  # SGD chromosome lengths
        "chrI": 230218, "chrII": 813184, "chrIII": 316620, "chrIV": 1531933,
        "chrV": 576874, "chrVI": 270161, "chrVII": 1090940, "chrVIII": 562643,
        "chrIX": 439888, "chrX": 745751, "chrXI": 666816, "chrXII": 1078177,
        "chrXIII": 924431, "chrXIV": 784333, "chrXV": 1091291, "chrXVI": 948066,
    },
    "arabidopsis": {  # TAIR10
        "chr1": 30427671, "chr2": 19698289, "chr3": 23459830, "chr4": 18585056, "chr5": 26975502,
    },
    "rice": {  # IRGSP-1.0 / RAP-DB
        "Chr1": 43270923, "Chr2": 35937250, "Chr3": 36413819, "Chr4": 35502694,
        "Chr5": 29958434, "Chr6": 31248787, "Chr7": 29697621, "Chr8": 28443022,
        "Chr9": 23012720, "Chr10": 23207287, "Chr11": 29021106, "Chr12": 27531856,
    },
    "human": {  # GRCh38, autosomes only
        "1": 248956422, "2": 242193529, "3": 198295559, "4": 190214555, "5": 181538259,
        "6": 170805979, "7": 159345973, "8": 145138636, "9": 138394717, "10": 133797422,
        "11": 135086622, "12": 133275309, "13": 114364328, "14": 107043718, "15": 101991189,
        "16": 90338345, "17": 83257441, "18": 80373285, "19": 58617616, "20": 64444167,
        "21": 46709983, "22": 50818468,
    },
}


def load_validation_catalogue(species: str) -> dict[str, list[int]]:
    """Return {chr -> [position, ...]} for catalogue features."""
    catalog: dict[str, list[int]] = {}

    if species == "rice":
        df = pd.read_csv(RICE_GT, sep="\t")
        # columns: gene_id, chr, start, end, gene_symbol, ...
        for _, r in df.iterrows():
            ch = str(r.get("chr") or r.get("chromosome") or "").strip()
            try:
                start = int(r.get("start") or r.get("position") or 0)
                end = int(r.get("end") or start)
                pos = (start + end) // 2
            except Exception:
                continue
            if ch:
                catalog.setdefault(ch, []).append(pos)

    elif species == "arabidopsis":
        try:
            df = pd.read_csv(ARABI_GT)
        except Exception:
            return catalog
        for _, r in df.iterrows():
            ch = f"chr{r.get('chrom') or r.get('chr') or r.get('chromosome')}"
            try:
                pos = int(r.get("pos") or r.get("position") or 0)
            except Exception:
                continue
            ch = ch.replace("chrchr", "chr")
            catalog.setdefault(ch, []).append(pos)

    elif species == "yeast":
        # Use the validation TSV's own nearest_gt_gene column to derive
        # candidate gene positions (we don't have a separate yeast catalogue file).
        # Instead, sample at chromosome-length scale and use the observed lead density
        # as the null reference. We'll handle yeast as a special case below.
        pass

    elif species == "human":
        # Human catalogue is the GWAS-Catalog + lipid-genetics canon — not stored as a single
        # file in the repo. Use the validation TSV's nearest_gt_gene column instead:
        # any catalogue feature with a coordinate is implicit in the "validated" column.
        # Special-case below.
        pass

    return catalog


def load_validation_table(species: str) -> pd.DataFrame:
    if species == "human":
        path = VALIDATION_DIR / "validation_human_gw.tsv"
    else:
        path = VALIDATION_DIR / f"validation_{species}.tsv"
    return pd.read_csv(path, sep="\t")


def overlap_in_window(chrom: str, pos: int, catalog_positions: list[int], window_bp: int) -> bool:
    """Return True iff any catalog position is within window_bp of pos."""
    if not catalog_positions:
        return False
    for cp in catalog_positions:
        if abs(cp - pos) <= window_bp:
            return True
    return False


def run_species(species: str, vt: pd.DataFrame) -> dict:
    """Run the null-overlap permutation for one species."""
    print(f"\n=== {species} ===", flush=True)

    n_obs = len(vt)
    chr_col = "lead_chr"
    pos_col = "lead_pos"

    # observed validation rate from the TSV
    n_validated_obs = int(vt["validated"].sum())
    rate_obs = n_validated_obs / n_obs

    # build a chromosome distribution from the observed leads
    obs_chr_counts = Counter(vt[chr_col].astype(str))

    if species in ("rice", "arabidopsis"):
        catalog = load_validation_catalogue(species)
        chr_lens = CHR_LENGTHS[species]

        # Permute: keep chromosome distribution, randomise positions
        null_rates = []
        for _ in range(N_PERMS):
            n_match = 0
            for ch, count in obs_chr_counts.items():
                # canonical chr name in catalogue
                canon_ch = ch
                if species == "arabidopsis" and not canon_ch.startswith("chr"):
                    canon_ch = f"chr{canon_ch}"
                if species == "rice" and not canon_ch.startswith("Chr"):
                    canon_ch = f"Chr{canon_ch.lstrip('chr')}"
                length = chr_lens.get(canon_ch) or chr_lens.get(ch)
                if not length:
                    continue
                cat_pos = catalog.get(canon_ch) or catalog.get(ch) or []
                if not cat_pos:
                    continue
                rand_pos = RNG.integers(1, length, size=count)
                for p in rand_pos:
                    if overlap_in_window(canon_ch, int(p), cat_pos, WINDOW_KB * 1000):
                        n_match += 1
            null_rates.append(n_match / n_obs)

        null_rates_arr = np.array(null_rates)
        # empirical p-value: P(null >= observed)
        p_emp = float((null_rates_arr >= rate_obs).mean())
        p_emp_adj = max(p_emp, 1 / N_PERMS)  # avoid p=0

        return {
            "species": species,
            "n_leads": int(n_obs),
            "n_validated_obs": int(n_validated_obs),
            "rate_obs": round(float(rate_obs), 4),
            "null_mean": round(float(null_rates_arr.mean()), 4),
            "null_p25": round(float(np.percentile(null_rates_arr, 25)), 4),
            "null_p50": round(float(np.percentile(null_rates_arr, 50)), 4),
            "null_p75": round(float(np.percentile(null_rates_arr, 75)), 4),
            "null_p95": round(float(np.percentile(null_rates_arr, 95)), 4),
            "n_perms": int(N_PERMS),
            "p_empirical": round(p_emp_adj, 4),
            "enrichment_obs_over_null_mean": round(rate_obs / max(null_rates_arr.mean(), 1e-9), 2),
            "catalogue_features": sum(len(v) for v in catalog.values()),
            "method": "permute positions; same chromosome distribution; same n_obs",
        }

    elif species in ("yeast", "human"):
        # No coordinate-level catalogue file in the repo.
        # Instead, use the validation TSV itself: each row carries
        # nearest_gt_dist (in kb to nearest catalogue feature, NaN if none).
        # The "validated" column is True iff nearest_gt_dist <= 250 kb.
        # For null-overlap, we permute the lead positions while keeping
        # nearest_gt_dist distribution as the empirical catalogue density.
        # This is conservative (it preserves the catalogue density observed
        # at OBSERVED leads, which slightly understates the null overlap
        # rate at a uniformly random position).
        #
        # Pragmatic approach: shuffle the validated/not-validated flag
        # across leads N_PERMS times. This gives the chance-overlap rate
        # under the null that the binary validated outcome is independent
        # of the observed clustering — i.e., the random-sample rate.
        # NOTE: this is an upper bound on the random-position rate, since
        # observed leads are GW-significant and thus enriched near genes.
        # We label this as "permuted-validated-flag" rather than "null-position".
        validated = vt["validated"].astype(int).to_numpy()
        null_rates = []
        for _ in range(N_PERMS):
            null_rates.append(RNG.permutation(validated).mean())
        null_rates_arr = np.array(null_rates)
        # Under simple permutation, the null mean = observed mean (by construction).
        # So p_empirical = 0.5 trivially. Need a different baseline.
        # Use the global catalogue density: per the validation TSV's nearest_gt_dist,
        # what fraction of leads have ANY catalogue feature within 250 kb in a random
        # position? This is approximately sum(catalogue_kb)/genome_kb * 2 (window).
        #
        # Without a coordinate-level catalogue, we report the chromosome-marginal null
        # by sampling random positions per chromosome and asking how many random
        # positions would fall within 250 kb of A POSITION CHOSEN UNIFORMLY AT RANDOM
        # FROM THE OBSERVED CATALOGUE DENSITY. This is approximated by:
        #   density_per_chr = (n_validated_on_chr / chr_length_kb) * 250 * 2
        # i.e., for each random lead, prob = (window_kb * 2) / chr_length_kb * gene_density_factor
        # We can approximate using observed validation distance distribution.
        # For simplicity in this round, we use the leads' own nearest_gt_dist column:
        # null overlap = fraction of leads where ANY of N_OBS random catalogue features
        # would fall within 250 kb. As a best-effort approximation, we report the
        # observed rate and the permuted-flag null (mean = obs by construction)
        # plus a "random-position estimate" derived from chromosome lengths.
        chr_lens = CHR_LENGTHS[species]
        if species == "human":
            # observed: 322 leads in the GW table; chromosome distribution from vt
            chr_distribution = obs_chr_counts.most_common()
            # density estimate: assume catalogue density = 1 feature per 100 kb (genome-wide)
            # This is a coarse estimate; precise null requires the actual catalogue file.
            #   Genome length ~3 Gb, ~20k genes => ~150 kb/gene
            #   Within 250 kb of any gene => ~almost everywhere
            # So the "rate" depends critically on the catalogue's specificity.
            # We report the permuted-flag null and label as approximate.
            random_position_estimate = "see notes (coordinate catalogue not in repo)"
        else:
            random_position_estimate = "see notes (coordinate catalogue not in repo)"

        return {
            "species": species,
            "n_leads": int(n_obs),
            "n_validated_obs": int(n_validated_obs),
            "rate_obs": round(float(rate_obs), 4),
            "null_method": "permuted-validated-flag (placeholder; coordinate catalogue not in repo)",
            "null_mean": round(float(null_rates_arr.mean()), 4),
            "null_p25": round(float(np.percentile(null_rates_arr, 25)), 4),
            "null_p50": round(float(np.percentile(null_rates_arr, 50)), 4),
            "null_p75": round(float(np.percentile(null_rates_arr, 75)), 4),
            "null_p95": round(float(np.percentile(null_rates_arr, 95)), 4),
            "n_perms": int(N_PERMS),
            "note": ("Without a coordinate-level catalogue file, we permute the validated flag "
                     "across leads. This trivially gives null_mean = obs (by construction). "
                     "A coordinate-level catalogue would let us compute the random-position "
                     "null. Treat the per-species rate as descriptive, not as evidence of "
                     "fine-mapping accuracy beyond chance."),
            "random_position_estimate": random_position_estimate,
        }


def main() -> None:
    out: dict = {"window_kb": WINDOW_KB, "n_perms": N_PERMS, "results": {}}
    for species in ("yeast", "arabidopsis", "rice", "human"):
        try:
            vt = load_validation_table(species)
        except Exception as e:
            print(f"  {species}: validation TSV missing ({e})")
            continue
        out["results"][species] = run_species(species, vt)
        print(f"  {species}: rate_obs={out['results'][species]['rate_obs']:.3f}, "
              f"null_mean={out['results'][species]['null_mean']:.3f}")

    with (OUT / "null_overlap.json").open("w") as f:
        json.dump(out, f, indent=2)
    print(f"\nWrote {OUT}/null_overlap.json")

    # TSV summary
    cols = ["species", "n_leads", "n_validated_obs", "rate_obs",
            "null_mean", "null_p95", "p_empirical", "enrichment_obs_over_null_mean"]
    with (OUT / "null_overlap.tsv").open("w") as f:
        f.write("\t".join(cols) + "\n")
        for sp, d in out["results"].items():
            f.write("\t".join(str(d.get(c, "")) for c in cols) + "\n")
    print(f"Wrote {OUT}/null_overlap.tsv")


if __name__ == "__main__":
    main()
