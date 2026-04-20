"""
Yeast Biological Discovery: HBP fine-mapping on real GWAS data.

Runs HBP and L1 on the top GWAS loci across 5 traits,
then queries gene annotations to identify biologically supported causal variants.

Usage:
    python tests/yeast_biological_discovery.py 2>/dev/null
"""

import csv, time, json, os
import numpy as np

import graphgwas.config as cfg
cfg._N_SAMPLES_DETECTED = False

from graphgwas.db import GraphGWASConnection
from graphgwas.config import detect_n_samples
from graphgwas.finemapping_v2 import (
    _load_locus_variants, _compute_association_stats, _compute_ld_matrix,
    _ld_deconvolve, _softmax, _build_credible_set,
    _cache_chr_graph_structure, fast_hbp_finemap,
)
from graphgwas.genotype import get_all_indices, get_phenotype_values

OUTDIR = "/mnt/data/GraphGWAS/results/yeast_discovery"
os.makedirs(OUTDIR, exist_ok=True)
RESULTS_DIR = "/mnt/data/GraphGWAS/results/grammar_corrected"

# Traits with known biology for validation
TRAITS = {
    "YPETHANOL": "Ethanol tolerance (ADH genes, SSU1, ALD6)",
    "YPDCUSO410MM": "Copper resistance (CUP1, CUP2, CTR1)",
    "YPD42": "Heat shock (HSP genes, trehalose pathway)",
    "YPDCAFEIN50": "Caffeine resistance (TOR pathway, PDR genes)",
    "YPGALACTOSE": "Galactose utilization (GAL pathway)",
}

WINDOW = 100000  # 100kb for yeast (compact genome)
TOP_LOCI_PER_TRAIT = 10


def load_gwas_results(trait):
    """Load GWAS results, return list sorted by p-value."""
    path = f"{RESULTS_DIR}/gwas_{trait}_grammar.tsv"
    hits = []
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            p = float(row["p_value"])
            if p < 1e-5:  # suggestive threshold
                hits.append({
                    "vid": row["variantId"], "chr": row["chr"],
                    "pos": int(row["pos"]), "p": p,
                    "logp": float(row["p_value_log10"]),
                    "beta": float(row["beta"]), "af": float(row["af"]),
                })
    return sorted(hits, key=lambda x: x["p"])


def cluster_into_loci(hits, window=100000):
    """Cluster hits into independent loci."""
    hits_sorted = sorted(hits, key=lambda x: (x["chr"], x["pos"]))
    loci = []
    for h in hits_sorted:
        merged = False
        for loc in loci:
            if h["chr"] == loc["chr"] and abs(h["pos"] - loc["lead_pos"]) < window:
                loc["hits"].append(h)
                if h["logp"] > loc["lead_logp"]:
                    loc["lead_pos"] = h["pos"]
                    loc["lead_logp"] = h["logp"]
                    loc["lead_vid"] = h["vid"]
                merged = True
                break
        if not merged:
            loci.append({
                "chr": h["chr"], "lead_pos": h["pos"],
                "lead_logp": h["logp"], "lead_vid": h["vid"],
                "hits": [h],
            })
    loci.sort(key=lambda x: -x["lead_logp"])
    return loci


def load_trait_phenotype(conn, trait, all_idx):
    """Load a trait's phenotype into gwas_value on Sample nodes."""
    prop_name = f"pheno_{trait}"

    # Read phenotype values from the trait-specific property
    result = conn.execute_read(f"""
        MATCH (s:Sample)
        WHERE s.`{prop_name}` IS NOT NULL
        RETURN s.packed_index AS idx, s.`{prop_name}` AS val
        ORDER BY s.packed_index
    """)

    pheno = np.full(len(all_idx), np.nan)
    idx_map = {int(all_idx[i]): i for i in range(len(all_idx))}
    for rec in result:
        i = idx_map.get(int(rec["idx"]))
        if i is not None and rec["val"] is not None:
            pheno[i] = float(rec["val"])

    # Store as gwas_value for genotype functions
    valid_data = [(int(all_idx[i]), float(pheno[i]))
                   for i in range(len(all_idx)) if not np.isnan(pheno[i])]
    for batch_start in range(0, len(valid_data), 500):
        batch = [{"idx": d[0], "val": d[1]} for d in valid_data[batch_start:batch_start+500]]
        conn.execute_write("""
            UNWIND $batch AS row
            MATCH (s:Sample {packed_index: row.idx})
            SET s.gwas_value = row.val, s.is_case = null, s.is_control = null
        """, {"batch": batch})

    n_valid = int((~np.isnan(pheno)).sum())
    return pheno, n_valid


def run_discovery():
    print("=" * 90)
    print("YEAST BIOLOGICAL DISCOVERY — HBP Fine-Mapping on Real GWAS Data")
    print("=" * 90)

    with GraphGWASConnection("bolt://localhost:7688", "neo4j", "graphgwas") as conn:
        detect_n_samples(conn)
        all_idx = get_all_indices(conn)
        print(f"N_SAMPLES={cfg.N_SAMPLES}, n_idx={len(all_idx)}")

        # Cache ALL chromosome graph structures (yeast has 16 chromosomes)
        print("\nCaching graph structure for all chromosomes...", flush=True)
        graph_caches = {}
        for i in range(1, 17):
            chr_name = f"chromosome{i}"
            t0 = time.time()
            gc = _cache_chr_graph_structure(conn, chr_name)
            graph_caches[chr_name] = gc
            if gc:
                print(f"  {chr_name}: {len(gc)} annotated variants ({time.time()-t0:.1f}s)")

        total_annotated = sum(len(v) for v in graph_caches.values())
        print(f"  Total: {total_annotated} annotated variants across 16 chromosomes\n")

        all_discoveries = {}

        for trait, description in TRAITS.items():
            print(f"\n{'='*90}")
            print(f"TRAIT: {trait} — {description}")
            print("="*90)

            # Load phenotype
            pheno, n_valid = load_trait_phenotype(conn, trait, all_idx)
            print(f"  Phenotype: {n_valid} samples with values")

            if n_valid < 100:
                print(f"  SKIP: too few samples")
                continue

            # Get GWAS hits and cluster into loci
            hits = load_gwas_results(trait)
            loci = cluster_into_loci(hits, WINDOW)
            print(f"  {len(hits)} suggestive hits in {len(loci)} independent loci")

            trait_results = []

            for loc_idx, loc in enumerate(loci[:TOP_LOCI_PER_TRAIT]):
                chr_name = loc["chr"]
                lead_pos = loc["lead_pos"]

                # Load variants in the locus
                variants = _load_locus_variants(
                    conn, chr_name, lead_pos, WINDOW, all_idx)
                if len(variants) < 5:
                    continue

                # L1 fine-mapping (no annotations)
                z_stats = _compute_association_stats(variants, pheno)
                R_sq = _compute_ld_matrix(variants)
                unique, n_nb = _ld_deconvolve(z_stats, R_sq, 0.3)
                l1_pip = _softmax(0.9 * unique)

                # HBP fine-mapping
                gc = graph_caches.get(chr_name, {})
                hbp_results = fast_hbp_finemap(
                    variants, pheno, gc,
                    chr_name=chr_name)

                hbp_pip = np.array([c.pip for c in hbp_results])
                hbp_vids = [c.variant_id for c in hbp_results]

                # Query gene annotations for top variants
                top_l1_idx = np.argsort(-l1_pip)[:5]
                top_hbp = hbp_results[:5]

                # Get gene info for all variants in top 10
                top_vids = set()
                for i in top_l1_idx:
                    top_vids.add(variants[i]["variantId"])
                for c in top_hbp:
                    top_vids.add(c.variant_id)

                gene_info = {}
                result = conn.execute_read("""
                    UNWIND $vids AS vid
                    MATCH (v:Variant {variantId: vid})-[:HAS_CONSEQUENCE]->(g:Gene)
                    OPTIONAL MATCH (g)-[:IN_PATHWAY]->(p:Pathway)
                    RETURN vid, g.symbol AS gene, g.description AS desc,
                           collect(DISTINCT p.name) AS pathways
                """, {"vids": list(top_vids)})
                for rec in result:
                    gene_info.setdefault(rec["vid"], []).append({
                        "gene": rec["gene"],
                        "desc": rec.get("desc", ""),
                        "pathways": [p for p in rec.get("pathways", []) if p],
                    })

                # Find cases where HBP and L1 disagree
                l1_top_vid = variants[top_l1_idx[0]]["variantId"]
                hbp_top_vid = hbp_results[0].variant_id
                disagree = l1_top_vid != hbp_top_vid

                locus_result = {
                    "locus": f"{chr_name}:{lead_pos}",
                    "n_variants": len(variants),
                    "lead_logp": loc["lead_logp"],
                    "l1_top": l1_top_vid,
                    "l1_top_pip": float(l1_pip[top_l1_idx[0]]),
                    "hbp_top": hbp_top_vid,
                    "hbp_top_pip": float(hbp_results[0].pip),
                    "disagree": disagree,
                    "l1_genes": gene_info.get(l1_top_vid, []),
                    "hbp_genes": gene_info.get(hbp_top_vid, []),
                }
                trait_results.append(locus_result)

                # Print locus summary
                tag = " *** DISAGREE ***" if disagree else ""
                print(f"\n  Locus {loc_idx+1}: {chr_name}:{lead_pos} "
                      f"({len(variants)} variants, -log10p={loc['lead_logp']:.1f}){tag}")

                print(f"    L1  top: {l1_top_vid} PIP={l1_pip[top_l1_idx[0]]:.4f}")
                for gi in gene_info.get(l1_top_vid, []):
                    pw = ", ".join(gi["pathways"][:3]) if gi["pathways"] else "none"
                    print(f"         → {gi['gene']} ({gi.get('desc','')[:60]}) [{pw}]")

                print(f"    HBP top: {hbp_top_vid} PIP={hbp_results[0].pip:.4f}")
                for gi in gene_info.get(hbp_top_vid, []):
                    pw = ", ".join(gi["pathways"][:3]) if gi["pathways"] else "none"
                    print(f"         → {gi['gene']} ({gi.get('desc','')[:60]}) [{pw}]")

                if disagree:
                    # Show both in context
                    l1_rank_of_hbp_top = next(
                        (i+1 for i in np.argsort(-l1_pip)
                         if variants[i]["variantId"] == hbp_top_vid), -1)
                    hbp_rank_of_l1_top = next(
                        (i+1 for i, c in enumerate(hbp_results)
                         if c.variant_id == l1_top_vid), -1)
                    print(f"    HBP's pick ranked #{l1_rank_of_hbp_top} by L1")
                    print(f"    L1's pick ranked #{hbp_rank_of_l1_top} by HBP")

                # Print top 5 HBP with gene annotations
                print(f"    {'Rank':<5s} {'Variant':<45s} {'PIP':>6s} {'Gene':<12s} {'Pathways'}")
                for rank, c in enumerate(hbp_results[:5], 1):
                    genes = gene_info.get(c.variant_id, [])
                    gname = genes[0]["gene"] if genes else "-"
                    pws = "; ".join(genes[0]["pathways"][:2]) if genes and genes[0]["pathways"] else "-"
                    print(f"    {rank:<5d} {c.variant_id:<45s} {c.pip:>6.4f} {gname:<12s} {pws}")

            all_discoveries[trait] = trait_results

        # Summary of disagreements
        print(f"\n{'='*90}")
        print("DISCOVERY SUMMARY — Loci Where HBP and L1 Disagree")
        print("="*90)

        total_loci = 0
        total_disagree = 0
        for trait, results in all_discoveries.items():
            n_disagree = sum(1 for r in results if r["disagree"])
            total_loci += len(results)
            total_disagree += n_disagree
            if n_disagree > 0:
                print(f"\n  {trait}: {n_disagree}/{len(results)} loci disagree")
                for r in results:
                    if r["disagree"]:
                        l1_gene = r["l1_genes"][0]["gene"] if r["l1_genes"] else "intergenic"
                        hbp_gene = r["hbp_genes"][0]["gene"] if r["hbp_genes"] else "intergenic"
                        print(f"    {r['locus']}: L1→{l1_gene} vs HBP→{hbp_gene}")

        print(f"\n  Total: {total_disagree}/{total_loci} loci with different top variants")

        # Save
        with open(f"{OUTDIR}/yeast_discovery.json", "w") as f:
            json.dump(all_discoveries, f, indent=2, default=str)
        print(f"\nSaved to {OUTDIR}/yeast_discovery.json")


if __name__ == "__main__":
    run_discovery()
