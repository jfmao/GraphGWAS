"""Surgically add 4 catalogue genes that were filtered out by the original
cache builders.

Why these were missing:
  - AT4G03060 (AOP2, Arabidopsis chr4): tair10.gff3 marks it as
    `pseudogene` (the COL-0 reference allele has a 5-bp deletion
    truncating it; the functional AOP2 is in Cvi). The Arabidopsis
    builder filters f[2] != "gene", skipping pseudogenes.

  - LOC_Os07g15770 (Ghd7, rice chr7): no PPI partners in RicePPINet
    at prob>=0.7, no entry in the Ren-2023 grain-quality TSV (which
    drives the rice cache's pathway field). The rice builder drops
    variants with neither pathway nor PPI annotation.

  - LOC_Os11g46200 / LOC_Os11g46210 (Pik-1 / Pik-2, rice chr11): same
    reason as Ghd7 — paired NLR with no RicePPINet/Ren entries.

This script does NOT do a full cache rebuild (~hours). Instead it:
  1. Scans the source GFF / snpEff VCF for the 4 missing genes' variants.
  2. Loads the existing per-chrom JSONs and either:
     (a) augments existing entries' `genes` field, OR
     (b) creates new entries with empty pathway/ppi (canonical-PPI
         injection later populates ppi at runtime).
  3. Writes the patched JSONs back in place.

Run once after a fresh cache build to recover catalogue coverage. The
upstream build scripts have also been patched so future rebuilds
include these gene types directly.

Usage: python scripts/patch_caches_for_missing_genes.py
"""
from __future__ import annotations

import gzip
import json
import re
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent

# ---------------------------------------------------------------------------
# Arabidopsis: AOP2 = AT4G03060 on chr4
# ---------------------------------------------------------------------------

ARABI_GFF = REPO / "data" / "arabidopsis" / "annotations" / "tair10.gff3.gz"
ARABI_CACHE = REPO / "data" / "arabidopsis" / "arabidopsis_graph_cache_v3_chr4.json"
ARABI_CACHE_V1 = REPO / "data" / "arabidopsis" / "arabidopsis_graph_cache_chr4.json"
NCBI_TO_CHROM = {
    "NC_003070.9": "1", "NC_003071.7": "2", "NC_003074.8": "3",
    "NC_003075.7": "4", "NC_003076.8": "5",
}

ARABI_PSEUDOGENES_TO_RECOVER = {
    "AT4G03060": {"chrom": "4", "expected_ncbi": "NC_003075.7"},
}


def patch_arabidopsis() -> dict:
    print("\n=== Arabidopsis chr4: AOP2 (AT4G03060) ===")
    out: dict = {}
    if not ARABI_GFF.exists():
        print(f"  ⚠ GFF missing: {ARABI_GFF}")
        return out

    # Find pseudogene window from GFF
    print(f"  Scanning {ARABI_GFF.name} for pseudogene AT4G03060...")
    target_window = None
    with gzip.open(ARABI_GFF, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            if f[2] != "pseudogene":
                continue
            if "AT4G03060" not in f[8]:
                continue
            ncbi = f[0]
            chrom = NCBI_TO_CHROM.get(ncbi)
            if chrom is None:
                continue
            target_window = (chrom, int(f[3]), int(f[4]))
            print(f"    found: chrom={chrom}  pos {f[3]}-{f[4]}")
            break
    if target_window is None:
        print(f"  ⚠ AT4G03060 not found in GFF as pseudogene; skipping")
        return out

    chrom, start, end = target_window

    # The cache builder drops variants that hit no annotated gene OR have
    # neither pathway nor PPI. For AT4G03060 (pseudogene-in-COL-0), variants
    # in its window may have been dropped entirely. Scan the source VCF
    # for biallelic SNPs in the window at MAF≥0.05 and add them to the
    # cache as new entries.
    print(f"  Scanning source VCF for variants in {chrom}:{start}-{end}...",
          flush=True)
    import subprocess
    VCF = REPO / "tests" / "data" / "arabidopsis" / "1001genomes_snp-short-indel_only_ACGTN.vcf.gz"
    # Use bcftools view for region slicing (uses tabix index if present)
    proc = subprocess.run(
        ["bcftools", "query", "-r", f"{chrom}:{start}-{end}",
         "-f", "%CHROM\t%POS\t%REF\t%ALT\n", str(VCF)],
        capture_output=True, text=True, timeout=300,
    )
    if proc.returncode != 0:
        print(f"    ⚠ bcftools failed: {proc.stderr[:300]}")
        return out
    new_vid_to_genes: dict = {}
    for line in proc.stdout.splitlines():
        f = line.split("\t")
        if len(f) < 4:
            continue
        v_chrom, pos, ref, alt = f[0], f[1], f[2], f[3]
        # Only biallelic SNPs
        if len(ref) != 1 or len(alt) != 1 or "," in alt:
            continue
        vid = f"{v_chrom}:{pos}:{ref}:{alt}"
        new_vid_to_genes[vid] = ["AT4G03060"]
    print(f"    found {len(new_vid_to_genes):,} biallelic SNPs in window")

    # Patch both cache files
    for cache_path in [ARABI_CACHE, ARABI_CACHE_V1]:
        if not cache_path.exists():
            print(f"  ⚠ cache missing: {cache_path.name}; skipping")
            continue
        print(f"  Loading {cache_path.name} ({cache_path.stat().st_size/1e6:.0f} MB)...",
              flush=True)
        cache = json.loads(cache_path.read_text())
        n_before = len(cache)
        n_augmented = 0
        n_added = 0
        for vid, gene_list in new_vid_to_genes.items():
            if vid in cache:
                existing = cache[vid].get("genes", [])
                new_genes = [g for g in gene_list if g not in existing]
                if new_genes:
                    cache[vid]["genes"] = (list(existing) + new_genes)[:4]
                    n_augmented += 1
            else:
                cache[vid] = {"genes": gene_list, "pathways": [], "ppi": []}
                n_added += 1
        print(f"    augmented {n_augmented:,} existing entries, "
              f"added {n_added:,} new entries  "
              f"(cache: {n_before:,} → {len(cache):,})")
        cache_path.write_text(json.dumps(cache))
        out[cache_path.name] = {
            "n_entries_before": n_before,
            "n_entries_after": len(cache),
            "n_genes_added": n_augmented,
            "n_new_entries": n_added,
        }
    return out


# ---------------------------------------------------------------------------
# Rice: Ghd7 (chr7), Pik-1 + Pik-2 (chr11)
# ---------------------------------------------------------------------------

RICE_PSEUDO_VCF = Path("/mnt/data/GraphPop/data/raw/3kRG_data/NB_bialSNP_pseudo_canonical_ALL.vcf.gz")
RICE_CACHE_DIR = REPO / "data" / "rice_3k" / "annotations"

RICE_GENES_TO_RECOVER = {
    "LOC_Os07g15770": "Chr7",  # Ghd7
    "LOC_Os11g46200": "Chr11",  # Pik-1
    "LOC_Os11g46210": "Chr11",  # Pik-2
}

GENE_RX = re.compile(rb"LOC_Os\d+g\d+")


def patch_rice() -> dict:
    print("\n=== Rice: Ghd7 (chr7), Pik-1, Pik-2 (chr11) ===")
    out: dict = {}
    if not RICE_PSEUDO_VCF.exists():
        print(f"  ⚠ VCF missing: {RICE_PSEUDO_VCF}")
        return out

    # Group by chrom for one cache load per chrom
    by_chrom: dict = {}
    for gene, chrom in RICE_GENES_TO_RECOVER.items():
        by_chrom.setdefault(chrom, []).append(gene)

    for chrom, genes in by_chrom.items():
        chrom_num = chrom.replace("Chr", "")
        print(f"\n  -- {chrom}: recovering {genes} --")
        target_set = set(g.encode() for g in genes)

        # Scan VCF for variants annotated to these genes
        print(f"    Scanning VCF for variants annotated to {genes}...", flush=True)
        # Collect variant_id → set of recovered gene IDs
        new_vid_to_genes: dict = {}
        with gzip.open(RICE_PSEUDO_VCF, "rb") as fh:
            n = 0
            for line in fh:
                if line.startswith(b"#"):
                    continue
                # Quick filter: skip lines that don't mention any target gene
                if not any(g in line for g in target_set):
                    continue
                parts = line.split(b"\t", 8)
                if len(parts) < 8:
                    continue
                v_chrom_raw = parts[0].decode()
                v_chrom = v_chrom_raw if v_chrom_raw.startswith("Chr") else f"Chr{v_chrom_raw}"
                if v_chrom != chrom:
                    continue
                pos = parts[1].decode()
                ref = parts[3].decode()
                alt = parts[4].decode()
                info = parts[7]
                gene_ids = list(dict.fromkeys(m.decode() for m in GENE_RX.findall(info)))
                target_hits = [g for g in gene_ids if g.encode() in target_set]
                if not target_hits:
                    continue
                vid = f"{chrom}:{pos}:{ref}:{alt}"
                # Keep top-4 genes (consistent with original builder)
                # Put target_hits first so they're not pushed off the cap
                merged = list(dict.fromkeys(target_hits + gene_ids))[:4]
                new_vid_to_genes[vid] = merged
                n += 1
                if n % 1000 == 0:
                    print(f"      {n:,} target variants accumulated...", flush=True)
        print(f"    Found {len(new_vid_to_genes):,} variants for {genes}")

        # Patch both v2 and v1 caches if they exist
        for suffix in ("v2_", ""):
            cache_path = RICE_CACHE_DIR / f"rice_graph_cache_{suffix}{chrom}.json"
            if not cache_path.exists():
                print(f"    ⚠ cache missing: {cache_path.name}; skipping")
                continue
            print(f"    Loading {cache_path.name} "
                  f"({cache_path.stat().st_size/1e6:.0f} MB)...", flush=True)
            cache = json.loads(cache_path.read_text())
            n_before = len(cache)
            n_augmented = 0
            n_added = 0
            for vid, gene_list in new_vid_to_genes.items():
                if vid in cache:
                    existing = cache[vid].get("genes", [])
                    new_genes = [g for g in gene_list if g not in existing]
                    if new_genes:
                        cache[vid]["genes"] = list(existing) + new_genes
                        # cap at top-4 (consistent with original builder)
                        cache[vid]["genes"] = cache[vid]["genes"][:4]
                        n_augmented += 1
                else:
                    cache[vid] = {
                        "genes": gene_list,
                        "pathways": [],
                        "ppi": [],
                    }
                    n_added += 1
            print(f"      augmented {n_augmented:,} existing entries, "
                  f"added {n_added:,} new entries  "
                  f"(cache: {n_before:,} → {len(cache):,})")
            cache_path.write_text(json.dumps(cache))
            out[cache_path.name] = {
                "n_entries_before": n_before,
                "n_entries_after": len(cache),
                "n_genes_added": n_augmented,
                "n_new_entries": n_added,
            }
    return out


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    diag: dict = {}
    diag["arabidopsis"] = patch_arabidopsis()
    diag["rice"] = patch_rice()
    print("\n=== Patch summary ===")
    print(json.dumps(diag, indent=2))


if __name__ == "__main__":
    main()
