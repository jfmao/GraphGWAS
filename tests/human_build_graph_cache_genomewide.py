"""Build a genome-wide human multi-omics graph_cache from GTEx + STRING.

Inputs:
  - GTEx_Analysis_v8_eQTL.tar (49 tissues × signif_variant_gene_pairs)
  - 9606.protein.links.v12.0.txt.gz (STRING v12 PPI)
  - 9606.protein.info.v12.0.txt.gz (ENSP→symbol)
  - gencode.v47.basic.annotation.gtf.gz (ENSG→symbol)

Schema:
  data/annotations/human_graph_cache_chr<N>.json (per chromosome)
  {variant_id_b38: {"genes": [symbol,...], "pathways": [], "ppi": [symbol,...]}}

The variant_id format matches the 1KG GRCh38 NYGC VCF: chr<N>:pos:ref:alt.
"""
from __future__ import annotations

import gzip
import json
import re
import tarfile
import time
from collections import defaultdict
from pathlib import Path

ANN = Path("/mnt/data/GraphGWAS/data/annotations")
GTEX_TAR = ANN / "GTEx_Analysis_v8_eQTL.tar"
STRING_LINKS = ANN / "9606.protein.links.v12.0.txt.gz"
STRING_INFO = ANN / "9606.protein.info.v12.0.txt.gz"
GENCODE = ANN / "gencode" / "gencode.v47.basic.annotation.gtf.gz"
OUT_DIR = ANN
STRING_THRESHOLD = 700


def load_ensg_to_symbol() -> dict[str, str]:
    """Parse GENCODE v47 GTF to build ENSG → symbol map."""
    print("[1/4] Parsing GENCODE v47 GTF for ENSG→symbol ...", flush=True)
    out: dict[str, str] = {}
    with gzip.open(GENCODE, "rt") as fh:
        for line in fh:
            if line.startswith("#"): continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene": continue
            attrs = f[8]
            ensg = sym = None
            for kv in attrs.split(";"):
                kv = kv.strip()
                if kv.startswith("gene_id "):
                    ensg = kv.split('"')[1].split(".")[0]
                elif kv.startswith("gene_name "):
                    sym = kv.split('"')[1]
            if ensg and sym:
                out[ensg] = sym
    print(f"  {len(out):,} ENSG→symbol mappings")
    return out


def load_ensp_to_symbol() -> dict[str, str]:
    """STRING info: ENSP → preferred_name (gene symbol)."""
    print("\n[2/4] Parsing STRING info ...", flush=True)
    out: dict[str, str] = {}
    with gzip.open(STRING_INFO, "rt") as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 2: continue
            ensp = f[0].split(".", 1)[1] if "." in f[0] else f[0]
            out[ensp] = f[1]
    print(f"  {len(out):,} ENSP→symbol mappings")
    return out


def load_string_ppi(ensp_to_sym: dict) -> dict[str, set[str]]:
    """STRING links → symbol↔symbol PPI graph at score≥THRESHOLD."""
    print(f"\n[3/4] Parsing STRING links (≥{STRING_THRESHOLD}) ...", flush=True)
    ppi: dict[str, set[str]] = defaultdict(set)
    n_total = 0; n_kept = 0
    with gzip.open(STRING_LINKS, "rt") as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split(" ")
            if len(f) < 3: continue
            n_total += 1
            try:
                s = int(f[2])
            except ValueError:
                continue
            if s < STRING_THRESHOLD: continue
            n_kept += 1
            a = f[0].split(".", 1)[1] if "." in f[0] else f[0]
            b = f[1].split(".", 1)[1] if "." in f[1] else f[1]
            sa, sb = ensp_to_sym.get(a), ensp_to_sym.get(b)
            if not sa or not sb: continue
            ppi[sa].add(sb); ppi[sb].add(sa)
    print(f"  {n_total:,} edges total; {n_kept:,} kept; {len(ppi):,} genes with ≥1 partner")
    return dict(ppi)


def stream_gtex_variant_genes(ensg_to_sym: dict) -> dict[str, dict]:
    """Stream GTEx tar; build {variant_id: {genes:[symbols], tissues:set}}.
    variant_id is GTEx-format chr<N>_pos_ref_alt_b38; we convert to chr<N>:pos:ref:alt.
    """
    print("\n[4/4] Streaming GTEx eQTL signif files ...", flush=True)
    var_genes: dict[str, set] = defaultdict(set)
    var_tissues: dict[str, set] = defaultdict(set)
    n_files = 0
    with tarfile.open(GTEX_TAR, "r") as tar:
        for member in tar.getmembers():
            if not member.name.endswith(".signif_variant_gene_pairs.txt.gz"):
                continue
            tissue = Path(member.name).name.replace(".v8.signif_variant_gene_pairs.txt.gz", "")
            n_files += 1
            t0 = time.time()
            f = tar.extractfile(member)
            if f is None: continue
            with gzip.open(f, "rt") as gz:
                hdr = gz.readline()  # header: variant_id\tgene_id\ttss_distance\t...
                cols = hdr.rstrip("\n").split("\t")
                v_idx = cols.index("variant_id")
                g_idx = cols.index("gene_id")
                n_rows = 0
                for line in gz:
                    f2 = line.rstrip("\n").split("\t")
                    if len(f2) <= max(v_idx, g_idx): continue
                    vraw = f2[v_idx]   # chr1_15211_T_C_b38
                    gid = f2[g_idx].split(".")[0]  # ENSG00000xxxxx
                    sym = ensg_to_sym.get(gid)
                    if not sym: continue
                    # Convert chr1_15211_T_C_b38 → chr1:15211:T:C
                    parts = vraw.rsplit("_", 4)
                    if len(parts) != 5: continue
                    chrom, pos, ref, alt, _ = parts
                    vid = f"{chrom}:{pos}:{ref}:{alt}"
                    var_genes[vid].add(sym)
                    var_tissues[vid].add(tissue)
                    n_rows += 1
            print(f"  {tissue:<35s} {n_rows:>9,} rows  "
                  f"({time.time()-t0:.0f}s, total cached: {len(var_genes):,})",
                  flush=True)
    print(f"\n  Total tissues processed: {n_files}")
    print(f"  Total annotated variants: {len(var_genes):,}")
    return var_genes, var_tissues


def main():
    t0 = time.time()
    ensg_to_sym = load_ensg_to_symbol()
    ensp_to_sym = load_ensp_to_symbol()
    ppi = load_string_ppi(ensp_to_sym)
    var_genes, var_tissues = stream_gtex_variant_genes(ensg_to_sym)

    print(f"\nBuilding per-chromosome caches ...", flush=True)
    by_chr: dict[str, dict] = defaultdict(dict)
    for vid, genes in var_genes.items():
        chrom = vid.split(":", 1)[0]  # chr1, chr2, ..., chrX, chrY
        if not re.match(r"^chr([1-9]|1[0-9]|2[0-2]|X|Y|M)$", chrom): continue
        # PPI partners (cap fan-out per variant)
        ppi_partners: list[str] = []
        for g in list(genes)[:4]:
            if g in ppi:
                for p in list(ppi[g])[:20]:
                    if p not in ppi_partners:
                        ppi_partners.append(p)
        gene_list = sorted(genes)[:4]
        if not gene_list and not ppi_partners:
            continue
        by_chr[chrom][vid] = {
            "genes": gene_list,
            "pathways": [],   # pathway layer not yet wired (Reactome download deferred)
            "ppi": ppi_partners[:30],  # cap
        }

    for chrom, cache in by_chr.items():
        out = OUT_DIR / f"human_graph_cache_{chrom}.json"
        with out.open("w") as fh:
            json.dump(cache, fh)
        print(f"  Wrote {out.name}: {len(cache):,} variants, "
              f"{out.stat().st_size/1e6:.1f} MB")
    print(f"\nTotal time: {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
