"""Genome annotation loaders for GraphGWAS.

Loads genes (GENCODE), pathways, eQTL (GTEx), PPI (STRING), regulatory elements
(ENCODE), and conservation (PhyloP) into the Neo4j graph. Designed to be
re-runnable: uses MERGE semantics with stable keys.
"""
from __future__ import annotations

import gzip
import multiprocessing as mp
from pathlib import Path
from typing import Iterable

import pandas as pd

from . import config
from .db import GraphGWASConnection

GENCODE_V47_BASIC = (
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/"
    "release_47/gencode.v47.basic.annotation.gtf.gz"
)


def parse_gencode_genes(gtf_path: Path, biotypes: Iterable[str] = ("protein_coding",)) -> pd.DataFrame:
    """Parse a GENCODE GTF and return one row per gene.

    Returns columns: geneId, symbol, chr, start, end, strand, biotype.
    geneId is stripped of ENSG version suffix to match cross-version queries.
    """
    biotype_set = set(biotypes)
    rows: list[dict] = []
    with gzip.open(gtf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            attrs = {}
            for kv in parts[8].split(";"):
                kv = kv.strip()
                if not kv:
                    continue
                key, _, val = kv.partition(" ")
                attrs[key] = val.strip('"')
            if attrs.get("gene_type") not in biotype_set:
                continue
            gene_id = attrs["gene_id"].split(".")[0]
            rows.append(
                dict(
                    geneId=gene_id,
                    symbol=attrs.get("gene_name", gene_id),
                    chr=parts[0],
                    start=int(parts[3]),
                    end=int(parts[4]),
                    strand=parts[6],
                    biotype=attrs["gene_type"],
                )
            )
    return pd.DataFrame(rows)


def load_genes_to_neo4j(genes: pd.DataFrame, batch_size: int = 2000) -> int:
    """MERGE Gene nodes keyed by geneId. Returns count loaded."""
    records = genes.to_dict(orient="records")
    cypher = """
    UNWIND $rows AS row
    MERGE (g:Gene {geneId: row.geneId})
    SET g.symbol = row.symbol,
        g.chr = row.chr,
        g.start = row.start,
        g.end = row.end,
        g.strand = row.strand,
        g.biotype = row.biotype
    """
    with GraphGWASConnection(config.NEO4J_URI, config.NEO4J_USER, config.NEO4J_PASSWORD) as conn:
        for i in range(0, len(records), batch_size):
            conn.execute_write(cypher, {"rows": records[i : i + batch_size]})
    return len(records)


def _build_edges_for_chr(chr_name: str, flanking: int, batch_size: int) -> tuple[str, int]:
    """Build HAS_CONSEQUENCE edges for all genes on one chromosome.

    Position-based: any variant within [gene.start - flanking, gene.end + flanking]
    gets a HAS_CONSEQUENCE edge with consequence="proximity". Existing VEP edges
    (e.g., consequence="missense_variant") are preserved; MERGE only dedupes on
    (variant)-[:HAS_CONSEQUENCE {consequence:"proximity"}]->(gene).
    """
    cypher = """
    UNWIND $batch AS row
    MATCH (g:Gene {geneId: row.geneId})
    WITH g, row
    MATCH (v:Variant {chr: $chr})
    WHERE v.pos >= row.lo AND v.pos <= row.hi
    MERGE (v)-[r:HAS_CONSEQUENCE {consequence: 'proximity'}]->(g)
      ON CREATE SET r.distance_bp =
          CASE
            WHEN v.pos < g.start THEN g.start - v.pos
            WHEN v.pos > g.end   THEN v.pos - g.end
            ELSE 0
          END
    RETURN count(r) AS n
    """
    total_edges = 0
    with GraphGWASConnection(config.NEO4J_URI, config.NEO4J_USER, config.NEO4J_PASSWORD) as conn:
        genes = conn.execute_read(
            "MATCH (g:Gene {chr: $chr}) "
            "RETURN g.geneId AS geneId, g.start AS start, g.end AS end",
            {"chr": chr_name},
        ).data()
        if not genes:
            return chr_name, 0
        for i in range(0, len(genes), batch_size):
            batch = [
                dict(
                    geneId=g["geneId"],
                    lo=g["start"] - flanking,
                    hi=g["end"] + flanking,
                )
                for g in genes[i : i + batch_size]
            ]
            result = conn.execute_write(cypher, {"chr": chr_name, "batch": batch}).single()
            if result:
                total_edges += int(result["n"])
    return chr_name, total_edges


def build_proximity_edges(
    chrs: Iterable[str] | None = None,
    flanking_kb: int = 5,
    genes_per_batch: int = 20,
    parallelism: int = 6,
) -> dict[str, int]:
    """Build HAS_CONSEQUENCE (proximity) edges genome-wide via parallel workers."""
    if chrs is None:
        chrs = [f"chr{i}" for i in range(1, 23)]
    flanking = flanking_kb * 1000
    with mp.get_context("spawn").Pool(parallelism) as pool:
        results = pool.starmap(
            _build_edges_for_chr,
            [(c, flanking, genes_per_batch) for c in chrs],
        )
    return dict(results)


def load_gtex_tissue(
    signif_gz: Path, tissue: str, batch_size: int = 5000
) -> tuple[int, int]:
    """Load one GTEx tissue's significant eQTL pairs into Neo4j.

    Creates (:Variant)-[:eQTL {tissue, beta, pval, slope_se, tss_distance}]->(:Gene)
    edges. Variants not in the 1KG database are silently skipped (1KG is a subset
    of all variants GTEx tested). Returns (rows_read, edges_created).
    """
    cypher = """
    UNWIND $batch AS row
    MATCH (v:Variant {chr: row.chr, pos: row.pos, ref: row.ref, alt: row.alt})
    MATCH (g:Gene {geneId: row.gene_id})
    CREATE (v)-[r:eQTL {
        tissue: $tissue,
        beta: row.beta,
        pval: row.pval,
        slope_se: row.slope_se,
        tss_distance: row.tss_distance
    }]->(g)
    RETURN count(r) AS n
    """
    rows_read = 0
    edges_created = 0
    batch: list[dict] = []

    def _flush(conn):
        nonlocal edges_created, batch
        if not batch:
            return
        result = conn.execute_write(cypher, {"tissue": tissue, "batch": batch}).single()
        if result:
            edges_created += int(result["n"])
        batch = []

    with GraphGWASConnection(config.NEO4J_URI, config.NEO4J_USER, config.NEO4J_PASSWORD) as conn:
        with gzip.open(signif_gz, "rt") as fh:
            header = fh.readline().rstrip("\n").split("\t")
            idx = {c: header.index(c) for c in header}
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                variant_id = parts[idx["variant_id"]]
                vid_parts = variant_id.split("_")
                if len(vid_parts) < 4:
                    continue
                chr_, pos, ref, alt = vid_parts[0], int(vid_parts[1]), vid_parts[2], vid_parts[3]
                gene_id = parts[idx["gene_id"]].split(".")[0]
                batch.append(
                    dict(
                        chr=chr_,
                        pos=pos,
                        ref=ref,
                        alt=alt,
                        gene_id=gene_id,
                        beta=float(parts[idx["slope"]]),
                        pval=float(parts[idx["pval_nominal"]]),
                        slope_se=float(parts[idx["slope_se"]]),
                        tss_distance=int(parts[idx["tss_distance"]]),
                    )
                )
                rows_read += 1
                if len(batch) >= batch_size:
                    _flush(conn)
            _flush(conn)
    return rows_read, edges_created


def _load_gtex_tissue_task(args):
    signif_gz, tissue, batch_size = args
    try:
        rows, edges = load_gtex_tissue(signif_gz, tissue, batch_size)
        return tissue, rows, edges, None
    except Exception as exc:  # noqa: BLE001
        return tissue, 0, 0, str(exc)


def load_string_ppi(
    links_gz: Path,
    info_gz: Path,
    score_threshold: int = 700,
    batch_size: int = 5000,
    replace: bool = True,
) -> tuple[int, int]:
    """Load STRING PPI as (:Gene)-[:INTERACTS_WITH {score}]->(:Gene).

    Maps ENSP protein IDs to gene symbols via the info file, then filters links
    to combined_score >= score_threshold. De-duplicates directionally (stores one
    edge per unordered pair). Requires Gene nodes with populated `symbol` field.
    """
    ensp_to_symbol: dict[str, str] = {}
    with gzip.open(info_gz, "rt") as fh:
        header = fh.readline()
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                ensp_to_symbol[parts[0]] = parts[1]

    if replace:
        with GraphGWASConnection(config.NEO4J_URI, config.NEO4J_USER, config.NEO4J_PASSWORD) as conn:
            conn.execute_write(
                "CALL { MATCH ()-[r:INTERACTS_WITH]->() DELETE r } IN TRANSACTIONS OF 50000 ROWS;"
            )

    cypher = """
    UNWIND $batch AS row
    MATCH (g1:Gene {symbol: row.a})
    MATCH (g2:Gene {symbol: row.b})
    MERGE (g1)-[r:INTERACTS_WITH]->(g2)
      ON CREATE SET r.score = row.score
      ON MATCH  SET r.score = row.score
    RETURN count(r) AS n
    """
    rows_read = 0
    edges_created = 0
    batch: list[dict] = []
    seen: set[tuple[str, str]] = set()

    def _flush(conn):
        nonlocal edges_created, batch
        if not batch:
            return
        result = conn.execute_write(cypher, {"batch": batch}).single()
        if result:
            edges_created += int(result["n"])
        batch = []

    with GraphGWASConnection(config.NEO4J_URI, config.NEO4J_USER, config.NEO4J_PASSWORD) as conn:
        with gzip.open(links_gz, "rt") as fh:
            fh.readline()  # skip header
            for line in fh:
                parts = line.rstrip("\n").split(" ")
                if len(parts) != 3:
                    continue
                try:
                    score = int(parts[2])
                except ValueError:
                    continue
                if score < score_threshold:
                    continue
                a = ensp_to_symbol.get(parts[0])
                b = ensp_to_symbol.get(parts[1])
                if not a or not b or a == b:
                    continue
                key = (a, b) if a < b else (b, a)
                if key in seen:
                    continue
                seen.add(key)
                batch.append(dict(a=key[0], b=key[1], score=score))
                rows_read += 1
                if len(batch) >= batch_size:
                    _flush(conn)
            _flush(conn)
    return rows_read, edges_created


def load_encode_ccre(
    bed_path: Path,
    batch_size: int = 10000,
    replace: bool = True,
) -> tuple[int, int]:
    """Load ENCODE cCRE regulatory elements.

    Creates (:RegulatoryElement {id, chr, start, end, classification, category}) nodes.
    Expects UCSC bigBedToBed output columns: chr start end ccreId score . ... classification ... category ...
    """
    if replace:
        with GraphGWASConnection(config.NEO4J_URI, config.NEO4J_USER, config.NEO4J_PASSWORD) as conn:
            conn.execute_write(
                "CALL () { MATCH (re:RegulatoryElement) DETACH DELETE re } IN TRANSACTIONS OF 10000 ROWS;"
            )
    cypher = """
    UNWIND $batch AS row
    MERGE (re:RegulatoryElement {id: row.id})
    SET re.chr = row.chr,
        re.start = row.start,
        re.end = row.end,
        re.classification = row.classification,
        re.category = row.category
    RETURN count(re) AS n
    """
    rows = 0
    created = 0
    batch: list[dict] = []

    def _flush(conn):
        nonlocal created, batch
        if not batch:
            return
        result = conn.execute_write(cypher, {"batch": batch}).single()
        if result:
            created += int(result["n"])
        batch = []

    with GraphGWASConnection(config.NEO4J_URI, config.NEO4J_USER, config.NEO4J_PASSWORD) as conn:
        with open(bed_path) as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 13:
                    continue
                batch.append(
                    dict(
                        chr=parts[0],
                        start=int(parts[1]),
                        end=int(parts[2]),
                        id=parts[3],
                        classification=parts[9],
                        category=parts[12],
                    )
                )
                rows += 1
                if len(batch) >= batch_size:
                    _flush(conn)
            _flush(conn)
    return rows, created


def _build_ccre_edges_for_chr(chr_name: str, batch_size: int) -> tuple[str, int]:
    """Create (:Variant)-[:IN_REGULATORY]->(:RegulatoryElement) edges for one chr."""
    cypher = """
    UNWIND $batch AS row
    MATCH (re:RegulatoryElement {id: row.id})
    WITH re, row
    MATCH (v:Variant {chr: $chr})
    WHERE v.pos >= row.start AND v.pos <= row.end
    MERGE (v)-[r:IN_REGULATORY]->(re)
    RETURN count(r) AS n
    """
    total = 0
    with GraphGWASConnection(config.NEO4J_URI, config.NEO4J_USER, config.NEO4J_PASSWORD) as conn:
        ccres = conn.execute_read(
            "MATCH (re:RegulatoryElement {chr: $chr}) RETURN re.id AS id, re.start AS start, re.end AS end",
            {"chr": chr_name},
        ).data()
        if not ccres:
            return chr_name, 0
        for i in range(0, len(ccres), batch_size):
            result = conn.execute_write(cypher, {"chr": chr_name, "batch": ccres[i : i + batch_size]}).single()
            if result:
                total += int(result["n"])
    return chr_name, total


def build_ccre_variant_edges(chrs: Iterable[str] | None = None, batch_size: int = 100, parallelism: int = 6) -> dict[str, int]:
    """Build IN_REGULATORY edges in parallel across chromosomes."""
    if chrs is None:
        chrs = [f"chr{i}" for i in range(1, 23)]
    with mp.get_context("spawn").Pool(parallelism) as pool:
        return dict(pool.starmap(_build_ccre_edges_for_chr, [(c, batch_size) for c in chrs]))


def load_gtex_all_tissues(
    gtex_dir: Path, parallelism: int = 4, batch_size: int = 5000
) -> dict[str, tuple[int, int]]:
    """Load all 49 tissues' significant eQTL pairs from a GTEx v8 directory."""
    tasks = []
    for fp in sorted(gtex_dir.glob("*.signif_variant_gene_pairs.txt.gz")):
        tissue = fp.name.replace(".v8.signif_variant_gene_pairs.txt.gz", "")
        tasks.append((fp, tissue, batch_size))
    print(f"Loading {len(tasks)} tissues with parallelism={parallelism}")
    with mp.get_context("spawn").Pool(parallelism) as pool:
        out = {}
        for tissue, rows, edges, err in pool.imap_unordered(_load_gtex_tissue_task, tasks):
            if err:
                print(f"  [FAIL] {tissue}: {err}")
            else:
                print(f"  [OK]   {tissue}: {rows:,} rows, {edges:,} edges")
            out[tissue] = (rows, edges)
    return out


def cli() -> None:
    import argparse
    ap = argparse.ArgumentParser(description="GraphGWAS annotation loaders")
    sub = ap.add_subparsers(dest="cmd", required=True)

    g = sub.add_parser("load-gencode", help="Load GENCODE genes")
    g.add_argument("--gtf", required=True, type=Path)
    g.add_argument("--biotypes", nargs="+", default=["protein_coding"])

    e = sub.add_parser("build-proximity-edges", help="Build HAS_CONSEQUENCE proximity edges")
    e.add_argument("--chrs", nargs="+", default=None)
    e.add_argument("--flanking-kb", type=int, default=5)
    e.add_argument("--parallelism", type=int, default=6)
    e.add_argument("--genes-per-batch", type=int, default=20)

    gt = sub.add_parser("load-gtex", help="Load one GTEx v8 tissue")
    gt.add_argument("--signif-gz", required=True, type=Path)
    gt.add_argument("--tissue", required=True)
    gt.add_argument("--batch-size", type=int, default=5000)

    gall = sub.add_parser("load-gtex-all", help="Load all 49 GTEx v8 tissues from a directory")
    gall.add_argument("--gtex-dir", required=True, type=Path)
    gall.add_argument("--parallelism", type=int, default=4)
    gall.add_argument("--batch-size", type=int, default=5000)

    enc = sub.add_parser("load-encode-ccre", help="Load ENCODE cCRE regulatory elements")
    enc.add_argument("--bed", required=True, type=Path)
    enc.add_argument("--no-replace", action="store_true")
    enc.add_argument("--batch-size", type=int, default=10000)

    encE = sub.add_parser("build-ccre-edges", help="Build IN_REGULATORY variant edges")
    encE.add_argument("--chrs", nargs="+", default=None)
    encE.add_argument("--batch-size", type=int, default=100)
    encE.add_argument("--parallelism", type=int, default=6)

    s = sub.add_parser("load-string", help="Load STRING PPI high-confidence edges")
    s.add_argument("--links-gz", required=True, type=Path)
    s.add_argument("--info-gz", required=True, type=Path)
    s.add_argument("--score-threshold", type=int, default=700)
    s.add_argument("--no-replace", action="store_true", help="Don't delete existing INTERACTS_WITH first")
    s.add_argument("--batch-size", type=int, default=5000)

    args = ap.parse_args()
    if args.cmd == "load-gencode":
        df = parse_gencode_genes(args.gtf, biotypes=args.biotypes)
        print(f"parsed {len(df)} genes from {args.gtf}")
        n = load_genes_to_neo4j(df)
        print(f"loaded {n} Gene nodes")
    elif args.cmd == "build-proximity-edges":
        counts = build_proximity_edges(
            chrs=args.chrs,
            flanking_kb=args.flanking_kb,
            genes_per_batch=args.genes_per_batch,
            parallelism=args.parallelism,
        )
        total = sum(counts.values())
        print(f"built proximity edges: {total:,} total")
        for c, n in sorted(counts.items()):
            print(f"  {c}: {n:,}")
    elif args.cmd == "load-gtex":
        rows, edges = load_gtex_tissue(args.signif_gz, args.tissue, args.batch_size)
        print(f"{args.tissue}: {rows:,} rows read, {edges:,} eQTL edges created")
    elif args.cmd == "load-gtex-all":
        out = load_gtex_all_tissues(args.gtex_dir, args.parallelism, args.batch_size)
        total_edges = sum(v[1] for v in out.values())
        print(f"\nTotal: {len(out)} tissues, {total_edges:,} eQTL edges")
    elif args.cmd == "load-encode-ccre":
        rows, created = load_encode_ccre(args.bed, args.batch_size, replace=not args.no_replace)
        print(f"ENCODE cCRE: {rows:,} rows, {created:,} RegulatoryElement nodes MERGEd")
    elif args.cmd == "build-ccre-edges":
        counts = build_ccre_variant_edges(args.chrs, args.batch_size, args.parallelism)
        total = sum(counts.values())
        print(f"built IN_REGULATORY edges: {total:,} total")
        for c, n in sorted(counts.items()):
            print(f"  {c}: {n:,}")
    elif args.cmd == "load-string":
        rows, edges = load_string_ppi(
            args.links_gz, args.info_gz,
            score_threshold=args.score_threshold,
            batch_size=args.batch_size,
            replace=not args.no_replace,
        )
        print(f"STRING: {rows:,} pairs above threshold, {edges:,} edges MERGEd")


if __name__ == "__main__":
    cli()
