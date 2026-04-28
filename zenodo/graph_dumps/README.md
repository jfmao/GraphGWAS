# Graph-database dumps

Pre-built Neo4j 5.26 Community Edition dump files. The schema and
counts match those reported in the accompanying paper (Supplementary
Table S1 for human 1KG).

## Files

### `yeast_1011_v0.1.dump` (0.5 GB)
1,011 *Saccharomyces cerevisiae* isolates from the 1011 Yeast Genomes
Project (Peter *et al.* 2018):

- 1,919,547 biallelic variants (SNPs + indels)
- 34 populations, 16 chromosomes
- 6,613 SGD ORF gene annotations
- 679 GO pathway terms
- 1.18 M variant→gene HAS_CONSEQUENCE edges
- 50 K gene→pathway IN_PATHWAY edges
- 35 growth-trait phenotypes loaded (971/1011 strains matched)
- 6,517 gene expression + 630 protein-abundance phenotypes

### `human_1kg_multiomics_v0.1.dump` (17 GB)
1000 Genomes Phase 3 high-coverage (3,202 samples) with full
multi-omics annotation layer:

- 70,691,875 biallelic variants across 22 autosomes
- 20,092 GENCODE v47 protein-coding genes
- 43,170,154 GTEx v8 tissue-specific eQTL edges (49 tissues)
- 230,850 STRING v12 INTERACTS_WITH edges (combined score ≥ 700)
- 370,000 ENCODE cCRE v4 regulatory-element nodes
- 38,940,085 HAS_CONSEQUENCE edges (VEP-annotated)

## Restore instructions

```bash
# Prerequisites: Neo4j 5.26 Community Edition installed
# (https://neo4j.com/download-center/#community)

# Allocate enough memory in neo4j.conf BEFORE restoring the human dump:
echo "server.memory.heap.initial_size=16g" >> conf/neo4j.conf
echo "server.memory.heap.max_size=16g"      >> conf/neo4j.conf
echo "server.memory.pagecache.size=16g"     >> conf/neo4j.conf
echo "server.bolt.listen_address=:7688"     >> conf/neo4j.conf  # avoid default port collision

# Restore (the instance must be stopped first)
bin/neo4j stop
bin/neo4j-admin database load --from-path=/path/to/dump.dir neo4j --overwrite-destination

# Start
bin/neo4j start

# Verify
cypher-shell -a bolt://localhost:7688 -u neo4j -p graphgwas \
    'CALL apoc.meta.stats() YIELD labels, relTypes RETURN labels, relTypes'
```

Expected output for the human dump: `labels.Variant ≈ 70.7M`, 
`relTypes.HAS_CONSEQUENCE ≈ 39M`, `relTypes.eQTL ≈ 43M`, etc.

## SHA-256 checksums

(To be computed during the post-acceptance upload; see
`../MANIFEST.md`.)

## Licence

MIT for schema and derived annotations. Upstream data subject to
original licences: 1000 Genomes Project (IGSR unrestricted with
citation), 1011 Yeast Genomes (unrestricted), GENCODE (open), GTEx
(open-access sumstats), STRING (CC-BY 4.0), ENCODE (CC-BY 4.0).
