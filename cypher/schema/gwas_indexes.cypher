// GraphGWAS Index Definitions
// Run via: cypher-shell -u neo4j -p graphpop < cypher/schema/gwas_indexes.cypher

// Pattern A: Region scan (Layer 2 core loop)
CREATE RANGE INDEX variant_chr_pos IF NOT EXISTS FOR (v:Variant) ON (v.chr, v.pos);

// Pattern B: Case/control partitioning
CREATE INDEX sample_is_case IF NOT EXISTS FOR (s:Sample) ON (s.is_case);
CREATE INDEX sample_is_control IF NOT EXISTS FOR (s:Sample) ON (s.is_control);

// Pattern C: Allele frequency filtering (burden, rare variant)
CREATE RANGE INDEX variant_af_total IF NOT EXISTS FOR (v:Variant) ON (v.af_total);

// Pattern D: Functional annotation traversal
CREATE INDEX gene_symbol IF NOT EXISTS FOR (g:Gene) ON (g.symbol);
CREATE INDEX pathway_name IF NOT EXISTS FOR (p:Pathway) ON (p.name);

// Pattern E: Association result queries
CREATE INDEX assoc_run_id IF NOT EXISTS FOR (ar:AssociationResult) ON (ar.run_id);
CREATE RANGE INDEX assoc_pvalue IF NOT EXISTS FOR (ar:AssociationResult) ON (ar.p_value);
CREATE INDEX assoc_phenotype IF NOT EXISTS FOR (ar:AssociationResult) ON (ar.phenotype_key);
CREATE CONSTRAINT gwas_study_id IF NOT EXISTS FOR (gs:GWASStudy) REQUIRE gs.id IS UNIQUE;
