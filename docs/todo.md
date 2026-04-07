# GraphGWAS — Task List

## CURRENT PHASE: Phase 0 — Schema Audit & Environment Setup

### [ ] 0.1 Connect to Neo4j and run full schema audit
- [ ] Run `CALL db.labels()` — record all node labels
- [ ] Run `CALL db.relationshipTypes()` — record all relationship types
- [ ] Run `MATCH (v:Variant) RETURN keys(v) LIMIT 1` — record Variant properties
- [ ] Run `MATCH (s:Sample) RETURN keys(s) LIMIT 1` — record Sample properties
- [ ] Run `MATCH ()-[r:CARRIES]->() RETURN count(r) LIMIT 1` — CARRIES status (CRITICAL)
- [ ] Run `MATCH (g:Gene) RETURN count(g)` — functional annotation status
- [ ] Run `MATCH ()-[r:LD]->() RETURN count(r) LIMIT 1` — LD edge status

### [ ] 0.2 Create docs/CARRIES_STATUS.md with audit results
- [ ] Record all findings from 0.1
- [ ] Select genotype representation strategy (Options 1–4 from CLAUDE.md)
- [ ] Document reason for strategy choice

### [ ] 0.3 Verify environment
- [ ] Java 21+ version check
- [ ] Maven available
- [ ] Python 3.10+ available
- [ ] Neo4j version (`CALL dbms.components()`)
- [ ] GDS plugin (`CALL gds.version()`)

### [ ] 0.4 Create test data
- [ ] Generate or locate small test VCF (500–5000 samples)
- [ ] Create synthetic phenotype CSV (sample_id, case_control, covariates)
- [ ] Verify test data loadable

### [ ] 0.5 Answer all 10 unresolved questions from IMPLEMENTATION_PLAN.md
- [ ] Document answers in docs/CARRIES_STATUS.md

---

## UPCOMING: Phase 1 — Data Foundation
(Start only after Phase 0 is complete and CARRIES_STATUS.md is written)

### [ ] 1.1 Implement phenotype import pipeline
### [ ] 1.2 Implement genotype layer (strategy TBD from Phase 0)
### [ ] 1.3 Verify functional annotation layer
### [ ] 1.4 Create all GWAS indexes
### [ ] 1.5 Implement QC procedure
