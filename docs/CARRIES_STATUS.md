# docs/CARRIES_STATUS.md — Genotype Representation Audit

**STATUS: COMPLETED — 2026-03-21**

This file documents the schema audit results from the live Neo4j database.
It is the single most important architectural document in GraphGWAS.

---

## Schema Audit Results

### Node Labels Present
```
Chromosome
GOTerm
Gene
GenomicWindow
Pathway
Population
Sample
Variant
```

### Relationship Types Present
```
HAS_CONSEQUENCE
HAS_GO_TERM
IN_PATHWAY
IN_POPULATION
NEXT
ON_CHROMOSOME
```

### Variant Node Properties (sample)
```
chr, pos, ref, alt, variant_type,
pop_ids, ac, an, af, het_count, hom_alt_count,
ac_total, an_total, af_total, call_rate,
het_exp, is_polarized,
gt_packed, phase_packed,        ← PACKED GENOTYPES (key finding)
variantId
```

**Sample values:**
- `variantId = "chr1:10390:CCCCTAA...:C"`
- `chr = "chr1"`, `pos = 10390`
- `af_total = 0.0075`, `call_rate = 1.0`
- `gt_packed` size = 801 bytes (2 bits × 3202 samples / 8 = 800.5 → 801)
- `phase_packed` size = 401 bytes (1 bit × 3202 samples / 8 = 400.25 → 401)

### Sample Node Properties (sample)
```
sampleId, population, sex, packed_index,
roh_chr{1..22}_n_roh, roh_chr{1..22}_froh, roh_chr{1..22}_total_kb
```

**Sample values:**
- `sampleId = "HG00096"`, `population = "GBR"`, `sex = 1`
- `packed_index = 0..3201` (contiguous, maps to byte position in gt_packed)

### CARRIES Edge Status
```
CARRIES edge count: 0
```

**Interpretation:**
- [ ] CARRIES edges exist (count > 0)
- [x] CARRIES edges absent — **genotypes stored as packed byte arrays on Variant nodes instead**

### Functional Annotation Status
```
Gene node count:              91,973
Pathway node count:            2,279
GOTerm node count:            18,858
HAS_CONSEQUENCE edge count:  891,689
HAS_GO_TERM edge count:       present (not counted)
IN_PATHWAY edge count:         present (not counted)
```
**Functional annotation layer is COMPLETE.**

### LD Edge Status
```
LD edge count: 0
```
**LD edges not pre-computed.** Will need to compute if needed for Phase 3.

### Sample Phenotype Status
```
Sample node count:                    3,202
Samples with phenotypes property:     0
```
**No phenotype data loaded yet.** Phenotype import is a Phase 1 task.

### Additional Counts
```
Variant node count:          70,691,875
Population node count:       26
Neo4j version:               2026.01.4 (Cypher 5/25)
GDS plugin:                  NOT INSTALLED
APOC plugin:                 NOT TESTED
```

---

## The gt_packed Encoding — Key Discovery

CARRIES edges were never created. Instead, GraphPop stores **full individual-level
genotypes as packed byte arrays directly on Variant nodes**. This is a 5th option
not anticipated in the original project design.

### Encoding Specification

**`gt_packed`** — 2 bits per sample, 4 samples per byte, LSB-first:

| Bits | Value | Meaning        | VCF |
|------|-------|----------------|-----|
| `00` |   0   | HOM_REF        | 0/0 |
| `01` |   1   | HET            | 0/1 |
| `10` |   2   | HOM_ALT        | 1/1 |
| `11` |   3   | MISSING        | ./. |

**Decode formula (Java):**
```java
int gt = (gtPacked[sampleIdx >> 2] >> ((sampleIdx & 3) << 1)) & 0x03;
```

**`phase_packed`** — 1 bit per sample, 8 samples per byte:
- `0` = ALT on haplotype 0 (1|0)
- `1` = ALT on haplotype 1 (0|1)

**`packed_index`** on Sample nodes = VCF column index (0-based).
Maps sample → byte position in gt_packed.

### Source Files (in /mnt/data/GraphPop)
- `graphpop-procedures/src/main/java/org/graphpop/procedures/PackedGenotypeReader.java`
- `graphpop-procedures/src/main/java/org/graphpop/procedures/GenotypeLoader.java`
- `graphpop-import/src/graphpop_import/vcf_parser.py`
- `graphpop-import/src/graphpop_import/csv_emitter.py`

---

## Known Reason for CARRIES Absence

CARRIES edges were deliberately not created in GraphPop. The packed genotype
arrays were chosen instead for population genomics statistics (pi, Tajima's D,
F_ST, SFS) where the access pattern is "read all genotypes for one variant
at a time" — perfectly served by a byte array on the Variant node.

At 70.7M variants × 3,202 samples, full CARRIES edges would be:
~70.7M × ~1,600 avg non-ref samples = ~113 billion edges — infeasible.

The packed array representation stores the same information in
70.7M × 801 bytes = ~56.7 GB of property data, which is manageable.

---

## Selected Genotype Representation Strategy

- [ ] **Option 1 — Full CARRIES edges** (infeasible at 70.7M variants × 3,202 samples)
- [ ] **Option 2 — Scoped CARRIES** (possible for candidate variants only)
- [ ] **Option 3 — Carrier Arrays on Variant nodes** (case_carrier_ids[], control_carrier_ids[])
- [ ] **Option 4 — External Parquet + graph index**
- [x] **Option 5 — Packed Genotype Arrays (EXISTING)** + Scoped CARRIES for traversal-heavy queries

**Selected option:** Option 5 — Hybrid: Use existing `gt_packed` as primary genotype
source for all GWAS computation, with optional Scoped CARRIES import for Phase 3
traversal operations (epistasis, cosegregation) on candidate loci only.

**Rationale:**
1. `gt_packed` already contains full individual-level genotypes for all 70.7M variants.
   No data import needed for genotype access.
2. For single-locus GWAS (Phase 2): decode gt_packed per variant, split by case/control
   using packed_index, run association test. This is the natural access pattern.
3. For multi-locus/epistasis (Phase 3): decode gt_packed to build carrier sets on-the-fly,
   or selectively create Scoped CARRIES edges for small candidate regions only.
4. GraphPop's existing Java procedures (PackedGenotypeReader, GenotypeLoader) already
   implement efficient decode — we can reuse or extend them.

**Estimated Scoped CARRIES (if needed for Phase 3):**
```
N phenotyped samples: 3,202 (all, initially)
N candidate variants: ~1,000–10,000 (post-GWAS significant + functional)
Estimated CARRIES edges: 3.2M–32M (feasible)
Estimated storage: 0.5–5 GB
Feasible on current hardware: YES
```

---

## Answers to 10 Unresolved Questions

1. **Genotype representation:** Packed byte arrays (`gt_packed`, `phase_packed`) on
   Variant nodes. 2 bits per sample genotype, 1 bit per sample phase. Full
   individual-level data for all 3,202 samples across all 70.7M variants.

2. **Sample nodes implemented? Phenotype properties?** YES, 3,202 Sample nodes exist
   with `sampleId`, `population`, `sex`, `packed_index`, and ROH statistics.
   NO phenotype data yet — `phenotypes` MAP property not set on any sample.

3. **Functional annotation (Gene/Pathway) implemented?** YES, fully implemented:
   91,973 Gene nodes, 2,279 Pathway nodes, 18,858 GOTerm nodes,
   891,689 HAS_CONSEQUENCE edges with `impact`, `consequence`, `feature_type`, `feature`.

4. **LD edges computed?** NO. Zero LD edges in the database. Will need computation
   if required for Phase 3 phenotype diffusion or fine-mapping.

5. **Reason CARRIES was dropped:** Deliberate architectural decision. Packed arrays
   on Variant nodes are more space-efficient for the "all genotypes for one variant"
   access pattern used by population genetics statistics. Full CARRIES would be
   ~113 billion edges — infeasible.

6. **Extending existing GraphPop DB or starting fresh?** EXTENDING the existing
   GraphPop database. All GraphPop data (Variant, Sample, Gene, Pathway, Population,
   GenomicWindow, population statistics) is already in place.

7. **Target cohort (rice/human/crop)?** HUMAN — 1000 Genomes Project Phase 3.
   3,202 samples, 26 populations (GBR, ASW, BEB, CDX, CEU, etc.).
   Population nodes have pre-computed summary statistics (pi, Tajima's D, Fay-Wu's H).

8. **Genome build and annotation version:** GRCh38 (inferred from chr1 variant positions
   and 1000 Genomes Phase 3). Annotation via VEP (consequence, impact on HAS_CONSEQUENCE edges).

9. **VEP output available?** YES — VEP annotations are already loaded as HAS_CONSEQUENCE
   edges with `consequence`, `impact`, `feature_type`, `feature` properties.
   891,689 variant-gene consequence edges exist.

10. **Hardware tier (laptop/cluster)?** Tier 1.5 — Workstation/Native Ubuntu. Neo4j allocated
    16 GB heap (`-Xmx16777216k`). 70.7M variants loaded. Sufficient for GWAS on
    3,202 samples. Not suitable for full CARRIES import at this scale.

---

## GWAS Implications Summary

| GWAS Capability | Feasibility | Approach |
|---|---|---|
| Single-locus association | READY | Decode gt_packed per variant, 2×2 table |
| Firth regression for rare variants | READY | Same decode, feed to Firth procedure |
| Population stratification (PCA) | NEEDS WORK | Compute from gt_packed or import from PLINK |
| Gene burden test | READY | Decode gt_packed for variants in gene, aggregate per sample |
| Pairwise epistasis (candidate loci) | READY | Decode gt_packed for both variants, co-carrier count |
| Genome-wide epistasis | INFEASIBLE | O(V^2) even with packed decode |
| Phenotype diffusion | NEEDS CARRIES | Build Scoped CARRIES for candidate region, then diffuse |
| Disease subgraph extraction | READY | Decode + enrichment scoring, no traversal needed |
| Kinship-aware cosegregation | NEEDS KINSHIP | No KINSHIP edges yet; compute or import |
| GNN export | READY | Decode gt_packed to build PyG edge_index |

---

## Sign-off

- [x] Schema audit complete
- [x] Strategy selected and documented
- [x] All 10 questions answered
- [ ] Reviewed with user before Phase 1 begins

**Completed:** Schema audit and CARRIES status analysis
**Date:** 2026-03-21
**Reviewed by:** User (pending)
