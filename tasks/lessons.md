# GraphGWAS — Lessons Learned

## How to Use This File
After ANY correction from the user, add an entry here describing:
- What went wrong
- Why it went wrong
- The rule to prevent it next time

Review this file at the start of every session.

---

## Session 1 Lessons (Pre-populated from Design Phase)

### LESSON 001: Always audit schema before writing graph code
**Pattern:** Attempted to write Cypher assuming CARRIES edges exist.
**Cause:** Design documents specify CARRIES but implementation dropped it for unknown reasons.
**Rule:** ALWAYS run the schema discovery protocol from GRAPH_SCHEMA_SKILL.md section 1
before writing any Cypher or Java stored procedure. Never assume a relationship type
or property exists without verifying it in the actual database.

### LESSON 002: CARRIES is the architecturally critical dependency
**Pattern:** Proceeded with multi-locus GWAS implementation without confirming genotype representation.
**Cause:** Assumed the design-specified schema was implemented as designed.
**Rule:** The CARRIES edge status determines which GWAS capabilities are possible.
Single-locus GWAS via fast-path (allele counts only) works without CARRIES.
All individual-level operations (epistasis, cosegregation, disease subgraph) require CARRIES
or a confirmed equivalent. Do not implement Phase 3 features before Phase 1 confirms
the genotype layer is in place.

### LESSON 003: Rare disease ≠ just "add Firth correction"
**Pattern:** Treating rare disease as "same GWAS with small N, use Firth."
**Cause:** Classical GWAS framing applied to graph context.
**Rule:** Rare disease in GraphGWAS requires a fundamentally different approach:
kinship-aware cosegregation (exploit family structure, don't correct it away),
pathway diffusion (aggregate signal across many rare variants), and disease subgraph
extraction (structural rather than statistical signal). Implement all three before
declaring rare disease support complete.

### LESSON 004: AssociationResult nodes must have run_id
**Pattern:** Stored GWAS results without run_id, making it impossible to distinguish runs.
**Cause:** Did not include versioning metadata.
**Rule:** Every AssociationResult node MUST have run_id, timestamp, method, n_cases,
n_controls, and covariates. Without these, graph accumulates results from different
runs that cannot be compared or cleaned up.

### LESSON 005: CO_OCCURS edges must be bounded
**Pattern:** Attempted to build co-occurrence graph for entire chromosome.
**Cause:** Did not appreciate O(V²) complexity.
**Rule:** Always bound CO_OCCURS edge computation to:
(1) minimum co-carrier count ≥ 3,
(2) maximum genomic window ≤ 1 Mb,
(3) only variants above minimum case carrier count.
Full-genome epistasis screen requires pre-filtering to candidate loci.

---

## Session 2 Lessons (Phase 1 Implementation)

### LESSON 006: Neo4j cannot store MAP properties
**Pattern:** Tried to SET s.phenotypes = {key: value, ...} as a MAP property.
**Cause:** Neo4j property values can only be primitives or arrays, NOT nested maps.
**Rule:** Store each phenotype trait as a separate property: `s.pheno_case_control`,
`s.pheno_age`, etc. Use a `phenotypes_loaded: BOOLEAN` flag instead of checking
for MAP existence. GraphMana schema.md says `phenotypes: MAP` but this is aspirational,
not implementable as a native Neo4j property.

### LESSON 007: gt_packed decode is the foundation — validate it first
**Pattern:** Built entire association pipeline before validating that genotype decode matches stored allele counts.
**Cause:** Assumed decode was correct without verification.
**Rule:** Before any statistical computation, verify that unpack_genotypes() on gt_packed
produces AC/AN that exactly matches the stored ac_total/an_total on the same Variant node.
This single check validates the entire genotype access layer. Do it first.

---

## Session 3 Lessons (Phase 2 Implementation)

### LESSON 008: Firth beta matches logistic for common variants
**Pattern:** Firth and standard logistic give nearly identical beta estimates for common
variants (diff < 0.001), confirming the implementation is correct.
**Rule:** Always validate a new regression method by comparing to the existing method on
common variants where both should agree. Divergence on common variants = bug.

### LESSON 009: MPAT directed test misses mixed-direction effects within genes
**Pattern:** APOE gene shows p=0.07 with directed MPAT despite strong single-variant signal.
**Cause:** Variants within APOE have mixed effect directions; the directed sum cancels.
**Rule:** Always run BOTH directed (preserves direction, better for pathways) and undirected
(captures mixed effects, better for genes) MPAT. Report the more significant of the two
for gene-level tests.

---

## Session 4 Lessons (Phase 3-5 Implementation)

### LESSON 010: Neo4j may crash under heavy query load
**Pattern:** Neo4j went down during flow network construction on chr19.
**Cause:** Large result sets from queries returning gt_packed for many variants.
**Rule:** Always use SKIP/LIMIT batching for variant queries. Monitor Neo4j process
during heavy scans. Implement connection retry logic in production code.

### LESSON 011: Graph-native methods require careful scoping
**Pattern:** Co-occurrence graph on 100kb window produced 205K edges — manageable.
**Cause:** All-pairs carrier-set AND is O(V²) but bounded by window + min_cocarrier filter.
**Rule:** Always enforce: (1) window ≤ 1Mb, (2) min_cocarriers ≥ 3, (3) max_variants cap.
Without bounds, a chromosome-wide scan would OOM.

---

## Session 5 Lessons (Yeast 1011 Validation)

### LESSON 012: PLINK .bed/.bim/.fam loses nucleotide identity
**Pattern:** PLINK GWAS matrix encoded alleles as "1"/"2" instead of A/C/G/T.
**Cause:** PLINK binary format only stores allele indices, not actual bases.
**Rule:** Always use the original VCF/gVCF when nucleotide identity matters.
Never rely on PLINK→VCF conversion for data that will be loaded into Neo4j.
bcftools norm can set REF from FASTA but cannot recover ALT from numeric codes.

### LESSON 013: Linear regression for quantitative traits needs different sample selection
**Pattern:** Parallel GWAS hung because quantitative traits set is_case=null,
and get_phenotype_indices() returned empty arrays.
**Cause:** single_locus_scan always called get_phenotype_indices() even for linear method.
**Rule:** For method="linear", use get_all_indices() + get_phenotype_values() instead
of the case/control pathway. Fixed in assoc.py.

### LESSON 014: Population structure inflation is the #1 priority for GWAS calibration
**Pattern:** Many yeast traits showed lambda > 2.0 (caffeine=2.45, benomyl=2.43).
**Cause:** 1011 yeast strains span 34 clades with strong population structure. No PCA covariates.
**Rule:** Always check genomic inflation (lambda) before interpreting GWAS results.
Traits with lambda > 1.1 need population structure correction. GraphGWAS should use
graph Laplacian eigenvectors (not imported PCA) for graph-native correction.

### LESSON 015: Spectral h² captures more variance than SNP-h²
**Pattern:** Spectral h²=0.70 vs published SNP-h²=0.18-0.29 for YPETHANOL.
**Cause:** Laplacian eigendecomposition on rare-variant similarity graph captures
structure that common-variant GRM cannot represent.
**Rule:** This gap is GraphGWAS's strongest scientific argument. Validate across
multiple traits before claiming as a methodological advance.

### LESSON 016: Neo4j result storage is too slow for genome-wide writes
**Pattern:** Storing 1.9M AssociationResult nodes took > 30 min and was killed.
**Cause:** Each batch creates nodes + MATCH for variant links + MATCH for study links.
**Rule:** Default to TSV output (like conventional GWAS tools). Neo4j storage is
optional (--store flag). For graph-native queries on results, load TSV post-hoc.

### LESSON 017: FINEMAP needs signed correlation, not r²
**Pattern:** FINEMAP produced empty .snp files (0.001s runtime) on 200+ variant loci.
**Cause:** LD matrix was r² (squared correlation) but FINEMAP expects signed Pearson
correlation. With r² (all positive), the Bayesian model becomes degenerate.
**Rule:** Always check what LD matrix format a tool expects. FINEMAP: signed correlation.
SuSiE: dosage matrix (computes LD internally). Different tools, different formats.

### LESSON 018: Annotation prior weight must be gentle (α≈0.9)
**Pattern:** L1 with α=0.5 (50% stat + 50% annotation) HURT ranking performance.
Annotations boosted many non-causal variants equally, diluting the statistical signal.
**Cause:** eQTL annotations cover ~26% of variants in a locus. Non-causal eQTL variants
get the same annotation boost as the true causal variant.
**Rule:** Set α=0.85-0.95. Annotations should be a gentle nudge, not equal weight.
Only reduce α (more annotation weight) when annotations strongly discriminate (high
variance across variants in the locus).

### LESSON 020: LASSO/ADMM over-shrinks in fine-mapping; use ridge + Bayes factors
**Pattern:** GRSD with L1+graph-TV penalty via ADMM produced rank 151-517 for causal variants.
After switching to ridge + graph Laplacian, still rank 4-53 (better but not competitive).
**Cause:** LASSO aggressively shrinks most betas to zero, including the causal variant among
hundreds of correlated proxies. Ridge preserves more signal but over-smooths.
**Rule:** For fine-mapping, proper Bayesian variable selection (SuSiE/FINEMAP) beats
penalized regression. Graph structure should augment priors, not replace the likelihood.
HBP works because it uses graph structure for the prior only, keeping statistical evidence
as the primary signal.

### LESSON 021: Pre-cache graph structure to avoid Neo4j bottleneck
**Pattern:** HBP took 223s/locus because each call queried Neo4j for variant→gene→pathway.
After pre-caching chr22 graph structure (9.4s once), HBP runs at 0.08s/locus.
**Rule:** For benchmarks or batch operations, cache the full chromosome's graph structure
in a Python dict once, then pass it to per-locus functions. The cache for chr22 has 514K
entries and loads in ~9s.

### LESSON 019: FAME is not a fine-mapping tool
**Pattern:** Tried to compare FAME to SuSiE/FINEMAP/L1 for per-variant ranking.
**Cause:** FAME estimates epistatic variance components (σ²) partitioned by annotation
categories. It produces no per-variant PIPs.
**Rule:** Match tools by output type. Fine-mapping → PIPs/ranks (SuSiE, FINEMAP, L1).
Heritability partitioning → variance components (FAME, LDSC, S-LDSC).

### LESSON 022: Verify array order when composing sorted-output functions
**Pattern:** v0.1.5 GAFM-MX/HBP-MX produced 0/30 rank-1 on chr22 single-causal F1
sims when first benchmarked, vs 16/30 for base GAFM/HBP. Looked like a method failure;
turned out to be two compounding bugs: (a) `_build_candidates` returns the
FinemapCandidate list sorted by combined_score, not in input variants order, but the
wrappers were calling `apply_mixture_posterior(base_pips, z_in, ...)` where base_pips
was in sorted order while z_in was in variants order — element-wise multiplication
misaligned the arrays. (b) Even with correct alignment, multiplying GAFM's
LD-deconvolved PIPs by ABF(raw z) re-introduces LD spread because LD-correlated
noise variants have similar |z| as the causal.
**Cause:** Neither bug surfaced on the rice 3kRG inflated panels where v0.1.5 was
originally validated, because (a) the inflation regime has different LD structure
and (b) the per-locus PIP gain was assumed to come from the mixture prior. The
chr22 head-to-head (clean weak-signal regime) exposed both bugs.
**Rule:** When composing functions that return sorted output, never assume the array
order matches the input — explicitly map by id (`{c.variant_id: c.pip for c in out}`)
and re-index back to input order before any element-wise operation. When applying
a Bayes-factor reweighting to LD-aware PIPs, use LD-deconvolved z (unique_stats
from `_ld_deconvolve`) so LD-aware ranking is preserved during reweighting.
Always run a head-to-head benchmark at clean weak signal AND inflation regime
before claiming a v0.x.y enhancement.

### LESSON 023: Run cross-species benchmarks BEFORE celebrating a methodological gain
**Pattern:** v0.1.5 was first validated on rice 3kRG inflated grain-weight panels,
showed Tier-1 recovery 38-48% → 62-67%, was committed as a paper-ready release.
Later cross-species reruns under the corrected library showed: yeast 0% → 24.5%
CS=1, Arabidopsis 3.7% → 64.8% CS=1, dramatic and consistent gains. The single-
species result was directionally correct but arguably insufficient for a Nature
Genetics methods paper.
**Cause:** A v0.x.y release motivated by one species' headline number can hide
both bugs and regime-specific failures. Cross-species replication is the cheapest
way to falsify or confirm.
**Rule:** Before marking a paper-ready release, run the core benchmark on AT LEAST
3 species + 1 clean simulation (chr22 single-causal F1 here). If the gain holds
in 4 of 4 settings, ship; otherwise hold and diagnose.

### LESSON 024: Per-h² (or per-config) benchmark outputs need per-config filenames
**Pattern:** The chr22 head-to-head benchmark wrote `v15_chr22.json` regardless of
H2 env var. Running h²=0.05, h²=0.10, h²=0.02 in sequence overwrote each prior
result; only the last survived. The supplementary table covering all three
heritability regimes had to be reconstructed from /tmp log files.
**Cause:** Single hard-coded output filename in a script that's parameterised by
env vars. Easy to miss because the script "just works" the first time you run it.
**Rule:** Any benchmark script parameterised by env vars or CLI args must encode
the parameter into the output filename: `v15_chr22_h{int(round(H2*100)):02d}.json`,
not `v15_chr22.json`. Do this when you first add the parameter, not after losing
results to overwrite.

### LESSON 025: Make long batch jobs resumable via per-locus JSON saves
**Pattern:** The IRRI 72-lead fine-mapping job died at lead 19/72 (OOM under swap
exhaustion — the script holds chromosome graph caches and per-trait sumstats in
memory cumulatively). The script wrote no intermediate state, so a restart
re-did all 19 completed leads from scratch.
**Cause:** Default batch design — accumulate results in memory, write a single
summary TSV at the end. Works on a 5-minute job; fails on a 30-min job that
might OOM.
**Rule:** Any per-locus / per-rep batch script should write a per-locus JSON to
disk after each locus completes (`{trait}_{chrom}_{pos}.json`), and check for
that file at the start of each loop iteration to skip already-done work. The
final summary then loads all per-locus JSONs. Same pattern used in the rice
grain script and rice IRRI script after this lesson; should be the default.

### LESSON 026: Back up old results before regenerating; never overwrite a paper-cited number in place
**Pattern:** The mixture-prior bug fix (LESSON 022) required regenerating yeast,
Arabidopsis, IRRI, and rice grain fine-mapping outputs. Running the new scripts
would have overwritten the buggy results that previous paper drafts cited.
**Cause:** Most fine-mapping scripts write to a fixed output directory; running
them with corrected library code obliterates the prior numbers.
**Rule:** Before regenerating, `mv` the previous output directory to a
`*_pre_<bugfix-tag>` (or `*_pre_<commit-sha>`) sibling. Keeps the audit trail
intact in case the new numbers turn out to be wrong, and enables quick
before/after diffs if a reviewer asks. Never destroy paper-cited numbers in
place.

### LESSON 027: No internal version tags or codenames in the paper text
**Pattern:** The manuscript carried "v0.1.4" / "v0.1.5" tags and codename
identifiers ("L1", "L4", "M1", "Phase 1.B") inherited from internal refactors.
A reader proofread the rendered PDF and flagged 14 occurrences of v0.1.x and
multiple "L4" references that read like private dev notes.
**Cause:** Internal release labels are useful for pinning a benchmark JSON to a
software state, but they have no place in a paper that will outlive any single
release. Codenames like "L1" / "L4" are project-internal shorthand (Layer 1, 4
of an architecture stack) that readers reliably misread (L1 → "L1-regularised
lasso").
**Rule:** Paper-facing names must describe the method, not its position in the
implementation. Established renames in this project: L1 → GAFM, L4 → GLEM,
M1 → LPCE. Keep the historical Python prefixes (`l1_*`, `l4_*`, `m1_*`) for
backwards compatibility but document the mapping in `Code Availability`. No
"v0.x.y" anywhere in the paper. Pre-submission grep:
`grep -nE "\bL[0-9]\b|\bM[0-9]\b|v0\.[0-9]\.[0-9]|Phase\s+[0-9]" *.tex`.

### LESSON 028: Abstract must be self-contained, focused, and free of attribution detail
**Pattern:** A reader proofread the rendered abstract and flagged three
problems: (i) author-year phrase "Niu et al. 2021" appearing in the abstract;
(ii) a sentence describing SBayesRC's recovery rate (a competitor's number that
does not advance the paper's own contribution); (iii) emphasis on "augmented
with a SBayesRC-style 4-component mixture-prior posterior reweight" rather
than what the new methods *do*.
**Cause:** Abstracts written incrementally inherit detail that belongs in
Methods/Results/Sup.
**Rule:** Three audit items for any abstract before submission:
(1) No author-year phrases — use numerical citations (`\cite{}`) or drop
attribution and let Methods carry it. Pre-submission grep
`grep -E "et al\.\s*[12][0-9]{3}" abstract.tex` should return nothing.
(2) Every sentence's subject is a thing the paper contributes, not a thing
a competitor contributes. If the abstract describes a competitor's number,
ask whether it is needed for the paper's central claim; if not, drop it.
(3) Describe what new methods *implement* (e.g. "mixture-prior posterior
reweighting"), not what they "are derived from" (e.g. "SBayesRC-style").
The latter reads like a private design note.

### LESSON 029: No revision-history meta-commentary in published prose
**Pattern:** Discussion contained "We retain the Table 1 intervention result
but reinterpret it in this light" — a phrase that betrays the manuscript's
revision history. Same class as "earlier draft", "previously stated", "we have
now revised", "in this light".
**Cause:** Phrases that work in a response-to-reviewers letter leak into the
manuscript when the response is folded into the prose.
**Rule:** A scientific paper should read as a single forward-looking statement,
not as a diff against an earlier draft. The substance of the
"reinterpret it in this light" sentence — that the typed factor graph is the
substrate for incorporating non-uniform priors, and that informative-vs-permuted
placebo on real biobank leads is open — should be stated directly.
Pre-submission grep:
`grep -nE "we retain|reinterpret|earlier draft|previously stated|in this light|we have now revised" *.tex`.

### LESSON 030: A paper that depends on internal repo paths is not self-contained
**Pattern:** Methods text said "the full nine-method comparison is in
\texttt{data/rice\_3k/results/grain\_finemap/grain\_finemap\_summary.tsv}". A
reader (correctly) flagged this — readers cannot follow paths inside a private
working directory.
**Cause:** Comfortable shortcut while drafting Methods: defer the full table to
"the repo" rather than fitting all data in a Sup Table.
**Rule:** Every quantitative claim that references a comparison readers would
want to inspect must point to a Sup Table (or Sup File) carried *with* the
paper. A 41-row × 9-method credible-set comparison fits as a `sidewaystable`
(landscape orientation, scriptsize fontsize, 2.5pt tabcolsep). Pre-submission
grep:
`grep -nE "data/|results/|tests/|src/|/mnt/" abstract.tex introduction.tex results.tex discussion.tex`
should be empty (Methods/Sup may reference paths for code and external data).

### LESSON 031: First-use citations for every named tool
**Pattern:** PLINK2 was cited (`chang2015second`) and tabix was cited
(`li2011tabix`), but BGEN had no citation despite being mentioned 5+ times in
Methods. Caught only on pre-submission grep.
**Cause:** When a tool name appears in many places in Methods, it's natural to
remember adding a `\cite{}` somewhere — but "somewhere" might be the third
mention, not the first.
**Rule:** Every named external resource (database, tool, data format,
package) must have a `\cite{}` on its FIRST appearance in Methods. Maintain
a checklist of common tools and grep each one's first mention before
submission: PLINK2, BGEN, tabix, BCFtools, Hail, gctb, susieR, susieinf,
finemapinf, SBayesRC, FINEMAP, SuSiE, GENCODE, GTEx, STRING, ENCODE,
RegulomeDB, etc.

### LESSON 032: Define every acronym on first use, even the ones that "everyone knows"
**Pattern:** Discussion used "S-LDSC" with no expansion. The acronym is
common in human-statistical-genetics circles but our paper targets a broader
fine-mapping audience (crop genomicists, AI-method developers) who may not
have seen it.
**Cause:** Author calibration drifted toward the human-GWAS subfield's
shorthand vocabulary.
**Rule:** Pre-submission, grep every all-caps token of length ≥ 2 in the
manuscript and confirm it is spelled out at first use. Common offenders in
this project: S-LDSC, GTEx, eQTL, cCRE, MAF, LD, PIP, CS, GWAS, F_ST, FPR,
TDR, BF. The bar for "everyone knows it" should be a non-domain-expert
reviewer, not the author's daily reading list.

### LESSON 033: Caption-figure consistency must be enforced by regenerating both together
**Pattern:** Sup Fig S4 caption described seven sub-panels (a–g) but the
rendered `fig4_calibration_null_fpr.pdf` still showed only the original four
(a–d). The text-edit pass had updated the caption with mixture-prior panels
e–g without regenerating the figure to match.
**Cause:** Asymmetric tooling — caption is a `.tex` string the author types,
figure is a `.pdf` produced by a Python script. Updating one without the
other is a one-line edit; updating both takes a paper-figure regeneration.
**Rule:** When changing a figure caption to add panels (or change axis
labels, units, method names), the same commit must update the Python
figure-generator and regenerate the artifact. Pre-submission verification:
`pdftotext main.pdf` the rendered paper and grep for the panel labels
(`(\textbf{e})`, `(\textbf{f})`) to confirm both the caption and the figure
render them. If only the caption renders the label, the figure has not been
regenerated.

### LESSON 034: Method names appear in code, paper, CLI, docs — keep a memory mapping current
**Pattern:** When renaming L4 → GLEM in the paper, the CLI tree figure
(`tests/generate_cli_tree.py`) still emitted "GAFM / HBP / CLGF / L4
fine-mapping on a locus" because the rename only swept `.tex` files. The
regenerated Sup Fig S3 PDF leaked the old name back into the rendered
manuscript.
**Cause:** A rename touches at minimum five surfaces: source code (Python
classes, function names, CLI choice strings), CLI documentation
(`docs/manual/commands/*.md`), README, paper `.tex` files, and any
auto-generated figure that hard-codes the name as a string.
**Rule:** When renaming a paper-facing method, sweep all six surfaces in the
same commit. Maintain a memory note (`project_graphgwas_method_names.md`)
that lists the current paper-facing names, the historical Python prefixes,
and the mapping; update it together with any rename. Pre-submission
verification: rebuild the paper and `pdftotext main.pdf` to confirm the old
name does not appear anywhere in the rendered output (including inside
auto-generated figures).

### LESSON 035: Build artefacts (PDFs, generated tables) belong in .gitignore; only source is committed
**Pattern:** Multiple paper-figure PDFs and the rendered `main.pdf` were
intermittently tracked by git, leading to confusing "modified" lines in
`git status` after every figure regeneration. Some were already gitignored
(`*.pdf`); some weren't.
**Cause:** Mixed history of how figures entered the repo — some by `git add`,
some only ever in the working tree.
**Rule:** Treat the whole `paper/finemapping_v1/figures/` directory and
`paper/finemapping_v1/main.pdf`, `cover_letter.pdf` as build artefacts:
gitignore them globally, commit only the source `.tex`, `.bib`, and the
Python/R generators in `tests/`. The reproducibility command in
`docs/REPRODUCIBILITY.md` rebuilds them. Authoritative paper PDF for sharing
goes to a separate release artefact (Zenodo / GitHub release), not into
git history.

### LESSON 036: Clean up dead inline code paths after porting to library wrappers
**Pattern:** After porting `apply_mixture_posterior` from
`tests/rice3k_grain_shape_finemap.py` (inline) into the
`graphgwas.finemapping_v2` library and switching the script to call the
library wrappers, the inline definition + an aliased import
(`apply_mixture_posterior as _lib_apply_mixture_posterior`) lingered as dead
code. Future maintainers would have hit it and assumed it was authoritative.
**Cause:** Port + cleanup are two separate edits; tempting to leave the
inline code "as a smoke-test reference" but nothing actually calls it.
**Rule:** When a script switches from inline math to a library call, delete
the dead inline definition in the same commit. Defensive aliasing (`X as _lib_X`)
is a smell when nothing calls `_lib_X`. Verify with a quick grep that no other
file imports the dead symbol.

### LESSON 037: Composite figures that need parallel info should be paired in the figure, not in two separate figures
**Pattern:** Original Figure 2 showed only HBP message-passing on the
factor graph. GAFM was described in prose without a corresponding diagram, even
though both methods are core contributions of the paper.
**Cause:** The figure was added when HBP was the only graph-native method;
when GAFM was added later, no one redesigned Figure 2.
**Rule:** When the paper introduces method A as the core contribution, then
later adds method B as a "complement", revisit the early figures —
specifically, can the architecture/data-flow figure be made into a paired
two-panel figure showing A and B side-by-side? In this paper the redesign
was: Fig 2(a) HBP factor graph; Fig 2(b) GAFM data flow (z → LD-deconvolve →
α-blend with graph functional score → softmax → PIPs). Two panels in one
figure beats one figure for A and a hand-wave for B.

### LESSON 038: Auto-mode does not authorise external publish actions
**Pattern:** While auto mode was active, I committed source changes locally
and pushed only when explicitly told to. PyPI publish, GitHub release, Zenodo
refresh remained pending throughout because "publish" is an external action
that requires explicit per-action authorisation.
**Cause:** Auto mode is a permission to execute autonomously on local,
reversible work — it is not blanket authorisation for actions visible to the
outside world.
**Rule:** Even in auto mode, external publish actions (PyPI upload, package
registry pushes, Zenodo / Figshare uploads, social posts, anything that
creates a public artefact tied to the user's identity) require an explicit
prompt from the user for that specific action. `git push` to an existing
remote is borderline — push only after explicit confirmation in the same
session.

### LESSON 039: Hand-typed subsection numbers diverge silently from LaTeX auto-numbered captions
**Pattern:** A reader of the rendered PDF noticed pages 60–65 contained
subsection running headers like "Supplementary Table S5b — 3kRG grain
weight + shape …" while the body rendered the caption as
"Table S6 3kRG grain weight + shape …". The mismatch then cascaded:
subsection "S6" rendered as "Table S7", "S7" as "Table S8", and so on.
In-text cross-references like "Supplementary Table~S7" pointed to the
wrong table, and at least one back-reference inside a caption pointed to
a completely different table than the author intended.
**Cause:** Two parallel numbering systems running independently:
(1) hand-typed `\subsection*{Supplementary Table S…}` headers, where the
author had inserted "S5b" as an addendum to S5 without renumbering the
rest, and (2) LaTeX's built-in `table` counter, which advances by one
each time `\caption{}` is called and ignores the hand labels. Once
"S5b" was inserted, every table after it had its rendered caption number
one ahead of its subsection header. The error was invisible while
drafting (subsection headers always looked correct) but obvious in the
rendered PDF.
**Rule:** Never hand-write supplementary table or figure numbers in
subsection headers if the body uses `\caption{}` — pick one source of
truth. Two practical strategies:
(a) Delete the number from the subsection header (`\subsection*{IRRI
18-trait whole-genome scan summary}` with no "Table S5") and let
the auto-numbered caption be the only labelled occurrence.
(b) Keep the redundant subsection number but ensure subsections appear
in the same order as `\caption{}` calls and that no number is skipped
or duplicated; verify with a render pass.
Pre-submission verification: `pdftotext -layout main.pdf` and grep
`Supplementary (Table|Figure) S[0-9]+` against `(Table|Fig\.) S[0-9]+`
captions; every adjacent pair should share the same number.
Cross-references in the prose are a second layer of pain — they point
to whichever number was current when the author typed them, so a
mid-paper insertion ("S5b") leaves the prose pointing at stale numbers
on either side of the insertion. After any insertion, sweep all
`Supplementary Table~S\d+` and `Supplementary Figure~S\d+` references
in `*.tex` and reconcile by content (the table the author *meant*),
not by number.

### LESSON 040: Listing every supplementary item in main.tex but providing no `\includegraphics` leaves an orphan
**Pattern:** Main.tex listed `Figure S7 — 3kRG grain weight + shape pass:
4-trait Manhattan, Q–Q, and 5-method recovery scorecard` in the
Supplementary Figures section, but `supplementary.tex` had no
`\subsection*{Supplementary Figure S7 …}` and no `\includegraphics`
calling the rice_grain figure files. The PDFs (rice_grain_manhattan.pdf,
rice_grain_qq.pdf, rice_grain_recovery_scorecard.pdf) sat in
`figures/`, were generated by a working script, and were never
rendered into the paper — orphans on both sides of the listing.
**Cause:** Sup-figures planned in advance and listed in main.tex's
front-matter, then forgotten when adding the `\subsection*` /
`\begin{figure}` / `\includegraphics` block in supplementary.tex. The
inconsistency is invisible in the rendered PDF unless the reader cross-
checks the front-matter listing against actual figure pages.
**Rule:** A pre-submission audit must confirm bidirectional consistency:
every `\item[Figure S\d+]` in main.tex's Sup-Figures listing must have
a matching `\subsection*{Supplementary Figure S\d+ …}` in
supplementary.tex with at least one `\includegraphics`; conversely
every `\includegraphics` of a figure under `figures/` must be referred
to by exactly one main- or sup-figure entry. Cheap orphan check:
```
for f in figures/*.pdf; do
  bn=$(basename "$f" .pdf)
  grep -q "$bn" *.tex || echo "orphan figure file: $bn"
done
```
Symmetrically for stale listings:
```
grep "item\[Figure S" main.tex | sed 's/.*Figure S\([0-9]\+\).*/\1/' \
  | while read n; do
      grep -q "Supplementary Figure S$n " supplementary.tex \
        || echo "main.tex lists Figure S$n but supplementary.tex has no subsection"
    done
```

### LESSON 041: Every supplementary figure or table must be cited from the main text at least once
**Pattern:** A reader spotted that Sup Fig S7 (rice grain Manhattan / Q–Q /
recovery scorecard) was rendered in the supplementary, listed in
main.tex's Sup-Figures front-matter, and even cross-referenced inside
the supplementary text — but the **main text never cited it**. A reader
who reads main results never learns that a supporting figure exists.
Same pattern can quietly affect any sup-table: add a sup item, list it
up front, never link from the main results paragraph that it supports.
**Cause:** Sup figures and tables are usually drafted alongside the main
text but the cross-reference (`Supplementary Figure~SN`) is added by
hand in the body. If the body sentence that "this would naturally cite
S7" is split off into a different drafting session, the cite never goes
in. The misuse is invisible: LaTeX builds without warnings, the figure
renders correctly, the front-matter listing looks complete.
**Rule:** A pre-submission audit must confirm that every
`\subsection*{Supplementary (Figure|Table) S\d+}` in supplementary.tex
is cited by at least one `Supplementary (Figure|Table)~S\d+` reference
in `abstract.tex`, `introduction.tex`, `results.tex`, `discussion.tex`,
or `methods.tex`. Cheap check:
```
for n in $(grep -oE "Supplementary (Figure|Table) S[0-9]+" supplementary.tex \
            | grep -oE "S[0-9]+" | sort -u); do
  uses=$(grep -hcE "Supplementary (Figure|Table)~?\s*$n\b" \
            abstract.tex introduction.tex results.tex discussion.tex methods.tex \
            2>/dev/null | paste -sd+ | bc)
  [[ "$uses" == "0" ]] && echo "$n is uncited from main text"
done
```
Same audit symmetrically for sup-tables.

### LESSON 042: Headline numbers in the body must point to the supporting evidence
**Pattern:** The rice 3kRG section asserted several headline numbers
("24/72 IRRI–Ren-2023 hits", "11/72 GAFM CS=1", "20/21 Niu QTNs") with
*one* sup-table reference for the entire ~600-line subsection. A reader
sees a wall of claims with no per-claim trail; this is exactly the
configuration that erodes confidence in a reviewing session.
**Cause:** The first draft of a section often inlines numbers because
they are fresh in the author's head. The supporting figures and tables
are added later, but the in-line claims aren't retrofitted with
`(Supplementary Table~S{N})` parentheticals after each number.
**Rule:** Every numerical claim in the body that is large enough to
warrant a reader's attention ("24/72", "20/21", "57\%", "47.6\%") must
end (or be parenthesised) with a pointer to the supporting figure,
table, or supplementary note. The minimum is one sup-pointer per
paragraph; preferably one per number. Pre-submission grep:
```
grep -nE "[0-9]+/[0-9]+\s*\([0-9]+%\)|[0-9]+\.[0-9]+%" results.tex \
  | grep -v "Supplementary\|Figure~\|Table~"
```
returns numerical claims with no nearby cross-reference; review each.

### LESSON 043: The paper's own new methods must appear in at least one main-text figure or table
**Pattern:** The v0.1.5 mixture-prior contributions (GAFM-MX, HBP-MX,
ENS) were the headline methodological advance of the paper but appeared
in **no main-text figure or table** — only in the supplementary. A
reader skimming the main figures would not have known the paper's own
new methods existed. The contribution was effectively buried by its own
authors.
**Cause:** Mixture-prior validation was added later in the manuscript
revision cycle, slotted into the supplementary because that is where
benchmark tables conventionally go, and never promoted to a main figure
when its primacy in the paper's contribution became clear. Easy to miss
during normal review because the main-figure list still looked
plausible (one figure per major section).
**Rule:** Before submission, audit every claimed methodological
contribution against the main-text display items. If the abstract or
introduction promises method X (here: GAFM-MX, HBP-MX, ENS), the body
should refer at least once to a main-text figure or table that *shows*
method X numerically — not merely a sentence "see Supplementary
Table S{N} for the cross-species evaluation". A new headline figure is
acceptable; promoting an existing supplementary table to the main text
is acceptable. Burying the paper's own contribution in the
supplementary is not. Pre-submission grep:
```
abstract_methods=$(grep -oE "GAFM-MX|HBP-MX|ENS|<your-new-method>" abstract.tex)
main_text_appearance=$(grep -oE "GAFM-MX|HBP-MX|ENS" \
   results.tex methods.tex discussion.tex | grep -v "Supplementary")
echo "abstract names: $abstract_methods"
echo "main-text references that point at a main figure or table: $main_text_appearance"
```
The two should agree.

### LESSON 044: Inserting a new main figure mid-paper renumbers the remainder; sweep the front-matter Figure list
**Pattern:** Adding a new Figure 6 (mixture-prior headline) between the
existing Pan-UKB Figure 5 and Multi-omics Coverage Figure 6 caused
LaTeX to renumber the multi-omics figure to Figure 7 in body captions
and `\ref` resolutions. The hand-typed `\item[Figure 6]` /
`\item[Figure 7]` entries in main.tex's front-matter Figure list
remained at the *old* numbers, so the listing said
"Figure 6 = Multi-omics coverage; Figure 7 = Mixture-prior" while the
PDF rendered them in the opposite order.
**Cause:** Two parallel numbering systems — LaTeX auto-numbering of
`\caption{}` calls inside `\begin{figure}` (correct) versus hand-typed
ordinal labels in the front-matter description list (stale). Same root
cause as LESSON 039 (sup-table S5b drift), now applied to the main
Figure list. `\ref` cross-references inside the body resolve correctly
because they go through the LaTeX label, but the front-matter list is
not part of the label graph.
**Rule:** When inserting a new main figure mid-paper, the same commit
must reorder the `\item[Figure N]` entries in main.tex's
"Main-text Figures" description list to match the new numerical order.
Cheap verification: `pdftotext main.pdf | grep -E "^Figure [0-9]+"` and
diff against `pdftotext main.pdf | grep -E "^Fig\. [0-9]+"`. The two
sequences must agree.

### LESSON 045: Wide tables — pick portrait + abbreviated headers OR sidewaystable + tiny font
**Pattern:** Two tables truncated on the right edge in the rendered PDF:
Sup Table S9 (9-method × 41-loci credible-set comparison; 13 columns)
was already a sidewaystable but rendered the rightmost SBayesRC column
truncated; Sup Table S10 (cross-species CS=1 sharpening; 7 columns)
was a portrait table with long per-cell strings ("60/245 (24.5%)") and
ran past the right margin.
**Cause:** LaTeX silently overruns the text width when a table's
content exceeds the available column budget. There is no warning by
default; the rendered PDF clips the cell content visually but compiles
without error.
**Rule:** Two corner-case tactics for wide tables:
(a) **Sidewaystable** + `\tiny` + `\setlength{\tabcolsep}{2pt}` for
13-column tables. `\scriptsize` with 2.5pt padding is usually too wide
for ≥ 12 columns at the typical sn-jnl page size. Verify by
`pdftotext -layout` and grepping for the rightmost column header to
confirm it is not truncated.
(b) **Portrait + abbreviated headers** + `\footnotesize` +
explicit `\hspace{}`-padded column specs (`l@{\hspace{6pt}}r@{\hspace{4pt}}r…`)
for 7-column tables with long cells. Two-line super-headers
(species name on row 1, "n=N" on row 2) save horizontal space without
losing information.
Pre-submission verification: `pdftotext -layout main.pdf` for every
sup-table page and confirm the rendered last column header matches the
.tex source last column header character-for-character.

### LESSON 046: 1×N strip layouts waste page space; use √N × √N grids for repeated panels
**Pattern:** The rice grain Q–Q figure was originally laid out as four
panels in a 1×4 row at figsize 14×3.6 in. Each Q–Q occupied
~3.5 × 3.6 in of paper, with most of the paper width wasted between
panels. The reader complained that the Q–Q plots were "very short and
small".
**Cause:** Default-grade scientific Python plotting defaults to one row
per category; when the page allows, this produces a wide, short strip
that fights the natural aspect ratio of a Q–Q plot (square-ish).
**Rule:** When N panels share the same axis type and are intended to be
read against a common scale (Q–Q, Manhattan strips, calibration curves,
ranked-rank plots), arrange them as a √N × √N grid (2×2 for 4 panels;
3×3 for 9 panels). The aspect ratio per panel improves from 1×4
horizontal to roughly 1×1, and the total figure size shrinks while
each panel doubles in apparent area. For Q–Q plots specifically, a
square aspect ratio is also methodologically correct because both axes
are on the same scale.

### LESSON 047: Q–Q sub-sampling must be uniform on the visualisation axis, not on rank
**Pattern:** A reader said the four Q–Q plots in Sup Fig S7 looked
"weird, not common" — instead of the canonical "diagonal in the bulk,
deflect upward at the tail" shape, the rendered curves rose at a
near-uniform steep angle from the origin, with the bulk looking
under-populated.
**Cause:** Naive uniform-rank sub-sampling. The script had:
`idx = unique(linspace(0, n-1, 5000) ∪ arange(min(5000, n)))`. With
n ≈ 27 M, the most-significant 5000 ranks span expected −log₁₀(p) ≈
3.7 to 7.4 (the upper tail). The `linspace` over [0, n-1] covers the
bulk's −log₁₀(p) ∈ [0, 3] range with only ~1500 points, separated by
~5400-th-rank gaps. Visual effect: dense cloud on the right half of
the panel + sparse trail on the left half = a ramp-like shape that
does not look like a canonical Q–Q.
**Rule:** Sub-sample on the visualisation axis, not on the rank axis.
For Q–Q plots specifically:
```
top_k    = 2000  # all top-K most-significant
log_idx  = unique(round(logspace(log10(top_k), log10(n-1), 5000)))
idx      = unique(arange(top_k) ∪ log_idx)
```
This places points uniformly along expected −log₁₀(p) (the natural
axis for visual interpretation), giving dense bulk coverage AND
faithful tail coverage. The same principle applies to any plot where
the rank axis is logarithmically related to the natural reading axis:
sub-sample on the reading axis, not on rank. For Manhattan plots, the
analogous rule is "sub-sample uniformly on chromosome position, not
on rank-by-p-value".

### LESSON 048: Overlay genomic-control reference lines on every Q–Q so inflation vs signal is readable at a glance
**Pattern:** Even after fixing the sub-sampling (LESSON 047), the
Q–Q plots showed only the y = x null reference line. A reader
familiar with classical GWAS knows the bulk should sit along a slope
of √λ_GC under chi-square inflation alone, but a fine-mapping /
crop-genomics audience may not have that intuition cached. Without
both reference lines, "the bulk is between the lines = inflation;
the tail deflects above both = signal" is not visible.
**Cause:** Default Python plotting libraries draw the y = x reference
but not the √λ_GC slope. Authors often forget to add the inflation
reference because they internalise λ_GC as a number (1.41, 1.78, etc)
rather than as a visual slope.
**Rule:** Every Q–Q plot in a stratified-cohort or inflated-panel
GWAS paper must have **two** reference lines:
```
ax.plot([0, x_max], [0, x_max], "k--", lw=0.8, label="y = x (null)")
ax.plot([0, x_max], [0, x_max * sqrt(λ_GC)], "gray", ls=":",
        label=f"slope = √λ_GC = {sqrt(λ_GC):.2f}")
```
Combined with a clear legend, the reader can decompose the curve
into "inflation contribution" (between the two lines) and "signal
contribution" (above both lines) without further explanation.

### LESSON 049: Long figure captions go on a separate page; figures use `[p]` placement
**Pattern:** A fully-annotated Sup Fig S7 caption (4 panel-level
sub-sections × ~10 lines each) overflowed the available page space
when the caption tried to live in the same `\begin{figure}[htbp]` as
the image. Pushing the figure to `[h]` placement squeezed the caption
to unreadable density; pushing the caption to `[t]` orphaned it from
the figure.
**Cause:** A `\caption{}` block lives inside the float — LaTeX
flows it directly under the included graphic. There is no
"continued on next page" mechanism for caption text in the standard
`figure` environment.
**Rule:** When a figure caption would exceed half a page:
(a) Set the figure to full-page float: `\begin{figure}[p]` (not
    `[htbp]`). This guarantees the figure gets its own page.
(b) Inside `\caption{}`, write only a 3–5-sentence short summary.
    End with "Extended caption with full annotation on the following
    page." or similar pointer.
(c) After `\end{figure}` and `\clearpage`, write the full annotation
    as ordinary text, opened by a bold paragraph header
    `\noindent\textbf{Supplementary Figure SN (extended caption).}`
    Use further `\noindent\textbf{...}` paragraph headers for each
    panel or each interpretive sub-section. The full annotation is
    then on a clean page, not crammed under the figure.
This is the same idea Nature Genetics's "Extended Data" mechanism
uses for figures whose interpretation needs more space than a
display-item caption can hold.
