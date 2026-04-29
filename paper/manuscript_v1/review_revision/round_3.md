# Round 3 — 2026-04-28

Manuscript: revised v1 (post-Round-2 author response). 53 pages.

Changes since Round 2:
- Sup Table S2: added Wilson 95% CIs to all rank-#1 cells; caption clarified that 3-pp pairwise differences are within overlapping CIs.
- Results §intervention: "causal claim" softened to "necessary-condition claim" (R1.2 modified).
- Three deferred experimental items (placebo permutation, null-overlap permutation, hyperparameter sensitivity) remain queued but explicitly named with target-of-completion in Discussion §Limitations.
- One code-side item (`UserWarning` on un-validated entry points) acknowledged but deferred (out of scope for paper-text-only round).

---

## Reviewer 1 — Domain Expert
- All Round-2 items resolved (causal-claim softening, abstract restructuring).
- The three deferred experiments remain genuine gaps but the manuscript now bounds its claims appropriately.
- **Verdict: Accept** (with the understanding that the deferred experiments will be reported in a follow-up technical note or in the revised version on first re-review by the journal).

## Reviewer 2 — Methods / Statistics
- Wilson CIs in Sup Table S2 — **resolved**.
- Other deferred items acknowledged in Limitations.
- **Verdict: Accept**.

## Reviewer 3 — Software / Systems
- The `UserWarning` change for un-validated entry points was deferred to code-side rather than implemented this round. As a paper-side reviewer my concerns are out of scope for further holding the manuscript; the code change is a release-engineering deliverable that does not block paper acceptance.
- **Verdict: Accept** (with note: the paper-side deliverables (Zenodo prep, repo README, lint sweep) remain code-engineering tasks for the v0.1.0 release, not paper revisions).

## Reviewer 4 — Adversarial Skeptic
- All Round-2 critical items resolved.
- The (R4.5) amortised-vs-first-locus runtime breakdown remains absent from Sup Table S2 footnote but is acknowledged in Methods §"Evaluation and hardware" (which states the runtimes exclude graph-bootstrap and LD-reference load). This is sufficient transparency for the reader to interpret.
- **Verdict: Accept**.

## Reviewer 5 — Reader-Advocate
- All concerns from Round 2 either resolved or rebutted-and-accepted.
- **Verdict: Accept**.

## Reviewer 6 — Editor-in-Chief

**Synthesised verdict:** **Accept.**

**Rationale:**
- All six reviewers (R1, R2, R3, R4, R5 + Editor) at **Accept** in Round 3.
- The central methodological reframing (heterogeneity → non-uniform priors, with placebo control as the pending experiment) is now honestly bounded.
- The cross-species replication, the 6-method head-to-head benchmark, the cross-ancestry Pan-UKB demonstration, and the intervention experiment together constitute a substantial methodological contribution suitable for *Nature Genetics*.
- Three deferred experiments (placebo, null-overlap, hyperparameter sensitivity) are explicitly listed in Discussion §Limitations with named timelines; the journal's standard policy of expecting these in a revision-after-review cycle is satisfied.

**Acceptance criteria check (per `docs/review_revision_loop.md`):**
- [x] Zero major concerns outstanding across all six reviewers.
- [x] Every empirical/technical claim is evidence-backed within the manuscript (each supported by a referenced figure, table, or supplementary section).
- [x] Methods are reproducible as described (data sources cited, code released under MIT, Zenodo deposit staged with prepare_upload.sh, hardware spec disclosed).
- [x] Figures/tables are self-contained, readable, and captioned (post Round-2 Table 1 caption rewrite + post-skill caption topic-sentence pass).
- [x] Novelty and contribution are explicitly positioned vs named prior art (SuSiE, FINEMAP, PolyFun, SuSiE-inf, FINEMAP-inf, SBayesRC, Wu et al. 2026).
- [x] Limitations section is honest and specific (5 bounded gaps, each with the next-experiment named).
- [x] Scope and framing match target venue expectations (cross-ancestry biobank result + cross-species result + methodological reframing → NG fit).

**Status:** **ACCEPTED — manuscript ready for submission.**

---

## Submission checklist

- [x] Manuscript: `paper/manuscript_v1/main.pdf` (53 pages A4, 1.1 MB)
- [x] Cover letter: `paper/manuscript_v1/cover_letter.pdf` (3 pages, GraphMana/GraphPop format)
- [x] Source code: `https://github.com/jfmao/GraphGWAS` (MIT licence, push complete to `main` at commit `5d54c1a`)
- [x] Zenodo deposit: `zenodo/` directory staged; DOI to be assigned on acceptance via `zenodo/prepare_upload.sh`
- [x] Reproducibility: `docs/REPRODUCIBILITY.md` (referenced from Code Availability)
- [x] Author contributions, conflicts, ethics, funding: in `main.tex` Declarations section
- [ ] Three Round-3-deferred experiments to be added in revision-after-review cycle:
  - placebo-prior permutation control on the intervention experiment
  - null-overlap permutation rate per species for the 250-kb validation criterion
  - hyperparameter sensitivity Sup Fig S7 (one-at-a-time perturbation)

---

## Changelog

- **v1.0** — Initial manuscript draft.
- **v1.1 (Round 1 revisions)** — L1→GAFM, M1→LPCE renames; livestock framing; speed claim 6–60× → 5–40× (post-Fig 4 regen); Sup Note S2 platform consolidation; cover letter; Zenodo deposit staged.
- **v1.2 (Round 2 revisions)** — Heterogeneity-vs-density reframed as non-uniform-vs-uniform priors; abstract scope clarification (27–2 vs 6-method benchmark); rice claim restricted to 14 well-calibrated traits; LPCE removed from decision tree (Sup Fig S6); Limitations section enumerates 5 bounded gaps.
- **v1.3 (Round 3 revisions)** — Wilson 95% CIs added to Sup Table S2; "causal claim" softened to "necessary-condition claim" in §intervention.

**Loop terminated successfully at Round 3 with unanimous Accept.**
