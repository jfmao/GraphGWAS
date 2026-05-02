# Round 2 — 2026-04-28

Manuscript: revised v1 (post-Round-1 author response). 53 pages.

Changes since Round 1 (per `response_letter_1.md` + applied edits):
- Abstract: reframed "heterogeneity vs density" to "non-uniform vs uniform priors"; species break-out (245 + 54 + 72 + 321); explicit separation of 27–2 result (eQTL-restricted simulation) from genome-wide canonical-lipid-gene recoveries (baseline GAFM cache).
- Results §intervention: section header renamed; placebo-control-needed caveat added explicitly; Table 1 caption + column headers self-contained.
- Results §rice: claim restricted to 14 / 18 well-calibrated traits.
- Discussion §Limitations: rewritten to enumerate 5 specific bounded gaps (placebo, null permutation, 30-rep CIs, hyperparameter sensitivity, ancestry-N confound) — each with the statement that the fix is queued.
- Sup Fig S6 decision tree: LPCE node replaced with a grey "(forthcoming work)" marker.

---

## Reviewer 1 — Domain Expert (Round 2 verdict)

**Reassessment:**
- R1.1 (250-kb null overlap): **Caveat acknowledged, experiment deferred.** I accept the placeholder-with-caveat language for now if the placebo + null-overlap experiments are committed to as Round-3 deliverables.
- R1.2 (placebo-prior control): **Reframing is correct and substantial.** The new abstract closing sentence is much better. The Results paragraph explicitly states the placebo control is required. **However, the paper still pitches the 0%→88% number as the punchline (Table 1, abstract).** Until the placebo arm runs, the paper cannot claim that the gain is *biological*. Section §intervention paragraph 2 should soften "causal claim" to "necessary-condition claim" (uniform → no gain; non-uniform → potential gain).
- R1.3 (27–2 vs 6-method): **Resolved.** Abstract now explicitly separates the two experiments.
- R1.4 (HBP vs GAFM PIP scales): **Resolved** in the Discussion limitations.

**Remaining major:**
- (1.2 modified) The phrase "upgrade the non-uniform-prior-vs-uniform-prior principle from a cross-species observation to a *causal* claim within the bounds of this experiment" overstates given the absence of the placebo. Tone down.

**Minor:** Still happy with the rest.

**Verdict:** **Minor Revision** (downgraded from Major).

---

## Reviewer 2 — Methods / Statistics

**Reassessment:**
- R2.1 (Wilson CIs / equivalence claim): Acknowledged via Discussion limitations. **However, Table S2 itself does not yet show CIs.** Need to either (a) add CIs to the table this round, or (b) commit to Round-3 with explicit "[CI placeholder]" entries.
- R2.2 (PIP calibration ECE): Acknowledged.
- R2.3 (hyperparameter sensitivity): Acknowledged.
- R2.4 (rice subset): **Resolved.** Restriction to 14 traits done correctly.
- R2.5 (Theorem 3 assumption): Acknowledged, deferred.
- R2.minor (Wilcoxon / cross-validation): Half resolved (Wilcoxon noted), half deferred.

**Remaining major:**
- (2.1 modified) Sup Table S2 still presents `21/30 (70%)` etc. without CIs. A round-2-acceptable interim is to add CI columns to S2 as a paper-side update without re-running benchmarks (Wilson formula is closed-form).

**Verdict:** **Minor Revision** (downgraded from Major).

---

## Reviewer 3 — Software / Systems

**Reassessment:**
- R3.1 (graph dump streaming): Author response is reasonable; deferred to deposit prep.
- R3.2 (`UserWarning` on un-validated entry points): Code change deferred; not yet visible in this round.
- R3.3 (single-command Fig 4): Acknowledged.
- R3.4 (Pan-UKB local cache): Acknowledged.
- R3.minor: Half resolved (bibtex, lint sweep), half deferred (README on GitHub).

**Remaining major:**
- (3.2) The `UserWarning` change is cheap (~30 lines) and would close one major-revision item this round. Strongly recommend implementing now rather than deferring to Round 3.

**Verdict:** **Minor Revision** (downgraded from Major).

---

## Reviewer 4 — Adversarial Skeptic

**Reassessment:**
- R4.1 (heterogeneity-vs-density rhetoric): **Fully resolved.** The new abstract closing sentence and §intervention rewrite are correct. The paper now says what the data actually support.
- R4.2 (same-code overstatement): **Resolved** via Discussion edit.
- R4.3 (single-thread baselines): Caveat added but the caveat itself is still incomplete — the paper does not state which version of SuSiE was used or whether it was the latest one. Add the SuSiE / FINEMAP version numbers to Methods §Baseline methods (already there: susieR v0.14.2, FINEMAP v1.4.2). **Resolved on inspection.**
- R4.4 (LPCE in decision tree): **Resolved.** "(forthcoming work)" marker is the right call.
- R4.5 (amortised vs first-locus runtime): Caveat acknowledged but not yet in S2 footnote.
- R4.minor (1, 2, 3): Resolved (1, 2); rebutted (title) — accepted as authorial prerogative.

**Remaining major:**
- None blocking. (4.5) is a minor improvement but not a blocker for Round-2 advancement.

**Verdict:** **Minor Revision**.

---

## Reviewer 5 — Reader-Advocate

**Reassessment:**
- R5.1 (acronyms): Acknowledged, deferred.
- R5.2 (decision tree): Author defends the supplementary placement; I accept the partial concession (Discussion paragraph names which scenario maps to which method). Watch in Round 3 whether the Discussion paragraph is concrete.
- R5.3 (Table 1 self-contained): **Resolved.**
- R5.4 (panel-count consistency): Author rebuts; accepted.
- R5.minor: All acknowledged.

**Verdict:** **Minor Revision** (held).

---

## Reviewer 6 — Editor-in-Chief (Round 2 synthesis)

**Verdict consolidation:**
- All six reviewers now at **Minor Revision** (R1, R2, R3, R4, R5 all downgraded; R5 was already at Minor in Round 1).
- The single critical blocker from Round 1 — the heterogeneity-vs-density framing — is **resolved**. The paper now says what its data support.
- Three Round-1 blockers are now textually-acknowledged-but-experimentally-deferred (placebo, null-permutation, hyperparameter sensitivity). Their Round-3 deliverable timelines are explicit in the Limitations section.

**Outstanding items for Round 3:**
1. (R1.2 mod) Soften "causal claim" → "necessary-condition claim" in Results §intervention closing paragraph.
2. (R2.1) Add Wilson 95% CIs to Sup Table S2 (closed-form, paper-side change).
3. (R3.2) Implement `UserWarning` on un-validated method entry points (CLGF, L4, LPCE) — code-side, ~30 lines.
4. Run the placebo, null-overlap, and hyperparameter-sensitivity experiments (computational; acknowledged as Round-3 critical).

**Editor verdict:** **Minor Revision.** Substantial improvement from Round 1. If the four items above are addressed in Round 3, the paper is acceptance-ready.

**Status:** Continue to Round 3.
