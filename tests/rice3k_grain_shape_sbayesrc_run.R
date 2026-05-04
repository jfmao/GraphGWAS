#!/usr/bin/env Rscript
# SBayesRC LD-reference construction + per-trait sbayesrc() runs for the
# rice grain weight + shape rerun, restricted to the 23 targeted blocks
# covering all 41 fine-mapped leads.
#
# Steps:
#   LDstep1 → generate per-block gctb commands (writes ld.sh + ldm.info)
#   LDstep2 → run gctb per block (LD matrix in b{N}.ldm.bin)
#   LDstep3 → eigen decomposition per block (block{N}.eigen.bin)
#   LDstep4 → merge per-block info into final LD reference
#   sbayesrc → per-trait posterior effects + per-variant PIPs
#
# Usage:
#   Rscript tests/rice3k_grain_shape_sbayesrc_run.R 2>&1 | tee /tmp/sbayesrc_rice.log

suppressMessages(library(SBayesRC))
suppressMessages(library(data.table))

DATA      <- "/mnt/data/GraphGWAS/data/rice_3k"
SBAYES    <- file.path(DATA, "sbayesrc_rice")
GENO      <- file.path(SBAYES, "rice_3k_sbayes_subset")
BLOCK_REF <- file.path(SBAYES, "rice_targeted_blocks.txt")
SUMS      <- file.path(SBAYES, "sumstats")
OUT_LD    <- file.path(SBAYES, "ldref")
OUT_FM    <- file.path(SBAYES, "results")
GCTB      <- "~/bin/gctb"

dir.create(OUT_FM, showWarnings = FALSE, recursive = TRUE)
# OUT_LD must NOT pre-exist; LDstep1 creates it itself.
unlink(OUT_LD, recursive = TRUE)

cat("\n=== Step 1: generate per-block gctb commands ===\n")
LDstep1(mafile = file.path(SUMS, "grain_TGW.ma"),
        genoPrefix = GENO,
        outDir = OUT_LD,
        blockRef = BLOCK_REF,
        tool = GCTB)

cat("\n=== Step 2: run gctb per block ===\n")
info <- fread(file.path(OUT_LD, "ldm.info"))
cat(sprintf("  %d blocks total\n", nrow(info)))
for (b in info$Block) {
    cat(sprintf("  block %d / %d\n", b, nrow(info)))
    tryCatch(
        LDstep2(outDir = OUT_LD, blockIndex = b),
        error = function(e) cat(sprintf("    ERROR: %s\n", e$message))
    )
}

cat("\n=== Step 3: eigen decomposition per block ===\n")
for (b in info$Block) {
    cat(sprintf("  block %d / %d\n", b, nrow(info)))
    tryCatch(
        LDstep3(outDir = OUT_LD, blockIndex = b),
        error = function(e) cat(sprintf("    ERROR: %s\n", e$message))
    )
}

cat("\n=== Step 4: merge per-block info ===\n")
LDstep4(outDir = OUT_LD)

for (trait in c("TGW", "GL", "GW", "RLW")) {
    cat(sprintf("\n=== sbayesrc for %s ===\n", trait))
    out_prefix <- file.path(OUT_FM, sprintf("sbayesrc_%s", trait))
    tryCatch(
        sbayesrc(mafile = file.path(SUMS, sprintf("grain_%s.ma", trait)),
                 LDdir = OUT_LD,
                 outPrefix = out_prefix),
        error = function(e) cat(sprintf("  ERROR: %s\n", e$message))
    )
}

cat("\n=== Done ===\n")
