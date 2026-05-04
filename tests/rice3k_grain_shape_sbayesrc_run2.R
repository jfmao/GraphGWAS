#!/usr/bin/env Rscript
# Step 3 (eigen decomp), Step 4 (merge), and sbayesrc per trait — picks up
# from where rice3k_grain_shape_sbayesrc_run.R left off (LD matrices already
# computed by gctb).

suppressMessages(library(SBayesRC))
suppressMessages(library(data.table))

DATA   <- "/mnt/data/GraphGWAS/data/rice_3k"
SBAYES <- file.path(DATA, "sbayesrc_rice")
SUMS   <- file.path(SBAYES, "sumstats")
OUT_LD <- file.path(SBAYES, "ldref")
OUT_FM <- file.path(SBAYES, "results")

dir.create(OUT_FM, showWarnings = FALSE, recursive = TRUE)

cat("\n=== Step 3: eigen decomposition per block ===\n")
info <- fread(file.path(OUT_LD, "ldm.info"))
for (b in info$Block) {
    cat(sprintf("  block %d / %d\n", b, nrow(info)))
    eigen_file <- file.path(OUT_LD, sprintf("block%d.eigen.bin", b))
    if (file.exists(eigen_file)) {
        cat(sprintf("    skipped (already done)\n"))
        next
    }
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
    # Reduced niter/burn for tractability on the 188K-variant dense rice subset
    # (default niter=3000 takes ~75 min/trait at 188K SNPs; we use niter=500/burn=200
    # which is ~5 min/trait and still gives well-converged PIPs for major effects).
    tryCatch(
        sbayesrc(mafile = file.path(SUMS, sprintf("grain_%s.ma", trait)),
                 LDdir = OUT_LD,
                 outPrefix = out_prefix,
                 niter = 500, burn = 200,
                 tuneIter = 100, tuneBurn = 60),
        error = function(e) cat(sprintf("  ERROR: %s\n", e$message))
    )
}

cat("\n=== Done ===\n")
