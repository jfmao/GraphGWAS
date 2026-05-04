#!/usr/bin/env Rscript
# v0.1.5: rerun SBayesRC on λ_GC-deflated .ma files. Reuses the
# 23-block 3kRG LD reference built in v0.1.4 (data/rice_3k/sbayesrc_rice/ldref/).

suppressMessages(library(SBayesRC))
suppressMessages(library(data.table))

DATA   <- "/mnt/data/GraphGWAS/data/rice_3k"
SBAYES <- file.path(DATA, "sbayesrc_rice")
SUMS   <- file.path(SBAYES, "sumstats_v15")
OUT_LD <- file.path(SBAYES, "ldref")
OUT_FM <- file.path(SBAYES, "results_v15")

dir.create(OUT_FM, showWarnings = FALSE, recursive = TRUE)

for (trait in c("TGW", "GL", "GW", "RLW")) {
    cat(sprintf("\n=== sbayesrc v0.1.5 (λ_GC deflated) for %s ===\n", trait))
    out_prefix <- file.path(OUT_FM, sprintf("sbayesrc_%s", trait))
    tryCatch(
        sbayesrc(mafile = file.path(SUMS, sprintf("grain_%s.ma", trait)),
                 LDdir = OUT_LD,
                 outPrefix = out_prefix),
        error = function(e) cat(sprintf("  ERROR: %s\n", e$message))
    )
}
cat("\n=== Done ===\n")
