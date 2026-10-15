#!/usr/bin/env Rscript
# Recompute res_null in an existing scan checkpoint, in place.
#
# The scan scores do not depend on the null: the scan produces a value at every
# position, and the null only supplies the tau those values are thresholded at.
# res_null happens to be COMPUTED inside scan_cells and cached alongside the
# scores, but that is where it is stored, not what it depends on -- so a change
# to the null requires re-simulating it (~15 s) and re-running call (seconds),
# not re-scanning every cell in every chromosome (20-60 min per sample).
#
#   Usage: refresh_null.R <care_id> [suffix]     suffix: multi (default) or atac
suppressPackageStartupMessages(library(data.table))
RR <- "<MINTCNA_ROOT>/MintCNA_git"
for (f in list.files(file.path(RR, "R"), "\\.[rR]$", recursive = TRUE, full.names = TRUE))
  try(suppressWarnings(sys.source(f, globalenv())), silent = TRUE)
setwd("<MINTCNA_ROOT>/ATAC_RNA")
a <- commandArgs(TRUE); S <- a[1]; SUF <- if (length(a) > 1) a[2] else "multi"
f <- sprintf("Result/MintCNA/%s_%s/%s_scan_h15.rds", S, SUF, S)
stopifnot(file.exists(f))
sc <- readRDS(f)

# Rebuild exactly the input scan_cells used: the LAST cell, after qn_all.
# The checkpoint stores h / penalty / mode / simT / seed at the TOP level, not
# under a params list -- read them from there so a refresh cannot silently
# re-simulate at different settings than the scan was run with.
`%or%` <- function(a, b) if (is.null(a) || is.na(a)) b else a
cl   <- sc$cell_list
mode <- sc$mode %or% (if (is.matrix(cl[[1]]) && nrow(cl[[1]]) > 1) "multi" else "dna")
h    <- sc$h
simT <- sc$simT %or% 10
pen  <- sc$penalty %or% "l2"
stopifnot(!is.null(h))
old_tau <- sc$res_null$q[which.min(abs(sc$res_null$probs - 0.90))]

set.seed(sc$seed %or% 1)
ref_cell <- as_cell_matrix(cl[[length(cl)]], mode)
sc$res_null <- if (mode == "dna") {
  univarInverse_chr(as.numeric(ref_cell), sc$ref, h = h, simT = simT)
} else {
  multiInverse_chr(ref_cell, sc$ref, h = h, simT = simT, penalty = pen)
}
new_tau <- sc$res_null$q[which.min(abs(sc$res_null$probs - 0.90))]
saveRDS(sc, f)
cat(sprintf("REFRESH %-8s %-6s h=%s simT=%s  tau(0.10) %.4f -> %.4f  (%+.1f%%)\n",
            S, mode, h, simT, old_tau, new_tau, 100*(new_tau/old_tau - 1)))
