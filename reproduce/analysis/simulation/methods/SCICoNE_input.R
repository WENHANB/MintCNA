#!/usr/bin/env Rscript
# SCICoNE input for one simulation setting: cells x bins read counts, comma separated, no header.
#   Rscript reproduce/analysis/simulation/methods/SCICoNE_input.R <setting>    (called by SCICoNE_run.sh)
# Input : <setting>/DNA_RC.RData (Y)
# Output: <setting>/Calling/SCICoNE/input/Y_cells_by_bins.csv
# Run from the root of the MintCNA repository with MINTCNA_DATA set (see reproduce/README.md).
ROOT <- Sys.getenv("MINTCNA_DATA"); if (!nzchar(ROOT)) stop("set MINTCNA_DATA (see reproduce/README.md)")
args <- commandArgs(trailingOnly = TRUE)
BASE <- if (length(args) >= 1) args[1] else "c3_N100_100kb_l1e4"
CNA_DIR <- file.path(ROOT, "Simulation/cnasim", BASE)
load(file.path(CNA_DIR, "DNA_RC.RData"))

out_dir <- file.path(CNA_DIR, "Calling/SCICoNE/input"); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.table(t(Y), file = file.path(out_dir, "Y_cells_by_bins.csv"), sep = ",",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
cat(sprintf("n_cells=%d n_bins=%d -> %s\n", ncol(Y), nrow(Y), file.path(out_dir, "Y_cells_by_bins.csv")))
