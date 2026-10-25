#!/usr/bin/env Rscript
# SCICoNE output -> bins x cells copy-number matrix (inferred copy number + 1 = state code).
#   Rscript reproduce/analysis/simulation/methods/SCICoNE_output.R <setting>   (called by SCICoNE_run.sh)
# Input : <setting>/Calling/SCICoNE/run1_inferred_cnvs.csv (cells x bins)
# Output: <setting>/Calling/SCICoNE.RData  (cna_scicone)
# Run from the root of the MintCNA repository with MINTCNA_DATA set (see reproduce/README.md).
suppressPackageStartupMessages(library(data.table))
ROOT <- Sys.getenv("MINTCNA_DATA"); if (!nzchar(ROOT)) stop("set MINTCNA_DATA (see reproduce/README.md)")
args <- commandArgs(trailingOnly = TRUE)
BASE <- if (length(args) >= 1) args[1] else "c3_N100_100kb_l1e4"
CNA_DIR <- file.path(ROOT, "Simulation/cnasim", BASE)

cna_scicone <- t(as.matrix(fread(file.path(CNA_DIR, "Calling/SCICoNE/run1_inferred_cnvs.csv"), header = FALSE))) + 1
colnames(cna_scicone) <- paste0("cell", 1:ncol(cna_scicone))
save(cna_scicone, file = file.path(CNA_DIR, "Calling/SCICoNE.RData"))
print(table(cna_scicone))
