#!/usr/bin/env Rscript
# rcCAE output -> bins x cells state matrix (1 loss / 2 neutral / 3 gain relative to each cell's ploidy).
#   Rscript reproduce/analysis/simulation/methods/rcCAE_output.R <setting>     (called by rcCAE_run.sh)
# Input : <setting>/Calling/rcCAE/out/copynumber.csv (cells x bins), ploidy.csv (one value per cell)
# Output: <setting>/Calling/rcCAE.RData  (cna_rccae)
# Run from the root of the MintCNA repository with MINTCNA_DATA set (see reproduce/README.md).
suppressPackageStartupMessages(library(data.table))
ROOT <- Sys.getenv("MINTCNA_DATA"); if (!nzchar(ROOT)) stop("set MINTCNA_DATA (see reproduce/README.md)")
args <- commandArgs(trailingOnly = TRUE)
BASE <- if (length(args) >= 1) args[1] else "c3_N100_100kb_l1e4"
CNA_DIR <- file.path(ROOT, "Simulation/cnasim", BASE)

CN <- fread(file.path(CNA_DIR, "Calling/rcCAE/out", "copynumber.csv"))
P  <- fread(file.path(CNA_DIR, "Calling/rcCAE/out", "ploidy.csv"))
CN[, cell := paste0("cell", .I)]
bin_cols <- setdiff(names(CN), "cell")

base_ploidy <- as.integer(round(as.numeric(P[1, ])))
if (length(base_ploidy) != nrow(CN)) {
  warning("ploidy length differs from the number of cells; using each cell's modal copy number")
  mode_int <- function(x) { tx <- table(x); as.integer(names(tx)[which.max(tx)]) }
  base_ploidy <- CN[, mode_int(unlist(.SD)), .SDcols = bin_cols]
}

delta <- as.matrix(CN[, ..bin_cols]) - base_ploidy
code <- matrix(2L, nrow = nrow(delta), ncol = ncol(delta))
code[delta <= -1] <- 1L
code[delta >= 1]  <- 3L

cna_rccae <- t(code); colnames(cna_rccae) <- CN$cell
save(cna_rccae, file = file.path(CNA_DIR, "Calling/rcCAE.RData"))
