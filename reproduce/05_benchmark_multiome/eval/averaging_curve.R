#!/usr/bin/env Rscript
# Why does pseudobulk beat per-cell? Measure the whole curve, not just the ends.
#   Usage: Rscript averaging_curve.R [care_id]
#
# Each cell sees the SAME underlying CNA -- they are the same tumour -- but each
# carries its OWN independent measurement noise. Averaging n cells therefore
# leaves the signal untouched and shrinks the noise by sqrt(n):
#
#     mean_j = s_j + noise_j / sqrt(n)
#
# so the correlation with the truth should climb as 1/sqrt(1 + (sigma_e^2/n +
# sigma_b^2)/sigma_s^2) and level off at the systematic-bias ceiling, NOT at 1.
# Where it levels off is the interesting number: that is the part of the error
# that is shared by every cell and therefore cannot be averaged away.
#
# Pseudobulk is simply the right-hand end of this curve (n = all cells).

suppressPackageStartupMessages({library(data.table)})
ROOT <- "<MINTCNA_ROOT>/ATAC_RNA"
setwd(ROOT)
CHRS <- paste0("chr", 1:22)
source("<REPRODUCE_DIR>/lib/eval/pearson_helpers.R")

args <- commandArgs(trailingOnly = TRUE)
CARE <- if (length(args)) args[1] else "P110T1"
NS   <- c(1, 2, 5, 10, 25, 50, 100, 250, 500, 1000, 2000)
REPS <- 30
set.seed(1)

truth <- fread(file.path("Result/truth", paste0(CARE, "_truth_100kb.tsv.gz")))[chr %in% CHRS]
cells <- fread(file.path("Data/cell_lists", paste0(CARE, "_cells.tsv")))
keep  <- cells[cell_type == "Malignant", barcode]
M <- percell_matrices(CARE, truth, keep)

tv <- truth$log2
res <- rbindlist(lapply(names(M), function(nm) {
  A <- M[[nm]]; if (is.null(A) || nrow(A) < 2) return(NULL)
  ok <- is.finite(tv) & colSums(is.finite(A)) == nrow(A)
  A <- A[, ok, drop = FALSE]; t2 <- tv[ok]
  N <- nrow(A)
  rbindlist(lapply(NS[NS <= N], function(n) {
    r <- vapply(seq_len(REPS), function(k) {
      i <- sample.int(N, n)
      v <- if (n == 1) A[i, ] else colMeans(A[i, , drop = FALSE])
      if (sd(v) == 0) NA_real_ else cor(v, t2)
    }, numeric(1))
    data.table(method = nm, n = n, r = median(r, na.rm = TRUE))
  }))
}))
# the right-hand end: all cells, i.e. the pseudobulk value we report
allc <- rbindlist(lapply(names(M), function(nm) {
  A <- M[[nm]]; if (is.null(A) || nrow(A) < 2) return(NULL)
  ok <- is.finite(tv) & colSums(is.finite(A)) == nrow(A)
  v <- colMeans(A[, ok, drop = FALSE])
  data.table(method = nm, n = nrow(A), r = cor(v, tv[ok]))
}))
res <- rbind(res, allc)
fwrite(res, sprintf("Result/qc/averaging_curve_%s.tsv", CARE), sep = "\t")

w <- dcast(res, n ~ method, value.var = "r")
cat(sprintf("=== %s: Pearson vs truth as cells are averaged (median of %d draws) ===\n\n",
            CARE, REPS))
print(w, nrow = 20)
cat("\nn = 1 is the per-cell value; the last row is the pseudobulk value we report.\n")
