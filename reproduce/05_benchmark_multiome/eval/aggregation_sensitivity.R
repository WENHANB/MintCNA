#!/usr/bin/env Rscript
# How much does the cell-aggregation rule change the pseudobulk Pearson?
#   Usage: Rscript aggregation_sensitivity.R [care_id ...]
#
# Averaging malignant cells into one profile is the field convention, but the
# papers do not all aggregate the same way:
#   epiAneufinder  mean of states mapped to 1/2/3
#   AtaCNV         mean over tumour cells -> log2, OR the modal state across cells
#   Copy-scAT      neither -- it binarises both sides and counts an arm as altered
#                  if ANY subclone carries it, which is an OR, not an average
#
# An OR is systematically more permissive than a mean, so the choice is not
# cosmetic. This script measures how much it actually moves our numbers, so the
# Methods section can state the rule and show it is not doing the work.
#
# Three rules, all applied in each method's native units before the single log2:
#   mean    arithmetic mean over malignant cells      (what we report)
#   median  robust to a minority of outlier cells
#   trimmed 10% trimmed mean, between the two

suppressPackageStartupMessages({library(data.table)})
ROOT <- "<MINTCNA_ROOT>/ATAC_RNA"
setwd(ROOT)
CHRS <- paste0("chr", 1:22)
source("Code/eval/pearson_helpers.R")

args <- commandArgs(trailingOnly = TRUE)
SAMPLES <- if (length(args)) args else c("P103T2", "P105T1", "P59T2", "P110T1", "P104T1")

agg <- list(
  mean    = function(M) colMeans(M, na.rm = TRUE),
  median  = function(M) apply(M, 2, median, na.rm = TRUE),
  trimmed = function(M) apply(M, 2, mean, trim = 0.1, na.rm = TRUE))

man <- fread("Data/meta/sample_manifest.csv"); setkey(man, care_id)

res <- rbindlist(lapply(SAMPLES, function(s) {
  truth <- fread(file.path("Result/truth", paste0(s, "_truth_100kb.tsv.gz")))[chr %in% CHRS]
  cells <- fread(file.path("Data/cell_lists", paste0(s, "_cells.tsv")))
  keep  <- cells[cell_type == "Malignant", barcode]
  # percell_matrices returns log2 already, so aggregating here is a mean OF log2
  # rather than a log2 OF the mean. For this sensitivity check that is the point:
  # if the two orders gave the same answer the question would not arise.
  M <- percell_matrices(s, truth, keep)
  rbindlist(lapply(names(M), function(nm) {
    if (is.null(M[[nm]]) || !nrow(M[[nm]])) return(NULL)
    rbindlist(lapply(names(agg), function(a) {
      v  <- agg[[a]](M[[nm]])
      ok <- is.finite(v) & is.finite(truth$log2)
      data.table(care_id = s, method = nm, rule = a,
                 r = if (sum(ok) > 100 && sd(v[ok]) > 0) cor(v[ok], truth$log2[ok]) else NA_real_)
    }))
  }))
}))

fwrite(res, "Result/qc/aggregation_sensitivity.tsv", sep = "\t")
w <- dcast(res, care_id + method ~ rule, value.var = "r")
w[, spread := round(pmax(mean, median, trimmed, na.rm = TRUE) -
                    pmin(mean, median, trimmed, na.rm = TRUE), 3)]
cat("=== Pearson under three aggregation rules (log2 space) ===\n\n")
print(w[order(method, care_id)], nrow = 40)
cat("\n=== how much the rule matters, per method ===\n")
print(res[!is.na(r), .(median_r = round(median(r), 3)), by = .(method, rule)] |>
      dcast(method ~ rule, value.var = "median_r"))
cat("\nmax spread across rules within a sample:", round(max(w$spread, na.rm = TRUE), 3), "\n")
