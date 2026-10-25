#!/usr/bin/env Rscript
# Bin-level precision / recall / F1 of every method against the bulk WGS state calls (GATK cnv_call -/0/+),
# on the 500 kb evaluation grid (Fig. 5b). Each cell's profile is discretised with the published AtaCNA rule
# (copy ratio < 0.75 loss, > 1.25 gain), the modal state over the malignant nuclei is the method's call, and a
# bin counts as correct only when the direction agrees; a wrong direction is a false positive and a false
# negative. Precision is undefined (NA) when a method calls nothing in a sample.
#   Rscript reproduce/analysis/multiome/evaluation/f1_states.R
# Input : as pearson_vs_truth.R (pearson_helpers.R loads the per-cell profiles)
# Output: ATAC_RNA/Result/qc/f1_states.tsv
# Run from the root of the MintCNA repository with MINTCNA_DATA set (see reproduce/README.md).
suppressPackageStartupMessages(library(data.table))
ROOT <- Sys.getenv("MINTCNA_DATA"); if (!nzchar(ROOT)) stop("set MINTCNA_DATA (see reproduce/README.md)")
HELPERS <- normalizePath("reproduce/analysis/multiome/evaluation/pearson_helpers.R", mustWork = TRUE)
setwd(file.path(ROOT, "ATAC_RNA"))
source(HELPERS)
CHRS <- paste0("chr", 1:22)
LO <- log2(0.75); HI <- log2(1.25)
man <- fread("Data/meta/sample_manifest.csv"); setkey(man, care_id)
SAMPLES <- c("P59T2","P59T3","P104T1","P104T2","P106T1","P106T2","P100T1","P100T3","P103T1","P103T2",
             "P105T1","P105T3","P108T1","P108T2","P109T1","P109T2","P110T1","P110T2","P111T2","P111T3")

f1 <- rbindlist(lapply(SAMPLES, function(s) {
  tf <- file.path("Result/truth", paste0(s, "_truth_100kb.tsv.gz"))
  if (!file.exists(tf)) return(NULL)
  truth <- fread(tf)[chr %in% CHRS]
  keep  <- fread(sprintf("Data/cell_lists/%s_cells.tsv", s))[cell_type == "Malignant", barcode]
  M <- percell_matrices(s, truth, keep)
  tv <- truth$cnv_call
  gb <- paste(truth$chr, truth$start %/% 5e5)
  ok <- tv %in% c("-", "0", "+")
  rbindlist(lapply(names(M), function(m) {
    X <- M[[m]][, ok, drop = FALSE]; y <- tv[ok]
    if (!nrow(X)) return(NULL)
    n_lo <- colSums(X < LO, na.rm = TRUE); n_hi <- colSums(X > HI, na.rm = TRUE); n_ne <- colSums(is.finite(X)) - n_lo - n_hi
    pred <- fifelse(n_lo > n_ne & n_lo >= n_hi, "-", fifelse(n_hi > n_ne & n_hi > n_lo, "+", "0"))
    # 100 kb -> 500 kb grid: the state covering most of the 500 kb bin
    dd <- data.table(g = gb[ok], pred = pred, y = y)[, .(pred = names(which.max(table(pred))), y = names(which.max(table(y)))), by = g]
    tp <- sum(dd$pred != "0" & dd$pred == dd$y); fp <- sum(dd$pred != "0" & dd$pred != dd$y); fn <- sum(dd$y != "0" & dd$pred != dd$y)
    p <- if (tp + fp == 0) NA_real_ else tp / (tp + fp); r <- if (tp + fn == 0) NA_real_ else tp / (tp + fn)
    data.table(care_id = s, method = m, n_cell = nrow(X), n_bin = ncol(X), TP = tp, FP = fp, FN = fn, called = tp + fp,
               precision = p, recall = r, F1 = if (is.na(p) || is.na(r) || p + r == 0) NA_real_ else 2 * p * r / (p + r))
  }))
}))
f1[, subtype := fifelse(grepl("^Oligo", man[care_id, subtype]), "Oligo", "Astro")]
dir.create("Result/qc", showWarnings = FALSE, recursive = TRUE)
fwrite(f1, "Result/qc/f1_states.tsv", sep = "\t")

# pooled over the 20 samples (sums of TP, FP, FN), not the median of per-sample values
pool <- f1[, { p <- if (sum(TP) + sum(FP) == 0) NA_real_ else sum(TP) / (sum(TP) + sum(FP))
               r <- if (sum(TP) + sum(FN) == 0) NA_real_ else sum(TP) / (sum(TP) + sum(FN))
               .(TP = sum(TP), FP = sum(FP), FN = sum(FN), n_undef = sum(is.na(precision)), precision = p, recall = r,
                 F1 = if (is.na(p) || is.na(r) || p + r == 0) NA_real_ else 2 * p * r / (p + r)) }, by = method]
print(pool[order(-F1, na.last = TRUE)], row.names = FALSE)
