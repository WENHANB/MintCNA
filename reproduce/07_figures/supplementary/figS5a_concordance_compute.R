#!/usr/bin/env Rscript
# Supplementary Fig. 5a data: MintCNA event-level F1 (DEL / DUP) against the RNA
# concordance percentage (0, 25, 50, 75, 100 %) on the baseline simulation
# (c3_N100_100kb_l1e4: 0.1X, no noise).
#
# Re-implements Simulation/code/Evaulation/Concor_F1.R (2026-03-18) as a script that
# can actually be run: that file relied on an interactive session in which the
# segment/matching helpers from Simulation/code/Evaulation_all.R were already
# defined. The helpers are pulled from that file by parsing it and evaluating ONLY
# the `name <- function(...)` assignments, so none of its top-level benchmark code
# runs. Truth is built with the same idiom as RNA_benchmark/eval/s5_cn_pr_3cov.R.
#
# WHAT IS DELIBERATELY NOT REPRODUCED: Concor_F1.R contains
#     summary_all[3:4, 9] <- summary_all[3:4, 9] - 0.005
# i.e. the two F1 values at 25 % concordance (DEL and DUP) were lowered by 0.005 by
# hand before plotting. This script writes the values as computed. The published
# S5a therefore differs from the new one at the 25 % point by exactly 0.005.

.libPaths(c("<R_LIBS_USER>", .libPaths()))
suppressPackageStartupMessages({ library(data.table); library(stringr); library(dplyr); library(IRanges) })  # IRanges: findOverlaps inside cnvwise_counts_nostrata
ROOT <- "<MINTCNA_ROOT>"
CNA_DIR <- file.path(ROOT, "Simulation/cnasim/c3_N100_100kb_l1e4")
OUT_CSV <- file.path(ROOT, "figure4paper/figureS-supplementary/_work/data/figS5a_concordance_f1.csv")

## ---- helpers from Evaulation_all.R: function definitions only ----------------
ex <- parse(file.path("<REPRODUCE_DIR>/06_evaluation/Evaulation_all.R"), keep.source = FALSE)
n_fun <- 0L
for (e in ex) if (is.call(e) && identical(e[[1]], as.name("<-")) && is.call(e[[3]]) &&
                  identical(e[[3]][[1]], as.name("function"))) { eval(e, envir = globalenv()); n_fun <- n_fun + 1L }
cat("helpers loaded from Evaulation_all.R:", n_fun, "\n")
stopifnot(exists("summarize_cna_all_cells"), exists("cnvwise_counts_nostrata"), exists("coerce_123"),
          exists("align_to_truth_by_colname_fillNA"))

## ---- truth (same idiom as s5_cn_pr_3cov.R / Binsize_F1.R) ---------------------
build_CNmat_sub <- function(CNA_DIR, key_need) {
  profiles <- fread(file.path(CNA_DIR, "profiles.tsv"))
  colnames(profiles)[5] <- "cn_states"; setnames(profiles, tolower(names(profiles)))
  profiles[, cn_total := { x <- trimws(as.character(cn_states))
    ifelse(grepl("\\.", x), rowSums(matrix(as.numeric(unlist(tstrsplit(x,"\\.",perl=TRUE))), ncol=2, byrow=FALSE)),
           round(as.numeric(x))) |> as.integer() }]
  profiles[!grepl("^chr", chrom), chrom := paste0("chr", chrom)]
  profiles[, bin_id := paste0(chrom, ":", start, "-", end)]
  profiles[, cell_num := as.integer(sub("^cell","",cell))]; profiles[, chr_num := as.integer(sub("^chr","",chrom))]
  profiles <- profiles[grepl("^cell\\d+$", cell)]
  setorder(profiles, chr_num, start, end, cell_num)
  cells <- profiles[, unique(.SD[order(cell_num)]$cell), .SDcols=c("cell","cell_num")]
  CNdt <- dcast(profiles, bin_id ~ cell, value.var="cn_total")
  CNmat <- as.matrix(CNdt[, ..cells]); rownames(CNmat) <- CNdt$bin_id
  ord <- order(as.integer(sub("^chr","",sub(":.*$","",rownames(CNmat)))),
               as.integer(sub("^.*:(\\d+)-.*$","\\1",rownames(CNmat))))
  CNmat <- CNmat[ord,,drop=FALSE]
  cn_bins <- rownames(CNmat)
  key_CN <- paste(sub("^chr","",sub(":.*$","",cn_bins)),
                  as.integer(sub("-.*$","",sub("^.*:","",cn_bins))),
                  as.integer(sub("^.*-","",sub("^.*:","",cn_bins))), sep="_")
  idx <- match(key_need, key_CN); if (anyNA(idx)) stop("bins missing in CNmat")
  label <- fread(file.path(CNA_DIR,"cell_types.tsv")); nc <- label[which(label$normal=="normal"), root]
  nor_idx <- which(colnames(CNmat) %in% nc)
  CNmat_sub <- CNmat[idx,,drop=FALSE]; CNmat_sub[, -nor_idx, drop=FALSE]
}

## ---- the five calls: 0 % = DNA-only, 25/50/75 % = partial, 100 % = full --------
# FIGS5A_P_SUFFIX="_2026-10-09" switches the 25/50/75 % points to the re-called files
# written by mintcna_concordance_rerun.R (same pipeline version as the endpoints).
PSFX <- Sys.getenv("FIGS5A_P_SUFFIX", "")
pfile <- function(k) if (PSFX == "") sprintf("Calling/MintCNA_h5_p%d.Rdata", k) else sprintf("Calling/MintCNA_h5_p%d%s.RData", k, PSFX)
files <- c(`0`   = "Calling/MintCNA_h5_single.RData", `25` = pfile(1), `50` = pfile(2), `75` = pfile(3),
           `100` = "Calling/MintCNA_h5.RData")
cat("25/50/75 % files:", pfile(1), pfile(2), pfile(3), "\n")
calls <- lapply(files, function(f) { e <- new.env(); load(file.path(CNA_DIR, f), e)
  if (exists("cna_multicna_single", e)) e$cna_multicna_single else e$cna_multicna })
key_need <- rownames(calls[["100"]])
CNmat_sub <- build_CNmat_sub(CNA_DIR, key_need)
truth_mat <- coerce_123(CNmat_sub + 1)              # Concor_F1.R: cn_total 0/1 -> 1, 2 -> 2, >=3 -> 3
cat("truth:", nrow(truth_mat), "bins x", ncol(truth_mat), "cells\n")
truth_seg <- summarize_cna_all_cells(truth_mat, baseline_cn = 2, keep_neutral = FALSE)

res <- rbindlist(lapply(names(calls), function(p) {
  M <- coerce_123(calls[[p]]); rownames(M) <- rownames(truth_mat)
  M <- align_to_truth_by_colname_fillNA(M, truth_mat)
  pred_seg <- summarize_cna_all_cells(M, baseline_cn = 2, keep_neutral = FALSE)
  cnt <- cnvwise_counts_nostrata(truth_seg, pred_seg, overlap_frac = 0.20)
  cnt[, .(TP = sum(TP, na.rm = TRUE), FP = sum(FP, na.rm = TRUE), FN = sum(FN, na.rm = TRUE),
          n_cells = uniqueN(cell)), by = .(state)][, percentage := as.integer(p)]
}))
res[, precision := TP / (TP + FP)]; res[, recall := TP / (TP + FN)]
res[, F1 := 2 * precision * recall / (precision + recall)]
res[, state := factor(state, levels = c("del", "dup"))]
setorder(res, state, percentage)
res[, published_F1 := ifelse(percentage == 25L, F1 - 0.005, F1)]   # the hand edit in Concor_F1.R
print(res[, .(state, percentage, TP, FP, FN, precision = round(precision, 4), recall = round(recall, 4),
              F1 = round(F1, 4), published_F1 = round(published_F1, 4))], row.names = FALSE)
fwrite(res, OUT_CSV); cat("wrote", OUT_CSV, "\n")
