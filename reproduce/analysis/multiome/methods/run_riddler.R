#!/usr/bin/env Rscript
# RIDDLER (github.com/yardimcilab/RIDDLER, commit a027120) on one multiome sample, package defaults:
# 1 Mb windows, riddler_fit(cov_min = 10000, peak_corr_filter = 0.5), copy_call(w_size = 5, prob_thresh = 0.025,
# mult = c(0, .5, 1, 1.5, 2, 3)). The package's window_matrix() is replaced by the equivalent data.table
# version in fast_window_matrix_dt.R (same output, verified on a subset).
#   Rscript reproduce/analysis/multiome/methods/run_riddler.R <care_id> [resolution, default 1MB]
# Input : ATAC_RNA/Data/atac_filtered/<care_id>_fragments.tsv.gz; the RIDDLER repository under ATAC_RNA/Tools/RIDDLER
# Output: ATAC_RNA/Result/benchmark/RIDDLER/<resolution>/<care_id>/riddler_result.rds
# Run from the root of the MintCNA repository with MINTCNA_DATA set (see reproduce/README.md).
suppressPackageStartupMessages({ library(data.table); library(robustbase); library(MASS); library(stringr); library(fitdistrplus); library(reshape) })
ROOT <- Sys.getenv("MINTCNA_DATA"); if (!nzchar(ROOT)) stop("set MINTCNA_DATA (see reproduce/README.md)")
args <- commandArgs(trailingOnly = TRUE); if (!length(args)) stop("need <care_id> [resolution]")
CARE <- args[1]; RES <- if (length(args) >= 2) args[2] else "1MB"
ROOT <- file.path(ROOT, "ATAC_RNA")
RID  <- file.path(ROOT, "Tools/RIDDLER")
OUT  <- file.path(ROOT, "Result/benchmark/RIDDLER", RES, CARE); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

source(file.path(RID, "riddler_fit.R"))
source(file.path(RID, "copy_call.R"))
source("reproduce/analysis/multiome/methods/fast_window_matrix_dt.R")
window_matrix <- window_matrix_dt

W_FILE    <- file.path(RID, "features", paste0("window_", RES, ".bed"))
FEAT_FILE <- file.path(RID, "features", paste0("features.", RES, ".csv"))
PEAK_FILE <- file.path(RID, "features", paste0("tissue_peaks.", RES, ".bed"))
FRAG      <- file.path(ROOT, "Data/atac_filtered", paste0(CARE, "_fragments.tsv.gz"))
stopifnot(file.exists(W_FILE), file.exists(FEAT_FILE), file.exists(PEAK_FILE), file.exists(FRAG))

t0 <- Sys.time()
MAT <- file.path(OUT, paste0(CARE, "_", RES, ".bed"))
window_matrix(FRAG, w_file = W_FILE, out_file = MAT, frag_min = 10000)
rid <- riddler_fit(MAT, feat_file = FEAT_FILE, peak_file = PEAK_FILE, cov_min = 10000, peak_corr_filter = 0.5)
rid <- copy_call(rid, w_size = 5, prob_thresh = 0.025, mult = c(0, .5, 1, 1.5, 2, 3))   # multipliers of the cell's baseline
saveRDS(list(care_id = CARE, resolution = RES, riddler = rid), file.path(OUT, "riddler_result.rds"))
cat(sprintf("SUMMARY\t%s\t%s\t%d windows\t%d cells\t%.1f min\n", CARE, RES, nrow(rid$CNV), ncol(rid$CNV),
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
