#!/usr/bin/env Rscript
# RIDDLER on one sample.  Usage: Rscript run_riddler.R <care_id> [resolution]
#
# Moore TW et al., "Robust and generalizable CNV detection for single-cell
# sequencing assays", Nucleic Acids Research 2026;54(14):gkag636.
# Source: github.com/yardimcilab/RIDDLER (vendored under Tools/RIDDLER).
#
# Why it is in the benchmark. It is the only caller here that is REFERENCE-FREE:
# a robust Poisson GLM predicts each window's expected coverage from GC,
# mappability and multi-tissue peak covariates, and CNVs are the outliers. Every
# other method leans on a set of reference nuclei, so our choice to exclude
# oligodendrocytes and astrocytes from that set is an assumption a reviewer can
# question. RIDDLER answers that question without needing the assumption.
# It also consumes fragments directly and needs no allele data, unlike
# Numbat-multiome and muscadet, whose CNA calling both hard-require BAM-derived
# allele counts this dataset cannot provide.
#
# Parameters are the package's own function defaults throughout:
#   window_matrix(frag_min = 5000)
#   riddler_fit(cov_min = 10000, peak_corr_filter = 0.5, maxit = 100,
#               huber = 1.345, acc = 1e-4)
#   copy_call(knn = 0, w_size = 5, prob_thresh = 0.025,
#             mult = c(0, .5, 1, 1.5, 2, 3))
# Resolution defaults to 1 MB, the one the authors demonstrate in run_riddler.R.
# The repo also ships 100KB/250KB/500KB feature sets; pass it as argument 2.

args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("need <care_id> [resolution]")
CARE <- args[1]
RES  <- if (length(args) >= 2) args[2] else "1MB"

suppressPackageStartupMessages({
  library(data.table); library(robustbase); library(MASS)
  library(stringr); library(fitdistrplus); library(reshape)
})

ROOT <- "<MINTCNA_ROOT>/ATAC_RNA"
RID  <- file.path(ROOT, "Tools/RIDDLER")
OUT  <- file.path(ROOT, "Result/benchmark/RIDDLER", RES, CARE)
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# RIDDLER's scripts source() each other by relative path, so give them absolute ones
source(file.path(RID, "riddler_fit.R"))
source(file.path(RID, "copy_call.R"))

# [MOD] window_matrix() replaced by a data.table equivalent. The shipped version
# does fread(data.table=FALSE) then frags[frags[,4] %in% barcodes,] on ~40M rows
# of a data.frame, then rescans per chromosome; one sample did not clear the
# barcode filter in 8 hours. The replacement copies the window-assignment rule
# line for line and was verified to reproduce the original's output exactly
# (identical counts, max abs difference 0, identical barcode columns) on a
# 400k-fragment subset where the original could still finish.
source(file.path(ROOT, "Code/benchmark/fast_window_matrix_dt.R"))
window_matrix <- window_matrix_dt

W_FILE    <- file.path(RID, "features", paste0("window_", RES, ".bed"))
FEAT_FILE <- file.path(RID, "features", paste0("features.", RES, ".csv"))
PEAK_FILE <- file.path(RID, "features", paste0("tissue_peaks.", RES, ".bed"))
FRAG      <- file.path(ROOT, "Data/atac_filtered", paste0(CARE, "_fragments.tsv.gz"))
stopifnot(file.exists(W_FILE), file.exists(FEAT_FILE), file.exists(PEAK_FILE), file.exists(FRAG))

cat(sprintf("[%s] %s at %s\n", format(Sys.time()), CARE, RES))
t0 <- Sys.time()

# ---- 1. bin fragments into windows -----------------------------------------
# window_matrix() reads columns 1-4 as chr/start/end/barcode; our fragments carry
# a fifth duplicate-count column, which it ignores. The barcodes are already the
# authors' ArchR-QC'd whitelist, so frag_min acts on top of that filter.
MAT <- file.path(OUT, paste0(CARE, "_", RES, ".bed"))
window_matrix(FRAG, w_file = W_FILE, out_file = MAT, frag_min = 10000)   # README value; 5000 was ours

# ---- 2. robust Poisson GLM --------------------------------------------------
rid <- riddler_fit(MAT, feat_file = FEAT_FILE, peak_file = PEAK_FILE,
                   cov_min = 10000, peak_corr_filter = 0.5)
cat(sprintf("fit: %d windows x %d cells\n", nrow(rid$obs_norm), ncol(rid$obs_norm)))

# ---- 3. call copy number ----------------------------------------------------
# RIDDLER expresses copy number as a MULTIPLIER of the cell's own baseline, where
# 1 = normal, not as absolute copies. mult = c(0,.5,1,1.5,2,3) therefore means
# homozygous loss / one-copy loss / neutral / one-copy gain / doubling / triple.
rid <- copy_call(rid, w_size = 5, prob_thresh = 0.025, mult = c(0, .5, 1, 1.5, 2, 3))
el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))

saveRDS(list(care_id = CARE, resolution = RES, riddler = rid),
        file.path(OUT, "riddler_result.rds"))

cn <- rid$CNV
cat(sprintf("SUMMARY\t%s\t%s\t%d\t%d\t%.1f\n",
            CARE, RES, nrow(cn), ncol(cn), el))
cat("  (care_id, resolution, n_windows, n_cells, minutes)\n")
cat("CNV value distribution:\n"); print(table(as.vector(cn)))
