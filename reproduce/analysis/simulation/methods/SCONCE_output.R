#!/usr/bin/env Rscript
# SCONCE output -> bins x cells copy-number state matrix (states 1-5, capped at 5).
#   Rscript reproduce/analysis/simulation/methods/SCONCE_output.R <setting>    (called by SCONCE_run.sh)
# Input : <setting>/Calling/SCONCE/out/<cell>_k5.bed (chr, start, end, copy number per bin)
# Output: <setting>/Calling/SCONCE.RData  (cna_sconce, ref_sconce)
# Run from the root of the MintCNA repository with MINTCNA_DATA set (see reproduce/README.md).
suppressPackageStartupMessages({ library(GenomicRanges); library(data.table) })
ROOT <- Sys.getenv("MINTCNA_DATA"); if (!nzchar(ROOT)) stop("set MINTCNA_DATA (see reproduce/README.md)")
args <- commandArgs(trailingOnly = TRUE)
BASE <- if (length(args) >= 1) args[1] else "c3_N100_100kb_l1e4"
CNA_DIR <- file.path(ROOT, "Simulation/cnasim", BASE)
load(file.path(CNA_DIR, "DNA_RC.RData"))
K <- 5                                                 # --maxKploid used in SCONCE_run.sh

bins_df <- data.frame(chr = as.character(seqnames(gr)), start = start(gr), end = end(gr), stringsAsFactors = FALSE)
bins_key <- paste0("chr", bins_df$chr, ":", bins_df$start + 1, "-", bins_df$end)

bed_files <- list.files(file.path(CNA_DIR, "Calling/SCONCE/out"), pattern = paste0("_k", K, "\\.bed$"), full.names = TRUE)
stopifnot(length(bed_files) > 0)
cell_ids <- sub(paste0("_k", K, "\\.bed$"), "", basename(bed_files))

CN <- matrix(NA_integer_, nrow = length(bins_key), ncol = length(bed_files), dimnames = list(bins_key, cell_ids))
for (i in seq_along(bed_files)) {
  dt <- fread(bed_files[i], header = FALSE)              # chr, start, end, copy state
  idx <- match(paste0(dt[[1]], ":", dt[[2]], "-", dt[[3]]), bins_key)
  if (anyNA(idx)) stop("bins in ", bed_files[i], " not found in the reference bin set")
  CN[idx, i] <- as.integer(dt[[4]])
}

ref_sconce <- bins_df
cna_sconce <- CN + 1L; cna_sconce[cna_sconce >= 5] <- 5L
save(cna_sconce, ref_sconce, file = file.path(CNA_DIR, "Calling/SCONCE.RData"))
