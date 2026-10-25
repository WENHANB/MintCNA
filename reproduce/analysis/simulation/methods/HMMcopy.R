#!/usr/bin/env Rscript
# HMMcopy (1.46.0) on one simulation setting.
#   Rscript reproduce/analysis/simulation/methods/HMMcopy.R <setting>
# Input : <setting>/DNA_RC.RData (Y, gr), <setting>/cell_types.tsv
# Output: <setting>/Calling/HMMcopy.RData  (cna_hmmcopy: bins x cells HMM states, ref_hmmcopy)
# Run from the root of the MintCNA repository with MINTCNA_DATA set (see reproduce/README.md).
suppressPackageStartupMessages({ library(GenomicRanges); library(data.table); library(HMMcopy) })
ROOT <- Sys.getenv("MINTCNA_DATA"); if (!nzchar(ROOT)) stop("set MINTCNA_DATA (see reproduce/README.md)")
args <- commandArgs(trailingOnly = TRUE)
BASE <- if (length(args) >= 1) args[1] else "c3_N100_100kb_l1e4"
CNA_DIR <- file.path(ROOT, "Simulation/cnasim", BASE)
load(file.path(CNA_DIR, "DNA_RC.RData"))
ref <- gr

copy.output <- NULL
for (i in 1:ncol(Y)) {
  reads <- data.frame(chr   = rep(ref@seqnames@values, ref@seqnames@lengths),
                      start = ref@ranges@start,
                      end   = ref@ranges@start + ref@ranges@width,
                      reads = Y[, i], gc = ref$gc, map = ref$mapp)
  copy <- correctReadcount(reads)
  param <- HMMsegment(copy, getparam = TRUE)
  param$e <- 0.9999999            # long-segment parameters, as in the HMMcopy vignette
  param$strength <- 1e7
  seg <- HMMsegment(copy, param, verbose = FALSE)
  copy.output <- cbind(copy.output, seg$state)
}
colnames(copy.output) <- colnames(Y)

cna_hmmcopy <- copy.output
ref_hmmcopy <- gr
dir.create(file.path(CNA_DIR, "Calling"), showWarnings = FALSE, recursive = TRUE)
save(cna_hmmcopy, ref_hmmcopy, file = file.path(CNA_DIR, "Calling/HMMcopy.RData"))
