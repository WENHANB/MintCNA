#!/usr/bin/env Rscript
# SCCNV on one simulation setting (R re-implementation of the published algorithm, SCCNV_functions.R).
#   Rscript reproduce/analysis/simulation/methods/SCCNV.R <setting>
# Input : <setting>/DNA_RC.RData (Y, gr)
# Output: <setting>/Calling/SCCNV.RData  (cna_sccnv: bins x cells copy-number state 1-5, ref_sccnv)
# Run from the root of the MintCNA repository with MINTCNA_DATA set (see reproduce/README.md).
suppressPackageStartupMessages({ library(GenomicRanges); library(data.table) })
ROOT <- Sys.getenv("MINTCNA_DATA"); if (!nzchar(ROOT)) stop("set MINTCNA_DATA (see reproduce/README.md)")
source("reproduce/analysis/simulation/methods/SCCNV_functions.R")
args <- commandArgs(trailingOnly = TRUE)
BASE <- if (length(args) >= 1) args[1] else "c3_N100_100kb_l1e4"
CNA_DIR <- file.path(ROOT, "Simulation/cnasim", BASE)
load(file.path(CNA_DIR, "DNA_RC.RData"))

seqlevelsStyle(gr) <- "UCSC"
ranges(gr) <- IRanges(start = start(gr) + 1, end = end(gr))
gr$gc <- as.numeric(mcols(gr)$gc) * 100

# SCCNV input: col 1 chr, col 2 mappability (filter), col 4 normalisation divisor, col 6 GC (%),
# col 7 N fraction, then one column of read counts per cell
dat1 <- data.frame(chr = as.character(seqnames(gr)), mapp = gr$mapp, start = start(gr),
                   map_norm = pmax(gr$mapp, 1e-6), end = end(gr), gc = gr$gc, N = 0,
                   stringsAsFactors = FALSE)
dat1 <- cbind(dat1, as.data.frame(Y, check.names = FALSE))
keep <- !is.na(dat1$gc) & !is.na(dat1$mapp) & dat1$gc >= 20 & dat1$gc <= 80 & dat1$mapp >= 0.90
dat1 <- dat1[keep, , drop = FALSE]

post <- SCCNV(dat1, minmap = 0.90, maxN = 0.1, ploidy = 2, windowsize = 20,
              binsize = as.integer(median(width(gr))))
cn <- sapply(post, function(pp) as.integer(max.col(as.matrix(pp), ties.method = "first") - 1))   # CN 0..20
cna_sccnv <- as.matrix(cn)
colnames(cna_sccnv) <- colnames(dat1)[8:ncol(dat1)]
cna_sccnv <- cna_sccnv + 1L; cna_sccnv[cna_sccnv >= 5] <- 5L                                   # states 1-5
ref_sccnv <- as.data.frame(gr[keep]); colnames(ref_sccnv)[1] <- "chr"

dir.create(file.path(CNA_DIR, "Calling"), showWarnings = FALSE, recursive = TRUE)
save(cna_sccnv, ref_sccnv, file = file.path(CNA_DIR, "Calling/SCCNV.RData"))
