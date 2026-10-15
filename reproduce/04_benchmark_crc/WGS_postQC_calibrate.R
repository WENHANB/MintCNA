#!/usr/bin/env Rscript
# Calibrate the minimum trustworthy CNA length in the bulk WGS truth.
#   Rscript WGS_postQC_calibrate.R [--nrep 10] [--seed 1]
# Run once; writes Result/CRC/Benchmark/wgs_fp_null.csv, which WGS_postQC.R reads.
#
# The truth in DNAcopy.RData is CBS output thresholded at |seg.mean| >= 0.15.
# CBS at alpha = 0.01 over 25,262 bins emits a predictable number of spurious
# short segments, and those are what show up as one- and two-bin spikes in the
# genome-wide tracks. To know where to cut, we need the length distribution of
# CBS's false positives on THIS data, not a rule of thumb.
#
# THE NULL
#   residual = smoothed log2R - fitted segment mean.
#   By construction the residual holds no CNA, but it keeps the real
#   autocorrelated technical noise (coverage waves, residual GC/mappability).
#   Re-running the identical CBS + threshold on it produces segments that are
#   false positives BY DEFINITION. A plain bin permutation would destroy the
#   autocorrelation and understate spurious segment length, so replicates come
#   from a block bootstrap (contiguous 200-bin blocks resampled within
#   chromosome) which preserves short-range structure.
#
# The calibration deliberately uses only the WGS data. Using the single-cell
# callers to decide which WGS segments are real would make the benchmark
# circular -- the truth would be defined by the methods it is meant to judge.

suppressPackageStartupMessages({
  library(data.table); library(DNAcopy); library(GenomicRanges)
})
ROOT <- "<MINTCNA_ROOT>"
OUT  <- file.path(ROOT, "Result/CRC/Benchmark")
args <- commandArgs(trailingOnly = TRUE)
NREP <- as.integer(if (is.na(match("--nrep", args))) 10 else args[match("--nrep", args) + 1])
SEED <- as.integer(if (is.na(match("--seed", args)))  1 else args[match("--seed", args) + 1])
set.seed(SEED)

# ---- CBS settings: must stay identical to WGS.R or the null is not the null ---
ALPHA <- 0.01; UNDO <- "sdundo"; UNDO_SD <- 3.0; MIN_W <- 5; THR <- 0.15
BLK   <- 200                       # block length for the bootstrap

# ------------------------------------- 1. rebuild the log2R exactly as WGS.R ---
load(file.path(ROOT, "Data/CRC/DNAmeth/dnameth_chr_all_cells_100k_gc_mapp.RData"))  # -> gr
covdir  <- file.path(ROOT, "Data/CRC/WGS/Coverage100k")
samples <- c("PT", "LN1", "ML1", "MP1", "NB")
read_cov <- function(s) {
  dt <- fread(file.path(covdir, paste0(s, ".regions.bed.gz")),
              col.names = c("chr", "start", "end", "mean"))
  cl <- paste0("chr", 1:22); dt <- dt[chr %in% cl]
  dt[, chr := factor(chr, levels = cl)]; setorder(dt, chr, start); dt
}
cov_list <- setNames(lapply(samples, read_cov), samples)
stopifnot(all(sapply(cov_list, \(x) identical(x[, 1:3], cov_list[[1]][, 1:3]))))
depth <- do.call(cbind, lapply(cov_list, \(x) x$mean)); colnames(depth) <- samples

gr$gc <- gr$gc / 100
idx <- which(gr$gc > 0.2 & gr$gc < 0.8 & gr$mapp > 0.9)
source(file.path(ROOT, "Code/MultiCNA/Normalization.functions.R"))
rcn <- Seurat::NormalizeData(
  GC.MAP.normalization(Y = depth[idx, ], ref_qc = as.data.frame(gr)[idx, ])$RC_norm,
  normalization.method = "RC", scale.factor = 1e6)
nc   <- apply(as.matrix(rcn[, "NB"]), 1, median)
tum  <- setdiff(colnames(rcn), "NB")
l2r  <- apply(as.matrix(rcn[, tum]), 2, function(x) log2((x + 0.01) / (nc + 0.01)))
rownames(l2r) <- rownames(rcn); colnames(l2r) <- tum
l2r  <- apply(l2r, 2, function(c) c - median(c, na.rm = TRUE))
l2rs <- modSaRa2::smooth(l2r, R = 10); colnames(l2rs) <- tum

ref    <- as.data.frame(gr)[idx, ]; colnames(ref)[1] <- "chr"
gb     <- GRanges(ref$chr, IRanges(ref$start, ref$end))
chrom  <- as.character(GenomicRanges::seqnames(gb))
maploc <- start(gb) + (width(gb) %/% 2)

seg_run <- function(y, sid)
  as.data.table(segment(CNA(genomdat = y, chrom = chrom, maploc = maploc,
                            data.type = "logratio", sampleid = sid),
                        alpha = ALPHA, undo.splits = UNDO, undo.SD = UNDO_SD,
                        min.width = MIN_W, verbose = 0)$output)[
    , .(chrom, loc.start, loc.end, num.mark, seg.mean)]

bin_from_seg <- function(sd) {
  h <- as.data.frame(findOverlaps(gb, GRanges(sd$chrom, IRanges(sd$loc.start, sd$loc.end))))
  m <- rep(NA_real_, length(gb)); m[h$queryHits] <- sd$seg.mean[h$subjectHits]
  m[is.na(m)] <- 0; m
}

segmean <- matrix(0, nrow(l2rs), ncol(l2rs), dimnames = dimnames(l2rs))
segs    <- list()
for (s in tum) { segs[[s]] <- seg_run(as.numeric(l2rs[, s]), s)
                 segmean[, s] <- bin_from_seg(segs[[s]]) }
state <- matrix(2L, nrow(l2rs), ncol(l2rs), dimnames = dimnames(l2rs))
state[segmean <= -THR] <- 1L; state[segmean >= THR] <- 3L

# ---- The rebuild must land on the stored truth, or it is not the pipeline that ---
# ---- produced it and the null below would describe something else. It is checked ---
# ---- as agreement, not identity: DNAcopy's default hybrid p-value uses random ---
# ---- permutation, so a handful of borderline breakpoints move with the seed. ---
# ---- Observed across seeds 1,2,3,7,42: 99.873%-100% per sample (worst case ML1, ---
# ---- 32 of 25262 bins). Bit-exactness is not required -- the residual is taken ---
# ---- from this run's own fit -- but a drop below 99.5% would mean the ---
# ---- normalisation, not the segmentation, has drifted. ---
e <- new.env(); load(file.path(OUT, "DNAcopy.RData"), envir = e)
ag <- vapply(tum, function(s) mean(state[, s] == e$cna_dnacopy[, s]), numeric(1))
cat("agreement with stored cna_dnacopy:\n"); print(round(ag, 5))
stopifnot(min(ag) >= 0.995)

# --------------------------------------------------- 2. empirical FP lengths ---
resid <- l2rs - segmean
cat("\nresidual sd:", paste(sprintf("%s=%.4f", tum, apply(resid, 2, sd)), collapse = "  "), "\n")

block_boot <- function(v, ch) unlist(lapply(split(seq_along(v), ch), function(i) {
  n <- length(i); st <- sample(seq_len(max(1, n - BLK + 1)), ceiling(n / BLK), replace = TRUE)
  head(unlist(lapply(st, function(s) v[i][s:min(s + BLK - 1, n)])), n)
}))

null <- rbindlist(lapply(tum, function(s) rbindlist(lapply(0:NREP, function(k) {
  y <- if (k == 0) as.numeric(resid[, s]) else block_boot(as.numeric(resid[, s]), chrom)
  d <- seg_run(y, sprintf("%s_%d", s, k))[abs(seg.mean) >= THR]
  if (nrow(d)) d[, `:=`(sample = s, rep = k)] else NULL
}))))

nrun <- (NREP + 1) * length(tum)
cat(sprintf("\n%d false-positive segments over %d null runs (%.1f per run)\n",
            nrow(null), nrun, nrow(null) / nrun))
print(quantile(null$num.mark, c(.5, .75, .9, .95, .99, 1)))
fwrite(null, file.path(OUT, sprintf("wgs_fp_null%s.csv", if (SEED == 1) "" else paste0("_seed", SEED))))
cat("wrote null for seed", SEED, "\n")
