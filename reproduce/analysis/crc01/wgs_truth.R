#!/usr/bin/env Rscript
# CRC01 bulk WGS truth: copy-number states of the four tumour sites (PT, LN1, ML1, MP1) against the
# normal blood (NB) on the 100 kb bins of the single-cell data.
# Coverage per bin -> GC / mappability correction and CPM scaling (the DNA preprocessing of the package)
# -> log2 ratio to NB, median-centred, outlier-smoothed -> CBS (DNAcopy) -> |segment mean| >= 0.15 called
# loss (1) or gain (3), otherwise neutral (2).
#   Rscript reproduce/analysis/crc01/wgs_truth.R
# Input : Data/CRC/WGS/Coverage100k/<sample>.regions.bed.gz (mosdepth 100 kb mean depth),
#         Data/CRC/mintcna_input/DNA_RC.RData (gr: the single-cell bins with gc and mapp)
# Output: Result/CRC/Benchmark/DNAcopy.RData  (cna_dnacopy: bins x sites states, ref_dnacopy)
# Run from the root of the MintCNA repository with MINTCNA_DATA set (see reproduce/README.md).
suppressPackageStartupMessages({ library(data.table); library(GenomicRanges); library(DNAcopy) })
ROOT <- Sys.getenv("MINTCNA_DATA"); if (!nzchar(ROOT)) stop("set MINTCNA_DATA (see reproduce/README.md)")
source("MintCNA/preprocess.R")                                   # GC.MAP.normalization()
e <- new.env(); load(file.path(ROOT, "Data/CRC/mintcna_input/DNA_RC.RData"), e); gr <- e$gr
covdir  <- file.path(ROOT, "Data/CRC/WGS/Coverage100k")
samples <- c("PT", "LN1", "ML1", "MP1", "NB")

# ---- 1. mean depth per bin and sample -----------------------------------------------------------
read_cov <- function(s) {
  dt <- fread(file.path(covdir, paste0(s, ".regions.bed.gz")), col.names = c("chr", "start", "end", "mean"))
  dt <- dt[chr %in% paste0("chr", 1:22)]
  dt[, chr := factor(chr, levels = paste0("chr", 1:22))]
  setorder(dt, chr, start)
}
cov_list <- setNames(lapply(samples, read_cov), samples)
stopifnot(all(sapply(cov_list, function(x) identical(x[, 1:3], cov_list[[1]][, 1:3]))))
depth <- do.call(cbind, lapply(cov_list, function(x) x$mean)); colnames(depth) <- samples
rownames(depth) <- paste(cov_list[[1]]$chr, cov_list[[1]]$start, sep = ":")

# ---- 2. same bin QC as the single cells, GC/mappability correction, CPM, log2 ratio to NB -----------
idx <- which(gr$gc > 0.2 & gr$gc < 0.8 & gr$mapp > 0.9)
ref <- as.data.frame(gr)[idx, ]; colnames(ref)[1] <- "chr"
rcn <- Seurat::NormalizeData(GC.MAP.normalization(Y = depth[idx, ], ref_qc = ref)$RC_norm,
                             normalization.method = "RC", scale.factor = 1e6)
tum <- setdiff(samples, "NB")
l2r <- apply(as.matrix(rcn[, tum]), 2, function(x) log2((x + 0.01) / (rcn[, "NB"] + 0.01)))
l2r <- apply(l2r, 2, function(x) x - median(x, na.rm = TRUE))
l2r <- modSaRa2::smooth(l2r, R = 10); colnames(l2r) <- tum

# ---- 3. CBS per site, segment means mapped back to bins, three states -----------------------------
gr_bins <- GRanges(ref$chr, IRanges(ref$start, ref$end))
chrom <- as.character(seqnames(gr_bins)); maploc <- start(gr_bins) + (width(gr_bins) %/% 2)
THR <- 0.15
bin_state <- matrix(2L, nrow(l2r), ncol(l2r), dimnames = dimnames(l2r))
for (s in tum) {
  seg <- segment(CNA(genomdat = as.numeric(l2r[, s]), chrom = chrom, maploc = maploc, data.type = "logratio", sampleid = s),
                 alpha = 0.01, undo.splits = "sdundo", undo.SD = 3.0, min.width = 5, verbose = 0)
  segdf <- as.data.frame(seg$output)
  hits <- as.data.frame(findOverlaps(gr_bins, GRanges(segdf$chrom, IRanges(segdf$loc.start, segdf$loc.end))))
  segmean <- rep(0, length(gr_bins)); segmean[hits$queryHits] <- segdf$seg.mean[hits$subjectHits]
  st <- rep(2L, length(segmean)); st[segmean <= -THR] <- 1L; st[segmean >= THR] <- 3L
  bin_state[, s] <- st
}

cna_dnacopy <- bin_state; ref_dnacopy <- ref
dir.create(file.path(ROOT, "Result/CRC/Benchmark"), showWarnings = FALSE, recursive = TRUE)
save(cna_dnacopy, ref_dnacopy, file = file.path(ROOT, "Result/CRC/Benchmark/DNAcopy.RData"))
cat(sprintf("%d bins; altered fraction per site: %s\n", nrow(bin_state),
            paste(sprintf("%s %.1f%%", tum, 100 * colMeans(bin_state != 2L)), collapse = ", ")))
