load("<MINTCNA_ROOT>/Data/CRC/DNAmeth/dnameth_chr_all_cells_100k_gc_mapp.RData")

setwd("<MINTCNA_ROOT>/")
library(data.table)

# ==== Paths ====
covdir  <- "<MINTCNA_ROOT>/Data/CRC/WGS/Coverage100k"
samples <- c("PT","LN1","ML1","MP1","NB")   # NB is the reference

# ==== 1. Helper: read and enforce proper ordering ====
read_cov <- function(s) {
  dt <- fread(file.path(covdir, paste0(s, ".regions.bed.gz")),
              col.names = c("chr","start","end","mean"))
  # enforce canonical order chr1–chr22
  chr_levels <- paste0("chr", 1:22)
  dt <- dt[chr %in% chr_levels]
  dt[, chr := factor(chr, levels = chr_levels)]
  setorder(dt, chr, start)   # ensures sorted by chr then start
  dt
}

# ==== 2. Read all samples ====
cov_list <- setNames(lapply(samples, read_cov), samples)

# sanity: identical bins
stopifnot(all(sapply(cov_list, \(x) identical(x[,1:3], cov_list[[1]][,1:3]))))

# ==== 3. Build matrices ====
bins <- cov_list[[1]][, .(chr, start, end)]
depth_mat <- do.call(cbind, lapply(cov_list, \(x) x$mean))
colnames(depth_mat) <- samples
rownames(depth_mat) <- paste(bins$chr, bins$start, sep = ":")

gr$gc <- gr$gc/100

idx <- which(gr$gc > 0.2 & gr$gc < 0.8 & gr$mapp > 0.9)

source("<MINTCNA_ROOT>/Code/MultiCNA/Normalization.functions.R")
read_count_matrix_norm_1 <- GC.MAP.normalization(Y = depth_mat[idx,],
                                                 ref_qc = as.data.frame(gr)[idx,])
### Scale to the same RC
read_count_matrix_norm <-Seurat::NormalizeData(read_count_matrix_norm_1$RC_norm, 
                                               normalization.method = "RC",scale.factor = 1000000)
nc_samples <- c("NB")
nc_median <- apply(as.matrix(read_count_matrix_norm[, nc_samples]),1,median)
log2R_matrix <- apply(as.matrix(read_count_matrix_norm[, setdiff(colnames(read_count_matrix_norm), nc_samples)]), 
                      2, function(sample_counts) {log2(((sample_counts+0.01) / (nc_median+0.01)))})
rownames(log2R_matrix) <- rownames(read_count_matrix_norm)
colnames(log2R_matrix) <- setdiff(colnames(read_count_matrix_norm), nc_samples)

## Median shift normalization
log2R_matrix <- apply(log2R_matrix, 2, function(col) {col - median(col, na.rm = TRUE)})

## smooth() the l2r
log2R_matrix_smoothed <- modSaRa2::smooth(log2R_matrix,R = 10)
colnames(log2R_matrix_smoothed) <- colnames(log2R_matrix)


suppressPackageStartupMessages({
  library(DNAcopy)
  library(GenomicRanges)
  library(data.table)
})

# ---- prep: chromosomes & maploc from ref ----
ref <- as.data.frame(gr)[idx,];colnames(ref)[1] <- "chr"
gr_bins <- GRanges(ref$chr, IRanges(ref$start, ref$end))  # 1-based for DNAcopy
chrom   <- as.character(seqnames(gr_bins))
maploc  <- start(gr_bins) + (width(gr_bins) %/% 2)

samples <- colnames(log2R_matrix_smoothed)

# ---- parameters you can tweak ----
ALPHA   <- 0.01      # segmentation significance
UNDO    <- "sdundo"  # undo small splits by SD
UNDO_SD <- 3.0
MIN_W   <- 5         # min #bins per segment
THR     <- 0.15     # |seg.mean| threshold to call DEL/DUP

# ---- outputs ----
bin_state <- matrix(2L, nrow = nrow(log2R_matrix_smoothed), ncol = ncol(log2R_matrix_smoothed),
                    dimnames = dimnames(log2R_matrix_smoothed))
seg_list  <- vector("list", length(samples)); names(seg_list) <- samples

# ---- run DNAcopy per site and map segments back to bins ----
for (s in samples) {
  y <- as.numeric(log2R_matrix_smoothed[, s])
  
  cna <- CNA(genomdat = y, chrom = chrom, maploc = maploc,
             data.type = "logratio", sampleid = s)
  # (input is already smoothed; skip smooth.CNA)
  seg <- segment(cna, alpha = ALPHA, undo.splits = UNDO,
                 undo.SD = UNDO_SD, min.width = MIN_W, verbose = 0)
  
  segdf <- as.data.frame(seg$output)[, c("chrom","loc.start","loc.end","num.mark","seg.mean")]
  seg_list[[s]] <- segdf
  
  # assign segment mean back to each bin
  seg_gr <- GRanges(segdf$chrom, IRanges(segdf$loc.start, segdf$loc.end))
  hits   <- findOverlaps(gr_bins, seg_gr)
  segmean_per_bin <- rep(NA_real_, length(gr_bins))
  if (length(hits)) {
    h <- as.data.frame(hits)
    segmean_per_bin[h$queryHits] <- segdf$seg.mean[h$subjectHits]
  }
  # default any NA to neutral (rare)
  segmean_per_bin[is.na(segmean_per_bin)] <- 0
  
  # call states per bin: 1/2/3
  st <- rep(2L, length(segmean_per_bin))
  st[segmean_per_bin <= -THR] <- 1L
  st[segmean_per_bin >=  THR] <- 3L
  bin_state[, s] <- st
}

code_map3 <- c(del = 1, nor = 2, dup = 3)
source("<MINTCNA_ROOT>/Code/MultiCNA/CNAout.R")
plot_S_heatmap(bin_state, ref, code_map = code_map3, 
               file = "Result/CRC/Segmentation/CNV/WGS.pdf",
               width = 10, height = 6)
cna_dnacopy <- bin_state; ref_dnacopy <- ref
save(cna_dnacopy, ref_dnacopy, file = "<MINTCNA_ROOT>/Result/CRC/Benchmark/DNAcopy.RData")
