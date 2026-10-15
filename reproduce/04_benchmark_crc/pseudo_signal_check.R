## ===== Setup (assumes your working dir and files from your script) =====
suppressPackageStartupMessages({
  library(data.table); library(rtracklayer); library(ggplot2)
})

load("<MINTCNA_ROOT>/Result/CRC/Benchmark/SCOPE.RData")
load("<MINTCNA_ROOT>/Result/CRC/Benchmark/MultiCNA.RData")
load("<MINTCNA_ROOT>/Result/CRC/Benchmark/HMMcopy.RData")
load("<MINTCNA_ROOT>/Result/CRC/Benchmark/FLCNA.RData")

# Gene annotation & bins reference
hg19_gene_info <- readGFF("Data/CRC/hg19/gencode.v19.annotation.gtf.gz")
hg19_gene_info <- hg19_gene_info[hg19_gene_info$type == "gene",]

# Unify ref (bins) to data.frame with chr,start,end
if (exists("ref_gr")) {
  ref_scope <- as.data.frame(ref_gr)
  names(ref_scope)[1:3] <- c("chr","start","end")
}
if (exists("ref_flcna")) names(ref_flcna)[1:3] <- c("chr","start","end")
if (exists("ref"))       names(ref)[1:3]       <- c("chr","start","end")

# SCOPE integer caps and column alignment
if (exists("iCN")) {
  iCN <- iCN + 1L
  iCN[iCN >= 5] <- 5L
}

## ===== If needed, build RNA/DNA L2R matrices from your list =====
# Expect: filtered_multiomic_list_qn_all (one per cell), and a vector "cells"
if (!exists("rna_l2r_qn")) {
  stopifnot(exists("filtered_multiomic_list_qn_all"))
  cells <- names(filtered_multiomic_list_qn_all)
  rna_l2r_qn <- do.call(cbind, lapply(cells, \(id) as.numeric(filtered_multiomic_list_qn_all[[id]][, 2])))
  colnames(rna_l2r_qn) <- cells
}
if (!exists("dna_l2r_qn")) {
  dna_l2r_qn <- do.call(cbind, lapply(cells, \(id) as.numeric(filtered_multiomic_list_qn_all[[id]][, 3])))
  colnames(dna_l2r_qn) <- cells
}

############################################

suppressPackageStartupMessages({ library(data.table) })
out_dir <- "./Result/CRC/Segmentation/Evaluation/pseudo"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## ----- 0) Align columns across all objects -----
common_cells <- Reduce(intersect, list(colnames(S), colnames(iCN),
                                       colnames(rna_l2r_qn), colnames(dna_l2r_qn)))
S           <- S[,           common_cells, drop = FALSE]
iCN         <- iCN[,         common_cells, drop = FALSE]
rna_l2r_qn  <- rna_l2r_qn[,  common_cells, drop = FALSE]
dna_l2r_qn  <- dna_l2r_qn[,  common_cells, drop = FALSE]

## ----- helpers -----
ov_bins <- function(chr, s, e, ref_df) which(ref_df$chr==chr & ref_df$start<e & ref_df$end>s)
to_signs <- function(mat){
  z <- matrix(0L, nrow(mat), ncol(mat))
  if (max(mat, na.rm=TRUE) > 5) { z[mat < 2] <- -1L; z[mat > 2] <-  1L } else { z[mat <=2] <- -1L; z[mat >=4] <-  1L }
  z
}
# robust region label (vector named by cells; NA if no bins)
region_labels <- function(mat, ref_df, chr, s, e, alpha = 0.5){
  b <- ov_bins(chr, s, e, ref_df)
  labs <- rep(NA_character_, ncol(mat)); names(labs) <- colnames(mat)
  if (!length(b)) return(labs)
  sub <- to_signs(mat)[b,,drop=FALSE]
  fg  <- colMeans(sub== 1L); fl <- colMeans(sub==-1L); any_alt <- colSums(sub!=0)>0
  labs[fg >= alpha]        <- "dup"
  labs[fl >= alpha]        <- "del"
  labs[is.na(labs) & !any_alt] <- "nor"
  labs
}
# sets for rule MultiCNA dup & SCOPE nor
cell_sets_by_label <- function(mat, ref_df, chr, s, e, alpha = 0.5){
  labs <- region_labels(mat, ref_df, chr, s, e, alpha)
  list(dup = names(labs)[labs=="dup"], del = names(labs)[labs=="del"], nor = names(labs)[labs=="nor"])
}

## ----- 1) Plotter (two panels; cytoband ± shoulders; gene bins in red) -----
plot_cytoband_gene_two <- function(gene, cell, shoulders_bins = 5,
                                   same_ylim = TRUE, min_half_range = 0.5,
                                   x_extent = c("window","chrom")) {
  x_extent <- match.arg(x_extent)
  
  r <- subset(region_cnv, gene_name == gene)
  if (!nrow(r)) stop("gene not in region_cnv: ", gene)
  r <- r[1,]
  chr  <- as.character(r$chrom)
  cb_s <- r$cytoband_start; cb_e <- r$cytoband_end
  gs   <- r$gene_start;     ge   <- r$gene_end
  band <- as.character(r$name)
  
  # plotting indices on MultiCNA grid (ref)
  b_band <- ov_bins(chr, cb_s, cb_e, ref)
  if (!length(b_band)) return(invisible(NULL))
  idx <- seq(max(1, min(b_band)-shoulders_bins),
             min(nrow(ref), max(b_band)+shoulders_bins))
  b_gene <- ov_bins(chr, gs, ge, ref)
  is_gene <- idx %in% b_gene
  
  # x positions
  xpos <- (ref$start[idx] + ref$end[idx]) / 2
  if (x_extent == "chrom") {
    chr_idx <- which(ref$chr == chr)
    xlim <- c(min(ref$start[chr_idx]), max(ref$end[chr_idx]))
  } else {
    xlim <- range(xpos)  # window = cytoband ± shoulders
  }
  
  # labels for title (robust)
  labs_M <- region_labels(S,   ref,       chr, gs, ge, alpha = 0.5)
  labs_S <- region_labels(iCN, ref_scope, chr, gs, ge, alpha = 0.5)
  lab_M  <- if (cell %in% names(labs_M)) labs_M[[cell]] else NA_character_
  lab_S  <- if (cell %in% names(labs_S)) labs_S[[cell]] else NA_character_
  ttl <- sprintf("%s (%s):%s:%s–%s  |  MultiCNA@%s  SCOPE@%s",
                 gene, band, chr, format(cb_s, big.mark=","), format(cb_e, big.mark=","),
                 ifelse(is.na(lab_M),"NA",lab_M), ifelse(is.na(lab_S),"NA",lab_S))
  
  # y data + shared y-limits
  y_r <- as.numeric(rna_l2r_qn[idx, cell])
  y_d <- as.numeric(dna_l2r_qn[idx, cell])
  if (same_ylim) {
    h <- max(min_half_range, max(abs(range(c(y_r,y_d), finite=TRUE))))
    ylim_r <- ylim_d <- c(-h, h)
  } else {
    hr <- max(min_half_range, max(abs(range(y_r, finite=TRUE))))
    hd <- max(min_half_range, max(abs(range(y_d, finite=TRUE))))
    ylim_r <- c(-hr, hr); ylim_d <- c(-hd, hd)
  }
  
  op <- par(mfrow=c(2,1), mar=c(3.3,4,3,1), oma=c(0,0,1,0)); on.exit(par(op), add=TRUE)
  
  # RNA plot
  plot(xpos, y_r, pch=16, cex=0.9, col=ifelse(is_gene,"red","black"),
       ylim=ylim_r, xlim=xlim,
       xlab="", ylab="RNA L2R (QN)", main=paste(cell, ttl, sep=" | "))
  abline(h=0, col="red", lty=2); abline(v=c(cb_s, cb_e), col="blue", lty=2)
  
  # DNA plot
  plot(xpos, y_d, pch=16, cex=0.9, col=ifelse(is_gene,"red","black"),
       ylim=ylim_d, xlim=xlim,
       xlab="Genomic position (bp)", ylab="DNA L2R (QN)")
  abline(h=0, col="red", lty=2); abline(v=c(cb_s, cb_e), col="blue", lty=2)
}

## ----- 2) Driver: find cells (MultiCNA dup ∩ SCOPE nor) and plot only -----
plot_discordant_dup_vs_nor <- function(gene, alpha = 0.5, shoulders_bins = 5,
                                       width = 10, height = 7){
  r <- subset(region_cnv, gene_name == gene)[1,]
  chr <- as.character(r$chrom); gs <- r$gene_start; ge <- r$gene_end
  
  sets_M <- cell_sets_by_label(S,   ref,       chr, gs, ge, alpha)
  sets_S <- cell_sets_by_label(iCN, ref_scope, chr, gs, ge, alpha)
  cells  <- intersect(sets_M$dup, sets_S$nor)
  
  if (!length(cells)) { message(sprintf("[Info] %s: no MultiCNA dup & SCOPE nor cells.", gene)); return(invisible(NULL)) }
  
  pdf(file.path(out_dir, sprintf("%s_dup_vs_nor_plots.pdf", gene)), width = width, height = height)
  for (cc in cells) plot_cytoband_gene_two(gene, cc, shoulders_bins = shoulders_bins)
  dev.off()
  message(sprintf("[Saved] %s plots -> %s", gene, file.path(out_dir, sprintf("%s_dup_vs_nor_plots.pdf", gene))))
}

## ----- Run (example) -----
plot_discordant_dup_vs_nor("MYC",   alpha = 0.5, shoulders_bins = 10, width = 11, height = 10)
plot_discordant_dup_vs_nor("PLAG1", alpha = 0.5, shoulders_bins = 10, width = 11, height = 10)



################### Normal Region: MultiCNA vs HMMcopy ####################
## Assumes in session:
## ref (chr,start,end), region_cnv
## S (MultiCNA ints, bins x cells; on ref)
## cna_mat (HMMcopy ints/values, bins x cells; on ref)
## rna_l2r_qn, dna_l2r_qn (bins x cells; on ref)

suppressPackageStartupMessages({ library(data.table) })
out_dir <- "./Result/CRC/Segmentation/Evaluation/pseudo"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## ---- align cells across objects ----
common_cells <- Reduce(intersect, list(colnames(S), colnames(cna_mat),
                                       colnames(rna_l2r_qn), colnames(dna_l2r_qn)))
S          <- S[,          common_cells, drop = FALSE]
cna_mat    <- cna_mat[,    common_cells, drop = FALSE]
rna_l2r_qn <- rna_l2r_qn[, common_cells, drop = FALSE]
dna_l2r_qn <- dna_l2r_qn[, common_cells, drop = FALSE]

## ---- helpers (reuseable) ----
ov_bins <- function(chr, s, e, ref_df) which(ref_df$chr==chr & ref_df$start<e & ref_df$end>s)
to_signs <- function(mat){
  z <- matrix(0L, nrow(mat), ncol(mat))
  if (max(mat, na.rm=TRUE) > 5) { z[mat < 2] <- -1L; z[mat > 2] <-  1L } else { z[mat <=2] <- -1L; z[mat >=4] <-  1L }
  z
}
region_labels <- function(mat, ref_df, chr, s, e, alpha = 0.5){
  b <- ov_bins(chr, s, e, ref_df)
  labs <- rep(NA_character_, ncol(mat)); names(labs) <- colnames(mat)
  if (!length(b)) return(labs)
  sub <- to_signs(mat)[b,,drop=FALSE]
  fg  <- colMeans(sub== 1L); fl <- colMeans(sub==-1L); any_alt <- colSums(sub!=0)>0
  labs[fg >= alpha]            <- "dup"
  labs[fl >= alpha]            <- "del"
  labs[is.na(labs) & !any_alt] <- "nor"
  labs
}
cell_sets_by_label <- function(mat, ref_df, chr, s, e, alpha = 0.5){
  labs <- region_labels(mat, ref_df, chr, s, e, alpha)
  list(dup = names(labs)[labs=="dup"], del = names(labs)[labs=="del"], nor = names(labs)[labs=="nor"])
}

## ---- your existing plotter, unchanged name (shared x-axis) ----
plot_cytoband_gene_two <- function(gene, cell, shoulders_bins = 5,
                                   same_ylim = TRUE, min_half_range = 0.5,
                                   x_extent = c("window","chrom")) {
  x_extent <- match.arg(x_extent)
  r <- subset(region_cnv, gene_name == gene); if (!nrow(r)) stop("gene not in region_cnv: ", gene)
  r <- r[1,]; chr <- as.character(r$chrom); cb_s <- r$cytoband_start; cb_e <- r$cytoband_end
  gs <- r$gene_start; ge <- r$gene_end; band <- as.character(r$name)
  
  b_band <- ov_bins(chr, cb_s, cb_e, ref); if (!length(b_band)) return(invisible(NULL))
  idx <- seq(max(1, min(b_band)-shoulders_bins), min(nrow(ref), max(b_band)+shoulders_bins))
  b_gene <- ov_bins(chr, gs, ge, ref); is_gene <- idx %in% b_gene
  xpos <- (ref$start[idx] + ref$end[idx]) / 2
  xlim <- if (x_extent=="chrom") { chr_idx <- which(ref$chr==chr); c(min(ref$start[chr_idx]), max(ref$end[chr_idx])) } else range(xpos)
  
  # labels: MultiCNA vs HMMcopy (both on ref)
  labs_M <- region_labels(S,       ref, chr, gs, ge, alpha = 0.5)
  labs_H <- region_labels(cna_mat, ref, chr, gs, ge, alpha = 0.5)
  lab_M  <- if (cell %in% names(labs_M)) labs_M[[cell]] else NA_character_
  lab_H  <- if (cell %in% names(labs_H)) labs_H[[cell]] else NA_character_
  ttl <- sprintf("%s (%s):%s:%s–%s  |  MultiCNA@%s  HMMcopy@%s",
                 gene, band, chr, format(cb_s, big.mark=","), format(cb_e, big.mark=","),
                 ifelse(is.na(lab_M),"NA",lab_M), ifelse(is.na(lab_H),"NA",lab_H))
  
  y_r <- as.numeric(rna_l2r_qn[idx, cell]); y_d <- as.numeric(dna_l2r_qn[idx, cell])
  if (same_ylim) { h <- max(min_half_range, max(abs(range(c(y_r,y_d), finite=TRUE)))); ylim_r <- ylim_d <- c(-h,h) }
  else { hr <- max(min_half_range, max(abs(range(y_r, finite=TRUE)))); hd <- max(min_half_range, max(abs(range(y_d, finite=TRUE))))
  ylim_r <- c(-hr,hr); ylim_d <- c(-hd,hd) }
  
  op <- par(mfrow=c(2,1), mar=c(3.3,4,3,1), oma=c(0,0,1,0)); on.exit(par(op), add=TRUE)
  plot(xpos, y_r, pch=16, cex=0.85, col=ifelse(is_gene,"red","black"),
       ylim=ylim_r, xlim=xlim, xlab="", ylab="RNA L2R (QN)", main=paste(cell, ttl, sep=" | "))
  abline(h=0, col="red", lty=2); abline(v=c(cb_s, cb_e), col="blue", lty=2)
  plot(xpos, y_d, pch=16, cex=0.85, col=ifelse(is_gene,"red","black"),
       ylim=ylim_d, xlim=xlim, xlab="Genomic position (bp)", ylab="DNA L2R (QN)")
  abline(h=0, col="red", lty=2); abline(v=c(cb_s, cb_e), col="blue", lty=2)
}

## ---- driver: MultiCNA normal & HMMcopy non-normal (dup or del) ----
plot_discordant_MultiCNA_nor_vs_HMMcopy_alt <- function(gene, alpha = 0.5,
                                                        shoulders_bins = 5,
                                                        width = 11, height = 8,
                                                        x_extent = "window") {
  r <- subset(region_cnv, gene_name == gene)[1,]
  chr <- as.character(r$chrom); gs <- r$gene_start; ge <- r$gene_end
  
  sets_M <- cell_sets_by_label(S,       ref, chr, gs, ge, alpha)
  sets_H <- cell_sets_by_label(cna_mat, ref, chr, gs, ge, alpha)
  
  # MultiCNA normal = no alt; HMMcopy alt = dup ∪ del
  cells <- intersect(sets_M$nor, union(sets_H$dup, sets_H$del))
  if (!length(cells)) { message(sprintf("[Info] %s: no cells where MultiCNA=Nor & HMMcopy=Alt.", gene)); return(invisible(NULL)) }
  
  pdf(file.path(out_dir, sprintf("%s_MultiCNA_nor_vs_HMMcopy_alt.pdf", gene)), width = width, height = height)
  for (cc in cells) plot_cytoband_gene_two(gene, cc, shoulders_bins = shoulders_bins,
                                           same_ylim = TRUE, min_half_range = 0.5,
                                           x_extent = x_extent)
  dev.off()
  message(sprintf("[Saved] %s -> %s", gene,
                  file.path(out_dir, sprintf("%s_MultiCNA_nor_vs_HMMcopy_alt.pdf", gene))))
}

## === Run on examples ===
plot_discordant_MultiCNA_nor_vs_HMMcopy_alt("ALK",  alpha = 0.5, shoulders_bins = 10, x_extent = "window")
plot_discordant_MultiCNA_nor_vs_HMMcopy_alt("KLF4", alpha = 0.5, shoulders_bins = 10, x_extent = "window")

