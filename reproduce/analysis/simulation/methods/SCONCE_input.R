#!/usr/bin/env Rscript
# SCONCE input for one simulation setting: one BED of read depths per cell (chr, start, end, depth),
# a diploid (normal cell) file list and a tumour file list.
#   Rscript reproduce/analysis/simulation/methods/SCONCE_input.R <setting>     (called by SCONCE_run.sh)
# Input : <setting>/DNA_RC.RData (Y, gr), <setting>/cell_types.tsv
# Output: <setting>/Calling/SCONCE/{in/<cell>.bed, ref/diploidFileList.txt, tumorFileList.txt}
# Run from the root of the MintCNA repository with MINTCNA_DATA set (see reproduce/README.md).
suppressPackageStartupMessages({ library(GenomicRanges); library(data.table) })
ROOT <- Sys.getenv("MINTCNA_DATA"); if (!nzchar(ROOT)) stop("set MINTCNA_DATA (see reproduce/README.md)")
args <- commandArgs(trailingOnly = TRUE)
BASE <- if (length(args) >= 1) args[1] else "c3_N100_100kb_l1e4"
CNA_DIR  <- file.path(ROOT, "Simulation/cnasim", BASE)
base_dir <- file.path(CNA_DIR, "Calling/SCONCE")
load(file.path(CNA_DIR, "DNA_RC.RData"))
label <- fread(file.path(CNA_DIR, "cell_types.tsv"))
nc_samples <- label[normal == "normal", root]

seqlevelsStyle(gr) <- "UCSC"
ranges(gr) <- IRanges(start = start(gr) + 1, end = end(gr))
gr$gc <- as.numeric(mcols(gr)$gc) * 100
keep <- !is.na(gr$gc) & !is.na(gr$mapp) & gr$gc >= 20 & gr$gc <= 80 & gr$mapp >= 0.90
Y <- Y[keep, ]; ref_gr <- gr[keep]
stopifnot(all(nc_samples %in% colnames(Y)))

for (d in c("in", "ref", "out", "logs")) dir.create(file.path(base_dir, d), recursive = TRUE, showWarnings = FALSE)
bins <- data.frame(chr = as.character(seqnames(ref_gr)), start = start(ref_gr), end = end(ref_gr))
cell_ids <- colnames(Y)
bed_files <- file.path(base_dir, "in", paste0(cell_ids, ".bed"))
for (j in seq_along(cell_ids))
  write.table(cbind(bins, readDepth = as.numeric(Y[, j])), bed_files[j],
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

is_normal <- cell_ids %in% nc_samples
writeLines(bed_files[is_normal],  file.path(base_dir, "ref", "diploidFileList.txt"))
writeLines(bed_files[!is_normal], file.path(base_dir, "tumorFileList.txt"))
cat(sprintf("%d diploid and %d tumour BED files in %s\n", sum(is_normal), sum(!is_normal), base_dir))
