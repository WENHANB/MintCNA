#!/usr/bin/env Rscript
# CRC01 scTrio-seq2 (GSE97693): single-cell methylome reads -> 100 kb read counts with GC content and
# mappability, in the input format of MintCNA and of the compared DNA methods.
#   Rscript reproduce/analysis/crc01/bin_counts.R
# Input : Data/CRC/DNAmeth/raw/*<cell>.singleC.txt.gz  (per-cell CpG tables: #Chr, Pos, Total, Met, UnMet)
# Output: Data/CRC/mintcna_input/DNA_RC.RData  Y (bins x cells read counts) + gr (GRanges, hg19 autosomes,
#                                              gc = fraction, mapp = CODEX2 100-mer mappability)
#         Data/CRC/mintcna_input/cell_types.tsv  (CRC01_NC_* cells are the normal controls)
# Run from the root of the MintCNA repository with MINTCNA_DATA set (see reproduce/README.md).
suppressPackageStartupMessages({
  library(data.table); library(GenomicRanges); library(BSgenome.Hsapiens.UCSC.hg19); library(Biostrings); library(CODEX2)
})
ROOT <- Sys.getenv("MINTCNA_DATA"); if (!nzchar(ROOT)) stop("set MINTCNA_DATA (see reproduce/README.md)")
IN  <- file.path(ROOT, "Data/CRC/DNAmeth/raw")
OUT <- file.path(ROOT, "Data/CRC/mintcna_input"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
BIN <- 100e3
CHRS <- paste0("chr", 1:22)

# ---- 1. reads per 100 kb bin, one cell at a time --------------------------------------------
bin_cell <- function(file) {
  d <- fread(file, select = c("#Chr", "Pos", "Total"))
  setnames(d, c("chr", "pos", "total"))
  d <- d[chr %in% CHRS]
  d[, bin := (pos - 1) %/% BIN]
  d[, .(total = sum(total, na.rm = TRUE)), by = .(chr, bin)]
}
files <- list.files(IN, pattern = "\\.singleC\\.txt\\.gz$", full.names = TRUE)
cells <- sub(".*_(CRC01_[A-Z0-9]+_\\d+)\\.singleC\\.txt\\.gz$", "\\1", files)
counts <- rbindlist(Map(function(f, cell) bin_cell(f)[, cell := cell], files, cells))

# ---- 2. bins x cells matrix on the union of covered bins ---------------------------------------
bins <- unique(counts[, .(chr, bin)])
bins[, chr := factor(chr, levels = CHRS)]
setorder(bins, chr, bin)
bins[, `:=`(start = bin * BIN, end = (bin + 1) * BIN, chr = as.character(chr))]
Y <- matrix(0L, nrow(bins), length(cells), dimnames = list(NULL, cells))
counts[, `:=`(i = match(paste(chr, bin), paste(bins$chr, bins$bin)), j = match(cell, cells))]
Y[cbind(counts$i, counts$j)] <- as.integer(round(counts$total))

# ---- 3. bin coordinates (1-based, clipped at chromosome ends), GC content and mappability ---------
gr <- GRanges(bins$chr, IRanges(bins$start + 1, bins$end))
hg <- BSgenome.Hsapiens.UCSC.hg19
seqlengths(gr) <- seqlengths(hg)[seqlevels(gr)]
gr <- trim(gr)
keep <- width(gr) > 0
gr <- gr[keep]; Y <- Y[keep, , drop = FALSE]
gr$gc   <- rowSums(letterFrequency(getSeq(hg, gr), c("G", "C"), as.prob = TRUE))
gr$mapp <- CODEX2::getmapp(gr)

save(Y, gr, file = file.path(OUT, "DNA_RC.RData"))
fwrite(data.table(root = cells, normal = ifelse(grepl("_NC_", cells), "normal", "tumor")),
       file.path(OUT, "cell_types.tsv"), sep = "\t")
cat(sprintf("%d bins x %d cells (%d normal controls) -> %s\n", nrow(Y), ncol(Y), sum(grepl("_NC_", cells)), OUT))
