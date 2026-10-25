#!/usr/bin/env Rscript
# Binned scDNA read counts and MintCNA input files for one simulation setting.
#   Rscript reproduce/analysis/simulation/make_mintcna_input.R <setting> [sigma]
#     sigma = 0 (default), 0.05 or 0.1: multiplicative Gaussian noise on the binned counts
#             (the <setting>_n0.05 / _n0.1 settings; set.seed(1))
# Input : <setting>/readcounts.tsv, cell_types.tsv, rna_counts_gene_by_cell.tsv, genes_hg38.tsv
# Output: <setting>/DNA_RC.RData            Y (bins x cells counts) + gr (GRanges with gc, mapp);
#                                           the input of MintCNA and of every compared DNA method
#         <setting>/mintcna_input/{cell_types.tsv, rna_counts.tsv, genes.tsv}
# GC content and mappability come from CODEX2 (getgc / getmapp, BSgenome.Hsapiens.UCSC.hg38);
# bins with mappability <= 0.9 or GC outside (0.2, 0.8) are removed.
# Run from the root of the MintCNA repository with MINTCNA_DATA set (see reproduce/README.md).
suppressPackageStartupMessages({
  library(data.table); library(GenomicRanges); library(CODEX2); library(BSgenome.Hsapiens.UCSC.hg38)
})
ROOT <- Sys.getenv("MINTCNA_DATA"); if (!nzchar(ROOT)) stop("set MINTCNA_DATA (see reproduce/README.md)")
args  <- commandArgs(trailingOnly = TRUE)
BASE  <- if (length(args) >= 1) args[1] else "c3_N100_100kb_l1e4"
SIGMA <- if (length(args) >= 2) as.numeric(args[2]) else 0
CNA_DIR <- file.path(ROOT, "Simulation/cnasim", BASE)
OUT <- file.path(CNA_DIR, "mintcna_input"); dir.create(OUT, showWarnings = FALSE)

# ---- 1. CNAsim read counts (long table: cell, chrom, start, end, acount, bcount) -> bins x cells ----
rc <- fread(file.path(CNA_DIR, "readcounts.tsv"))
setnames(rc, tolower(names(rc)))
rc[, bin_id := paste0(chrom, ":", start, "-", end)]
rc[, total := acount + bcount]
RC_wide <- dcast(rc[, .(total = sum(total)), by = .(bin_id, cell)], bin_id ~ cell, value.var = "total", fill = 0)
RCmat <- as.matrix(RC_wide[, -1]); rownames(RCmat) <- RC_wide$bin_id

ref <- data.table(bin_id = rownames(RCmat))
ref[, c("chr", "pos") := tstrsplit(bin_id, ":", fixed = TRUE)]
ref[, c("start", "end") := tstrsplit(pos, "-", fixed = TRUE)]
ref[, `:=`(start = as.integer(start), end = as.integer(end), pos = NULL)]
ref[, chr := factor(chr, levels = c(paste0("chr", 1:22), "chrX", "chrY"))]
setorder(ref, chr, start, end)
RCmat <- RCmat[ref$bin_id, order(as.integer(sub("^cell", "", colnames(RCmat)))), drop = FALSE]
ref <- ref[, .(chr = as.character(chr), start, end)]

# ---- 2. optional noise -----------------------------------------------------------------------
if (SIGMA > 0) {
  set.seed(1)
  eps <- matrix(rnorm(length(RCmat), mean = 1, sd = SIGMA), nrow = nrow(RCmat))
  RCmat <- round(RCmat * eps); RCmat[RCmat < 0] <- 0
  storage.mode(RCmat) <- "integer"
}

# ---- 3. GC / mappability (CODEX2) and bin QC -------------------------------------------------
gr <- GRanges(seqnames = gsub("chr", "", ref$chr), ranges = IRanges(start = ref$start, end = ref$end))
gc   <- getgc(gr, genome = BSgenome.Hsapiens.UCSC.hg38) * 0.01
mapp <- getmapp(gr, genome = BSgenome.Hsapiens.UCSC.hg38)
keep <- which(!is.na(gc) & mapp > 0.9 & gc > 0.2 & gc < 0.8)

Y  <- RCmat[keep, ]
gr <- GRanges(data.frame(chr = as.character(seqnames(gr)), start = start(gr), end = end(gr),
                         gc = as.numeric(gc), mapp = as.numeric(mapp))[keep, ])
save(gr, Y, file = file.path(CNA_DIR, "DNA_RC.RData"))
cat(sprintf("%s: %d of %d bins kept, %d cells, sigma = %g -> DNA_RC.RData\n",
            BASE, length(keep), nrow(RCmat), ncol(Y), SIGMA))

# ---- 4. the other MintCNA inputs, in the package's format ------------------------------------
file.copy(file.path(CNA_DIR, "cell_types.tsv"), file.path(OUT, "cell_types.tsv"), overwrite = TRUE)
file.copy(file.path(CNA_DIR, "rna_counts_gene_by_cell.tsv"), file.path(OUT, "rna_counts.tsv"), overwrite = TRUE)
genes <- fread(file.path(CNA_DIR, "genes_hg38.tsv"))
setnames(genes, "chrom", "chr")
fwrite(genes, file.path(OUT, "genes.tsv"), sep = "\t")
cat("wrote", OUT, "\n")
