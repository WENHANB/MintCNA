#!/usr/bin/env Rscript
# Build the RNA-arm input for one multiome sample, on exactly the nuclei the ATAC
# arm kept.
#   Usage: Rscript make_rna_input.R P111T3 [bin_dir_suffix]
#
# The gene annotation comes from EnsDb.Hsapiens.v86, the same source the CONGAS+
# and inferCNV runners use, so all three sections sit on one coordinate system.

suppressPackageStartupMessages({
  library(data.table); library(Matrix); library(EnsDb.Hsapiens.v86)
})
ROOT <- "<MINTCNA_ROOT>/ATAC_RNA"; setwd(ROOT)
args <- commandArgs(trailingOnly = TRUE)
S <- if (length(args)) args[1] else "P111T3"
CHRS <- paste0("chr", 1:22)
OUT <- file.path("Data/mintcna_input", S)

# ---- the cells the ATAC arm ended up with -----------------------------------
# align_omics() needs the two modalities on the same nuclei; the ATAC side has
# already been through the >3,000-fragment depth QC, so start from its cell table.
ct <- fread(file.path(OUT, "cell_types.tsv"))
cat(sprintf("ATAC arm kept %d nuclei (%d reference)\n", nrow(ct), sum(ct$normal == "normal")))

# ---- GEX matrix --------------------------------------------------------------
man  <- fread("Data/meta/sample_manifest.csv")[care_id == S]
gfil <- strsplit(man$gex_files, ";")[[1]]
mtx  <- readMM(file.path("Data/gex_matrices", grep("matrix",   gfil, value = TRUE)))
bc   <- fread(file.path("Data/gex_matrices", grep("barcodes", gfil, value = TRUE)),
              header = FALSE)$V1
ft   <- fread(file.path("Data/gex_matrices", grep("features", gfil, value = TRUE)),
              header = FALSE)
setnames(ft, c("gene_id", "gene_name", "type")[seq_len(ncol(ft))])
stopifnot(ncol(mtx) == length(bc), nrow(mtx) == nrow(ft))

keep <- match(ct$cell, bc)
if (anyNA(keep)) {
  cat(sprintf("  %d ATAC nuclei absent from the GEX matrix; dropped from BOTH arms\n",
              sum(is.na(keep))))
  ct <- ct[!is.na(keep)]; keep <- keep[!is.na(keep)]
}
mtx <- mtx[, keep, drop = FALSE]
colnames(mtx) <- ct$cell
rownames(mtx) <- ft$gene_id

# ---- gene annotation, autosomes only, one row per gene ----------------------
gg <- as.data.table(ensembldb::genes(EnsDb.Hsapiens.v86,
        columns = c("gene_id", "gene_name", "seq_name", "gene_seq_start", "gene_seq_end")))
gg <- gg[, .(gene_id, gene_name, chr = paste0("chr", seqnames),
             start = start, end = end)][chr %in% CHRS]
gg <- unique(gg, by = "gene_id")

i <- match(rownames(mtx), gg$gene_id)
ok <- !is.na(i)
mtx <- mtx[ok, , drop = FALSE]
ann <- gg[i[ok]]
setorder(ann, chr, start)
mtx <- mtx[match(ann$gene_id, rownames(mtx)), , drop = FALSE]
cat(sprintf("%d genes on chr1-22 with coordinates, %d nuclei\n", nrow(mtx), ncol(mtx)))

# ---- write -------------------------------------------------------------------
# run_mintcna.R reads rna_counts with read.table(header=TRUE, row.names=1)
d <- as.data.table(as.matrix(mtx))
d <- cbind(gene = ann$gene_id, d)
fwrite(d, file.path(OUT, "rna_counts.tsv"), sep = "\t", quote = FALSE)
fwrite(ann[, .(gene_id, gene_name, chr, start, end)],
       file.path(OUT, "gene_annotation.tsv"), sep = "\t", quote = FALSE)
# keep the two arms' cell tables identical
fwrite(ct, file.path(OUT, "cell_types.tsv"), sep = "\t")
cat("wrote rna_counts.tsv, gene_annotation.tsv; cell_types.tsv re-synced\n")

cat(sprintf("\nUMI per nucleus: median %s   genes detected: median %d\n",
            format(round(median(Matrix::colSums(mtx))), big.mark = ","),
            round(median(Matrix::colSums(mtx > 0)))))
