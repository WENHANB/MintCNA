#!/usr/bin/env Rscript
# RNA-arm input of MintCNA for one multiome sample: the GEX count matrix restricted to the nuclei the ATAC
# arm kept, with gene coordinates from EnsDb.Hsapiens.v86 (autosomes), ordered along the genome.
#   Rscript reproduce/analysis/multiome/make_rna_input.R <care_id>         (called by run_mintcna.sh)
# Input : ATAC_RNA/Data/mintcna_input/<care_id>/cell_types.tsv (scripts/make_atac_input.R),
#         ATAC_RNA/Data/gex_matrices/ (download_geo.sh), ATAC_RNA/Data/meta/sample_manifest.csv
# Output: ATAC_RNA/Data/mintcna_input/<care_id>/{rna_counts.tsv, gene_annotation.tsv}; cell_types.tsv is
#         rewritten without the nuclei absent from the GEX matrix so both arms list the same cells
# Run from the root of the MintCNA repository with MINTCNA_DATA set (see reproduce/README.md).
suppressPackageStartupMessages({ library(data.table); library(Matrix); library(EnsDb.Hsapiens.v86) })
ROOT <- Sys.getenv("MINTCNA_DATA"); if (!nzchar(ROOT)) stop("set MINTCNA_DATA (see reproduce/README.md)")
setwd(file.path(ROOT, "ATAC_RNA"))
args <- commandArgs(trailingOnly = TRUE)
S <- if (length(args)) args[1] else "P111T2"
CHRS <- paste0("chr", 1:22)
OUT <- file.path("Data/mintcna_input", S)

ct <- fread(file.path(OUT, "cell_types.tsv"))
man  <- fread("Data/meta/sample_manifest.csv")[care_id == S]
gfil <- strsplit(man$gex_files, ";")[[1]]
mtx  <- readMM(file.path("Data/gex_matrices", grep("matrix",   gfil, value = TRUE)))
bc   <- fread(file.path("Data/gex_matrices", grep("barcodes", gfil, value = TRUE)), header = FALSE)$V1
ft   <- fread(file.path("Data/gex_matrices", grep("features", gfil, value = TRUE)), header = FALSE)
setnames(ft, c("gene_id", "gene_name", "type")[seq_len(ncol(ft))])
stopifnot(ncol(mtx) == length(bc), nrow(mtx) == nrow(ft))

keep <- match(ct$cell, bc)
if (anyNA(keep)) {
  cat(sprintf("%d ATAC nuclei absent from the GEX matrix are dropped from both arms\n", sum(is.na(keep))))
  ct <- ct[!is.na(keep)]; keep <- keep[!is.na(keep)]
}
mtx <- mtx[, keep, drop = FALSE]
colnames(mtx) <- ct$cell; rownames(mtx) <- ft$gene_id

gg <- as.data.table(ensembldb::genes(EnsDb.Hsapiens.v86,
        columns = c("gene_id", "gene_name", "seq_name", "gene_seq_start", "gene_seq_end")))
gg <- unique(gg[, .(gene_id, gene_name, chr = paste0("chr", seqnames), start, end)][chr %in% CHRS], by = "gene_id")
i <- match(rownames(mtx), gg$gene_id); ok <- !is.na(i)
mtx <- mtx[ok, , drop = FALSE]
ann <- gg[i[ok]]; setorder(ann, chr, start)
mtx <- mtx[match(ann$gene_id, rownames(mtx)), , drop = FALSE]

fwrite(cbind(gene = ann$gene_id, as.data.table(as.matrix(mtx))), file.path(OUT, "rna_counts.tsv"), sep = "\t", quote = FALSE)
fwrite(ann[, .(gene_id, gene_name, chr, start, end)], file.path(OUT, "gene_annotation.tsv"), sep = "\t", quote = FALSE)
fwrite(ct, file.path(OUT, "cell_types.tsv"), sep = "\t")
cat(sprintf("%s: %d genes x %d nuclei -> %s\n", S, nrow(mtx), ncol(mtx), OUT))
