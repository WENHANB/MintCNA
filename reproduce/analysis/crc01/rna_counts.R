#!/usr/bin/env Rscript
# CRC01 scTrio-seq2 RNA (GSE97693): RSEM expected counts per cell -> genes x cells count matrix and gene
# annotation (GENCODE v19, hg19 autosomes), the RNA input of MintCNA.
#   Rscript reproduce/analysis/crc01/rna_counts.R
# Input : Data/CRC/RNA/rsem/<SRR>.genes.results  (STAR + RSEM, one file per cell)
#         Data/CRC/SraRunTable.csv, Data/CRC/sample.csv  (SRR run -> GEO sample -> cell name)
#         Data/CRC/RNA/Reference/star/gencode.v19.annotation.gtf
# Output: Data/CRC/mintcna_input/rna_counts.tsv, Data/CRC/mintcna_input/genes.tsv
# Run from the root of the MintCNA repository with MINTCNA_DATA set (see reproduce/README.md).
suppressPackageStartupMessages({ library(data.table); library(rtracklayer) })
ROOT <- Sys.getenv("MINTCNA_DATA"); if (!nzchar(ROOT)) stop("set MINTCNA_DATA (see reproduce/README.md)")
OUT <- file.path(ROOT, "Data/CRC/mintcna_input"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ---- gene coordinates ----------------------------------------------------------------------------
gtf <- import(file.path(ROOT, "Data/CRC/RNA/Reference/star/gencode.v19.annotation.gtf"))
gg  <- gtf[gtf$type == "gene"]
ann <- data.table(gene_id = gg$gene_id, chr = as.character(seqnames(gg)), start = start(gg), end = end(gg))
ann <- ann[chr %in% paste0("chr", 1:22)][!duplicated(gene_id)]

# ---- expected counts per cell ------------------------------------------------------------------
files <- list.files(file.path(ROOT, "Data/CRC/RNA/rsem"), pattern = "genes\\.results$", recursive = TRUE, full.names = TRUE)
srr <- sub("\\..*$", "", basename(files))
counts <- lapply(files, function(f) fread(f, select = c("gene_id", "expected_count")))
genes <- counts[[1]]$gene_id
M <- sapply(counts, function(d) d$expected_count[match(genes, d$gene_id)])
rownames(M) <- genes

# ---- SRR run -> cell name (GEO sample title "scTrioSeq2Rna_<cell>") -------------------------------
sra <- fread(file.path(ROOT, "Data/CRC/SraRunTable.csv"))[, .(Run, Accession = `GEO_Accession (exp)`)]
smp <- fread(file.path(ROOT, "Data/CRC/sample.csv"))[, .(Accession, Title)]
map <- merge(sra, smp, by = "Accession")
colnames(M) <- sub("^scTrioSeq2Rna_", "", map$Title[match(srr, map$Run)])
stopifnot(!anyNA(colnames(M)))

# ---- autosomal genes with coordinates, ordered along the genome ----------------------------------
ann <- ann[gene_id %in% rownames(M)]
setorder(ann, chr, start)
M <- M[ann$gene_id, , drop = FALSE]
fwrite(data.table(gene_id = rownames(M), as.data.table(M)), file.path(OUT, "rna_counts.tsv"), sep = "\t")
fwrite(ann, file.path(OUT, "genes.tsv"), sep = "\t")
cat(sprintf("%d genes x %d cells -> %s\n", nrow(M), ncol(M), OUT))
