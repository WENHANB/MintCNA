#!/usr/bin/env Rscript
# Measure sequencing depth / coverage for the 20 benchmark samples and compare
# against what a standard 10x Multiome run delivers.
#
# ATAC side comes free: the authors' ArchR metadata already carries nFrags,
# TSSEnrichment, NucleosomeRatio. RNA side has to be computed from the mtx
# matrices -- their metadata carries cell types but no nCount/nFeature.
#
# Both sides are restricted to the paired, QC-passing nuclei we actually use
# (Data/cell_lists), not the raw CellRanger output, so the numbers are the ones
# that describe the data going into the benchmark.

suppressPackageStartupMessages({library(data.table); library(Matrix)})
ROOT <- "<MINTCNA_ROOT>/ATAC_RNA"
setwd(ROOT)

S <- c("P59T2","P59T3","P104T1","P104T2","P106T1","P106T2","P100T1","P100T3","P103T1","P103T2",
       "P105T1","P105T3","P108T1","P108T2","P109T1","P109T2","P110T1","P110T2","P111T2","P111T3")

man <- fread("Data/meta/sample_manifest.csv")[care_id %in% S]
setkey(man, care_id)

out <- rbindlist(lapply(S, function(samp) {
  gsm <- man[samp, gsm_gex]
  # the gex_files column names the exact prefix, e.g. GSM9625460_SJ03-0_matrix.mtx.gz
  mtxf <- file.path("Data/gex_matrices", grep("matrix", strsplit(man[samp, gex_files], ";")[[1]], value = TRUE))
  bcf  <- file.path("Data/gex_matrices", grep("barcodes", strsplit(man[samp, gex_files], ";")[[1]], value = TRUE))
  if (!file.exists(mtxf)) { message("MISSING gex: ", samp); return(NULL) }

  m  <- readMM(mtxf)
  bc <- fread(bcf, header = FALSE)$V1
  stopifnot(ncol(m) == length(bc))

  # keep only the paired nuclei that survived ATAC QC -- same barcode namespace
  keep_bc <- fread(file.path("Data/cell_lists", paste0(samp, "_cells.tsv")))$barcode
  idx <- which(bc %in% keep_bc)

  umi  <- Matrix::colSums(m)
  gene <- Matrix::colSums(m > 0)

  data.table(care_id = samp,
             n_gex_all   = ncol(m),
             n_paired    = length(idx),
             umi_med_all = median(umi),
             umi_med     = if (length(idx)) median(umi[idx])  else NA_real_,
             umi_q25     = if (length(idx)) quantile(umi[idx], .25) else NA_real_,
             umi_q75     = if (length(idx)) quantile(umi[idx], .75) else NA_real_,
             gene_med    = if (length(idx)) median(gene[idx]) else NA_real_,
             gene_q25    = if (length(idx)) quantile(gene[idx], .25) else NA_real_,
             n_genes_detected = sum(Matrix::rowSums(m[, idx, drop = FALSE]) > 0))
}), fill = TRUE)

fwrite(out, "Result/qc/rna_depth.tsv", sep = "\t")
print(out)
