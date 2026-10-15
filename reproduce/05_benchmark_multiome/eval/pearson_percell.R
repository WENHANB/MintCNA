#!/usr/bin/env Rscript
# Per-CELL Pearson against the bulk WGS truth, split by tumour subtype.
#   Usage: Rscript pearson_percell.R
#
# Pseudobulk Pearson gives one number per method per sample. Per-cell Pearson
# gives a distribution, which shows something pseudobulk cannot: whether a method
# is accurate because most of its cells are right, or because a noisy population
# averages out. Those are different claims and only the second is visible here.
#
# Two things this metric cannot do, both of which have to be stated rather than
# hidden:
#   - CONGAS+ emits ONE profile per clone, so every cell in a clone gets an
#     identical correlation. Its "distribution" has as many distinct values as
#     there are clones, usually 1-4. It is plotted, but the box is an artefact of
#     clone count, not of cell-to-cell variability.
#   - epiAneufinder emits three discrete states per cell, so a per-cell
#     correlation is computed on a 3-level vector. That is legitimate (Pearson on
#     ordinal data is a monotone summary) but it is coarser than the continuous
#     callers and will sit lower for that reason alone.
#
# Correlations are computed by matrix algebra rather than a loop: centre each
# cell's profile and the truth, then r = (Mc %*% tc) / (||Mc|| * ||tc||).

suppressPackageStartupMessages({library(data.table)})
ROOT <- "<MINTCNA_ROOT>/ATAC_RNA"
setwd(ROOT)
CHRS <- paste0("chr", 1:22)
source("Code/eval/pearson_helpers.R")

SAMPLES <- c("P59T2","P59T3","P104T1","P104T2","P106T1","P106T2","P100T1","P100T3",
             "P103T1","P103T2","P105T1","P105T3","P108T1","P108T2","P109T1","P109T2",
             "P110T1","P110T2","P111T2","P111T3")

# per-cell correlation of a cells x bins matrix against a truth vector
cor_rows <- function(M, tv) {
  ok <- is.finite(tv) & colSums(is.finite(M)) == nrow(M)
  if (sum(ok) < 100) return(rep(NA_real_, nrow(M)))
  M <- M[, ok, drop = FALSE]; tv <- tv[ok]
  Mc <- M - rowMeans(M)
  tc <- tv - mean(tv)
  den <- sqrt(rowSums(Mc^2)) * sqrt(sum(tc^2))
  r <- as.vector(Mc %*% tc) / den
  r[!is.finite(r)] <- NA_real_
  r
}

man <- fread("Data/meta/sample_manifest.csv"); setkey(man, care_id)

res <- rbindlist(lapply(SAMPLES, function(s) {
  tf <- file.path("Result/truth", paste0(s, "_truth_100kb.tsv.gz"))
  if (!file.exists(tf)) return(NULL)
  truth <- fread(tf)[chr %in% CHRS]
  cells <- fread(file.path("Data/cell_lists", paste0(s, "_cells.tsv")))
  keep  <- cells[cell_type == "Malignant", barcode]
  message(s, ": ", length(keep), " malignant nuclei")

  M <- percell_matrices(s, truth, keep)     # named list of cells x truthbins
  rbindlist(lapply(names(M), function(nm) {
    if (is.null(M[[nm]]) || !nrow(M[[nm]])) return(NULL)
    r <- cor_rows(M[[nm]], truth$log2)
    data.table(care_id = s, subtype = man[s, subtype], method = nm,
               barcode = rownames(M[[nm]]), r = r)
  }))
}))

res <- res[is.finite(r)]
fwrite(res, "Result/qc/pearson_percell.tsv.gz", sep = "\t")

cat("\n=== per-cell Pearson, by subtype and method ===\n")
print(res[, .(n_cells = .N, median = round(median(r), 3),
              q25 = round(quantile(r, .25), 3), q75 = round(quantile(r, .75), 3)),
          by = .(subtype, method)][order(subtype, -median)], nrow = 30)
