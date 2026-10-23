#!/usr/bin/env Rscript
# Is a sample's subclone split just sequencing depth? AUC of per-nucleus depth
# separating the two clones: 0.5 = depth carries no information about the split,
# 1.0 = depth alone reproduces it. Reported alongside the between-clone state
# difference, which is the thing the split is supposed to be about.
.libPaths(c("<R_LIBS_USER>", .libPaths()))
suppressPackageStartupMessages(library(data.table))
ROOT <- "<MINTCNA_ROOT>"
TAG  <- "h15a005e020"
auc <- function(x, g) { r <- rank(x); n1 <- sum(g); n0 <- sum(!g)
  a <- (sum(r[g]) - n1*(n1+1)/2)/(n1*n0); max(a, 1-a) }

out <- rbindlist(lapply(c("P111T2","P111T3","P100T1"), function(S) {
  f <- sprintf("%s/ATAC_RNA/Result/MintCNA/%s_multi/%s_results_%s.RData", ROOT, S, S, TAG)
  if (!file.exists(f)) return(NULL)
  e <- new.env(); load(f, e)
  cna <- get("cna_state_matrix", e); sub <- get("subclones", e)
  keep <- fread(sprintf("%s/ATAC_RNA/Data/cell_lists/%s_cells.tsv", ROOT, S))[
            cell_type == "Malignant", barcode]
  i <- colnames(cna) %in% keep; cna <- cna[, i, drop = FALSE]; sub <- sub[i]
  sz <- names(sort(table(sub), decreasing = TRUE))
  g  <- sub == as.integer(sz[1])
  # depth proxy: the per-nucleus total of the input count matrix, if stored
  d <- NULL
  for (nm in c("depth","total_counts","lib_size","cell_depth")) if (!is.null(e[[nm]])) d <- e[[nm]]
  if (is.null(d)) { Y <- e[["Y"]]; if (!is.null(Y)) d <- colSums(Y) }
  a <- if (is.null(d)) NA_real_ else auc(as.numeric(d)[i], g)
  m1 <- rowMeans(cna[, g, drop = FALSE]); m2 <- rowMeans(cna[, !g, drop = FALSE])
  data.table(sample = S, n = ncol(cna), n_C1 = sum(g), n_C2 = sum(!g),
             depth_auc = a, state_diff_mean = mean(abs(m1-m2)),
             state_diff_p90 = quantile(abs(m1-m2), 0.90),
             frac_bins_diff_gt_0.5 = mean(abs(m1-m2) > 0.5))
}))
print(out, row.names = FALSE)
