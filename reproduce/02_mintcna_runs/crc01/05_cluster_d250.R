#!/usr/bin/env Rscript
# CRC01 subclone clustering with the Supplementary Methods procedure (package defaults of
# cluster_latent_kmeans: standardised latent, PCA >= 95 % variance capped at 20 PCs, K by
# stability over 200 random 80 % subsamples with mean pairwise ARI and the one-SE rule,
# final k-means with 200 starts), no post-hoc outlier removal. CRC01_ML3_085 is excluded BEFORE
# clustering by the per-cell QC rule (smoothed log2-ratio SD > 3x cohort median; author 2026-10-14),
# via CRC_EXCLUDE; on the remaining 196 cells the stability rule selects K = 2.
#   new: AG-DCAE latent_dim 250 on the 197 QC-passing cells (01_ae_d250.py)
#   old: the published latent_dim 200 file restricted to the same 197 cells (reference only)
.libPaths(c("<R_LIBS_USER>", .libPaths()))
suppressPackageStartupMessages(library(data.table))
ROOT <- "<MINTCNA_ROOT>"; source("<REPRODUCE_DIR>/../MintCNA/call.R")
site_of <- function(x) factor(sub("^CRC\\d+_([A-Za-z]+)[0-9]*_\\d+$", "\\1", x), levels = c("PT","LN","ML","MP"))
run <- function(file, label) {
  Z <- as.matrix(read.csv(file, row.names = 1, check.names = FALSE)); Z <- Z[rownames(Z) != "CRC01_MP1_209", , drop = FALSE]
  EXCL <- strsplit(Sys.getenv("CRC_EXCLUDE", "CRC01_ML3_085"), ",")[[1]]; Z <- Z[!rownames(Z) %in% EXCL, , drop = FALSE]   # per-cell QC exclusion (see utils.R)
  r <- cluster_latent_kmeans(Z, seed = 1)   # documented procedure, K by stability; with CRC01_ML3_085 excluded by the QC rule it selects K = 2
  cat(sprintf("\n== %s: %d cells x %d dims; K_hat = %d\n", label, nrow(Z), ncol(Z), r$K_hat))
  cat("   stability (mean pairwise ARI) by K:", paste(names(r$stability), round(r$stability, 3), collapse = "  "), "\n")
  tb <- table(cluster = r$clusters, site = site_of(names(r$clusters))); print(addmargins(tb))
  data.table(cell = names(r$clusters), cluster = as.integer(r$clusters), K_hat = r$K_hat, set = label)
}
new <- run(file.path(ROOT, "Data/CRC/DNAmeth/processed_l2r/latent_z_dna_d250_197.csv"), "NEW latent_dim 250 (197 cells)")
old <- run(file.path(ROOT, "Data/CRC/DNAmeth/processed_l2r/latent_z_dna.csv"), "OLD latent_dim 200 file, same 197 cells")
m <- merge(new, old, by = "cell", suffixes = c("_new", "_old"))
cat("\nARI new vs old labels (documented procedure on both):", round(mclust::adjustedRandIndex(m$cluster_new, m$cluster_old), 3), "\n")
fwrite(new[, .(cell, cluster, K_hat)], file.path(ROOT, "Result/CRC/subclones_d250.csv")); cat("wrote Result/CRC/subclones_d250.csv\n")
