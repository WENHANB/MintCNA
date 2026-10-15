#!/usr/bin/env Rscript
# What does eps = 0.20 actually cover, in each mode?
#
# eps is an ABSOLUTE band on the fused segment mean z. In multi mode z is a
# coverage-weighted average of the RNA and DNA segment means, which are quantile-
# normalised to sd_target 0.278 and 0.691 respectively; in dna mode z IS the DNA
# segment mean. Mixing in the lower-variance RNA arm shrinks z, so the same
# number is a different quantile of the two distributions and the two arms of the
# ablation would not be thresholded comparably.
suppressPackageStartupMessages(library(data.table))
RR <- "<MINTCNA_ROOT>/MintCNA_git"
for (f in list.files(file.path(RR,"R"), "\\.[rR]$", recursive=TRUE, full.names=TRUE))
  try(suppressWarnings(sys.source(f, globalenv())), silent=TRUE)
setwd("<MINTCNA_ROOT>/ATAC_RNA")
a <- commandArgs(TRUE); S <- if (length(a)) a[1] else "P111T3"; set.seed(1)

zdist <- function(suffix) {
  sc <- readRDS(sprintf("Result/MintCNA/%s_%s/%s_scan_h15.rds", S, suffix, S))
  mode <- sc$mode
  tau  <- sc$res_null$q[which.min(abs(sc$res_null$probs - 0.90))]
  cells <- sample(names(sc$cell_list), 200)
  # Mirror call_cells() exactly: cps_track[[cell]] IS the per-chromosome list
  # (no $scores wrapper), and per-chromosome positions become global indices
  # through reindex_cps, not through a which() on ref$chr.
  chr_levels  <- sc$chr_levels
  idx_by_chr  <- split(seq_len(nrow(sc$ref)), factor(sc$ref$chr, levels = chr_levels))
  chr_lengths <- vapply(chr_levels, function(c) length(idx_by_chr[[c]]), integer(1))
  z <- unlist(lapply(cells, function(cn) {
    X  <- as_cell_matrix(sc$cell_list[[cn]], mode)
    sb <- sc$cps_track[[cn]]
    cps_by_chr <- lapply(chr_levels, function(ch) {
      lm <- sb[[ch]]$local_maxima
      if (is.null(lm) || !nrow(lm)) integer(0) else sort(lm[abs(lm$score) > tau, "position"])
    })
    names(cps_by_chr) <- chr_levels
    g <- reindex_cps(cps_by_chr, chr_lengths)
    g <- g[g > 1 & g < ncol(X)]
    if (!length(g)) return(numeric(0))
    fusion_omics(X, g)$z
  }))
  z[is.finite(z)]
}
for (suf in c("multi","atac")) {
  f <- sprintf("Result/MintCNA/%s_%s/%s_scan_h15.rds", S, suf, S)
  if (!file.exists(f)) { cat(suf, ": no scan yet\n"); next }
  z <- zdist(suf)
  cat(sprintf("\n%s  (%s segments from 200 cells)\n", toupper(suf), format(length(z), big.mark=",")))
  cat(sprintf("  sd(z) = %.3f    IQR = %.3f\n", sd(z), IQR(z)))
  cat(sprintf("  %-24s %8s %8s\n", "eps", "|z|<=eps", "equivalent eps"))
  for (e in c(0.10, 0.15, 0.20, 0.25, 0.30))
    cat(sprintf("  eps = %.2f  covers %6.1f%% of segments\n", e, 100*mean(abs(z) <= e)))
  assign(paste0("z_", suf), z, envir = globalenv())
}
if (exists("z_multi") && exists("z_atac")) {
  cat(sprintf("SCALE %s  n_rna_na %.3f  sd_multi %.4f  sd_atac %.4f  ratio %.3f  eps_match %.3f\n",
      S, mean(!is.finite(readRDS(sprintf("Result/MintCNA/%s_multi/%s_scan_h15.rds", S, S))$cell_list[[1]]["RNA",])),
      sd(z_multi), sd(z_atac), sd(z_atac)/sd(z_multi),
      quantile(abs(z_atac), mean(abs(z_multi) <= 0.20))))
  cat("\n=== what eps on the ATAC-only arm matches eps = 0.20 on multi? ===\n")
  target <- mean(abs(z_multi) <= 0.20)
  eq <- quantile(abs(z_atac), target)
  cat(sprintf("  multi eps 0.20 leaves %.1f%% of segments neutral\n", 100*target))
  cat(sprintf("  the same %.1f%% on the ATAC-only arm needs eps = %.3f\n", 100*target, eq))
  cat(sprintf("  ratio sd(z_atac)/sd(z_multi) = %.2f  ->  scaled eps = %.3f\n",
              sd(z_atac)/sd(z_multi), 0.20*sd(z_atac)/sd(z_multi)))
}
