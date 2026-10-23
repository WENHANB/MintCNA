#!/usr/bin/env Rscript
# Round 3, the two remaining questions.
#
# (a) The DNA path ends with modSaRa2::smooth(R=10, t=2). That is DNAcopy-style
#     OUTLIER smoothing, not a running mean -- it rewrites only points more than
#     t SD from their neighbourhood and leaves everything else alone. At 7,000
#     reads per bin an outlier is an artefact. At 1.7 fragments per bin the upper
#     tail is real accessibility, so clipping it may be either a free win or a
#     distortion. Measure it.
#
# (b) Dividing by the reference assumes the tumour-to-reference ACCESSIBILITY
#     ratio is flat at 1. It is not: malignant glial cells and myeloid cells open
#     different chromatin. If the residual still tracks reference accessibility
#     after the division, a spline in log2(reference) removes what division could
#     not -- something no DNA pipeline needs and no ATAC caller currently does.

suppressPackageStartupMessages({library(data.table); library(matrixStats)})
ROOT <- "<MINTCNA_ROOT>/ATAC_RNA"; setwd(ROOT)
CHRS <- paste0("chr", 1:22); W <- 5e5

bins  <- fread("Data/ref/bins_500kb_hg38.tsv.gz")
cells <- fread("Data/cell_lists/P111T3_cells.tsv")
cnt   <- fread("Data/mintcna_input/P111T3_counts_500kb.tsv.gz")
tr <- fread("Result/truth/P111T3_truth_100kb.tsv.gz")[chr %in% CHRS]
tr <- tr[, .(truth = mean(log2, na.rm = TRUE)), by = .(chr, bin = start %/% W)]
bins <- merge(bins, tr, by = c("chr", "bin"), all.x = TRUE); bins[, key := paste(chr, bin)]

mal <- cells[cell_type == "Malignant", barcode]; ref <- cells[is_reference == TRUE, barcode]
cnt[, key := paste(chr, bin)]
bi <- match(cnt$key, bins$key); ci <- match(cnt$barcode, c(mal, ref)); ok <- !is.na(bi) & !is.na(ci)
Y <- matrix(0L, nrow(bins), length(c(mal, ref)), dimnames = list(bins$key, c(mal, ref)))
Y[cbind(bi[ok], ci[ok])] <- cnt$count[ok]
IM <- seq_along(mal); IR <- length(mal) + seq_along(ref)
keep <- bins$frac_valid >= .5 & bins$frac_blacklist <= .1 &
        bins$mapp >= .7 & bins$gc >= .2 & bins$gc <= .8 & !is.na(bins$truth)
Yk <- Y[keep, ]; bk <- bins[keep]; tv <- bk$truth

cor_rows <- function(M, t) {
  Mc <- M - rep(colMeans(M), each = nrow(M)); tc <- t - mean(t)
  as.vector(crossprod(Mc, tc)) / (sqrt(colSums(Mc^2)) * sqrt(sum(tc^2)))
}
score <- function(L, tag) {
  L <- sweep(L, 2, colMedians(L, na.rm = TRUE))
  data.table(variant = tag,
             per_cell   = round(median(cor_rows(L[, IM, drop=FALSE], tv), na.rm = TRUE), 3),
             pseudobulk = round(cor(rowMeans(L[, IM, drop=FALSE]), tv), 3))
}
lib  <- colSums(Yk)
p_in <- rowSums(Yk[, IR, drop=FALSE]) / sum(Yk[, IR, drop=FALSE])
L0   <- log2((Yk + 1) / (outer(p_in, lib) + 1))
res  <- list(score(L0, "baseline: pooled reference, a=1"))

# (a) outlier smoothing, exactly as the DNA path applies it
if (requireNamespace("modSaRa2", quietly = TRUE)) {
  for (tt in c(2, 3)) {
    S <- modSaRa2::smooth(L0, R = 10, t = tt)
    dimnames(S) <- dimnames(L0)
    res[[length(res)+1]] <- score(S, sprintf("+ modSaRa2::smooth(R=10, t=%d)", tt))
    cat(sprintf("  t=%d rewrote %.1f%% of entries\n", tt, 100*mean(abs(S - L0) > 1e-9)))
  }
}

# (b) is there accessibility structure LEFT after dividing by the reference?
pb  <- rowMeans(L0[, IM, drop = FALSE])
acc <- log2(p_in + 1e-8)
cat(sprintf("\ncorrelation of the pseudobulk RESIDUAL with reference accessibility: %.3f\n",
            cor(pb, acc)))
cat(sprintf("  (with the truth removed first: %.3f)\n",
            cor(residuals(lm(pb ~ tv)), acc)))
fit <- lm(pb ~ splines::ns(acc, df = 5))
cat(sprintf("  a 5-df spline in accessibility explains %.1f%% of the residual\n\n",
            100 * summary(fit)$r.squared))

# apply that spline correction per cell and see whether the truth survives it
Lc <- L0 - matrix(predict(fit), nrow(L0), ncol(L0))
res[[length(res)+1]] <- score(Lc, "+ spline correction on accessibility")

print(rbindlist(res))
