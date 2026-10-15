#!/usr/bin/env Rscript
# Does GC correction earn its place once there is a per-bin normal reference?
#
# epiAneufinder and AtaCNV both do cell-wise LOESS GC correction, and the earlier
# evidence against it here tested MintCNA's STRATIFIED MEDIAN SCALING, not LOESS.
# Those are different estimators, so the objection was not actually answered.
#
# The structural question is why each method needs it:
#
#   epiAneufinder  takes NO normal-cell reference -- it z-scores segment means
#                  within a chromosome. Nothing else removes GC bias, so LOESS
#                  is doing real work.
#   AtaCNV         takes a per-bin normal reference AND does GC correction.
#   MintCNA (here) takes a per-bin pooled reference.
#
# If GC bias is a property of the library chemistry, it is shared between tumour
# and reference nuclei and cancels in the ratio. If it varies cell to cell -- as
# it does in scDNA, where every cell gets its own whole-genome amplification --
# it does not cancel and must be corrected per cell. So measure the cell-to-cell
# spread first, then test the correction three ways.

suppressPackageStartupMessages({library(data.table); library(matrixStats)})
ROOT <- "<MINTCNA_ROOT>/ATAC_RNA"; setwd(ROOT)
CHRS <- paste0("chr", 1:22); W <- 5e5
set.seed(1)

bins  <- fread("Data/ref/bins_500kb_hg38.tsv.gz")
cells <- fread("Data/cell_lists/P111T3_cells.tsv")
cnt   <- fread("Data/mintcna_input/P111T3_counts_500kb.tsv.gz")
tr <- fread("Result/truth/P111T3_truth_100kb.tsv.gz")[chr %in% CHRS]
tr <- tr[, .(truth = mean(log2, na.rm = TRUE)), by = .(chr, bin = start %/% W)]
bins <- merge(bins, tr, by = c("chr", "bin"), all.x = TRUE)[, key := paste(chr, bin)]

mal <- cells[cell_type == "Malignant", barcode]; ref <- cells[is_reference == TRUE, barcode]
# subsample malignant nuclei: per-cell LOESS on 1,710 cells is slow and 400 is
# plenty to separate estimators that differ by 0.02 in median correlation
mal <- sample(mal, min(400, length(mal)))
use <- c(mal, ref)
cnt[, key := paste(chr, bin)]
bi <- match(cnt$key, bins$key); ci <- match(cnt$barcode, use); ok <- !is.na(bi) & !is.na(ci)
Y <- matrix(0L, nrow(bins), length(use), dimnames = list(bins$key, use))
Y[cbind(bi[ok], ci[ok])] <- cnt$count[ok]
IM <- seq_along(mal); IR <- length(mal) + seq_along(ref)

keep <- bins$frac_valid >= .5 & bins$frac_blacklist <= .1 & bins$mapp >= .5 &
        bins$gc >= .2 & bins$gc <= .8 & !is.na(bins$truth)
Yk <- Y[keep, ]; bk <- bins[keep]; tv <- bk$truth; gc <- bk$gc
lib <- colSums(Yk)
cat(sprintf("%d bins, %d malignant (subsampled), %d reference\n\n", nrow(Yk), length(mal), length(ref)))

cor_rows <- function(M, t) {
  Mc <- M - rep(colMeans(M), each = nrow(M)); tc <- t - mean(t)
  as.vector(crossprod(Mc, tc)) / (sqrt(colSums(Mc^2)) * sqrt(sum(tc^2)))
}
score <- function(Z, tag) {
  p <- rowSums(Z[, IR, drop = FALSE]); p <- p / sum(p)
  l <- log2((Z + 1) / (outer(p, colSums(Z)) + 1))
  l <- sweep(l, 2, colMedians(l, na.rm = TRUE))
  data.table(variant = tag,
             per_cell   = round(median(cor_rows(l[, IM, drop=FALSE], tv), na.rm = TRUE), 3),
             pseudobulk = round(cor(rowMeans(l[, IM, drop=FALSE]), tv), 3))
}

# ---- 1. how much does GC bias actually vary between nuclei? -----------------
cpm <- sweep(Yk, 2, lib, "/") * 1e6
rgc <- apply(cpm, 2, function(v) cor(v, gc))
cat("correlation of a nucleus's coverage with bin GC:\n")
cat(sprintf("  malignant  median %+.3f   sd across nuclei %.3f   range %+.3f to %+.3f\n",
            median(rgc[IM]), sd(rgc[IM]), min(rgc[IM]), max(rgc[IM])))
cat(sprintf("  reference  median %+.3f   sd across nuclei %.3f\n",
            median(rgc[IR]), sd(rgc[IR])))
cat(sprintf("  pooled tumour %+.3f  vs pooled reference %+.3f  (difference %+.3f)\n\n",
            cor(rowSums(Yk[, IM, drop=FALSE]), gc), cor(rowSums(Yk[, IR, drop=FALSE]), gc),
            cor(rowSums(Yk[, IM, drop=FALSE]), gc) - cor(rowSums(Yk[, IR, drop=FALSE]), gc)))

res <- list(score(Yk, "no GC correction (current)"))

# ---- 2. cell-wise LOESS, the epiAneufinder / AtaCNV recipe ------------------
loess_correct <- function(M) {
  out <- M
  for (j in seq_len(ncol(M))) {
    v <- M[, j]; mv <- mean(v)
    fit <- tryCatch(predict(loess(v ~ gc, span = 0.75, degree = 2,
                                  family = "gaussian", surface = "direct")),
                    error = function(e) rep(mv, length(v)))
    fit[!is.finite(fit) | fit < 1e-8] <- 1e-8   # numerical floor, per the plan
    out[, j] <- v * mv / fit
  }
  out
}
cat("running cell-wise LOESS on", ncol(cpm), "nuclei...\n")
gcm <- loess_correct(cpm)
rgc2 <- apply(gcm, 2, function(v) cor(v, gc))
cat(sprintf("  |cor with GC| median: %.3f before -> %.3f after\n",
            median(abs(rgc)), median(abs(rgc2))))
# put it back on a count-like scale so the same pseudocount means the same thing
res[[2]] <- score(sweep(gcm, 2, lib / 1e6, "*"), "cell-wise LOESS (epiAneufinder recipe)")

# ---- 3. one LOESS curve estimated on pooled counts, applied to every cell ----
# If GC bias is a shared chemistry property, this captures all of it and adds
# none of the noise of fitting a curve per cell to ~2 counts per bin.
pv  <- rowSums(Yk); fit <- predict(loess(pv ~ gc, span = 0.75, degree = 2,
                                         family = "gaussian", surface = "direct"))
fit[!is.finite(fit) | fit < 1e-8] <- 1e-8
res[[3]] <- score(Yk * (mean(pv) / fit), "pooled LOESS, one curve for all cells")

print(rbindlist(res))
