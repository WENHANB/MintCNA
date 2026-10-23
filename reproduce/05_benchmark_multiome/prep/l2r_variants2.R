#!/usr/bin/env Rscript
# Round 2. The reference estimator and the GC question are settled; these are the
# three that are left: how good could a per-cell profile possibly be at this depth,
# does a deeper cross-sample reference help, and how much does smoothing buy.

suppressPackageStartupMessages({library(data.table); library(matrixStats)})
ROOT <- "<MINTCNA_ROOT>/ATAC_RNA"
setwd(ROOT)
CHRS <- paste0("chr", 1:22); W <- 5e5

bins  <- fread("Data/ref/bins_500kb_hg38.tsv.gz")
cells <- fread("Data/cell_lists/P111T3_cells.tsv")
cnt   <- fread("Data/mintcna_input/P111T3_counts_500kb.tsv.gz")
tr <- fread("Result/truth/P111T3_truth_100kb.tsv.gz")[chr %in% CHRS]
tr <- tr[, .(truth = mean(log2, na.rm = TRUE)), by = .(chr, bin = start %/% W)]
bins <- merge(bins, tr, by = c("chr", "bin"), all.x = TRUE)
bins[, key := paste(chr, bin)]

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
lib <- colSums(Yk)
p_in <- rowSums(Yk[, IR, drop=FALSE]) / sum(Yk[, IR, drop=FALSE])

cat(sprintf("truth amplitude at 500 kb: sd = %.3f\n", sd(tv)))
lam <- median(lib) * p_in
cat(sprintf("median cell: %d fragments over %d bins, lambda = %.2f per bin\n",
            as.integer(median(lib)), length(tv), median(lam)))
# A per-cell profile is truth + Poisson noise. Correlation with the truth cannot
# exceed sd_signal / sqrt(sd_signal^2 + sd_noise^2), whatever the method.
sn <- median(1.443 / sqrt(lam))
cat(sprintf("Poisson log2 noise = %.2f  ->  per-cell Pearson CEILING = %.3f\n\n",
            sn, sd(tv) / sqrt(sd(tv)^2 + sn^2)))

res <- list()
res[[1]] <- score(log2((Yk + 1) / (outer(p_in, lib) + 1)), "within-sample reference (478 nuclei)")

# --- a deeper reference: reference nuclei pooled over all 20 samples ----------
# Myeloid open chromatin is largely shared across patients, so pooling multiplies
# the depth behind every bin. The risk is a patient-specific batch shift.
allref <- rbindlist(lapply(list.files("Data/cell_lists", "_cells.tsv$", full.names = TRUE), function(f) {
  s <- sub("_cells.tsv$", "", basename(f))
  d <- fread(f)[is_reference == TRUE, .(barcode)]; d[, care_id := s]; d
}))
pooled_file <- "Result/qc/pooled_reference_500kb.tsv.gz"
if (!file.exists(pooled_file)) {
  pr <- rbindlist(lapply(unique(allref$care_id), function(s) {
    f <- sprintf("Data/mintcna_input/%s_counts_500kb.tsv.gz", s)
    if (!file.exists(f) || file.size(f) < 1000) return(NULL)   # still being written
    d <- tryCatch(fread(f), error = function(e) NULL)
    if (is.null(d) || !nrow(d)) return(NULL)
    d[barcode %in% allref[care_id == s, barcode], .(count = sum(count)), by = .(chr, bin)]
  }))
  if (nrow(pr)) {
    cat(sprintf("[pooled reference built from %d samples]\n", uniqueN(allref$care_id)))
    fwrite(pr[, .(count = sum(count)), by = .(chr, bin)], pooled_file, sep = "\t")
  }
}
if (file.exists(pooled_file)) {
  pr <- fread(pooled_file); pr[, key := paste(chr, bin)]
  p_all <- pr$count[match(bk$key, pr$key)]; p_all[is.na(p_all)] <- 0
  p_all <- p_all / sum(p_all)
  res[[2]] <- score(log2((Yk + 1) / (outer(p_all, lib) + 1)), "cross-sample reference (all 20)")
} else cat("[cross-sample reference skipped: only P111T3 is binned so far]\n")

# --- smoothing ---------------------------------------------------------------
# The accessibility residual is high-frequency and the CNA signal is whole-arm,
# so a running mean should trade almost nothing away.
L0 <- log2((Yk + 1) / (outer(p_in, lib) + 1))
# rows are already in (chr, bin) order, so each chromosome is one contiguous
# block; smooth inside a block and never across a centromere or a chromosome end
blocks <- split(seq_len(nrow(L0)), factor(bk$chr, levels = unique(bk$chr)))
for (w in c(3, 5, 11, 21, 51)) {
  S <- L0
  for (ix in blocks) {
    sub <- L0[ix, , drop = FALSE]
    S[ix, ] <- apply(sub, 2, function(v)
      as.numeric(frollmean(v, min(w, length(v)), align = "center", na.rm = TRUE)))
  }
  fin <- rowSums(!is.finite(S)) == 0
  res[[length(res) + 1]] <- local({
    L <- sweep(S[fin, ], 2, colMedians(S[fin, ], na.rm = TRUE))
    data.table(variant = sprintf("+ running mean %2d bins (%.1f Mb)", w, w * .5),
               per_cell = round(median(cor_rows(L[, IM, drop=FALSE], tv[fin]), na.rm=TRUE), 3),
               pseudobulk = round(cor(rowMeans(L[, IM, drop=FALSE]), tv[fin]), 3))
  })
}
print(rbindlist(res))
