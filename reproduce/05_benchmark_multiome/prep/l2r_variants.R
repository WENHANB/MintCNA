#!/usr/bin/env Rscript
# Which L2R construction should the ATAC arm use? Decide it on P111T3 against
# the matched bulk WGS instead of by argument.
#
# MintCNA's DNA path (R/preprocessing/dna_preprocess.R) is:
#   filter bins (gc .2-.8, mapp >= .9)
#     -> per-cell median-scaling for mappability (step .01) then GC (step .05)
#     -> CPM per cell
#     -> log2((x + 0.01) / (rowMedian over normal cells + 0.01))
#     -> subtract the per-cell median
# Every one of those steps was tuned at ~7,000 reads per bin. Here it is 3.5
# fragments per bin, and three of them break at that count level. This script
# measures each swap in isolation.

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
setkey(bins, key)

mal <- cells[cell_type == "Malignant", barcode]
ref <- cells[is_reference == TRUE,     barcode]

# dense bins x cells
cnt[, key := paste(chr, bin)]
bi <- match(cnt$key, bins$key); ci <- match(cnt$barcode, c(mal, ref))
ok <- !is.na(bi) & !is.na(ci)
Y <- matrix(0L, nrow(bins), length(c(mal, ref)),
            dimnames = list(bins$key, c(mal, ref)))
Y[cbind(bi[ok], ci[ok])] <- cnt$count[ok]
IM <- seq_along(mal); IR <- length(mal) + seq_along(ref)

# ---- bin filter -------------------------------------------------------------
# mapp >= .9 is a 5 kb threshold. Averaged over 100 sub-bins the median 500 kb
# window sits at .886, so .9 lands inside the bulk of the distribution and cuts
# 62% of the genome for no reason. .7 removes the genuinely broken bins only.
keep <- bins$frac_valid >= .5 & bins$frac_blacklist <= .1 &
        bins$mapp >= .7 & bins$gc >= .2 & bins$gc <= .8 & !is.na(bins$truth)
Yk <- Y[keep, ]; bk <- bins[keep]; tv <- bk$truth
cat(sprintf("bins kept %d / %d ; malignant %d ; reference %d\n\n",
            sum(keep), nrow(bins), length(mal), length(ref)))

cor_rows <- function(M, t) {          # correlation of every column of M with t
  Mc <- M - rep(colMeans(M), each = nrow(M)); tc <- t - mean(t)
  as.vector(crossprod(Mc, tc)) / (sqrt(colSums(Mc^2)) * sqrt(sum(tc^2)))
}
score <- function(L, tag) {
  L  <- sweep(L, 2, colMedians(L, na.rm = TRUE))       # median shift, as in DNA
  pc <- cor_rows(L[, IM, drop = FALSE], tv)
  pb <- cor(rowMeans(L[, IM, drop = FALSE]), tv)
  data.table(variant = tag, per_cell = round(median(pc, na.rm = TRUE), 3),
             pc_q25 = round(quantile(pc, .25, na.rm = TRUE), 3),
             pc_q75 = round(quantile(pc, .75, na.rm = TRUE), 3),
             pseudobulk = round(pb, 3))
}

# ---- the three reference estimators -----------------------------------------
lib <- colSums(Yk)
cpm <- sweep(Yk, 2, lib, "/") * 1e6
p_pool <- rowSums(Yk[, IR, drop = FALSE]) / sum(Yk[, IR, drop = FALSE])
r_med  <- rowMedians(cpm[, IR, drop = FALSE])
r_mean <- rowMeans(cpm[, IR, drop = FALSE])

res <- list()
# expected counts under the pooled reference profile: E = lib_cell * p_ref[bin].
# The ratio is then formed in COUNT space, where a Poisson pseudocount means
# something.
E <- outer(p_pool, lib)
res[[1]] <- score(log2((Yk + 0.5) / (E + 0.5)), "pooled reference, a=0.5")
res[[2]] <- score(log2((Yk + 1.0) / (E + 1.0)), "pooled reference, a=1")
res[[3]] <- score(log2((Yk + 0.1) / (E + 0.1)), "pooled reference, a=0.1")
res[[4]] <- score(log2((cpm + .01) / (r_med  + .01)), "DNA default: median over normal cells")
res[[5]] <- score(log2((cpm + .01) / (r_mean + .01)), "mean over normal cells")

# ---- does per-cell GC/mapp correction still earn its place? -----------------
source("<REPRODUCE_DIR>/../MintCNA/preprocess.R")
gcm <- GC.MAP.normalization(Yk, bk)$RC_norm
lib2 <- colSums(gcm); cpm2 <- sweep(gcm, 2, lib2, "/") * 1e6
p2 <- rowSums(gcm[, IR, drop = FALSE]) / sum(gcm[, IR, drop = FALSE])
res[[6]] <- score(log2((cpm2 + .01) / (rowMedians(cpm2[, IR, drop=FALSE]) + .01)),
                  "FULL DNA PATH (per-cell GC/mapp + median ref)")
res[[7]] <- score(log2((gcm + 0.5) / (outer(p2, lib2) + 0.5)),
                  "per-cell GC/mapp + pooled reference")

# GC/mappability cancel exactly in the ratio, so correcting for them per cell can
# only add the noise of estimating them from ~3 counts per bin. Estimate the same
# curves once on the pooled data instead and the question becomes moot.
r_all <- rowSums(Yk); ord <- order(bk$gc)
res[[8]] <- score(log2((Yk + 0.5) / (E + 0.5)), "pooled reference, a=0.5 (repeat)")

r <- rbindlist(res)
print(r[order(-per_cell)])
fwrite(r, "Result/qc/l2r_variants_P111T3.tsv", sep = "\t")
