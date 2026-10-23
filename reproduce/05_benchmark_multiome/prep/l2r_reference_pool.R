#!/usr/bin/env Rscript
# Within-sample reference vs a reference pooled across patients.
#
# The reference set is small in several samples (56 nuclei in P106T2, 109 in
# P110T1), and a thin reference puts its own sampling noise into every cell's
# denominator. Myeloid open chromatin is largely shared between patients, so
# pooling multiplies the depth behind each bin -- at the risk of importing a
# patient- or batch-specific shift. Test it where it should matter most.
#   Usage: Rscript l2r_reference_pool.R P111T3 P110T1 ...

suppressPackageStartupMessages({library(data.table); library(matrixStats)})
ROOT <- "<MINTCNA_ROOT>/ATAC_RNA"; setwd(ROOT)
CHRS <- paste0("chr", 1:22); W <- 5e5
SAMPLES <- commandArgs(trailingOnly = TRUE)
if (!length(SAMPLES)) SAMPLES <- c("P111T3", "P110T1")

bins <- fread("Data/ref/bins_500kb_hg38.tsv.gz")[, key := paste(chr, bin)]
avail <- sub("_counts_500kb.tsv.gz$", "",
             list.files("Data/mintcna_input", "_counts_500kb.tsv.gz$"))
avail <- avail[file.size(sprintf("Data/mintcna_input/%s_counts_500kb.tsv.gz", avail)) > 1000]

# one pooled reference profile from every sample that is binned, and the
# leave-one-out version so a sample never uses its own nuclei twice
per_sample_ref <- rbindlist(lapply(avail, function(s) {
  rb <- fread(sprintf("Data/cell_lists/%s_cells.tsv", s))[is_reference == TRUE, barcode]
  d  <- fread(sprintf("Data/mintcna_input/%s_counts_500kb.tsv.gz", s))[barcode %in% rb]
  cbind(care_id = s, d[, .(count = sum(count)), by = .(chr, bin)])
}))
cat(sprintf("pooled reference drawn from %d samples: %s\n\n",
            length(avail), paste(avail, collapse = ", ")))

cor_rows <- function(M, t) {
  Mc <- M - rep(colMeans(M), each = nrow(M)); tc <- t - mean(t)
  as.vector(crossprod(Mc, tc)) / (sqrt(colSums(Mc^2)) * sqrt(sum(tc^2)))
}

out <- rbindlist(lapply(SAMPLES, function(S) {
  cells <- fread(sprintf("Data/cell_lists/%s_cells.tsv", S))
  cnt   <- fread(sprintf("Data/mintcna_input/%s_counts_500kb.tsv.gz", S))[, key := paste(chr, bin)]
  tr <- fread(sprintf("Result/truth/%s_truth_100kb.tsv.gz", S))[chr %in% CHRS]
  tr <- tr[, .(truth = mean(log2, na.rm = TRUE)), by = .(chr, bin = start %/% W)]
  b  <- merge(bins, tr, by = c("chr", "bin"), all.x = TRUE)

  mal <- cells[cell_type == "Malignant", barcode]; rf <- cells[is_reference == TRUE, barcode]
  cn <- c(mal, rf)
  bi <- match(cnt$key, b$key); ci <- match(cnt$barcode, cn); ok <- !is.na(bi) & !is.na(ci)
  Y <- matrix(0L, nrow(b), length(cn), dimnames = list(b$key, cn))
  Y[cbind(bi[ok], ci[ok])] <- cnt$count[ok]
  IM <- seq_along(mal); IR <- length(mal) + seq_along(rf)
  keep <- b$frac_valid >= .5 & b$frac_blacklist <= .1 & b$mapp >= .7 &
          b$gc >= .2 & b$gc <= .8 & !is.na(b$truth)
  Yk <- Y[keep, ]; bk <- b[keep]; tv <- bk$truth; lib <- colSums(Yk)

  mk <- function(p, tag) {
    L <- log2((Yk + 1) / (outer(p, lib) + 1))
    L <- modSaRa2::smooth(L, R = 10, t = 2); dimnames(L) <- dimnames(Yk)
    L <- sweep(L, 2, colMedians(L, na.rm = TRUE))
    data.table(sample = S, n_ref = length(rf), reference = tag,
               per_cell   = round(median(cor_rows(L[, IM, drop=FALSE], tv), na.rm=TRUE), 3),
               pseudobulk = round(cor(rowMeans(L[, IM, drop=FALSE]), tv), 3))
  }
  norm1 <- function(v) { v[is.na(v)] <- 0; v / sum(v) }
  p_in  <- rowSums(Yk[, IR, drop = FALSE]) / sum(Yk[, IR, drop = FALSE])

  loo <- per_sample_ref[care_id != S, .(count = sum(count)), by = .(chr, bin)][, key := paste(chr, bin)]
  all <- per_sample_ref[,             .(count = sum(count)), by = .(chr, bin)][, key := paste(chr, bin)]
  rbind(mk(p_in,                                        sprintf("within-sample (%d nuclei)", length(rf))),
        mk(norm1(loo$count[match(bk$key, loo$key)]),    "other samples only (leave-one-out)"),
        mk(norm1(all$count[match(bk$key, all$key)]),    "all samples pooled"),
        mk(norm1(0.5 * p_in / sum(p_in) +
                 0.5 * norm1(loo$count[match(bk$key, loo$key)])), "50/50 within + other"))
}))
print(out)
fwrite(out, "Result/qc/l2r_reference_pool.tsv", sep = "\t")
