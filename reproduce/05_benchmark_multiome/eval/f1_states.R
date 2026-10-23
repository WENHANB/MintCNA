#!/usr/bin/env Rscript
# Bin-level F1 against the bulk WGS state calls.
#
#   TP  truth is altered, the method calls altered, and the DIRECTION agrees
#   FP  the method calls altered where the truth is neutral, or in the wrong
#       direction
#   FN  the truth is altered and the method calls neutral, or the wrong direction
#
# A bin called in the wrong direction counts as both FP and FN: it invents an
# event that is not there and misses the one that is. That is the conservative
# reading and it is what "the direction has to agree" means when gains and losses
# are pooled into one score.
#
#   precision = TP/(TP+FP)   recall = TP/(TP+FN)   F1 = 2PR/(P+R)
#
# TRUTH: GATK's own `cnv_call` (-/0/+), not a threshold of ours. It carries the
# segment-level significance test already -- where it disagrees with a plain
# |log2| > 0.1 rule, almost every case is GATK calling neutral because the
# segment rests on too few WGS bins to support its amplitude.
#
# DISCRETISATION happens PER CELL, then the consensus is the modal state over
# malignant nuclei. The rule is AtaCNV's own published one, "elements larger than
# 1.25 were defined as gains, those less than 0.75 were defined as losses" (Cell
# Rep Methods 2025), applied there to Copy-scAT, inferCNV and CopyKAT. Using a
# competitor's published threshold removes any suggestion that we chose a cut
# that suits us.
#
# Per CELL, not per pseudobulk track. A first version thresholded the
# cell-averaged track and gave epiAneufinder F1 = 0.000: averaging 0/1/2 state
# codes across nuclei puts its track between -0.11 and 0, which never reaches
# log2(0.75) = -0.415, so nothing was ever called. The threshold was designed for
# a single cell's copy ratio and only means anything there.
#
# The same rule is LOSSLESS for the methods that already emit states, so there is
# no special case: log2(1/2) = -1.000 < -0.415 -> loss, log2(3/2) = +0.585 >
# +0.322 -> gain, log2(2/2) = 0 -> neutral. Every method is handled identically.
suppressPackageStartupMessages(library(data.table))
setwd("<MINTCNA_ROOT>/ATAC_RNA")
CHRS <- paste0("chr", 1:22)
source("<REPRODUCE_DIR>/lib/eval/pearson_helpers.R")     # percell_matrices(): every method as
                                          # cells x truth-bins, in log2
LO <- log2(0.75); HI <- log2(1.25)
man <- fread("Data/meta/sample_manifest.csv"); setkey(man, care_id)
SAMPLES <- c("P59T2","P59T3","P104T1","P104T2","P106T1","P106T2","P100T1","P100T3",
             "P103T1","P103T2","P105T1","P105T3","P108T1","P108T2","P109T1","P109T2",
             "P110T1","P110T2","P111T2","P111T3")

f1 <- rbindlist(lapply(SAMPLES, function(s) {
  tf <- file.path("Result/truth", paste0(s, "_truth_100kb.tsv.gz"))
  if (!file.exists(tf)) return(NULL)
  truth <- fread(tf)[chr %in% CHRS]
  keep  <- fread(sprintf("Data/cell_lists/%s_cells.tsv", s))[cell_type == "Malignant", barcode]
  M <- percell_matrices(s, truth, keep)
  tv <- truth$cnv_call
  gb <- paste(truth$chr, truth$start %/% 5e5)
  ok <- tv %in% c("-", "0", "+")
  rbindlist(lapply(names(M), function(m) {
    X <- M[[m]][, ok, drop = FALSE]; y <- tv[ok]
    if (!nrow(X)) return(NULL)
    # AtaCNV's procedure exactly: threshold PER CELL, then take the mode for each
    # bin over all tumour cells, then count bins. Their words: "we first converted
    # the copy ratios into copy number states using the fixed threshold, then
    # obtained the common result across all tumor cells by taking the mode for
    # each bin", and the counting unit is the bin.
    n_lo <- colSums(X < LO, na.rm = TRUE)
    n_hi <- colSums(X > HI, na.rm = TRUE)
    n_ne <- colSums(is.finite(X)) - n_lo - n_hi
    pred <- fifelse(n_lo > n_ne & n_lo >= n_hi, "-",
             fifelse(n_hi > n_ne & n_hi > n_lo, "+", "0"))
    # Aggregate to the SAME 500 kb grid the Pearson row uses. The truth's native
    # grid is 100 kb; leaving F1 there would put the two rows of one panel on two
    # different grids. The bin's state is the one covering most of it.
    g  <- gb[ok]
    dd <- data.table(g = g, pred = pred, y = y)[
            , .(pred = names(which.max(table(pred))), y = names(which.max(table(y)))), by = g]
    pred <- dd$pred; y <- dd$y
    tp <- sum(pred != "0" & pred == y)
    fp <- sum(pred != "0" & pred != y)     # neutral truth, or wrong direction
    fn <- sum(y != "0" & pred != y)        # missed, or wrong direction
    # precision is UNDEFINED when a method calls nothing, not zero. epiAneufinder
    # calls nothing at all in 17 of 20 samples; writing 0 there says "everything
    # it called was wrong" when it called nothing, and it drags the summary to 0
    # while the three samples where it does call are 1.000, 0.750 and 0.000.
    p <- if (tp + fp == 0) NA_real_ else tp/(tp + fp)
    r <- if (tp + fn == 0) NA_real_ else tp/(tp + fn)
    data.table(care_id = s, method = m, n_cell = nrow(X), n_bin = ncol(X),
               TP = tp, FP = fp, FN = fn, called = tp + fp, precision = p, recall = r,
               F1 = if (is.na(p) || is.na(r) || p + r == 0) NA_real_ else 2*p*r/(p + r))
  }))
}))
f1[, subtype := fifelse(grepl("^Oligo", man[care_id, subtype]), "Oligo", "Astro")]
fwrite(f1, "Result/qc/f1_states.tsv", sep = "\t")

# POOLED across samples, not the median of per-sample values. With 17 of 20
# samples undefined for one method, a median is a statement about how often the
# method is undefined, not about how well it does when it is not.
pool <- function(d) d[, {
  p <- if (sum(TP)+sum(FP) == 0) NA_real_ else sum(TP)/(sum(TP)+sum(FP))
  r <- if (sum(TP)+sum(FN) == 0) NA_real_ else sum(TP)/(sum(TP)+sum(FN))
  .(TP = sum(TP), FP = sum(FP), FN = sum(FN),
    n_undef = sum(is.na(precision)), precision = p, recall = r,
    F1 = if (is.na(p)||is.na(r)||p+r == 0) NA_real_ else 2*p*r/(p+r))
}, by = method]

cat("=== bin-level, POOLED over 20 samples; direction must agree ===\n\n")
o <- pool(f1)[order(-F1, na.last = TRUE)]
for (c in c("precision","recall","F1")) o[[c]] <- round(o[[c]], 3)
print(o, row.names = FALSE)
cat("\n  n_undef = samples in which the method called nothing (precision undefined)\n")
cat("\n=== by subtype (pooled) ===\n\n")
print(dcast(pool(f1[, .(TP,FP,FN,precision,method,subtype)][, .(TP,FP,FN,precision), by=.(method,subtype)]),
            method ~ ., value.var = "F1")[0])   # placeholder; per-subtype below
for (st in c("Oligo","Astro")) {
  cat(sprintf("  %s:\n", st)); q <- pool(f1[subtype == st])[order(-F1, na.last = TRUE)]
  for (c in c("precision","recall","F1")) q[[c]] <- round(q[[c]], 3)
  print(q[, .(method, precision, recall, F1)], row.names = FALSE)
}
W <- dcast(f1, care_id ~ method, value.var = "F1")
p <- sapply(setdiff(names(W), c("care_id", "MintCNA")), function(m)
       suppressWarnings(wilcox.test(W$MintCNA, W[[m]], paired = TRUE)$p.value))
cat("\npaired Wilcoxon vs MintCNA (Holm):\n")
print(data.table(method = names(p), p_holm = signif(p.adjust(p, "holm"), 3))[order(p_holm)])
