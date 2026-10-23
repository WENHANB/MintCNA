#!/usr/bin/env Rscript
# Step 9 of the plan: look at the segment-mean distribution BEFORE choosing the
# TGMM neutral band. Segment means do not depend on eps -- eps only enters the
# state assignment -- so the whole distribution is available from CP-A.
#
# The decisive view is the last one: where do segments that the bulk WGS calls
# LOST actually sit, versus segments it calls NEUTRAL. The band has to separate
# those two, and nothing else.
suppressPackageStartupMessages({library(data.table); library(matrixStats)})
setwd("<MINTCNA_ROOT>/ATAC_RNA")
for (f in list.files("<REPRODUCE_DIR>/../MintCNA","\\.[rR]$",
                     recursive=TRUE, full.names=TRUE)) try(suppressWarnings(sys.source(f, globalenv())), silent=TRUE)

sc  <- readRDS("Result/MintCNA/P111T3_multi/P111T3_scan_h5.rds")
ref <- sc$ref
cat(sprintf("cells %d   bins %d   h %d   simT %s\n", length(sc$cps_track), nrow(ref), sc$h, sc$simT))
sdt <- median(sapply(sc$cell_list, function(m) sd(m["DNA",], na.rm=TRUE)))
cat(sprintf("QN sd_target (DNA row): %.3f   <- the scale eps acts on\n\n", sdt))

tr <- fread("Result/truth/P111T3_truth_100kb.tsv.gz")[chr %in% paste0("chr",1:22)]
tr <- tr[, .(t=mean(log2,na.rm=TRUE)), by=.(chr,bin=start%/%5e5)][, k:=paste(chr,bin)]
tv <- tr$t[match(paste(ref$chr, ref$start%/%5e5), tr$k)]

chr_lv <- sc$chr_levels
idx <- split(seq_len(nrow(ref)), factor(ref$chr, levels=chr_lv))
clen <- vapply(chr_lv, function(c) length(idx[[c]]), integer(1))

for (alpha in c(0.05, 0.10)) {
  tau <- sc$res_null$q[sc$res_null$probs == (1-alpha)]
  out <- rbindlist(lapply(names(sc$cps_track), function(cn) {
    cps <- lapply(chr_lv, function(ch) { lm <- sc$cps_track[[cn]][[ch]]$local_maxima
                                          sort(lm[abs(lm$score) > tau, "position"]) })
    names(cps) <- chr_lv
    g <- reindex_cps(cps, clen); M <- sc$cell_list[[cn]]
    g <- g[g > 1 & g < ncol(M)]
    if (!length(g)) return(NULL)
    fo <- fusion_omics(M, g)
    data.table(cell=cn, z=fo$z, len=fo$length_bins, s=fo$start_bin, e=fo$end_bin)
  }))
  out <- out[len > 1]
  # what the truth says about each segment
  out[, truth := vapply(seq_len(.N), function(i) mean(tv[s[i]:e[i]], na.rm=TRUE), numeric(1))]
  out[, cls := ifelse(truth < -0.15, "LOST", ifelse(abs(truth) < 0.05, "NEUTRAL", "other"))]

  cat(sprintf("=== alpha = %.2f   (tau = %.3f) ===\n", alpha, tau))
  cat(sprintf("  change-points per cell: median %.0f  (IQR %.0f-%.0f)\n",
      median(out[, .N, by=cell]$N), quantile(out[, .N, by=cell]$N,.25), quantile(out[, .N, by=cell]$N,.75)))
  cat(sprintf("  segment length (bins):  median %.0f   (%.1f Mb)\n", median(out$len), median(out$len)*0.5))
  cat(sprintf("  segment mean z:         %s\n",
      paste(sprintf("%s=%.2f", c("p05","p25","p50","p75","p95"),
                    quantile(out$z, c(.05,.25,.5,.75,.95))), collapse="  ")))
  cat("  by what the WGS says:\n")
  print(out[cls != "other", .(n=.N, med_z=round(median(z),3),
                              q25=round(quantile(z,.25),3), q75=round(quantile(z,.75),3),
                              med_len=as.numeric(median(len))), by=cls])
  # how well would each candidate band separate them?
  cat("  a segment is called non-neutral if |z| > eps:\n")
  for (eps in c(0.05,0.10,0.15,0.20,0.30,0.50)) {
    tpr <- mean(out[cls=="LOST", z < -eps]); fpr <- mean(out[cls=="NEUTRAL", abs(z) > eps])
    cat(sprintf("    eps=%.2f   recall on LOST %.2f   false-call on NEUTRAL %.2f\n", eps, tpr, fpr))
  }
  cat("\n")
  fwrite(out, sprintf("Result/MintCNA/P111T3_multi/segmeans_a%.2f.tsv", alpha), sep="\t")
}
