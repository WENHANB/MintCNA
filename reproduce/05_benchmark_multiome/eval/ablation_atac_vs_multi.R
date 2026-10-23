#!/usr/bin/env Rscript
# The ablation: what does the RNA arm contribute?
#
# The two arms share the binning, the L2R, the AG-DCAE output (by symlink, so
# even the training randomness is the same), h, alpha and st. They differ in
# `mode` -- whether RNA enters the fusion -- and in eps, which is rescaled 0.20 ->
# 0.28 because the fused segment mean is 1.39x wider without RNA averaged into
# it. Running both arms at eps = 0.20 would confound fusion with threshold
# strictness, so that variant is scored too and reported alongside.
suppressPackageStartupMessages(library(data.table))
setwd("<MINTCNA_ROOT>/ATAC_RNA")
CHRS <- paste0("chr",1:22); GRID <- 5e5
src <- readLines("<REPRODUCE_DIR>/05_benchmark_multiome/eval/pearson_vs_truth.R")
i0 <- grep("^project <- function", src); i1 <- grep("^# --- per-method pseudobulk", src)[1]
eval(parse(text = paste(src[i0:(i1-1)], collapse="\n")))
man <- fread("Data/meta/sample_manifest.csv"); setkey(man, care_id)
S <- c("P59T2","P59T3","P100T1","P100T3","P103T1","P103T2","P104T1","P104T2","P105T1","P105T3",
       "P106T1","P106T2","P108T1","P108T2","P109T1","P109T2","P110T1","P110T2","P111T2","P111T3")
# eps 0.25 is the arm the ablation is judged on: it is the value that makes the
# two arms call the same share of the genome, so the comparison is not confounded
# by one of them being thresholded harder. 0.20 and 0.28 are kept to show how
# sensitive the conclusion is to that choice.
ARMS <- list(`ATAC+RNA`      = c("multi","h15a005e020"),
             `ATAC only`     = c("atac", "h15a005e025"),
             `ATAC only e28` = c("atac", "h15a005e028"),
             `ATAC only e20` = c("atac", "h15a005e020"))

res <- rbindlist(lapply(S, function(s) {
  truth <- fread(sprintf("Result/truth/%s_truth_100kb.tsv.gz", s))[chr %in% CHRS]
  keep  <- fread(sprintf("Data/cell_lists/%s_cells.tsv", s))[cell_type=="Malignant", barcode]
  rbindlist(lapply(names(ARMS), function(nm) {
    a <- ARMS[[nm]]
    f <- sprintf("Result/MintCNA/%s_%s/%s_results_%s.RData", s, a[1], s, a[2])
    if (!file.exists(f)) return(NULL)
    e <- new.env(); if (inherits(try(load(f,e), silent=TRUE), "try-error")) return(NULL)
    Sm <- base::get("cna_state_matrix", e); rf <- base::get("ref", e)
    i <- colnames(Sm) %in% keep; if (!any(i)) return(NULL)
    v  <- project(truth, rf$chr, rf$start, rf$end, log2(rowMeans(Sm[,i,drop=FALSE])/2))
    fg <- project(truth, rf$chr, rf$start, rf$end, rowMeans(Sm[,i,drop=FALSE]==3))
    fl <- project(truth, rf$chr, rf$start, rf$end, rowMeans(Sm[,i,drop=FALSE]==1))
    d  <- data.table(chr=truth$chr, gb=truth$start %/% GRID, v=v, fg=fg, fl=fl, t=truth$log2)
    d  <- d[, lapply(.SD, mean, na.rm=TRUE), by=.(chr,gb)][is.finite(v)&is.finite(t)]
    pr <- function(c_,t_){tp<-sum(c_&t_);fp<-sum(c_&!t_);fn<-sum(!c_&t_)
      c(rec=tp/max(tp+fn,1), prec=tp/max(tp+fp,1))}
    g <- pr(d$fg>0.5, d$t>0.1); l <- pr(d$fl>0.5, d$t< -0.1)
    data.table(care_id=s, arm=nm, subtype=substr(man[s,subtype],1,5), n_cell=sum(i),
               r=cor(d$v,d$t), pct_alt=mean(Sm[,i]!=2),
               rec=(g["rec"]+l["rec"])/2, prec=(g["prec"]+l["prec"])/2)
  }))
}))
fwrite(res, "Result/qc/ablation_atac_vs_multi.tsv", sep="\t")

cat("=== per-arm summary ===\n\n")
print(res[, .(n=.N, r_median=round(median(r),3), r_mean=round(mean(r),3),
              pct_alt=round(100*mean(pct_alt),1), recall=round(mean(rec),3),
              precision=round(mean(prec),3)), by=arm][match(names(ARMS), arm)])

cat("\n=== does the eps rescale actually align the two arms? ===\n")
cat("   (if 0.28 is right, ATAC-only pct_alt should sit near ATAC+RNA, and 0.20 should overshoot)\n\n")
p <- dcast(res, care_id + subtype ~ arm, value.var="pct_alt")
for (nm in names(p)[-(1:2)])
  cat(sprintf("   mean pct_alt  %-16s %.1f%%\n", nm, 100*mean(p[[nm]], na.rm=TRUE)))

cat("\n=== the ablation, paired by sample ===\n\n")
w <- dcast(res[arm %in% c("ATAC+RNA","ATAC only")], care_id + subtype ~ arm, value.var="r")
w[, delta := round(`ATAC+RNA` - `ATAC only`, 3)]
for (c in c("ATAC+RNA","ATAC only")) w[[c]] <- round(w[[c]],3)
print(w[order(-delta)], nrows=25)
if (nrow(w) > 3) {
  cat(sprintf("\nRNA helps in %d of %d samples;  mean delta %+.4f, median %+.4f\n",
      sum(w$delta>0), nrow(w), mean(w$delta), median(w$delta)))
  cat(sprintf("paired Wilcoxon p = %.4g;  paired t p = %.4g\n",
      wilcox.test(w$`ATAC+RNA`, w$`ATAC only`, paired=TRUE)$p.value,
      t.test(w$`ATAC+RNA`, w$`ATAC only`, paired=TRUE)$p.value))
  cat("\nby subtype:\n"); print(w[, .(n=.N, mean_delta=round(mean(delta),4),
                                      wins=sum(delta>0)), by=subtype])
}
