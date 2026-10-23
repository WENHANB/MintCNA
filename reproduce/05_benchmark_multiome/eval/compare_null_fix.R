#!/usr/bin/env Rscript
# Three configurations, same scans, same eps, differing only in the null and alpha:
#   h15a010e020        alpha 0.10, the null as shipped   (everything reported so far)
#   fixnull_a005e020   alpha 0.05, corrected null        (predicted to reproduce it)
#   fixnull_a010e020   alpha 0.10, corrected null        (genuinely looser)
suppressPackageStartupMessages(library(data.table))
setwd("<MINTCNA_ROOT>/ATAC_RNA")
CHRS <- paste0("chr",1:22); GRID <- 5e5
src <- readLines("<REPRODUCE_DIR>/05_benchmark_multiome/eval/pearson_vs_truth.R")
i0 <- grep("^project <- function", src); i1 <- grep("^# --- per-method pseudobulk", src)[1]
eval(parse(text = paste(src[i0:(i1-1)], collapse="\n")))
man <- fread("Data/meta/sample_manifest.csv"); setkey(man, care_id)
S <- c("P59T2","P59T3","P100T1","P100T3","P103T1","P103T2","P104T1","P104T2","P105T1","P105T3",
       "P106T1","P106T2","P108T1","P108T2","P109T1","P109T2","P110T1","P110T2","P111T2","P111T3")
TAGS <- c("h15a010e020","fixnull_a005e020","fixnull_a010e020")

res <- rbindlist(lapply(S, function(s) {
  truth <- fread(sprintf("Result/truth/%s_truth_100kb.tsv.gz", s))[chr %in% CHRS]
  keep  <- fread(sprintf("Data/cell_lists/%s_cells.tsv", s))[cell_type=="Malignant", barcode]
  rbindlist(lapply(TAGS, function(tg) {
    f <- sprintf("Result/MintCNA/%s_multi/%s_results_%s.RData", s, s, tg)
    if (!file.exists(f)) return(NULL)
    e <- new.env(); if (inherits(try(load(f,e), silent=TRUE), "try-error")) return(NULL)
    Sm <- base::get("cna_state_matrix", e); rf <- base::get("ref", e)
    i <- colnames(Sm) %in% keep
    v  <- project(truth, rf$chr, rf$start, rf$end, log2(rowMeans(Sm[,i,drop=FALSE])/2))
    fg <- project(truth, rf$chr, rf$start, rf$end, rowMeans(Sm[,i,drop=FALSE]==3))
    fl <- project(truth, rf$chr, rf$start, rf$end, rowMeans(Sm[,i,drop=FALSE]==1))
    d  <- data.table(chr=truth$chr, gb=truth$start %/% GRID, v=v, fg=fg, fl=fl, t=truth$log2)
    d  <- d[, lapply(.SD, mean, na.rm=TRUE), by=.(chr,gb)][is.finite(v)&is.finite(t)]
    pr <- function(c_,t_){tp<-sum(c_&t_);fp<-sum(c_&!t_);fn<-sum(!c_&t_)
      c(rec=tp/max(tp+fn,1), prec=tp/max(tp+fp,1))}
    g <- pr(d$fg>0.5, d$t>0.1); l <- pr(d$fl>0.5, d$t< -0.1)
    data.table(care_id=s, tag=tg, subtype=substr(man[s,subtype],1,5),
               r=cor(d$v,d$t), pct_alt=mean(Sm[,i]!=2),
               rec=(g["rec"]+l["rec"])/2, prec=(g["prec"]+l["prec"])/2)
  }))
}))
fwrite(res, "Result/qc/null_fix_comparison.tsv", sep="\t")

cat("=== Pearson, all 20 samples ===\n\n")
w <- dcast(res, care_id + subtype ~ tag, value.var="r")
for (c in TAGS) if (c %in% names(w)) w[[c]] <- round(w[[c]],3)
w[, d_a05 := round(fixnull_a005e020 - h15a010e020, 3)]
w[, d_a10 := round(fixnull_a010e020 - h15a010e020, 3)]
print(w[order(subtype, care_id)], nrows=25)
cat("\n=== summary by tag ===\n\n")
print(res[, .(n=.N, r_median=round(median(r),3), r_mean=round(mean(r),3),
              pct_alt=round(100*mean(pct_alt),1), recall=round(mean(rec),3),
              precision=round(mean(prec),3)), by=.(tag)][match(TAGS, tag)])
cat("\n=== by subtype ===\n\n")
print(res[, .(r=round(mean(r),3), recall=round(mean(rec),3), precision=round(mean(prec),3),
              pct_alt=round(100*mean(pct_alt),1)), by=.(subtype, tag)][order(subtype)])
cat(sprintf("\nalpha 0.05 on the fixed null vs the old alpha 0.10: median |dr| = %.4f, max %.4f\n",
    median(abs(w$d_a05)), max(abs(w$d_a05))))
