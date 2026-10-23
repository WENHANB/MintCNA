#!/usr/bin/env Rscript
# Every MintCNA configuration on disk, scored SEPARATELY BY SUBTYPE, with
# gain and loss kept apart.
#
# The astrocytomas are being tuned and the oligodendrogliomas are not, so the
# question is never "is this better on the astrocytomas" alone -- a setting that
# lifts them and drops the fourteen oligodendrogliomas is a worse global setting,
# and the only honest alternative to that is a per-subtype parameter, which needs
# its own argument in the Methods. Both columns are therefore always printed.
suppressPackageStartupMessages(library(data.table))
setwd("<MINTCNA_ROOT>/ATAC_RNA")
CHRS <- paste0("chr",1:22); GRID <- 5e5
src <- readLines("<REPRODUCE_DIR>/05_benchmark_multiome/eval/pearson_vs_truth.R")
i0 <- grep("^project <- function", src); i1 <- grep("^# --- per-method pseudobulk", src)[1]
eval(parse(text = paste(src[i0:(i1-1)], collapse = "\n")))
man <- fread("Data/meta/sample_manifest.csv"); setkey(man, care_id)

samples <- sub("_multi$", "", basename(list.dirs("Result/MintCNA", recursive = FALSE)))
res <- rbindlist(lapply(samples, function(S) {
  tf <- sprintf("Result/truth/%s_truth_100kb.tsv.gz", S)
  if (!file.exists(tf)) return(NULL)
  truth <- fread(tf)[chr %in% CHRS]
  keep  <- fread(sprintf("Data/cell_lists/%s_cells.tsv", S))[cell_type=="Malignant", barcode]
  fs <- list.files(sprintf("Result/MintCNA/%s_multi", S),
                   sprintf("^%s_results_h[0-9]+a[0-9]+e[0-9]+\\.RData$", S), full.names=TRUE)
  rbindlist(lapply(fs, function(f) {
    tag <- sub(".*_results_(.*)\\.RData$", "\\1", f)
    e <- new.env()
    if (inherits(try(load(f, e), silent=TRUE), "try-error")) return(NULL)
    Sm <- get("cna_state_matrix", e); rf <- get("ref", e)
    i <- colnames(Sm) %in% keep; if (!any(i)) return(NULL)
    fg <- project(truth, rf$chr, rf$start, rf$end, rowMeans(Sm[,i,drop=FALSE]==3))
    fl <- project(truth, rf$chr, rf$start, rf$end, rowMeans(Sm[,i,drop=FALSE]==1))
    mv <- project(truth, rf$chr, rf$start, rf$end, log2(rowMeans(Sm[,i,drop=FALSE])/2))
    d <- data.table(chr=truth$chr, gb=truth$start %/% GRID, t=truth$log2, fg=fg, fl=fl, mv=mv)
    d <- d[, lapply(.SD, mean, na.rm=TRUE), by=.(chr,gb)][is.finite(t) & is.finite(mv)]
    pr <- function(c_, t_) { tp<-sum(c_&t_); fp<-sum(c_&!t_); fn<-sum(!c_&t_)
      list(TP=tp, FP=fp, FN=fn, prec=tp/max(tp+fp,1), rec=tp/max(tp+fn,1)) }
    g <- pr(d$fg>0.5, d$t> 0.1); l <- pr(d$fl>0.5, d$t< -0.1)
    data.table(care_id=S, tag=tag, r=cor(d$mv, d$t),
               g_TP=g$TP, g_FP=g$FP, g_FN=g$FN, g_prec=g$prec, g_rec=g$rec,
               l_TP=l$TP, l_FP=l$FP, l_FN=l$FN, l_prec=l$prec, l_rec=l$rec)
  }))
}))
res[, subtype := fifelse(grepl("^Oligo", man[care_id, subtype]), "Oligo", "Astro")]
res[, `:=`(h     = as.integer(sub("^h(\\d+)a.*", "\\1", tag)),
           alpha = as.numeric(sub(".*a(\\d{3})e.*", "\\1", tag))/100,
           eps   = as.numeric(sub(".*e(\\d{3})$", "\\1", tag))/100)]
fwrite(res, "Result/qc/mintcna_subtype_scores.tsv", sep="\t")

agg <- res[, .(n=.N, r=round(mean(r),3),
               g_rec=round(sum(g_TP)/max(sum(g_TP)+sum(g_FN),1),3),
               g_prec=round(sum(g_TP)/max(sum(g_TP)+sum(g_FP),1),3),
               l_rec=round(sum(l_TP)/max(sum(l_TP)+sum(l_FN),1),3),
               l_prec=round(sum(l_TP)/max(sum(l_TP)+sum(l_FP),1),3)),
           by=.(subtype, h, alpha, eps)]
agg[, F1 := round(2*(g_rec+l_rec)/2*(g_prec+l_prec)/2 /
                  max(((g_rec+l_rec)/2 + (g_prec+l_prec)/2), 1e-9), 3)]
cat("=== ASTRO (the ones being tuned) ===\n\n")
print(agg[subtype=="Astro"][order(-r)], nrows=40)
cat("\n=== OLIGO at the same settings — the guard ===\n\n")
print(agg[subtype=="Oligo"][order(-r)], nrows=40)
cat("\n=== settings present for BOTH subtypes, ranked by the worse of the two ===\n\n")
b <- dcast(agg, h+alpha+eps ~ subtype, value.var=c("r","g_rec","l_rec"))
b <- b[!is.na(r_Astro) & !is.na(r_Oligo)]
b[, worse := pmin(r_Astro, r_Oligo)]
print(b[order(-worse)], nrows=30)
