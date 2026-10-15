#!/usr/bin/env Rscript
# Gain and loss scored separately, per subtype.
#
# MintCNA emits discrete states, so its calls need no threshold: a bin is called
# gain when more than half the malignant nuclei put it in state 3. AtaCNV emits a
# continuous copy ratio, so a threshold is unavoidable there -- it is set per
# sample to make AtaCNV call the SAME fraction of the genome altered as the truth
# does, which removes the arbitrary choice and makes precision comparable rather
# than a function of how liberal the threshold was.
suppressPackageStartupMessages(library(data.table))
setwd("<MINTCNA_ROOT>/ATAC_RNA")
CHRS <- paste0("chr",1:22); GRID <- 5e5; TAG <- "h15a010e020"
src <- readLines("Code/eval/pearson_vs_truth.R")
i0 <- grep("^project <- function", src); i1 <- grep("^# --- per-method pseudobulk", src)[1]
eval(parse(text = paste(src[i0:(i1-1)], collapse = "\n")))

man <- fread("Data/meta/sample_manifest.csv"); setkey(man, care_id)
SAMPLES <- sub("_multi$","", basename(dirname(
  list.files("Result/MintCNA", sprintf("_results_%s.RData", TAG), recursive=TRUE, full.names=TRUE))))
SAMPLES <- sub(sprintf("_results_%s.RData", TAG), "", SAMPLES)
SAMPLES <- unique(sub("_multi","", basename(list.dirs("Result/MintCNA", recursive=FALSE))))

res <- rbindlist(lapply(SAMPLES, function(S) {
  f <- sprintf("Result/MintCNA/%s_multi/%s_results_%s.RData", S, S, TAG)
  tf <- sprintf("Result/truth/%s_truth_100kb.tsv.gz", S)
  if (!file.exists(f) || !file.exists(tf)) return(NULL)
  truth <- fread(tf)[chr %in% CHRS]
  keep  <- fread(sprintf("Data/cell_lists/%s_cells.tsv", S))[cell_type=="Malignant", barcode]

  e <- new.env(); load(f, e); Sm <- get("cna_state_matrix", e); rf <- get("ref", e)
  i <- colnames(Sm) %in% keep
  fg <- project(truth, rf$chr, rf$start, rf$end, rowMeans(Sm[, i, drop=FALSE] == 3))
  fl <- project(truth, rf$chr, rf$start, rf$end, rowMeans(Sm[, i, drop=FALSE] == 1))
  mv <- project(truth, rf$chr, rf$start, rf$end, log2(rowMeans(Sm[, i, drop=FALSE])/2))

  af <- sprintf("Result/benchmark/AtaCNV/%s/atacnv_result.rds", S)
  av <- rep(NA_real_, nrow(truth))
  if (file.exists(af)) {
    A <- readRDS(af)$cnv$copy_ratio; j <- rownames(A) %in% keep
    cn <- do.call(rbind, strsplit(colnames(A), "_"))
    av <- project(truth, cn[,1], as.numeric(cn[,2]), as.numeric(cn[,3]),
                  log2(colMeans(A[j,,drop=FALSE])))
  }
  d <- data.table(chr=truth$chr, gb=truth$start %/% GRID, t=truth$log2,
                  fg=fg, fl=fl, mv=mv, av=av)
  d <- d[, lapply(.SD, mean, na.rm=TRUE), by=.(chr,gb)][is.finite(t)]
  d[, care_id := S][]
}))
res[, subtype := fifelse(grepl("^Oligo", man[care_id, subtype]), "Oligo", "Astro")]

conf <- function(x, truth_dir) {
  tp <- sum(x$called & x$truth); fp <- sum(x$called & !x$truth); fn <- sum(!x$called & x$truth)
  prec <- tp/max(tp+fp,1); rec <- tp/max(tp+fn,1)
  data.table(direction=truth_dir, TP=tp, FP=fp, FN=fn,
             precision=round(prec,3), recall=round(rec,3),
             F1=round(2*prec*rec/max(prec+rec,1e-9),3))
}
cat("=== MintCNA, native discrete calls (bin called gain if >50% of nuclei say dup) ===\n\n")
for (st in c("Astro","Oligo")) {
  d <- res[subtype==st & is.finite(mv)]
  g <- conf(data.table(called=d$fg>0.5, truth=d$t> 0.1), "gain")
  l <- conf(data.table(called=d$fl>0.5, truth=d$t< -0.1), "loss")
  cat(sprintf("%s  (%d samples, %d bins)\n", st, uniqueN(d$care_id), nrow(d)))
  print(rbind(g,l)); cat("\n")
}

cat("=== amplitude actually recovered (threshold-free) ===\n")
cat("median called log2 / median truth log2, in bins the truth calls gain or loss\n\n")
cat(sprintf("%-7s %-6s %10s %12s %12s\n","subtype","dir","truth","MintCNA","AtaCNV"))
for (st in c("Astro","Oligo")) for (dir in c("gain","loss")) {
  d <- res[subtype==st & if (dir=="gain") t>0.1 else t< -0.1]
  tt <- median(d$t, na.rm=TRUE)
  cat(sprintf("%-7s %-6s %10.3f %6.3f (%3.0f%%) %6.3f (%3.0f%%)\n", st, dir, tt,
      median(d$mv,na.rm=TRUE), 100*median(d$mv,na.rm=TRUE)/tt,
      median(d$av,na.rm=TRUE), 100*median(d$av,na.rm=TRUE)/tt))
}
