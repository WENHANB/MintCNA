#!/usr/bin/env Rscript
# What EVERY method says at a named locus, beside the bulk WGS.
#
# The point is to separate two explanations for a miss: the method is wrong, or
# the locus is below what this data and this evaluation grid can resolve. If no
# method sees an event, the second is the honest reading; if only one method
# misses it, the first is.
#
# Everything is read from Result/qc/profiles/<S>_profiles_100kb.tsv.gz, the
# projected pseudobulk tracks written by pearson_vs_truth.R, so the numbers are
# on the same grid and the same projection as the benchmark table.
suppressPackageStartupMessages(library(data.table))
setwd("<MINTCNA_ROOT>/ATAC_RNA")
a <- commandArgs(TRUE)
LOCI <- data.table(
  gene = c("CDKN2A", "CCND2", "1p", "19q"),
  chr  = c("chr9", "chr12", "chr1", "chr19"),
  from = c(21967752, 4273761,        1,  26200000),   # hg38; 1p and 19q are arms
  to   = c(21995301, 4305353, 123400000,  58600000))
WANT <- if (length(a)) a else LOCI$gene
S <- sub("_profiles_100kb.tsv.gz$", "", basename(
       list.files("Result/qc/profiles", "_profiles_100kb.tsv.gz$")))
man <- fread("Data/meta/sample_manifest.csv"); setkey(man, care_id)

out <- rbindlist(lapply(S, function(s) {
  pr <- fread(sprintf("Result/qc/profiles/%s_profiles_100kb.tsv.gz", s))
  rbindlist(lapply(WANT, function(g) {
    L <- LOCI[gene == g]
    x <- pr[chr == L$chr & end >= L$from & start <= L$to & is.finite(log2)]
    if (!nrow(x)) return(NULL)
    x[, .(care_id = s, gene = g, subtype = substr(man[s, subtype], 1, 5),
          tp = man[s, timepoint], n_bin = .N, med = median(log2)), by = method]
  }))
}))
setnames(out, "med", "log2")
fwrite(out, "Result/qc/locus_all_methods.tsv", sep = "\t")

for (g in WANT) {
  d <- out[gene == g]
  w <- dcast(d, care_id + subtype + tp ~ method, value.var = "log2")
  mm <- setdiff(names(w), c("care_id","subtype","tp"))
  for (c in mm) w[[c]] <- round(w[[c]], 3)
  setcolorder(w, c("care_id","subtype","tp","bulk WGS"))
  setorder(w, `bulk WGS`)
  cat(sprintf("\n=== %s — median pseudobulk log2 over the locus ===\n\n", g))
  print(w, nrows = 25)
  # does each method track the truth across samples?
  tv <- w[["bulk WGS"]]
  cat("\n  Spearman with the bulk WGS across samples:\n")
  for (m in setdiff(mm, "bulk WGS"))
    cat(sprintf("    %-14s %6.2f\n", m,
                suppressWarnings(cor(tv, w[[m]], method = "spearman", use = "complete.obs"))))
}
