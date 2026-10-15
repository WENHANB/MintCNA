#!/usr/bin/env Rscript
# tau is the ONLY quantity call_cells() takes from the null, so two settings with
# the same tau produce identical calls. Does one alpha recover the old tau on
# every sample, given the bug shifted tau by 7.8-17.0% depending on the sample?
suppressPackageStartupMessages(library(data.table))
setwd("<MINTCNA_ROOT>/ATAC_RNA")
o <- fread("Data/meta/oldtau.txt", header = FALSE, col.names = c("care_id","old_a10"))
r <- rbindlist(lapply(o$care_id, function(s) {
  f <- sprintf("Result/MintCNA/%s_multi/%s_scan_h15.rds", s, s)
  if (!file.exists(f)) return(NULL)
  n <- readRDS(f)$res_null
  g <- function(a) n$q[which.min(abs(n$probs - (1 - a)))]
  data.table(care_id = s, new_a10 = g(0.10), new_a05 = g(0.05))
}))
d <- merge(r, o, by = "care_id")[, ratio := new_a05 / old_a10][order(ratio)]
cat("Does alpha = 0.05 on the FIXED null reproduce the old alpha = 0.10 tau?\n\n")
print(d[, .(care_id, old_tau_a10 = round(old_a10,4), new_tau_a05 = round(new_a05,4),
            ratio = round(ratio,3), new_tau_a10 = round(new_a10,4))], nrows=25)
cat(sprintf("\nratio: median %.3f   range %.3f - %.3f\n", median(d$ratio), min(d$ratio), max(d$ratio)))
cat(sprintf("within 5%% of the old tau on %d of %d samples\n", sum(abs(d$ratio-1)<=0.05), nrow(d)))
