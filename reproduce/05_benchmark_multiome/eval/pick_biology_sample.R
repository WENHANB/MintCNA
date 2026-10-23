#!/usr/bin/env Rscript
# Which sample carries the biology panel?
#
# Three criteria, all stated in advance: the method must be working there
# (Pearson and F1), the subclone split must not be explained by sequencing depth
# (a split that depth alone predicts is a technical artefact, not a clone), and
# the minor clone must be big enough to compare against.
suppressPackageStartupMessages(library(data.table))
setwd("<MINTCNA_ROOT>/ATAC_RNA")
TAG <- "h15a005e020"
pe <- fread("Result/qc/pearson_vs_truth_own.tsv")[method == "MintCNA", .(care_id, pearson = r)]
f1 <- fread("Result/qc/f1_states.tsv")[method == "MintCNA", .(care_id, F1, precision, recall)]
man <- fread("Data/meta/sample_manifest.csv"); setkey(man, care_id)

q <- rbindlist(lapply(pe$care_id, function(s) {
  f <- sprintf("Result/MintCNA/%s_multi/%s_results_%s.RData", s, s, TAG)
  if (!file.exists(f)) return(NULL)
  e <- new.env(); load(f, e)
  sub <- base::get("subclones", e); Sm <- base::get("cna_state_matrix", e)
  cl <- fread(sprintf("Data/cell_lists/%s_cells.tsv", s))
  d <- cl$nFrags[match(colnames(Sm), cl$barcode)]
  tb <- sort(table(sub), decreasing = TRUE)
  if (length(tb) < 2) return(data.table(care_id = s, k = length(tb), minor_n = 0L,
                                        minor_frac = 0, clone_sep = NA_real_, depth_auc = NA_real_))
  o <- as.integer(names(tb)[1:2]); i1 <- sub == o[1]; i2 <- sub == o[2]
  x <- d[i1]; y <- d[i2]; x <- x[is.finite(x)]; y <- y[is.finite(y)]
  auc <- mean(outer(x, y, ">") + 0.5*outer(x, y, "=="))
  data.table(care_id = s, k = length(tb), minor_n = as.integer(min(tb)),
             minor_frac = min(tb)/sum(tb),
             clone_sep = mean(abs(rowMeans(Sm[, i1, drop=FALSE]) -
                                  rowMeans(Sm[, i2, drop=FALSE])), na.rm = TRUE),
             depth_auc = max(auc, 1 - auc))
}))
d <- merge(merge(pe, f1, by = "care_id"), q, by = "care_id")
d[, subtype := substr(man[care_id, subtype], 1, 5)]
d[, usable := depth_auc < 0.70 & minor_n >= 200 & pearson > 0.75 & F1 > 0.55]
setorder(d, -usable, -pearson)
cat("=== all 20, ranked ===\n\n")
print(d[, .(care_id, subtype, pearson = round(pearson,3), F1 = round(F1,3),
            k, minor_n, clone_sep = round(clone_sep,3),
            depth_auc = round(depth_auc,2), usable)], nrows = 25)
cat("\n=== candidates meeting all three criteria ===\n\n")
print(d[usable == TRUE, .(care_id, subtype, pearson = round(pearson,3), F1 = round(F1,3),
                          minor_n, clone_sep = round(clone_sep,3), depth_auc = round(depth_auc,2))])
