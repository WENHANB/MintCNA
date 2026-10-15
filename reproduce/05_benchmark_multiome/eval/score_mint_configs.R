#!/usr/bin/env Rscript
# Score every MintCNA configuration on disk, through the SAME projection the
# headline benchmark table uses. project() is sourced out of pearson_vs_truth.R
# rather than reimplemented, so a config number and a benchmark number are
# comparable rather than merely similar.
suppressPackageStartupMessages({library(data.table)})
setwd("<MINTCNA_ROOT>/ATAC_RNA")
CHRS <- paste0("chr", 1:22); GRID <- 5e5

src <- readLines("<REPRODUCE_DIR>/05_benchmark_multiome/eval/pearson_vs_truth.R")
i0 <- grep("^project <- function", src); i1 <- grep("^# --- per-method pseudobulk", src)[1]
eval(parse(text = paste(src[i0:(i1-1)], collapse = "\n")))

res <- rbindlist(lapply(c("P111T3","P110T1","P104T1"), function(S) {
  truth <- fread(sprintf("Result/truth/%s_truth_100kb.tsv.gz", S))[chr %in% CHRS]
  keep  <- fread(sprintf("Data/cell_lists/%s_cells.tsv", S))[cell_type=="Malignant", barcode]
  files <- list.files(sprintf("Result/MintCNA/%s_multi", S),
                      sprintf("^%s_results_.*\\.RData$", S), full.names = TRUE)
  rbindlist(lapply(files, function(f) {
    tag <- sub(".*_results_(.*)\\.RData$", "\\1", f)
    # One config file came back as 1.1 MB of nul bytes -- a write that failed on
    # the parallel filesystem, not a logic error. A single such file must not
    # take down the scoring of the other 82, so name it and carry on.
    e <- new.env()
    if (inherits(try(load(f, e), silent = TRUE), "try-error")) {
      message("UNREADABLE, skipped: ", basename(f)); return(NULL)
    }
    Sm <- get("cna_state_matrix", e); rf <- get("ref", e)
    i <- colnames(Sm) %in% keep; if (!any(i)) return(NULL)
    neutral <- if (grepl("^st5", tag)) 3 else 2
    v <- project(truth, rf$chr, rf$start, rf$end,
                 log2(rowMeans(Sm[, i, drop=FALSE]) / neutral))
    d <- data.table(chr = truth$chr, gb = truth$start %/% GRID, v = v, t = truth$log2)
    d <- d[, .(v = mean(v, na.rm=TRUE), t = mean(t, na.rm=TRUE)), by = .(chr, gb)]
    d <- d[is.finite(v) & is.finite(t)]
    alt <- abs(d$t) > 0.1
    data.table(care_id = S, tag = tag, n_cell = sum(i), n_bin = nrow(d),
               r = cor(d$v, d$t),
               sign_alt = mean(sign(d$v[alt]) == sign(d$t[alt])),
               mae = mean(abs(d$v - d$t)),
               pct_alt = mean(Sm[, i] != neutral))
  }))
}))
fwrite(res, "Result/qc/mintcna_config_scores.tsv", sep = "\t")

# parse the tag into its factors so the grid can be read as a grid
res[, `:=`(h  = fifelse(grepl("^h",  tag), as.integer(sub("^h(\\d+).*","\\1",tag)), 5L),
           sm = fifelse(grepl("^sm", tag), as.integer(sub("^sm(\\d+).*","\\1",tag)), 0L))]
res[, base := sub("^(h\\d+|sm\\d+)", "", tag)]
g <- res[grepl("^a0", base) & !grepl("^st5|asym", tag)]
w <- dcast(g, h + sm + base ~ care_id, value.var = "r")
w[, mean := round(rowMeans(.SD, na.rm=TRUE), 3), .SDcols = c("P104T1","P110T1","P111T3")]
for (c in c("P104T1","P110T1","P111T3")) w[[c]] <- round(w[[c]], 3)
cat("=== Pearson, every (h, smoothing, alpha/eps) on the three samples ===\n\n")
print(w[order(-mean)], nrows = 60)
cat("\n=== same configs, sign agreement on truth-altered bins (mean of 3) ===\n\n")
s <- dcast(g, h + sm + base ~ care_id, value.var = "sign_alt")
s[, mean := round(rowMeans(.SD, na.rm=TRUE), 3), .SDcols = c("P104T1","P110T1","P111T3")]
print(s[order(-mean), .(h, sm, base, mean)], nrows = 60)
