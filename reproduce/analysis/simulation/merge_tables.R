#!/usr/bin/env Rscript
# Merge the per-setting evaluation tables of the nine simulation settings into one table each,
# adding the coverage and noise level parsed from the setting name (0.1x and 0 when absent).
#   Rscript reproduce/analysis/simulation/merge_tables.R
# Input : Simulation/results/tables/<setting>/*.csv  (evaluate.R)
# Output: Simulation/results/tables/merged/{summary_del,summary_dup,summary_mss,summary_one,
#         bb_summary,bb_percell,RMSE_percell}.csv  (the source data of Figs 2-3)
# Run from the root of the MintCNA repository with MINTCNA_DATA set (see reproduce/README.md).
suppressPackageStartupMessages(library(data.table))
ROOT <- Sys.getenv("MINTCNA_DATA"); if (!nzchar(ROOT)) stop("set MINTCNA_DATA (see reproduce/README.md)")
TAB_ROOT <- file.path(ROOT, "Simulation/results/tables")
OUT_DIR  <- file.path(TAB_ROOT, "merged"); dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

settings <- c("c3_N100_100kb_l1e4", "c3_N100_100kb_l1e4_c0.05", "c3_N100_100kb_l1e4_c0.05_n0.05",
              "c3_N100_100kb_l1e4_c0.05_n0.1", "c3_N100_100kb_l1e4_c0.2", "c3_N100_100kb_l1e4_c0.2_n0.05",
              "c3_N100_100kb_l1e4_c0.2_n0.1", "c3_N100_100kb_l1e4_n0.05", "c3_N100_100kb_l1e4_n0.1")
parse_params <- function(s) {
  cov <- regmatches(s, regexpr("_c(0\\.[0-9]+)", s)); noi <- regmatches(s, regexpr("_n(0\\.[0-9]+)", s))
  list(coverage = if (length(cov) == 1) as.numeric(sub("^_c", "", cov)) else 0.1,
       noise    = if (length(noi) == 1) as.numeric(sub("^_n", "", noi)) else 0)
}

for (tbl in c("summary_del", "summary_dup", "summary_mss", "summary_one", "bb_summary", "bb_percell", "RMSE_percell")) {
  merged <- rbindlist(lapply(settings, function(s) {
    p <- parse_params(s)
    dt <- fread(file.path(TAB_ROOT, s, paste0(tbl, ".csv")))
    dt[, `:=`(coverage = p$coverage, noise = p$noise)]
  }), use.names = TRUE, fill = TRUE)
  setcolorder(merged, c("coverage", "noise", setdiff(names(merged), c("coverage", "noise"))))
  fwrite(merged, file.path(OUT_DIR, paste0(tbl, ".csv")))
}
message("wrote ", OUT_DIR)
