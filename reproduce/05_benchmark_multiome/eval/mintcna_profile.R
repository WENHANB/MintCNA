#!/usr/bin/env Rscript
# Turn a MintCNA result into the same pseudobulk track the six benchmark methods
# get, so it can sit in the same figure and the same Pearson table.
#   Usage: Rscript mintcna_profile.R <care_id> [tag ...]
#
# The convention is epiAneufinder's, which the whole comparison already uses:
# a discrete per-cell caller becomes continuous by AVERAGING THE STATE CODES
# ACROSS CELLS, and the average is converted once to log2 through the method's
# own neutral code. Nothing about MintCNA's amplitude is invented.
#   state 1 = del, 2 = nor, 3 = dup   ->  log2(mean_over_cells / 2)

suppressPackageStartupMessages({library(data.table)})
setwd("<MINTCNA_ROOT>/ATAC_RNA")
args <- commandArgs(trailingOnly = TRUE)
S <- if (length(args)) args[1] else "P111T3"
TAGS <- if (length(args) > 1) args[-1] else NULL
D <- sprintf("Result/MintCNA/%s_multi", S)
GRID <- 5e5; CHRS <- paste0("chr", 1:22)

files <- if (is.null(TAGS)) list.files(D, "_results_.*\\.RData$", full.names = TRUE) else
  file.path(D, sprintf("%s_results_%s.RData", S, TAGS))
files <- files[file.exists(files)]
if (!length(files)) stop("no MintCNA result files in ", D)

tr <- fread(sprintf("Result/truth/%s_truth_100kb.tsv.gz", S))[chr %in% CHRS]
tr <- tr[, .(truth = mean(log2, na.rm = TRUE)), by = .(chr, gb = start %/% GRID)]

out <- rbindlist(lapply(files, function(f) {
  tag <- sub(".*_results_(.*)\\.RData$", "\\1", f)
  e <- new.env(); load(f, e)
  Sm <- get("cna_state_matrix", e); rf <- get("ref", e)
  # neutral code is 2 for st = 3 (del/nor/dup) and 3 for st = 5
  neutral <- if (grepl("^st5", tag)) 3 else 2
  v <- log2(rowMeans(Sm) / neutral)                 # states -> continuous, then one log2
  d <- data.table(chr = rf$chr, gb = rf$start %/% GRID, v = v)[chr %in% CHRS]
  d <- d[, .(v = mean(v, na.rm = TRUE)), by = .(chr, gb)]
  m <- merge(d, tr, by = c("chr", "gb"))
  ok <- is.finite(m$v) & is.finite(m$truth)
  cat(sprintf("%-12s  bins %4d   coverage %.3f   Pearson %.3f   del %.1f%%  dup %.1f%%\n",
              tag, sum(ok), mean(is.finite(d$v)), cor(m$v[ok], m$truth[ok]),
              100*mean(Sm < neutral), 100*mean(Sm > neutral)))
  data.table(care_id = S, method = sprintf("MintCNA (%s)", tag),
             chr = d$chr, gb = d$gb, log2 = d$v)
}))
fwrite(out, sprintf("Result/qc/profiles/%s_mintcna_500kb.tsv.gz", S), sep = "\t")
cat("\nwrote", sprintf("Result/qc/profiles/%s_mintcna_500kb.tsv.gz", S), "\n")
