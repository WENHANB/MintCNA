#!/usr/bin/env Rscript
# Pseudobulk copy-number track of a MintCNA result on the 500 kb evaluation grid, in log2 ratio, the same
# representation the compared methods receive: the state codes (1 loss, 2 neutral, 3 gain) are averaged over
# the malignant nuclei and divided by the neutral code.
#   Rscript reproduce/analysis/multiome/evaluation/mintcna_profile.R <care_id> [tag ...]   (default: every result file)
# Input : ATAC_RNA/Result/MintCNA/<care_id>_multi/<care_id>_results_<tag>.RData, ATAC_RNA/Result/truth/<care_id>_truth_100kb.tsv.gz
# Output: ATAC_RNA/Result/qc/profiles/<care_id>_mintcna_500kb.tsv.gz
# Run from the root of the MintCNA repository with MINTCNA_DATA set (see reproduce/README.md).
suppressPackageStartupMessages(library(data.table))
ROOT <- Sys.getenv("MINTCNA_DATA"); if (!nzchar(ROOT)) stop("set MINTCNA_DATA (see reproduce/README.md)")
setwd(file.path(ROOT, "ATAC_RNA"))
args <- commandArgs(trailingOnly = TRUE)
S <- if (length(args)) args[1] else "P111T2"
TAGS <- if (length(args) > 1) args[-1] else NULL
D <- sprintf("Result/MintCNA/%s_multi", S)
GRID <- 5e5; CHRS <- paste0("chr", 1:22)

files <- if (is.null(TAGS)) list.files(D, "_results_.*\\.RData$", full.names = TRUE) else file.path(D, sprintf("%s_results_%s.RData", S, TAGS))
files <- files[file.exists(files)]
if (!length(files)) stop("no MintCNA result files in ", D)
tr <- fread(sprintf("Result/truth/%s_truth_100kb.tsv.gz", S))[chr %in% CHRS]
tr <- tr[, .(truth = mean(log2, na.rm = TRUE)), by = .(chr, gb = start %/% GRID)]

out <- rbindlist(lapply(files, function(f) {
  tag <- sub(".*_results_(.*)\\.RData$", "\\1", f)
  e <- new.env(); load(f, e)
  Sm <- get("cna_state_matrix", e); rf <- get("ref", e)
  neutral <- if (grepl("^st5", tag)) 3 else 2
  d <- data.table(chr = rf$chr, gb = rf$start %/% GRID, v = log2(rowMeans(Sm) / neutral))[chr %in% CHRS]
  d <- d[, .(v = mean(v, na.rm = TRUE)), by = .(chr, gb)]
  m <- merge(d, tr, by = c("chr", "gb")); ok <- is.finite(m$v) & is.finite(m$truth)
  cat(sprintf("%-12s bins %4d  Pearson %.3f\n", tag, sum(ok), cor(m$v[ok], m$truth[ok])))
  data.table(care_id = S, method = sprintf("MintCNA (%s)", tag), chr = d$chr, gb = d$gb, log2 = d$v)
}))
dir.create("Result/qc/profiles", showWarnings = FALSE, recursive = TRUE)
fwrite(out, sprintf("Result/qc/profiles/%s_mintcna_500kb.tsv.gz", S), sep = "\t")
