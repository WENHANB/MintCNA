#!/usr/bin/env Rscript
# ATAC-only counterpart of a sample's multi config: the SAME parameters with
# mode switched to "dna", writing to <sample>_atac so the multi results are not
# overwritten. The AG-DCAE output is shared by symlink rather than retrained --
# it depends only on the ATAC L2R, so retraining would inject run-to-run variance
# into the one comparison the ablation is supposed to isolate.
#   Usage: make_config_atac.R <care_id>
args <- commandArgs(trailingOnly = TRUE); S <- args[1]
setwd("<MINTCNA_ROOT>/ATAC_RNA")
tpl <- readLines(sprintf("<REPRODUCE_DIR>/05_benchmark_multiome/mintcna/config_%s_multi.yaml", S))
tpl <- sub('^mode: "multi".*$', 'mode: "dna"                       # ATAC only: the ablation arm', tpl)
tpl <- sub(sprintf('  dir:    "Result/MintCNA/%s_multi"', S),
           sprintf('  dir:    "Result/MintCNA/%s_atac"',  S), tpl, fixed = TRUE)
stopifnot(any(grepl('^mode: "dna"', tpl)), any(grepl(sprintf("%s_atac", S), tpl)))
f <- sprintf("<REPRODUCE_DIR>/05_benchmark_multiome/mintcna/config_%s_atac.yaml", S)
writeLines(tpl, f)

src <- sprintf("Result/MintCNA/%s_multi", S); dst <- sprintf("Result/MintCNA/%s_atac", S)
dir.create(dst, showWarnings = FALSE, recursive = TRUE)
# Replace, do not skip. P111T3_atac already held byte-identical copies from an
# exploratory run in August; they happened to match, but a stale real file that
# did NOT match would silently make the ablation compare two different
# preprocessings. A symlink cannot go stale.
for (n in c("ae_dna.csv", "ref_dna.csv", "l2r_dna.csv", "latent_z.csv")) {
  from <- normalizePath(file.path(src, sprintf("%s_%s", S, n)), mustWork = TRUE)
  to   <- file.path(dst, sprintf("%s_%s", S, n))
  if (file.exists(to) && !nzchar(Sys.readlink(to))) {
    if (tools::md5sum(to) != tools::md5sum(from))
      stop("stale non-symlink differs from the multi run: ", to)
    file.remove(to)
  }
  if (!file.exists(to)) file.symlink(from, to)
}
cat("wrote", f, "and linked the shared AG-DCAE output into", dst, "\n")
