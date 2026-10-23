#!/usr/bin/env Rscript
# Write a per-sample MintCNA config by substituting the sample name into the
# P111T3 one, so every sample runs the identical parameter set and a change to
# the template cannot silently apply to only some of them.
#   Usage: Rscript make_config.R <care_id>
args <- commandArgs(trailingOnly = TRUE)
S <- args[1]
setwd("<MINTCNA_ROOT>/ATAC_RNA")
tpl <- readLines("<REPRODUCE_DIR>/05_benchmark_multiome/mintcna/config_P111T3_multi.yaml")
out <- gsub("P111T3", S, tpl, fixed = TRUE)
f <- sprintf("<REPRODUCE_DIR>/05_benchmark_multiome/mintcna/config_%s_multi.yaml", S)
writeLines(out, f)
cat("wrote", f, "\n")
