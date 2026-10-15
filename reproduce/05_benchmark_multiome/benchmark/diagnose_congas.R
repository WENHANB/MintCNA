#!/usr/bin/env Rscript
# Why do all four CONGAS+ fits die with
#   ValueError: Expected parameter concentration (Tensor of shape (41,)) of
#   distribution Gamma(...) to satisfy GreaterThan(lower_bound=0.0)
#
# The Gamma priors come from auto_config_run(). Some entry is <= 0 or NaN.
#
# The hyperparameters cannot be inspected from the saved checkpoint: they are
# reticulate handles to Python tensors, and an external pointer does not survive
# serialisation -- reading them back gives "attempt to apply non-function". So
# auto_config_run() has to be re-run in the same session as the inspection.
#
#   Usage: Rscript diagnose_congas.R <care_id>

args <- commandArgs(trailingOnly = TRUE)
CARE <- if (length(args)) args[1] else "P106T2"

suppressPackageStartupMessages({library(Rcongas); library(reticulate); library(dplyr)})
ROOT <- "<MINTCNA_ROOT>/ATAC_RNA"
setwd(ROOT)

r <- readRDS(file.path("Result/benchmark/CONGASplus", CARE, "congas_input.rds"))
x <- r$object                       # pure R tibbles, survives the round-trip
cat("object segments:", nrow(Rcongas::get_input(x, what = "segmentation")), "\n")
print(Rcongas::stat(x)[c("nsegments", "ncells_RNA", "ncells_ATAC")])

hyper <- Rcongas::auto_config_run(x, 1:4)

tonum <- function(v) {
  if (inherits(v, "python.builtin.object")) as.numeric(py_to_r(v$detach()$numpy()))
  else as.numeric(v)
}

cat("\n================ hyperparameter audit ================\n")
segs <- Rcongas::get_input(x, what = "segmentation")
for (n in names(hyper)) {
  v <- tryCatch(tonum(hyper[[n]]), error = function(e) NULL)
  if (is.null(v) || !length(v)) next
  bad <- which(!is.finite(v) | v <= 0)
  cat(sprintf("%-20s len=%-4d min=%-12.5g max=%-12.5g  BAD=%d\n",
              n, length(v), min(v, na.rm = TRUE), max(v, na.rm = TRUE), length(bad)))
  if (length(bad)) {
    cat("    values:", paste(signif(v[head(bad, 8)], 4), collapse = ", "), "\n")
    if (length(v) == nrow(segs)) {
      cat("    segments:\n")
      print(segs[head(bad, 8), c("chr", "from", "to")])
    }
  }
}

# The per-segment counts are what auto_config_run derives the priors from, so
# look for the segments that carry no signal in one modality.
cat("\n================ per-segment counts ================\n")
for (mod in c("RNA", "ATAC")) {
  d <- Rcongas::get_input(x, what = "data") %>% dplyr::filter(modality == mod)
  s <- d %>% group_by(segment_id) %>%
       summarise(total = sum(value), nz_cells = sum(value > 0), n = n(), .groups = "drop") %>%
       arrange(total)
  cat("--", mod, ": segments with the lowest totals --\n"); print(head(s, 8))
  cat("   segments with total == 0:", sum(s$total == 0), "\n")
}
