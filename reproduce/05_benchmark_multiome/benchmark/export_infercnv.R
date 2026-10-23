#!/usr/bin/env Rscript
# Export inferCNV's final matrix to a plain RDS so downstream code does not need
# the infercnv package.  Usage: export_infercnv.R [care_id ...]   (default: all)
#
# The result of infercnv::run() is an S4 object of class "infercnv". Reading it
# back requires the class definition, which lives only in the conda env
# (Tools/infercnv_env) -- the figure and metrics scripts run under the system R
# and would fail on `@` access. Exporting the two slots we actually use keeps
# that dependency confined to this script.
#
# Run it with the conda env's Rscript:
#   Tools/infercnv_env/bin/Rscript Code/benchmark/export_infercnv.R

suppressPackageStartupMessages({library(infercnv)})
ROOT <- "<MINTCNA_ROOT>/ATAC_RNA"
setwd(ROOT)
# Results are written per CONFIGURATION, not straight into the sample directory:
#   Result/benchmark/inferCNV/<tag>/<care_id>/
# This script used to point at the un-tagged parent, so after the 2026-08-18 move
# to package defaults it exported nothing for the default run and inferCNV simply
# vanished from the Pearson table -- a missing file makes the getter return NULL,
# which is indistinguishable from "not run yet".
TAG  <- Sys.getenv("INFERCNV_TAG", "default")
BASE <- file.path("Result/benchmark/inferCNV", TAG)
if (!dir.exists(BASE)) stop("no such inferCNV configuration: ", BASE)

args <- commandArgs(trailingOnly = TRUE)
samples <- if (length(args)) args else list.dirs(BASE, recursive = FALSE, full.names = FALSE)

for (s in samples) {
  d <- file.path(BASE, s)
  # `run.final.infercnv_obj` is what infercnv::run() writes when the whole
  # pipeline has finished, whatever options it ran under. The previous version
  # matched `2*_denoise*` instead, which only exists when denoise = TRUE -- so at
  # the package default (denoise = FALSE) it found nothing and skipped every
  # sample while reporting "not finished yet".
  #
  # It is still the completion marker: it appears only at the end of the run, so
  # a job still in flight is skipped rather than half-exported. Verified to be the
  # same object as 22_denoise.* for the cut01_denoiseT configuration.
  f <- file.path(d, "run.final.infercnv_obj")
  if (!file.exists(f)) { message(s, ": not finished yet, skipping"); next }
  o <- tryCatch(readRDS(f[1]), error = function(e) NULL)
  if (is.null(o)) { message(s, ": object not readable yet, skipping"); next }

  ref_idx <- unlist(o@reference_grouped_cell_indices)
  obs_idx <- unlist(o@observation_grouped_cell_indices)

  saveRDS(list(care_id     = s, tag = TAG,
               expr        = o@expr.data,             # genes x cells, centred on 1
               gene_order  = o@gene_order,            # chr / start / stop
               ref_cells   = colnames(o@expr.data)[ref_idx],
               obs_cells   = colnames(o@expr.data)[obs_idx],
               source_file = basename(f[1])),
          file.path(d, "infercnv_export.rds"))
  message(sprintf("%s: %d genes x %d cells (%d ref, %d obs) from %s",
                  s, nrow(o@expr.data), ncol(o@expr.data),
                  length(ref_idx), length(obs_idx), basename(f[1])))
}
