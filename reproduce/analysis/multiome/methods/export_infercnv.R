#!/usr/bin/env Rscript
# Export inferCNV's final expression matrix and gene order to a plain RDS, so that the evaluation scripts do
# not need the infercnv package (the run.final.infercnv_obj file is an S4 object of its class).
#   Rscript reproduce/analysis/multiome/methods/export_infercnv.R [care_id ...]      (default: every sample found)
#   INFERCNV_TAG selects the configuration folder (default "default").
# Input : ATAC_RNA/Result/benchmark/inferCNV/<tag>/<care_id>/run.final.infercnv_obj
# Output: ATAC_RNA/Result/benchmark/inferCNV/<tag>/<care_id>/infercnv_export.rds
# Run from the root of the MintCNA repository with MINTCNA_DATA set (see reproduce/README.md).
suppressPackageStartupMessages(library(infercnv))
ROOT <- Sys.getenv("MINTCNA_DATA"); if (!nzchar(ROOT)) stop("set MINTCNA_DATA (see reproduce/README.md)")
setwd(file.path(ROOT, "ATAC_RNA"))
TAG  <- Sys.getenv("INFERCNV_TAG", "default")
BASE <- file.path("Result/benchmark/inferCNV", TAG)
if (!dir.exists(BASE)) stop("no such inferCNV configuration: ", BASE)
args <- commandArgs(trailingOnly = TRUE)
samples <- if (length(args)) args else list.dirs(BASE, recursive = FALSE, full.names = FALSE)

for (s in samples) {
  f <- file.path(BASE, s, "run.final.infercnv_obj")        # written when the whole run has finished
  if (!file.exists(f)) { message(s, ": not finished, skipped"); next }
  o <- readRDS(f)
  ref_idx <- unlist(o@reference_grouped_cell_indices); obs_idx <- unlist(o@observation_grouped_cell_indices)
  saveRDS(list(care_id = s, tag = TAG, expr = o@expr.data, gene_order = o@gene_order,
               ref_cells = colnames(o@expr.data)[ref_idx], obs_cells = colnames(o@expr.data)[obs_idx]),
          file.path(BASE, s, "infercnv_export.rds"))
  message(sprintf("%s: %d genes x %d cells (%d reference, %d observed)", s, nrow(o@expr.data), ncol(o@expr.data),
                  length(ref_idx), length(obs_idx)))
}
