#!/usr/bin/env Rscript
# SCOPE (1.16.0) on one simulation setting: CODEX2 normalisation with the normal cells, ploidy
# initialisation, cross-sample CBS segmentation, integer copy numbers.
#   Rscript reproduce/analysis/simulation/methods/SCOPE.R <setting>
# Input : <setting>/DNA_RC.RData (Y, gr), <setting>/cell_types.tsv
# Output: <setting>/Calling/SCOPE.RData  (cna_scope: bins x cells integer copy number, ref_scope)
# Run from the root of the MintCNA repository with MINTCNA_DATA set (see reproduce/README.md).
suppressPackageStartupMessages({
  library(GenomicRanges); library(SCOPE); library(WGSmapp); library(BSgenome.Hsapiens.UCSC.hg19)
  library(data.table); library(parallel)
})
ROOT <- Sys.getenv("MINTCNA_DATA"); if (!nzchar(ROOT)) stop("set MINTCNA_DATA (see reproduce/README.md)")
args <- commandArgs(trailingOnly = TRUE)
BASE <- if (length(args) >= 1) args[1] else "c3_N100_100kb_l1e4"
CNA_DIR <- file.path(ROOT, "Simulation/cnasim", BASE)
load(file.path(CNA_DIR, "DNA_RC.RData"))
label <- fread(file.path(CNA_DIR, "cell_types.tsv"))
nc_samples <- label[normal == "normal", root]

seqlevelsStyle(gr) <- "UCSC"
ranges(gr) <- IRanges(start = start(gr) + 1, end = end(gr))
gr$gc <- as.numeric(mcols(gr)$gc) * 100           # SCOPE expects GC in percent

norm_index <- match(nc_samples, colnames(Y)); norm_index <- norm_index[!is.na(norm_index)]
if (length(norm_index) < 3) warning("fewer than 3 normal cells found in Y")

gc_qc <- as.numeric(mcols(gr)$gc)
norm0 <- normalize_codex2_ns_noK(Y_qc = Y, gc_qc = gc_qc, norm_index = norm_index)
ploidy <- initialize_ploidy(Y = Y, Yhat = norm0$Yhat, ref = gr)
normObj <- normalize_scope_foreach(Y_qc = Y, gc_qc = gc_qc, K = 1, ploidyInt = ploidy,
                                   norm_index = norm_index, T = 1:5, beta0 = norm0$beta.hat, nCores = 2)
Yhat <- normObj$Yhat[[which.max(normObj$BIC)]]

ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
chr_levels <- paste0("chr", 1:22)
seglist <- mclapply(chr_levels, function(ch)
  segment_CBScs(Y = Y, Yhat = Yhat, sampname = colnames(Y), ref = gr, chr = ch,
                mode = "integer", max.ns = 1),
  mc.cores = ncores, mc.preschedule = FALSE)
iCN <- do.call(rbind, lapply(seglist, function(z) z[["iCN"]]))

cna_scope <- iCN; ref_scope <- as.data.frame(gr)
dir.create(file.path(CNA_DIR, "Calling"), showWarnings = FALSE, recursive = TRUE)
save(cna_scope, ref_scope, file = file.path(CNA_DIR, "Calling/SCOPE.RData"))
