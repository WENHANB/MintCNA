#!/usr/bin/env Rscript
# Pseudobulk Pearson correlation of every method with the bulk WGS truth, per sample (Fig. 5a).
# Each method's output is averaged over the malignant nuclei it retained, in its native units, converted
# once to a log2 ratio through the method's own neutral value, projected onto the truth's 100 kb grid
# (midpoint containment) and scored on the 500 kb evaluation grid against the truth's log2 ratio.
# The configurations scored are the methods' published defaults and MintCNA's reported setting
# (h 15, alpha 0.05, eps 0.20 for ATAC + RNA; eps 0.25 for the ATAC-only arm).
#   Rscript reproduce/analysis/multiome/evaluation/pearson_vs_truth.R [care_id ...]   (default: the 20 samples)
# Input : ATAC_RNA/Result/truth/, Result/benchmark/<method>/, Result/MintCNA/, Data/cell_lists/, Data/meta/sample_manifest.csv
# Output: ATAC_RNA/Result/qc/pearson_vs_truth_own.tsv, Result/qc/profiles/<care_id>_profiles_100kb.tsv.gz
# Run from the root of the MintCNA repository with MINTCNA_DATA set (see reproduce/README.md).
suppressPackageStartupMessages(library(data.table))
ROOT <- Sys.getenv("MINTCNA_DATA"); if (!nzchar(ROOT)) stop("set MINTCNA_DATA (see reproduce/README.md)")
setwd(file.path(ROOT, "ATAC_RNA"))
CHRS <- paste0("chr", 1:22)
GRID <- 5e5
EPI_TAG       <- Sys.getenv("EPI_TAG",       "w100kb_mf20000")
RIDDLER_DIR   <- Sys.getenv("RIDDLER_DIR",   "1MB")
INFERCNV_TAG  <- Sys.getenv("INFERCNV_TAG",  "default")
MINT_TAG      <- Sys.getenv("MINT_TAG",      "h15a005e020")
MINT_ATAC_TAG <- Sys.getenv("MINT_ATAC_TAG", "h15a005e025")
args <- commandArgs(trailingOnly = TRUE)
SAMPLES <- if (length(args)) args else c("P59T2","P59T3","P104T1","P104T2","P106T1","P106T2","P100T1","P100T3","P103T1","P103T2",
                                         "P105T1","P105T3","P108T1","P108T2","P109T1","P109T2","P110T1","P110T2","P111T2","P111T3")

# value of the method bin containing each truth bin's midpoint (NA where no bin does)
project <- function(truth, chr, from, to, val) {
  ok <- is.finite(val) & chr %in% CHRS
  if (!any(ok)) return(rep(NA_real_, nrow(truth)))
  d <- data.table(chr = as.character(chr[ok]), from = as.numeric(from[ok]), bin_from = as.numeric(from[ok]),
                  bin_to = as.numeric(to[ok]), val = as.numeric(val[ok]))
  setkey(d, chr, from)
  q <- data.table(chr = as.character(truth$chr), from = (as.numeric(truth$start) + as.numeric(truth$end)) / 2, i = seq_len(nrow(truth)))
  setkey(q, chr, from)
  j <- d[q, roll = TRUE]
  j[is.na(bin_to) | from < bin_from | from > bin_to, val := NA_real_]
  out <- rep(NA_real_, nrow(truth)); out[j$i] <- j$val
  out
}

# ---- one pseudobulk log2-ratio track per method ---------------------------------------------------
get_epi <- function(s, truth, keep) {          # states 0/1/2 -> copy number 1/2/3, mean over cells
  f <- file.path("Result/benchmark/epiAneufinder", EPI_TAG, s, "epiAneufinder_results/results_table.tsv")
  if (!file.exists(f)) return(NULL)
  ef <- fread(f); cc <- grep("^cell-", names(ef), value = TRUE); use <- cc[sub("^cell-", "", cc) %in% keep]
  if (!length(use)) return(NULL)
  project(truth, ef$seq, ef$start, ef$end, log2(rowMeans(as.matrix(ef[, ..use]) + 1) / 2))
}
get_atacnv <- function(s, truth, keep) {       # copy ratio, neutral = 1
  f <- file.path("Result/benchmark/AtaCNV", s, "atacnv_result.rds")
  if (!file.exists(f)) return(NULL)
  A <- readRDS(f)$cnv$copy_ratio; i <- rownames(A) %in% keep; if (!any(i)) return(NULL)
  cn <- do.call(rbind, strsplit(colnames(A), "_"))
  project(truth, cn[, 1], as.numeric(cn[, 2]), as.numeric(cn[, 3]), log2(colMeans(A[i, , drop = FALSE])))
}
get_copyscat <- function(s, truth, keep, arms) {   # arm-level pseudo copy number, neutral = 2; untested arms absent, uncalled = 2
  f <- file.path("Result/benchmark/CopyscAT", s, paste0(s, "clean_cnv_cnv_scores.csv"))
  if (!file.exists(f)) return(NULL)
  cs <- read.csv(f, row.names = 1, check.names = FALSE); i <- rownames(cs) %in% keep; if (!any(i)) return(NULL)
  tested <- as.character(readRDS(file.path("Result/benchmark/CopyscAT", s, "copyscat_result.rds"))$collapse[[1]])
  m <- setNames(rep(2, length(tested)), tested); cm <- colMeans(cs[i, , drop = FALSE]); m[names(cm)] <- cm
  a <- arms[label %in% names(m)]
  project(truth, a$chr, a$from, a$to, log2(m[a$label] / 2))
}
get_riddler <- function(s, truth, keep) {      # multiplier of the cell's baseline, neutral = 1
  f <- file.path("Result/benchmark/RIDDLER", RIDDLER_DIR, s, "riddler_result.rds")
  if (!file.exists(f)) return(NULL)
  rd <- readRDS(f)$riddler
  i <- sub("\\.(\\d+)$", "-\\1", colnames(rd$CNV)) %in% keep; if (!any(i)) return(NULL)
  w <- rd$window; mi <- match(rownames(rd$CNV), as.character(w[[4]]))
  project(truth, as.character(w[[1]])[mi], w[[2]][mi], w[[3]][mi], log2(rowMeans(rd$CNV[, i, drop = FALSE])))
}
get_congas <- function(s, truth, keep) {       # integer copy number per clone and segment, cell-weighted mean, neutral = 2
  f <- file.path("Result/benchmark/CONGASplus", s, "congas_result.rds")
  if (!file.exists(f)) return(NULL)
  bf <- readRDS(f)$fit$best_fit
  cna <- as.data.table(if (!is.null(bf$CNA_real)) bf$CNA_real else bf$CNA)
  ass <- as.data.table(bf$cluster_assignments); ass[, cell := sub("-(RNA|ATAC)$", "", cell)]
  ass <- unique(ass[, .(cell, cluster)])[cell %in% keep]; if (!nrow(ass)) return(NULL)
  p <- merge(cna, ass[, .(n = .N), by = cluster], by = "cluster")[, .(cn = sum(value * n) / sum(n)), by = segment_id]
  p[, c("c", "f", "t") := tstrsplit(segment_id, ":", fixed = TRUE)]
  project(truth, p$c, as.numeric(p$f), as.numeric(p$t), log2(p$cn / 2))
}
get_infercnv <- function(s, truth, keep) {     # residual expression, neutral = 1
  f <- file.path("Result/benchmark/inferCNV", INFERCNV_TAG, s, "infercnv_export.rds")
  if (!file.exists(f)) return(NULL)
  iv <- readRDS(f); i <- colnames(iv$expr) %in% keep; if (!any(i)) return(NULL)
  go <- iv$gene_order
  project(truth, as.character(go$chr), go$start, go$stop, log2(rowMeans(iv$expr[, i, drop = FALSE])))
}
mint_track <- function(f, truth, keep) {       # states averaged over cells, divided by the neutral code
  if (!file.exists(f)) return(NULL)
  e <- new.env(); load(f, e)
  Sm <- get("cna_state_matrix", e); rf <- get("ref", e)
  i <- colnames(Sm) %in% keep; if (!any(i)) return(NULL)
  neutral <- if (grepl("st5", basename(f))) 3 else 2
  project(truth, rf$chr, rf$start, rf$end, log2(rowMeans(Sm[, i, drop = FALSE]) / neutral))
}

arms <- as.data.table(readRDS("Result/benchmark/AtaCNV/P105T1/atacnv_result.rds")$bins)
arms <- arms[chr %in% CHRS & arm %in% c("p", "q"), .(from = min(start), to = max(end)), by = .(chr, arm)]
arms[, label := paste0(chr, arm)]
ATAC_METHODS <- c("epiAneufinder", "AtaCNV", "Copy-scAT", "RIDDLER", "CONGAS+", "MintCNA", "MintCNA(ATAC)")

res <- rbindlist(lapply(SAMPLES, function(s) {
  tf <- file.path("Result/truth", paste0(s, "_truth_100kb.tsv.gz"))
  if (!file.exists(tf)) return(NULL)
  truth <- fread(tf)[chr %in% CHRS]
  keep  <- fread(file.path("Data/cell_lists", paste0(s, "_cells.tsv")))[cell_type == "Malignant", barcode]
  P <- list(epiAneufinder = get_epi(s, truth, keep), AtaCNV = get_atacnv(s, truth, keep),
            `Copy-scAT` = get_copyscat(s, truth, keep, arms), RIDDLER = get_riddler(s, truth, keep),
            `CONGAS+` = get_congas(s, truth, keep), inferCNV = get_infercnv(s, truth, keep),
            MintCNA = mint_track(sprintf("Result/MintCNA/%s_multi/%s_results_%s.RData", s, s, MINT_TAG), truth, keep),
            `MintCNA(ATAC)` = mint_track(sprintf("Result/MintCNA/%s_atac/%s_results_%s.RData", s, s, MINT_ATAC_TAG), truth, keep))

  # the projected tracks, with the truth, for the genome-wide profile figures
  prof <- rbindlist(lapply(names(P), function(nm) if (is.null(P[[nm]])) NULL else
    data.table(care_id = s, method = nm, chr = truth$chr, start = truth$start, end = truth$end, log2 = P[[nm]])))
  if (nrow(prof)) {
    prof <- rbind(prof, data.table(care_id = s, method = "bulk WGS", chr = truth$chr, start = truth$start, end = truth$end, log2 = truth$log2))
    dir.create("Result/qc/profiles", showWarnings = FALSE, recursive = TRUE)
    fwrite(prof, sprintf("Result/qc/profiles/%s_profiles_100kb.tsv.gz", s), sep = "\t")
  }

  # score on the 500 kb grid: on the bins each method reports (r) and on the bins common to all scATAC methods (r_common)
  agg <- rbindlist(c(lapply(names(P), function(nm) if (is.null(P[[nm]])) NULL else
                       data.table(method = nm, chr = truth$chr, gb = truth$start %/% GRID, v = P[[nm]])),
                     list(data.table(method = "truth", chr = truth$chr, gb = truth$start %/% GRID, v = truth$log2))))
  agg <- agg[, .(v = mean(v, na.rm = TRUE)), by = .(method, chr, gb)]; agg[is.nan(v), v := NA_real_]
  W <- dcast(agg, chr + gb ~ method, value.var = "v")
  tv <- W[["truth"]]; mm <- setdiff(names(W), c("chr", "gb", "truth")); ok_atac <- intersect(ATAC_METHODS, mm)
  kA <- if (length(ok_atac)) Reduce(`&`, c(lapply(ok_atac, function(m) is.finite(W[[m]])), list(is.finite(tv)))) else rep(FALSE, nrow(W))
  rbindlist(lapply(mm, function(nm) {
    v <- W[[nm]]; o <- is.finite(v) & is.finite(tv)
    data.table(care_id = s, method = nm, cellset = "own", n_cell = length(keep), n_malignant = length(keep),
               n_bin = sum(o), n_common = sum(kA), grid_kb = GRID / 1e3, coverage = mean(is.finite(v)),
               r = if (sum(o) > 50 && sd(v[o]) > 0) cor(v[o], tv[o]) else NA_real_,
               r_common = if (nm %in% ok_atac && sum(kA) > 50 && sd(v[kA]) > 0) cor(v[kA], tv[kA]) else NA_real_)
  }))
}))

dir.create("Result/qc", showWarnings = FALSE, recursive = TRUE)
fwrite(res, "Result/qc/pearson_vs_truth_own.tsv", sep = "\t")
cat("wrote Result/qc/pearson_vs_truth_own.tsv\n\n")
print(res[!is.na(r), .(n = .N, median = round(median(r), 3), mean = round(mean(r), 3), min = round(min(r), 3), max = round(max(r), 3)),
          by = method][order(-median)])
