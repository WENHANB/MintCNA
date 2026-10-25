#!/usr/bin/env Rscript
# Precision / recall / F1 per copy-number state (CN 0, 1, 3, 4, 5) of the three MintCNA modes
# (DNA + RNA, DNA only, RNA only) at coverage 0.05x, 0.1x and 0.2x: the data of Supplementary Fig. 5b-d.
# An event is a run of bins with the same state; a truth event of state k is recovered when a called
# event of the same direction overlaps it; a called event is a false positive when it overlaps no truth
# event of its direction. Counts are pooled over cells. The RNA-only mode is scored on the bins with
# RNA coverage only.
#   Rscript reproduce/analysis/simulation/cn_state_precision_recall.R
# Input : <setting>/profiles.tsv, cell_types.tsv, Calling/MintCNA_h5{,_single,_rnaonly}.RData,
#         combined_multiomic_100kb_allcells.Rdata  (settings c3_N100_100kb_l1e4{,_c0.05,_c0.2})
# Output: Simulation/results/tables/s5_cn_pr_3cov.csv
# Run from the root of the MintCNA repository with MINTCNA_DATA set (see reproduce/README.md).
suppressPackageStartupMessages(library(data.table))
ROOT <- Sys.getenv("MINTCNA_DATA"); if (!nzchar(ROOT)) stop("set MINTCNA_DATA (see reproduce/README.md)")
CN <- file.path(ROOT, "Simulation/cnasim")

parse_bins <- function(bin_ids) {
  b <- gsub("^chr", "", bin_ids, ignore.case = TRUE)
  if (any(grepl(":", b))) {
    chr <- sub(":.*$", "", b); pos <- sub("^.*?:", "", b)
    start <- as.integer(sub("-.*$", "", pos)); end <- as.integer(sub("^.*?-", "", pos))
  } else {
    sp <- strsplit(b, "_", fixed = TRUE); chr <- vapply(sp, `[[`, "", 1)
    start <- as.integer(vapply(sp, `[[`, "", 2)); end <- as.integer(vapply(sp, `[[`, "", 3))
  }
  data.table(chr = chr, start = start, end = end)
}
make_segments <- function(states, bins_dt) {
  stopifnot(length(states) == nrow(bins_dt)); dt <- copy(bins_dt); dt[, state := states]
  dt[, { r <- rle(state); ends <- cumsum(r$lengths); starts <- ends - r$lengths + 1L
         data.table(state = r$values, start = start[starts], end = end[ends]) }, by = chr]
}
truth_dir <- function(cn) fifelse(cn %in% c(0L, 1L), "DEL", fifelse(cn == 2L, "NEU", "DUP"))
pred_dir  <- function(x)  fifelse(x == 1L, "DEL", fifelse(x == 2L, "NEU", "DUP"))

cnk_counts_one_cell <- function(truth_vec, pred_vec, bins_dt, cnk) {
  truth_vec <- as.integer(truth_vec); pred_vec <- as.integer(pred_vec)
  truth_segs_cn <- make_segments(truth_vec, bins_dt)
  pred_segs_dir <- make_segments(pred_dir(pred_vec), bins_dt)
  dir_k <- truth_dir(cnk)
  pred_k <- pred_segs_dir[state == dir_k]; pred_k[, pid := .I]
  truth_dir_pos <- truth_segs_cn[(dir_k == "DUP" & state %in% c(3L, 4L, 5L)) | (dir_k == "DEL" & state %in% c(0L, 1L))]
  FP <- 0L
  if (nrow(pred_k) > 0L) {
    if (nrow(truth_dir_pos) == 0L) FP <- nrow(pred_k)
    else { setkey(pred_k, chr, start, end); setkey(truth_dir_pos, chr, start, end)
      ov2 <- foverlaps(pred_k, truth_dir_pos, type = "any", nomatch = 0L)
      FP <- nrow(pred_k) - length(unique(ov2$pid)) }
  }
  truth_pos <- truth_segs_cn[state == cnk]
  if (nrow(truth_pos) == 0L) return(data.table(cn = paste0("CN", cnk), TP = 0L, FP = as.integer(FP), FN = 0L))
  truth_pos[, id := .I]
  if (nrow(pred_k) == 0L) { TP <- 0L; FN <- nrow(truth_pos) }
  else { setkey(truth_pos, chr, start, end); setkey(pred_k, chr, start, end)
    ov <- foverlaps(truth_pos, pred_k, type = "any", nomatch = 0L)
    TP <- length(unique(ov$id)); FN <- nrow(truth_pos) - TP }
  data.table(cn = paste0("CN", cnk), TP = as.integer(TP), FP = as.integer(FP), FN = as.integer(FN))
}
cnk_eval <- function(truth_mat, pred_mat, cnks = c(0L, 1L, 3L, 4L, 5L)) {
  stopifnot(all(dim(truth_mat) == dim(pred_mat)), identical(rownames(truth_mat), rownames(pred_mat)),
            identical(colnames(truth_mat), colnames(pred_mat)))
  bins_dt <- parse_bins(rownames(truth_mat))
  per_cell <- rbindlist(lapply(colnames(truth_mat), function(cell)
    rbindlist(lapply(cnks, function(k) cnk_counts_one_cell(truth_mat[, cell], pred_mat[, cell], bins_dt, k)))[, cell := cell]))
  out <- per_cell[, .(TP = sum(TP), FP = sum(FP), FN = sum(FN), n_cells = uniqueN(cell)), by = .(cn)]
  out[, precision := fifelse((TP + FP) == 0, NA_real_, TP / (TP + FP))]
  out[, recall    := fifelse((TP + FN) == 0, NA_real_, TP / (TP + FN))]
  out[, F1 := fifelse(is.na(precision) | is.na(recall) | (precision + recall) == 0, NA_real_, 2 * precision * recall / (precision + recall))]
  out[]
}
coerce_123 <- function(mat) { v <- as.integer(mat)
  if (all(v %in% c(1L, 2L, 3L, NA_integer_))) return(matrix(v, nrow(mat), ncol(mat), dimnames = dimnames(mat)))
  o <- rep(NA_integer_, length(v)); o[v %in% c(1L, 2L)] <- 1L; o[v == 3L] <- 2L; o[v >= 4L] <- 3L
  matrix(o, nrow(mat), ncol(mat), dimnames = dimnames(mat)) }

# truth: total copy number per bin and tumour cell, on MintCNA's bins
build_CNmat_sub <- function(CNA_DIR, key_need) {
  profiles <- fread(file.path(CNA_DIR, "profiles.tsv"))
  colnames(profiles)[5] <- "cn_states"; setnames(profiles, tolower(names(profiles)))
  profiles[, cn_total := { x <- trimws(as.character(cn_states))
    ifelse(grepl("\\.", x), rowSums(matrix(as.numeric(unlist(tstrsplit(x, "\\.", perl = TRUE))), ncol = 2, byrow = FALSE)),
           round(as.numeric(x))) |> as.integer() }]
  profiles[!grepl("^chr", chrom), chrom := paste0("chr", chrom)]
  profiles[, bin_id := paste0(chrom, ":", start, "-", end)]
  profiles[, cell_num := as.integer(sub("^cell", "", cell))]; profiles[, chr_num := as.integer(sub("^chr", "", chrom))]
  profiles <- profiles[grepl("^cell\\d+$", cell)]
  setorder(profiles, chr_num, start, end, cell_num)
  cells <- profiles[, unique(.SD[order(cell_num)]$cell), .SDcols = c("cell", "cell_num")]
  CNdt <- dcast(profiles, bin_id ~ cell, value.var = "cn_total")
  CNmat <- as.matrix(CNdt[, ..cells]); rownames(CNmat) <- CNdt$bin_id
  CNmat <- CNmat[order(as.integer(sub("^chr", "", sub(":.*$", "", rownames(CNmat)))),
                       as.integer(sub("^.*:(\\d+)-.*$", "\\1", rownames(CNmat)))), , drop = FALSE]
  cn_bins <- rownames(CNmat)
  key_CN <- paste(sub("^chr", "", sub(":.*$", "", cn_bins)), as.integer(sub("-.*$", "", sub("^.*:", "", cn_bins))),
                  as.integer(sub("^.*-", "", sub("^.*:", "", cn_bins))), sep = "_")
  idx <- match(key_need, key_CN); if (anyNA(idx)) stop("bins missing in the truth")
  nc <- fread(file.path(CNA_DIR, "cell_types.tsv"))[normal == "normal", root]
  CNmat[idx, !colnames(CNmat) %in% nc, drop = FALSE]
}

covs <- list(c("0.05X", "c3_N100_100kb_l1e4_c0.05"), c("0.1X", "c3_N100_100kb_l1e4"), c("0.2X", "c3_N100_100kb_l1e4_c0.2"))
res <- list()
for (cv in covs) {
  lab <- cv[1]; CNA_DIR <- file.path(CN, cv[2]); cat("\n===== coverage", lab, "(", cv[2], ") =====\n")
  e <- new.env()
  load(file.path(CNA_DIR, "Calling/MintCNA_h5.RData"), e)
  load(file.path(CNA_DIR, "Calling/MintCNA_h5_single.RData"), e)
  load(file.path(CNA_DIR, "Calling/MintCNA_h5_rnaonly.RData"), e)
  key_need <- rownames(e$cna_multicna)
  truth_mat <- build_CNmat_sub(CNA_DIR, key_need); truth_mat[truth_mat >= 5] <- 5L
  # bins with RNA coverage (the same set in every cell): the RNA-only mode is scored on these only
  fenv <- new.env(); load(file.path(CNA_DIR, "combined_multiomic_100kb_allcells.Rdata"), fenv)
  cmx <- fenv$combined_multiomic_list[[1]]
  keep_rna <- key_need %in% as.character(cmx[, 1])[!is.na(cmx[, 2])]
  cat(sprintf("  RNA-covered bins: %d of %d (%.1f%%)\n", sum(keep_rna), length(keep_rna), 100 * mean(keep_rna)))
  preds <- list(`DNA+RNA` = e$cna_multicna, `DNA-only` = e$cna_multicna_single, `RNA-only` = e$cna_multicna_rnaonly)
  for (nm in names(preds)) {
    M <- coerce_123(preds[[nm]])
    cc <- intersect(colnames(truth_mat), colnames(M))
    if (length(cc) != ncol(truth_mat)) cat("  NOTE:", nm, "cells", length(cc), "of", ncol(truth_mat), "\n")
    Tm <- truth_mat[, cc, drop = FALSE]; M <- M[, cc, drop = FALSE]
    rownames(M) <- rownames(Tm)
    if (nm == "RNA-only") { Tm <- Tm[keep_rna, , drop = FALSE]; M <- M[keep_rna, , drop = FALSE] }
    r <- cnk_eval(Tm, M)
    r[, `:=`(coverage = lab, method = nm, n_bins = nrow(Tm))]; res[[paste(lab, nm)]] <- r
  }
}
out <- rbindlist(res)
gain <- out[cn %in% c("CN3", "CN4", "CN5")][, .(coverage, method, cn, n_bins, precision = round(precision, 3),
                                                 recall = round(recall, 3), F1 = round(F1, 3), TP, FP, FN)]
setorder(gain, coverage, cn, method)
print(gain, row.names = FALSE)
dir.create(file.path(ROOT, "Simulation/results/tables"), recursive = TRUE, showWarnings = FALSE)
fwrite(gain, file.path(ROOT, "Simulation/results/tables/s5_cn_pr_3cov.csv"))
cat("\nwrote Simulation/results/tables/s5_cn_pr_3cov.csv\n")
