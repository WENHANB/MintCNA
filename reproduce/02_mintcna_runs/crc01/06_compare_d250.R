#!/usr/bin/env Rscript
# Before/after table for every CRC01 number in the paper: published files vs the 2026-10-14
# re-run (AG-DCAE latent_dim 250, 197 cells, DNA-only run with its own change points).
.libPaths(c("<R_LIBS_USER>", .libPaths()))
suppressPackageStartupMessages(library(data.table))
ROOT <- "<MINTCNA_ROOT>"; setwd(ROOT)
source("Result/CRC/figures/scripts/utils.R"); source("MintCNA_git/MintCNA/call.R")
B <- "Result/CRC/Benchmark"
site_of <- function(x) factor(sub("^CRC\\d+_([A-Za-z]+)[0-9]*_\\d+$", "\\1", x), levels = c("PT","LN","ML","MP"))
arm_loss <- function(cna, ref, chr, lo = -Inf, hi = Inf) {   # fraction of arm bins called loss (state 1) > 0.5
  idx <- ref$chr == chr & ref$start >= lo & ref$end <= hi; colMeans(cna[idx, , drop = FALSE] == 1, na.rm = TRUE) > 0.5 }
P17 <- 25248165; Q5 <- 49405641; APC <- c(112043195, 112181936)                 # hg19
load_set <- function(sfx) { e <- new.env(); load(file.path(B, paste0("MultiCNA_AE_197", sfx, ".RData")), e)
  cna <- get("cna_multicna_ae", e); ref <- get("ref_multicna_ae", e); ref$chr <- sub("^chr", "", as.character(ref$chr)); list(cna = cna, ref = ref) }
lines <- c("| quantity | before (published) | after (d = 250, 197 cells) |", "|---|---|---|")
add <- function(q, a, b) lines <<- c(lines, sprintf("| %s | %s | %s |", q, a, b))
S <- list(before = load_set(""), after = load_set("_d250"))
EXCL <- strsplit(Sys.getenv("CRC_EXCLUDE", "CRC01_ML3_085"), ",")[[1]]; EXCL <- EXCL[nzchar(EXCL)]
if (length(EXCL)) S$after$cna <- S$after$cna[, !colnames(S$after$cna) %in% EXCL, drop = FALSE]
for (nm in names(S)) S[[nm]]$site <- site_of(colnames(S[[nm]]$cna))
## 1. 17p / 17q loss by site -----------------------------------------------------------
for (arm in c("17p", "17q")) {
  v <- sapply(S, function(s) { l <- if (arm == "17p") arm_loss(s$cna, s$ref, "17", hi = P17) else arm_loss(s$cna, s$ref, "17", lo = P17 + 1)
    t <- tapply(l, s$site, function(x) sprintf("%.1f%% (%d/%d)", 100 * mean(x), sum(x), length(x))); paste(names(t), t, collapse = ", ") })
  add(paste(arm, "loss by site"), v[["before"]], v[["after"]]) }
## 2. subclones ------------------------------------------------------------------------
lat_old <- as.matrix(read.csv("Data/CRC/DNAmeth/processed_l2r/latent_z_dna.csv", row.names = 1, check.names = FALSE))
r0 <- cluster_latent_kmeans(lat_old, seed = 1, K_min = 2, K_max = 2, scale_Z = FALSE); sz <- table(r0$clusters)
outl <- names(r0$clusters)[r0$clusters == as.integer(names(sz)[which.min(sz)])]
r1 <- cluster_latent_kmeans(lat_old[setdiff(rownames(lat_old), outl), , drop = FALSE], seed = 1, K_min = 2, K_max = 2, scale_Z = FALSE)
cl_old <- r1$clusters[intersect(names(r1$clusters), colnames(S$before$cna))]
new_lab <- fread("Result/CRC/subclones_d250.csv"); new_lab <- new_lab[!cell %in% EXCL]; cl_new <- setNames(new_lab$cluster, new_lab$cell)
relabel <- function(cl) { o <- names(sort(table(cl), decreasing = TRUE)); setNames(paste0("C", match(as.character(cl), o)), names(cl)) }   # C1 = largest, C2, C3 ... by size
cl_old <- relabel(cl_old); cl_new <- relabel(cl_new)
desc <- function(cl, cna, ref) { s <- site_of(names(cl)); out <- c()
  for (k in sort(unique(cl))) { i <- names(cl)[cl == k]; t <- table(factor(s[cl == k], levels = levels(s)))
    q17 <- mean(arm_loss(cna[, i, drop = FALSE], ref, "17", lo = P17 + 1)); q5 <- mean(arm_loss(cna[, i, drop = FALSE], ref, "5", lo = Q5))
    apc <- ref$chr == "5" & ref$end >= APC[1] & ref$start <= APC[2]; apcl <- if (length(i) == 1) mean(cna[apc, i] == 1) > 0.5 else mean(colMeans(cna[apc, i, drop = FALSE] == 1) > 0.5)
    out <- c(out, sprintf("%s n=%d (%s); 17q loss %.0f%%, 5q loss %.0f%%, APC bins loss %.0f%%", k, length(i), paste(names(t), t, collapse = "/"), 100*q17, 100*q5, 100*apcl)) }
  paste(out, collapse = "<br>") }
add("chosen K", "2 (fixed, outlier CRC01_ML3_085 removed, unscaled latent)", sprintf("%d (documented procedure, no outlier removal)", new_lab$K_hat[1]))
add("clone sizes, site composition, defining events", desc(cl_old, S$before$cna, S$before$ref), desc(cl_new, S$after$cna, S$after$ref))
add("ARI(before labels, after labels) on shared cells", "", sprintf("%.3f (n = %d)", mclust::adjustedRandIndex(cl_old[intersect(names(cl_old), names(cl_new))], cl_new[intersect(names(cl_old), names(cl_new))]), length(intersect(names(cl_old), names(cl_new)))))
## 3. accuracy vs bulk WGS (Fig. 4b) ---------------------------------------------------
per_cell <- function(T, P) { P <- as.matrix(P); cm <- intersect(rownames(T), rownames(P)); if (length(cm)) { T <- T[cm, , drop = FALSE]; P <- P[cm, , drop = FALSE] }
  m <- cell_to_truthcol(colnames(P)); m <- m[!is.na(truth_col) & truth_col %in% colnames(T)]
  m[, rmse := vapply(seq_len(.N), function(i) sqrt(mean((T[, truth_col[i]] - P[, cell[i]])^2, na.rm = TRUE)), numeric(1))]
  m[, pearson := vapply(seq_len(.N), function(i) { x <- as.numeric(T[, truth_col[i]]); y <- as.numeric(P[, cell[i]]); k <- is.finite(x) & is.finite(y); if (sum(k) < 2 || sd(x[k]) == 0 || sd(y[k]) == 0) NA_real_ else cor(x[k], y[k]) }, numeric(1))]; m }
for (sfx in c("", "_d250")) { d <- load_all_methods(sfx = sfx, exclude = if (sfx == "") character(0) else EXCL); for (mt in c("MintCNA", "MintCNA(DNA)")) { pc <- per_cell(d$truth_mat, d$methods[[mt]])
  assign(paste0("acc", sfx, mt), sprintf("r %.3f / RMSE %.3f (n = %d)", median(pc$pearson, na.rm = TRUE), median(pc$rmse, na.rm = TRUE), nrow(pc))) } }
add("median Pearson r / RMSE vs bulk WGS, MintCNA", get("accMintCNA"), get("acc_d250MintCNA"))
add("median Pearson r / RMSE vs bulk WGS, MintCNA(DNA)", get("accMintCNA(DNA)"), get("acc_d250MintCNA(DNA)"))
## 4. runtime ------------------------------------------------------------------------
grab <- function(f, pat) { if (!file.exists(f)) return("pending"); x <- readLines(f, warn = FALSE); x <- x[grepl(pat, x)]; if (length(x)) paste(trimws(x), collapse = "; ") else "pending" }
L <- "<REPRODUCE_DIR>/02_mintcna_runs/crc01/logs"; TL <- "<REPRODUCE_DIR>/02_mintcna_runs/crc01/timing/logs"
ae_default <- grab(file.path(L, "ae_45883922.out"), "^AE rc=")
ae_turin   <- grab(Sys.glob(file.path(TL, "ae_turin_*.out"))[1], "^AE TIMING")
sc_turin   <- grab(file.path(TL, "single_cps_turin_R.log"), "^TIMING single cps")
add("runtime: AG-DCAE training (8 CPU cores)", "part of the 80 / 75 min totals (Computational Time.xlsx, composition not recorded)", sprintf("hpg-default: %s; hpg-turin: %s", ae_default, ae_turin))
add("runtime: change points + calling, multi-omics", "80 min total (Computational Time.xlsx, 8 cores)", paste("hpg-turin:", grab(file.path(L, "call_multi_R.log"), "^RUNTIME")))
add("runtime: change points + calling, DNA-only", "75 min total", sprintf("hpg-default: %s; hpg-turin: %s", grab(file.path(L, "call_single_R.log"), "^RUNTIME"), sc_turin))
writeLines(c("# CRC01 before/after (2026-10-14 re-run: AG-DCAE latent_dim 250, 197 QC-passing cells, raw RNA, DNA-only run with its own change points)", "", lines), "<REPRODUCE_DIR>/02_mintcna_runs/crc01/before_after.md")
cat(paste(lines, collapse = "\n"), "\n")
