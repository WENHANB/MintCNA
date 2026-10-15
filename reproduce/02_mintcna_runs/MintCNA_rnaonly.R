## MintCNA RNA-ONLY mode (per-cell). Mirror of MintCNA_single.R (DNA-only) but driven
## by the RNA channel: RNA (row 1 / col 2) drives change-points AND state-calling; the
## DNA channel is dropped. AG-DCAE (DNA denoiser) is not used. Reuses the existing
## combined_multiomic_100kb_allcells.Rdata (RNA values genuine). Output object
## cna_multicna_rnaonly -> Calling/MintCNA_h5_rnaonly.RData.
.libPaths(c("<R_LIBS_USER>", .libPaths()))  # use R4.4 lib (softImpute/InspectChangepoint)
args <- commandArgs(trailingOnly = TRUE)
BASE <- if (length(args) >= 1) args[1] else "c3_N100_100kb_l1e4"
print(BASE)
CNA_DIR <- file.path("<MINTCNA_ROOT>/Simulation/cnasim", BASE)

load(file.path(CNA_DIR,"combined_multiomic_100kb_allcells.Rdata"))
files <- list.files("<MINTCNA_ROOT>/Code/MultiCNA",
                    pattern = "\\.[rR]$", full.names = TRUE)
invisible(lapply(files, source))

filtered_multiomic_list <- combined_multiomic_list
ord <- with(ref, order(as.numeric(sub("chr", "", chr)), start))
ref <- ref[ord,]

## QN ALL (row 1 = RNA = df[,2]; row 2 = DNA = df[,3])
sd_2 <- median(as.numeric(lapply(filtered_multiomic_list, function(x)sd(x[,2], na.rm = T))))
sd_3 <- median(as.numeric(lapply(filtered_multiomic_list, function(x)sd(x[,3], na.rm = T))))
filtered_multiomic_list_qn_all <- lapply(filtered_multiomic_list, function(df) {
  X <- rbind(RNA = as.numeric(df[,2]), DNA = as.numeric(df[,3]))
  X_qn <- qn_modal(X, sd_target = c(sd_2, sd_3))
  df[,2] <- X_qn["RNA", ]; df[,3] <- X_qn["DNA", ]
  df
})

## Change-points detection (RNA-driven: zero the DNA row 2)
timing_cell <- data.frame(cell = character(), seconds = numeric(), stringsAsFactors = FALSE)
cps_track <- list(); h <- 5
for (cell_name in names(filtered_multiomic_list)) {
  t0_cell <- Sys.time(); print(cell_name)
  cell <- filtered_multiomic_list_qn_all[[cell_name]]
  rownames(cell) <- cell[,1]; cell <- cell[,-1]; storage.mode(cell) <- "numeric"
  cell <- t(cell)
  cell[2,] <- c(0.001, rep(0, ncol(cell)-1))          # zero DNA -> RNA drives change-points
  for (chr in unique(ref$chr)) {
    cell_chr <- cell[, grepl(paste0(chr,"_"), colnames(cell))]
    x.cusum.val_multi_l2 <- multiMissInspect_scores(cell_chr, h = h, penalty = "l2")
    cps_track[[cell_name]][[chr]] <- x.cusum.val_multi_l2
  }
  timing_cell <- rbind(timing_cell, data.frame(cell = cell_name,
    seconds = as.numeric(difftime(Sys.time(), t0_cell, units = "secs")), stringsAsFactors = FALSE))
}
save(timing_cell, cps_track, file = file.path(CNA_DIR,"MintCNA_cps_h5_rnaonly.RData"))
load(file.path(CNA_DIR,"MintCNA_cps_h5_rnaonly.RData"))

## cps thresholding
res_null <- multiInverse_chr(x = cell, ref = ref, h = h, precise = 1000, simT = 1)
alpha <- 0.30
cell_cps_list <- list()
for (cell_name in names(filtered_multiomic_list)) {
  chr_length <- c(); cps_multi_list <- list()
  for (chr in unique(ref$chr)) {
    ref_chr <- ref[which(ref$chr==chr),]
    x.cusum.val_multi_l2 <- cps_track[[cell_name]][[chr]]
    tau <- res_null$q[res_null$probs == (1 - alpha)]
    cps_multi <- sort(x.cusum.val_multi_l2$local_maxima[abs(x.cusum.val_multi_l2$local_maxima$score)>tau,1])
    chr_length <- c(chr_length, nrow(ref_chr)); cps_multi_list[[chr]] <- cps_multi
  }
  cell_cps_list[[cell_name]] <- cps_multi_list
}
summary(sapply(cell_cps_list, function(x) sum(sapply(x, length))))

## CNA calling (RNA-only: feed row 1 = RNA into CNVcluster)
cna_df_all <- c(); converge_cell <- c(); L <- 30
for (cell_name in names(filtered_multiomic_list)) {
  print(cell_name)
  cell <- filtered_multiomic_list_qn_all[[cell_name]]
  rownames(cell) <- cell[,1]; cell <- cell[,-1]; storage.mode(cell) <- "numeric"
  cell <- t(cell)
  # states: keep BOTH channels fused (like MintCNA_single.R). RNA log2R (~+/-0.25) is
  # below the TGMM truncation eps0=0.30, so RNA-only calling would collapse to all-Neutral;
  # fusing DNA preserves the dynamic range the truncation needs. RNA-only enters via
  # segmentation only (DNA zeroed in the change-point step above).
  idx <- which(names(cell_cps_list) == cell_name)
  cps_cell <- reindex_cps(cell_cps_list[[idx]], chr_length)
  cps_cell <- cps_cell[-c(1,length(cps_cell))]
  cps_cell <- cps_cell[which(cps_cell < ncol(cell))]
  result_cnv <- CNVcluster(Y = cell, cp = cps_cell, L = L, st = 3, method = "count")
  if (result_cnv$iter < L) converge_cell <- c(converge_cell, cell_name)
  cna_df <- make_cna_table(result_cnv$cnv.state, result_cnv$cnv.start, result_cnv$cnv.end, ref)
  cna_df$mean <- mapply(function(s, e) mean(cell[, s:e], na.rm = TRUE), cna_df$start_idx, cna_df$end_idx)
  cna_df$cell <- cell_name; cna_df$bin_length <- cna_df$end_idx-cna_df$start_idx
  cna_df <- cna_df[which(cna_df$end_idx-cna_df$start_idx > 1),]
  cna_df_all <- rbind(cna_df_all, cna_df)
}
cna_df_all <- cna_df_all[which(cna_df_all$cell %in% converge_cell),]
code_map3 <- c(del = 1, nor = 2, dup = 3)
cna_multicna_rnaonly <- build_cna_matrix(cna_df_all, ref, code_map = code_map3)
print(table(cna_multicna_rnaonly))
ref_multicna <- ref
save(cna_multicna_rnaonly, ref_multicna, file = file.path(CNA_DIR,"Calling/MintCNA_h5_rnaonly.RData"))
cat("saved Calling/MintCNA_h5_rnaonly.RData\n")
