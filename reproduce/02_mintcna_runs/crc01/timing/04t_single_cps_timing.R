# CRC01 re-run 2026-10-14: copy of Code/CRC/CNA_calling_QNALL_AE_single.R on the d=250 input, DNA-only change points.
setwd("<MINTCNA_ROOT>/"); T0 <- Sys.time()
load("Data/CRC/combined_multiomic_100kb_allcells_AE_d250.Rdata")   # d=250 re-run
load("Data/CRC/DNAmeth/cell_keep.RData")

files <- list.files(
  "<REPRODUCE_DIR>/lib/MultiCNA",
  pattern = "\\.[rR]$",
  full.names = TRUE
)
invisible(lapply(files, source))

## Filter out the outliners
# 1) per-cell indices of bins to drop (dna < -4)
# idx_per_df <- lapply(combined_multiomic_list, function(m) {which(as.numeric(m[, 3]) < -4)})
# 
# # 2) count dropouts per cell
# dropout_counts <- sapply(idx_per_df, length)

# 3) keep only cells with < 5000 dropouts
# cells_keep <- names(dropout_counts)[dropout_counts < 10]
#filtered_cells_list <- combined_multiomic_list[cell_keep]
#filtered_multiomic_list <- combined_multiomic_list[cell_keep]
filtered_multiomic_list <- combined_multiomic_list
# idx_per_df_keep <- idx_per_df[cells_keep]

# 4) union of bad rows across the kept cells
#drop_idx <- sort(unique(unlist(idx_per_df_keep)))

# 5) remove those rows from all retained cells
# filtered_multiomic_list <- lapply(filtered_cells_list, function(m) {
#   if (length(drop_idx)) {
#     keep <- setdiff(seq_len(nrow(m)), drop_idx)
#     m[keep, , drop = FALSE]
#   } else {
#     m
#   }
# })

#filtered_multiomic_list <- filtered_multiomic_list[!(names(filtered_multiomic_list) %in% c("CRC01_PT1_297", "CRC01_PT1_203")) ]

ord <- with(ref, order(as.numeric(sub("chr", "", chr)), start))
ref <- ref[ord,]

summary(as.numeric(lapply(filtered_multiomic_list, function(x)sd(x[,3], na.rm = T))))
summary(as.numeric(lapply(filtered_multiomic_list, function(x)sd(x[,2], na.rm = T))))

## QN ALL
filtered_multiomic_list_qn_all <- lapply(filtered_multiomic_list, function(df) {
  # 2×n matrix: DNA in row 1, RNA in row 2
  X <- rbind(RNA = as.numeric(df[,2]),
             DNA = as.numeric(df[,3]))
  
  # QN on all chr
  X_qn <- qn_modal(X,sd_target = c(0.25,0.5))
  
  # Replace original RNA and DNA columns
  df[,2] <- X_qn["RNA", ]
  df[,3] <- X_qn["DNA", ]
  
  df
})

summary(as.numeric(lapply(filtered_multiomic_list_qn_all, function(x)sd(x[,3], na.rm = T))))

## Change-points detection
#cell_name <- "CRC01_PT1_231"
h <- 15
if (file.exists("<REPRODUCE_DIR>/02_mintcna_runs/crc01/timing/h15_AE_197_single_d250_turin.RData")) {
  # rebuild the `cell` matrix the downstream code expects (last cell of the list, as in the loop)
  cell <- filtered_multiomic_list_qn_all[[length(filtered_multiomic_list_qn_all)]]
  rownames(cell) <- cell[,1]; cell <- cell[,-1]; storage.mode(cell) <- "numeric"; cell <- t(cell)
  cell[1,] <- c(0.001,rep(0,ncol(cell)-1))
  chr_length <- as.integer(table(factor(ref$chr, levels = unique(ref$chr))))
  load("<REPRODUCE_DIR>/02_mintcna_runs/crc01/timing/h15_AE_197_single_d250_turin.RData")   # cps_track + timing_cell
  cat("change points already computed -> loaded <REPRODUCE_DIR>/02_mintcna_runs/crc01/timing/h15_AE_197_single_d250_turin.RData\n")
} else {
timing_cell <- data.frame(cell = character(), seconds = numeric(), stringsAsFactors = FALSE)
cps_track <- list()
h <- 15 
for (cell_name in names(filtered_multiomic_list)) {
  t0_cell <- Sys.time() 
  print(cell_name)
  ## Select cell
  cell <- filtered_multiomic_list_qn_all[[cell_name]]
  rownames(cell) <- cell[,1];cell <- cell[,-1];storage.mode(cell) <- "numeric"
  cell <- t(cell)
  cell[1,] <- c(0.001,rep(0,ncol(cell)-1))
  
  chr_length <- c()
  cps_multi_list <- list()
  for (chr in unique(ref$chr)) {
    print(chr)
    cell_chr <- cell[, grepl(paste0(chr,"_"), colnames(cell))]
    ref_chr <- ref[which(ref$chr==chr),]
    
    h = h
    #x.cusum.val_multi_l1 <- multiMissInspect_scores(cell_chr, h = h, penalty = "l1")
    x.cusum.val_multi_l2 <- multiMissInspect_scores(cell_chr, h = h, penalty = "l2")
    cps_track[[cell_name]][[chr]] <- x.cusum.val_multi_l2
    
  }
  timing_cell <- rbind(timing_cell, data.frame(
    cell = cell_name,
    seconds = as.numeric(difftime(Sys.time(), t0_cell, units = "secs")),
    stringsAsFactors = FALSE
  ))
}
mean(timing_cell$seconds)
save(cps_track, timing_cell, file="<REPRODUCE_DIR>/02_mintcna_runs/crc01/timing/h15_AE_197_single_d250_turin.RData")
cat(sprintf("TIMING single cps %.1f min on %s\n", sum(timing_cell$seconds)/60, Sys.info()[["nodename"]])); quit(save = "no")
}


# (published script reloaded the multi-omics change points here; the DNA-only run now keeps its own)

## cps 
# (dead line removed 2026-10-14: multiInverse() no longer exists; res_null below is what is used)
alpha <- 0.001

res_null <- multiInverse_chr(x = cell, ref = ref, h = h, precise = 1000, simT = 1)
alpha <- 0.15

cell_cps_list <- list()
for (cell_name in names(filtered_multiomic_list)) {
  t0_cell <- Sys.time() 
  
  chr_length <- c()
  cps_multi_list <- list()
  for (chr in unique(ref$chr)) {
    #print(chr)
    cell_chr <- cell[, grepl(paste0(chr,"_"), colnames(cell))]
    ref_chr <- ref[which(ref$chr==chr),]
    
    x.cusum.val_multi_l2 <- cps_track[[cell_name]][[chr]]
    
    # D <- x.cusum.val_multi_l2$D
    # local_max <- x.cusum.val_multi_l2$local_maxima
    # v <- x.cusum.val_multi_l2$v_df
    
    #tau <- inv$q[inv$probs == (1 - alpha)]
    tau <- res_null$q[res_null$probs == (1 - alpha)]
    cps_multi <- sort(x.cusum.val_multi_l2$local_maxima[abs(x.cusum.val_multi_l2$local_maxima$score)>tau,1])
    
    ## Save result
    chr_length <- c(chr_length, nrow(ref_chr))
    cps_multi_list[[chr]] <- cps_multi
  }
  cell_cps_list[[cell_name]] <- cps_multi_list
}
summary(sapply(cell_cps_list, function(x) sum(sapply(x, length))))


## CNA clustering
cna_df_all <- c()
converge_cell <- c()
L=30
for (cell_name in names(filtered_multiomic_list)) {
  print(cell_name)
  ## Select cell
  cell <- filtered_multiomic_list_qn_all[[cell_name]]
  rownames(cell) <- cell[,1]
  cell <- cell[,-1]
  storage.mode(cell) <- "numeric"
  cell <- t(cell)
  
  idx <- which(names(cell_cps_list) == cell_name)
  cps_cell <- reindex_cps(cell_cps_list[[idx]],chr_length)
  cps_cell <- cps_cell[-c(1,length(cps_cell))]
  cps_cell <- cps_cell[which(cps_cell < ncol(cell))]
  result_cnv <- CNVcluster(Y = cell,cp = cps_cell, L = L, st = 3,method = "count")
  if(result_cnv$iter < L){
    converge_cell <- c(converge_cell, cell_name)
  }
  cna_df <- make_cna_table(result_cnv$cnv.state, result_cnv$cnv.start, result_cnv$cnv.end, ref)
  cna_df$mean <- mapply(function(s, e) mean(cell[, s:e], na.rm = TRUE),
                        cna_df$start_idx, cna_df$end_idx)
  cna_df$cell <- cell_name
  cna_df$bin_length <- cna_df$end_idx-cna_df$start_idx
  # keep only rows where dup has mean > 0 and del has mean < 0
  # cna_df <- cna_df[(cna_df$state == "dup" & cna_df$mean > 0) |
  #                    (cna_df$state == "del" & cna_df$mean < 0),]
  cna_df <- cna_df[which(cna_df$end_idx-cna_df$start_idx > 1),]
  cna_df_all <- rbind(cna_df_all, cna_df)
}
cna_df_all <- cna_df_all[which(cna_df_all$cell %in% converge_cell),]
## Calling result to output matrix
# 3 status
code_map3 <- c(del = 1, nor = 2, dup = 3)
cna_multicna_ae_single <- build_cna_matrix(cna_df_all, ref,code_map = code_map3)
# 5 status
# code_map5 <- c(del = 1, loss = 2, nor = 3, gain = 4, amp = 5)
# cna_multicna <- build_cna_matrix(cna_df_all, ref, code_map = code_map5, default_state = "nor")

ref_multicna_ae <- ref

## Cluster
# Compute cell-cell distance (use correlation distance)
# library(cluster)
# D <- as.dist(1 - cor(S, use = "pairwise.complete.obs"))
# hc <- hclust(D, method = "ward.D2")
# pick_k <- function(hc, D, ks = 2:6) {
#   sapply(ks, function(k) mean(silhouette(cutree(hc, k), D)[, 3]))
# }
# avg_s <- pick_k(hc, D, 2:6)
# if (max(avg_s, na.rm = TRUE) < 0.1) {
#   message("Likely only one major CNA cluster (no strong substructure).")
#   clusters <- rep(1, ncol(S))
# } else {
#   k_best <- (2:6)[which.max(avg_s)]
#   clusters <- cutree(hc, k = k_best)
# }
# pdf("Result/CRC/Segmentation/CNV/CNA_hierarchical_clustering.pdf", width = 7, height = 5)
# plot(hc, main = "Hierarchical clustering of CNA profiles",
#      xlab = "Cells", sub = "", cex = 0.6)
# dev.off()

## Heatmap
# S: bins × cells integer matrix (1..5)
# ref: data.frame with columns chr, start, end (row order = bin order)
plot_S_heatmap(cna_multicna_ae_single, ref_multicna_ae, code_map = code_map3, 
               file = "Result/CRC/Segmentation/CNV/MultiCNA_h15_AE_single_d250.pdf",
               width = 12, height = 20)

save(cna_multicna_ae_single, ref_multicna_ae, file = "<MINTCNA_ROOT>/Result/CRC/Benchmark/MultiCNA_AE_197_single_d250.RData"); cat("saved single, cells", ncol(cna_multicna_ae_single), "\n")
cat(sprintf("RUNTIME cps %.1f min, total %.1f min\n", sum(timing_cell$seconds)/60, as.numeric(difftime(Sys.time(), T0, units = "mins"))))
