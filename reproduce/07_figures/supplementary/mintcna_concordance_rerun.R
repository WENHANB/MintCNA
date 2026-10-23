# Re-run of the partial-concordance MintCNA calls (25 / 50 / 75 % of genes
# concordant) with the CURRENT pipeline. Verbatim copy of
# Simulation/code/MintCNA/MintCNA.R (2026-03-10) except for the three paths: the
# input is combined_multiomic_100kb_allcells_p<k>.Rdata and every output carries
# the suffix _p<k>_2026-10-09 so the 2026-02-16 files are never touched.
# Why: Calling/MintCNA_h5_p1..3.Rdata date from 2026-02-16, but the 0 % and 100 %
# endpoints (MintCNA_h5_single / MintCNA_h5, also used by Fig. 2) were re-called on
# 2026-03-19 after Code/MultiCNA changed (Multi_Inspect.R 03-10, GMM_CNstatus.R
# 03-18, cluster.R 03-17), so Supp Fig 5a mixed two pipeline versions.
#   Rscript mintcna_concordance_rerun.R <k>        k in 1, 2, 3
args <- commandArgs(trailingOnly = TRUE)
K <- as.integer(if (length(args) >= 1) args[1] else 1L); stopifnot(K %in% 1:3)
BASE <- "c3_N100_100kb_l1e4"; print(c(BASE, K))
CNA_DIR <- file.path("<MINTCNA_ROOT>/Simulation/cnasim", BASE)
SFX <- sprintf("_p%d_2026-10-09", K)
.libPaths(c("<R_LIBS_USER>", .libPaths()))

load(file.path(CNA_DIR, sprintf("combined_multiomic_100kb_allcells_p%d.Rdata", K)))
files <- list.files(
  "<REPRODUCE_DIR>/lib/MultiCNA",
  pattern = "\\.[rR]$",
  full.names = TRUE
)
invisible(lapply(files, source))


filtered_multiomic_list <- combined_multiomic_list

ord <- with(ref, order(as.numeric(sub("chr", "", chr)), start))
ref <- ref[ord,]

#summary(as.numeric(lapply(filtered_multiomic_list, function(x)sd(x[,3], na.rm = T))))

## QN ALL
sd_2 <- median(as.numeric(lapply(filtered_multiomic_list, function(x)sd(x[,2], na.rm = T))))
sd_3 <- median(as.numeric(lapply(filtered_multiomic_list, function(x)sd(x[,3], na.rm = T))))
filtered_multiomic_list_qn_all <- lapply(filtered_multiomic_list, function(df) {
  # 2×n matrix: DNA in row 1, RNA in row 2
  X <- rbind(RNA = as.numeric(df[,2]),
             DNA = as.numeric(df[,3]))
  
  # QN on all chr
  X_qn <- qn_modal(X,sd_target = c(sd_2,sd_3))
  
  # Replace original RNA and DNA columns
  df[,2] <- X_qn["RNA", ]
  df[,3] <- X_qn["DNA", ]
  
  df
})

#summary(as.numeric(lapply(filtered_multiomic_list_qn_all, function(x)sd(x[,2], na.rm = T))))

## Change-points detection
timing_cell <- data.frame(cell = character(), seconds = numeric(), stringsAsFactors = FALSE)
cps_track <- list()
h <- 5
for (cell_name in names(filtered_multiomic_list)) {
  t0_cell <- Sys.time() 
  print(cell_name)
  ## Select cell
  cell <- filtered_multiomic_list_qn_all[[cell_name]]
  rownames(cell) <- cell[,1];cell <- cell[,-1];storage.mode(cell) <- "numeric"
  cell <- t(cell)
  
  chr_length <- c() 
  cps_multi_list <- list()
  for (chr in unique(ref$chr)) {
    print(chr)
    cell_chr <- cell[, which(ref$chr==chr)]
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
#mean(paste0("Mean per cell in seconds:",timing_cell$seconds))
save(timing_cell,cps_track, file = file.path(CNA_DIR, paste0("MintCNA_cps_h5", SFX, ".RData")))

load(file.path(CNA_DIR, paste0("MintCNA_cps_h5", SFX, ".RData")))
## cps 
#inv <- multiInverse(cell, h = h, simT = 1,lambda = NULL, penalty = "l2")
res_null <- multiInverse_chr(x = cell, ref = ref, h = h, precise = 1000, simT = 1)
alpha <- 0.30

cell_cps_list <- list()
for (cell_name in names(filtered_multiomic_list)) {
  print(cell_name)
  
  chr_length <- c()
  cps_multi_list <- list()
  for (chr in unique(ref$chr)) {
    #print(chr)
    cell_chr <- cell[, which(ref$chr==chr)]
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
  cna_df <- cna_df[which(cna_df$end_idx-cna_df$start_idx > 1),]
  cna_df_all <- rbind(cna_df_all, cna_df)
}
cna_df_all <- cna_df_all[which(cna_df_all$cell %in% converge_cell),]
## Calling result to output matrix
# 3 status
code_map3 <- c(del = 1, nor = 2, dup = 3)
cna_multicna <- build_cna_matrix(cna_df_all, ref,code_map = code_map3)
#cna_multicna_ae <- build_cna_matrix(cna_df_all, ref,code_map = code_map3)

# 5 status
#code_map5 <- c(del = 1, loss = 2, nor = 3, gain = 4, amp = 5)
#cna_multicna <- build_cna_matrix(cna_df_all, ref, code_map = code_map5, default_state = "nor")

ref_multicna <- ref

save(cna_multicna, ref_multicna, file = file.path(CNA_DIR, paste0("Calling/MintCNA_h5", SFX, ".RData")))
cat("wrote", file.path(CNA_DIR, paste0("Calling/MintCNA_h5", SFX, ".RData")), "\n")
#save(cna_multicna_ae, ref_multicna, file = file.path(CNA_DIR,"Calling/MintCNA_h10_AE.RData"))

