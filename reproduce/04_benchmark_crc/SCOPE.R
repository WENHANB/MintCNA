setwd("<MINTCNA_ROOT>/")
load("Data/CRC/DNAmeth/dnameth_chr_all_cells_100k_gc_mapp.RData")  
load("Data/CRC/DNAmeth/cell_keep.RData")

suppressPackageStartupMessages({
  library(dplyr); library(GenomicRanges)
  library(SCOPE); library(WGSmapp)
  library(BSgenome.Hsapiens.UCSC.hg19)
})

nc_samples <- c("CRC01_NC_507","CRC01_NC_509","CRC01_NC_517",
                "CRC01_NC_521","CRC01_NC_525","CRC01_NC_527","CRC01_NC_529")

## --- 3) Bin QC (GC/mapp) + rebuild filtered GRanges ---
keep <- !is.na(gr$gc) & !is.na(gr$mapp) & gr$gc >= 20 & gr$gc <= 80 & gr$mapp >= 0.90
#Y      <- Y[keep, c(cell_keep, nc_samples)]
Y      <- Y[keep, ]
ref_gr <- gr[keep]

## --- 4) Normal indices (your 7 NC cells) ---
norm_index <- match(nc_samples, colnames(Y))
norm_index <- norm_index[!is.na(norm_index)]
if (length(norm_index) < 3) warning("Fewer than 3 normals found in Y; check names.")

## --- 5) Normalization -> ploidy -> SCOPE ---
gc_qc <- as.numeric(mcols(ref_gr)$gc) 
norm0 <- normalize_codex2_ns_noK(Y_qc = Y, gc_qc = gc_qc, norm_index = norm_index)

# ploidy <- initialize_ploidy(Y = Y, Yhat = norm0$Yhat, ref = ref_gr)
# normObj <- normalize_scope_foreach(
#   Y_qc = Y, gc_qc = gc_qc, K = 1, ploidyInt = ploidy,
#   norm_index = norm_index, T = 1:5, beta0 = norm0$beta.hat, nCores = 2
# )
# 
# pick <- which.max(normObj$BIC)
# Yhat <- normObj$Yhat[[pick]]

Yhat <- norm0$Yhat

chr_levels <- paste0("chr", 1:22)

## --- 6) Cross-sample segmentation per chr -> integer CN ---
suppressPackageStartupMessages({ library(parallel) })
# use all CPUs Slurm gave this task (fallback=4 locally)
ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
options(mc.cores = ncores)

seglist <- mclapply(chr_levels, function(ch){
  message("\n== ", ch, " ==\n")
  segment_CBScs(
    Y = Y, Yhat = Yhat, sampname = colnames(Y),
    ref = ref_gr, chr = ch, mode = "integer", max.ns = 1
  )
}, mc.preschedule = F)  # good when work per chr is uneven
names(seglist) <- chr_levels

iCN <- do.call(rbind, lapply(seglist, function(z){z[["iCN"]]}))

cna_scope <- iCN; ref_scope <- as.data.frame(ref_gr)

save(cna_scope, ref_scope, file = "<MINTCNA_ROOT>/Result/CRC/Benchmark/SCOPE_allcells.RData")

source("<REPRODUCE_DIR>/lib/MultiCNA/CNAout.R")
code_map5 <- c(del = 1, loss = 2, nor = 3, gain = 4, amp = 5)
plot_S_heatmap(cna_scope, ref_scope, code_map = code_map5, 
               file = "Result/CRC/Segmentation/CNV/SCOPE_allcells.pdf",
               width = 12, height = 20)
