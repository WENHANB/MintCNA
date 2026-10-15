suppressPackageStartupMessages({
  library(dplyr); library(GenomicRanges)
  library(SCOPE); library(WGSmapp)
  library(BSgenome.Hsapiens.UCSC.hg19)
  library(data.table)
})


args <- commandArgs(trailingOnly = TRUE)
BASE <- if (length(args) >= 1) args[1] else "100kb_K3_N200"
print(BASE)
CNA_DIR <- file.path("<MINTCNA_ROOT>/Simulation/cnasim", BASE)  
load(file.path(CNA_DIR,"DNA_RC.RData"))
label <- fread(file.path(CNA_DIR,"cell_types.tsv"))

nc_samples <- label[which(label$normal == "normal"), root]

seqlevelsStyle(gr) <- "UCSC"
ranges(gr) <- IRanges(start = start(gr) + 1,
                      end   = end(gr))
gr$gc <- as.numeric(mcols(gr)$gc)*100 

## --- 4) Normal indices---
norm_index <- match(nc_samples, colnames(Y))
norm_index <- norm_index[!is.na(norm_index)]
if (length(norm_index) < 3) warning("Fewer than 3 normals found in Y; check names.")

## --- 5) Normalization -> ploidy -> SCOPE ---
gc_qc <- as.numeric(mcols(gr)$gc)
norm0 <- normalize_codex2_ns_noK(Y_qc = Y, gc_qc = gc_qc, norm_index = norm_index)

ploidy <- initialize_ploidy(Y = Y, Yhat = norm0$Yhat, ref = gr)
normObj <- normalize_scope_foreach(
   Y_qc = Y, gc_qc = gc_qc, K = 1, ploidyInt = ploidy,
   norm_index = norm_index, T = 1:5, beta0 = norm0$beta.hat, nCores = 2
)
 
pick <- which.max(normObj$BIC)
Yhat <- normObj$Yhat[[pick]]

#Yhat <- norm0$Yhat

chr_levels <- paste0("chr", 1:22)

## --- 6) Cross-sample segmentation per chr -> integer CN ---
suppressPackageStartupMessages({ library(parallel)})
# use all CPUs Slurm gave this task (fallback=4 locally)
ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
options(mc.cores = ncores)

seglist <- mclapply(chr_levels, function(ch){
  message("\n== ", ch, " ==\n")
  segment_CBScs(
    Y = Y, Yhat = Yhat, sampname = colnames(Y),
    ref = gr, chr = ch, mode = "integer", max.ns = 1
  )
}, mc.preschedule = F)  # good when work per chr is uneven
names(seglist) <- chr_levels

iCN <- do.call(rbind, lapply(seglist, function(z){z[["iCN"]]}))

cna_scope <- iCN; ref_scope <- as.data.frame(gr)

save(cna_scope, ref_scope, file = file.path(CNA_DIR,"Calling/SCOPE_test.RData"))

