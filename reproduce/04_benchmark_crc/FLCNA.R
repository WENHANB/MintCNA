library(FLCNA)
library(data.table)
library(dplyr)
library(GenomicRanges)

setwd("<MINTCNA_ROOT>/")
load("Data/CRC/DNAmeth/dnameth_chr_all_cells_100k_gc_mapp.RData") 
load("Data/CRC/DNAmeth/cell_keep.RData")

nc_samples <- c("CRC01_NC_507","CRC01_NC_509","CRC01_NC_517",
                "CRC01_NC_521","CRC01_NC_525","CRC01_NC_527","CRC01_NC_529")

Y <- Y[,c(cell_keep, nc_samples)]

# # Quality Control 
QCobject <- FLCNA_QC(Y_raw=Y, ref_raw=gr,
                     mapp_thresh = 0.9,
                     gc_thresh = c(20, 80))

# # Normalization
# log2Rdata <- FLCNA_normalization_ref(Y=QCobject$Y, gc=QCobject$ref$gc, map=QCobject$ref$mapp, 
#                                      ref_id = which(colnames(Y) %in% nc_samples))
# 
# # Simultaneous CNA detection and subclone clustering
# tuning <- cbind(rep(c(2,3,4),each=4), rep(c(1.5, 2.5, 5, 10),3))
# output_FLCNA <- FLCNA(tuning, K=NULL, lambda=NULL, y=t(log2Rdata), ref=QCobject$ref)
# 
# print(paste0("best K:",output_FLCNA$K.best,"; best lambda:",output_FLCNA$lambda.best))
# save(output_FLCNA, file = "<MINTCNA_ROOT>/Result/CRC/Benchmark/FLCNA_all.RData")
#load("<MINTCNA_ROOT>/Result/CRC/Benchmark/FLCNA_all.RData")

## New try
ref <- read.csv("./Data/CRC/DNAmeth/processed_l2r/ref_all_cells_100kb.csv",row.names = 1)
l2r <- read.csv("./Data/CRC/DNAmeth/processed_l2r/l2r_smoothed_all_cells_100kb.csv",row.names = 1)
l2r <- l2r[,cell_keep]

ord <- with(ref, order(as.numeric(sub("chr", "", chr)), start))
l2r <- l2r[ord,]
ref <- ref[ord,]

tuning <- cbind(rep(c(2,3,4),each=4), rep(c(1.5, 2.5, 5, 10),3))
output_FLCNA <- FLCNA(tuning, K=NULL, lambda=NULL, y=t(l2r), ref=QCobject$ref)

print(paste0("best K:",output_FLCNA$K.best,"; best lambda:",output_FLCNA$lambda.best))
save(output_FLCNA, file = "<MINTCNA_ROOT>/Result/CRC/Benchmark/FLCNA_l2r.RData")

# CNA clustering
CNA.output <- CNA.out(mean.matrix = output_FLCNA$mu.hat.best, Clusters=output_FLCNA$s.hat.best,
                      LRR=log2Rdata, QC_ref=QCobject$ref, cutoff=0.95, L=100)

CNA.matrix.flcna <- matrix(rep(3, nrow(log2Rdata)*ncol(log2Rdata)), nrow(log2Rdata), ncol(log2Rdata))
colnames(CNA.matrix.flcna) <- colnames(log2Rdata)
rownames(CNA.matrix.flcna) <- rownames(log2Rdata)
for (i in 1:ncol(log2Rdata)){
  # clu_index <- unique(CNA.outdata$Cluster[CNA.outdata$sampleID==i])
  if (sum(CNA.output$sampleID==i)==0) {next}
  CNA <- CNA.output[CNA.output$sampleID==i,]
  for (j in 1:nrow(CNA)){
    start=CNA$start[j]
    end=CNA$end[j]
    state=CNA$state[j]
    if (state=="Del.d") {state=1}
    if (state=="Del.s") {state=2}
    if (state=="Dup.s") {state=4}
    if (state=="Dup.d") {state=5}
    CNA.matrix.flcna[start:end,i] <- state
  }
}

## Save CNA output and Heatmap
cna_flcna <- CNA.matrix.flcna; ref_flcna <- as.data.frame(QCobject$ref)
save(cna_flcna, ref_flcna, file = "<MINTCNA_ROOT>/Result/CRC/Benchmark/FLCNA.RData")

source("<REPRODUCE_DIR>/lib/MultiCNA/CNAout.R")
code_map5 <- c(del = 1, loss = 2, nor = 3, gain = 4, amp = 5); colnames(ref_flcna)[1] <- "chr"
plot_S_heatmap(CNA.matrix.flcna, ref_flcna, code_map = code_map5,
               file = "Result/CRC/Segmentation/CNV/FLCNA_l2r.pdf",
               width = 12, height = 20)
