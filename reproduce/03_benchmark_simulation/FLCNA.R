library(FLCNA)
library(data.table)
library(dplyr)
library(GenomicRanges)

args <- commandArgs(trailingOnly = TRUE)
BASE <- if (length(args) >= 1) args[1] else "100kb_c2_N50_2"
print(BASE)
CNA_DIR <- file.path("<MINTCNA_ROOT>/Simulation/cnasim", BASE)  
load(file.path(CNA_DIR,"DNA_RC.RData"))
label <- fread(file.path(CNA_DIR,"cell_types.tsv"))

nc_samples <- label[which(label$normal == "normal"), root]

seqlevelsStyle(gr) <- "UCSC"
ranges(gr) <- IRanges(start = start(gr) + 1,
                      end   = end(gr))
gr$gc <- as.numeric(mcols(gr)$gc)*100 

# # Quality Control 
QCobject <- FLCNA_QC(Y_raw=Y, ref_raw=gr,
                     mapp_thresh = 0.9,
                     gc_thresh = c(20, 80))

# Normalization
log2Rdata <- FLCNA_normalization_ref(Y=QCobject$Y, gc=QCobject$ref$gc, map=QCobject$ref$mapp,
                                     ref_id = which(colnames(Y) %in% nc_samples))

# Simultaneous CNA detection and subclone clustering
tuning <- cbind(rep(c(3,4,5),each=4), rep(c(1.5, 2.5, 5, 10),3))
output_FLCNA <- FLCNA(tuning, K=NULL, lambda=NULL, y=t(log2Rdata), ref=QCobject$ref)

print(paste0("best K:",output_FLCNA$K.best,"; best lambda:",output_FLCNA$lambda.best))
save(output_FLCNA, file = file.path(CNA_DIR,"Calling/FLCNA_all.RData"))

#load("<MINTCNA_ROOT>/Result/CRC/Benchmark/FLCNA_all.RData")

# CNA clustering
CNA.output <- CNA.out(mean.matrix = output_FLCNA$mu.hat.best, Clusters=output_FLCNA$s.hat.best,
                      LRR=log2Rdata, QC_ref=QCobject$ref, cutoff=0.80, L=100)

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

cna_flcna <- CNA.matrix.flcna; ref_flcna <- as.data.frame(QCobject$ref)
save(cna_flcna, ref_flcna, file = file.path(CNA_DIR,"Calling/FLCNA.RData"))

