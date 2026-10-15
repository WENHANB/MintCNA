library(tidyverse)
library(GenomicRanges)
library(data.table)
library(HMMcopy)

args <- commandArgs(trailingOnly = TRUE)
BASE <- if (length(args) >= 1) args[1] else "c3_N100_100kb_l1e4"
print(BASE)
CNA_DIR <- file.path("<MINTCNA_ROOT>/Simulation/cnasim", BASE)  
load(file.path(CNA_DIR,"DNA_RC.RData"))
label <- fread(file.path(CNA_DIR,"cell_types.tsv"))
ref <- gr

copy.output <- NULL
for (i in 1:ncol(Y)){
  normal_reads1 <- data.frame(chr=rep(ref@seqnames@values,ref@seqnames@lengths),
                              start=ref@ranges@start,
                              end=ref@ranges@start+ref@ranges@width,
                              reads=Y[,i],
                              gc=ref$gc,
                              map=ref$mapp)
  #normal_reads1$chr <- substr(normal_reads1$chr,4,nchar(as.character(normal_reads1$chr)))
  normal_copy <- correctReadcount(normal_reads1)
  
  default_param <- HMMsegment(normal_copy, getparam = TRUE)
  longseg_param <- default_param
  longseg_param$e <- 0.9999999
  longseg_param$strength <- 1e7 # strength higher, CNA less 
  longseg_segments <- HMMsegment(normal_copy, longseg_param, verbose = FALSE)
  table(longseg_segments$state)
  copy.output <- cbind(copy.output, longseg_segments$state)
}
colnames(copy.output) <- colnames(Y)

## convert segments -> bin×cell integer CN matrix
cna_hmmcopy <- copy.output
ref_hmmcopy <- gr

dir.create(file.path(CNA_DIR, "Calling"), showWarnings = FALSE, recursive = TRUE)
save(cna_hmmcopy, ref_hmmcopy, file = file.path(CNA_DIR, "Calling/HMMcopy.RData"))
