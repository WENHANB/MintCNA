library(data.table)

setwd("<MINTCNA_ROOT>/")
load("Data/CRC/DNAmeth/dnameth_chr_all_cells_100k_gc_mapp.RData")  

cna_scicone <- t(fread("<MINTCNA_ROOT>/Result/CRC/Benchmark/SCICoNE/run1_inferred_cnvs.csv", header = F))
cna_scicone <- cna_scicone + 1
colnames(cna_scicone) <- colnames(Y)
save(cna_scicone, file = "<MINTCNA_ROOT>/Result/CRC/Benchmark/SCICoNE.RData")
table(cna_scicone)