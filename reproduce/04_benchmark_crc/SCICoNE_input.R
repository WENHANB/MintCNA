setwd("<MINTCNA_ROOT>/")
load("Data/CRC/DNAmeth/dnameth_chr_all_cells_100k_gc_mapp.RData")  

library(data.table); library(GenomicRanges)

keep <- !is.na(gr$gc) & !is.na(gr$mapp) & gr$gc >= 20 & gr$gc <= 80 & gr$mapp >= 0.90
#Y      <- Y[keep, c(cell_keep, nc_samples)]
Y      <- Y[keep, ]
gr <- gr[keep]

nbins  <- nrow(Y)
ncells <- ncol(Y)
X <- t(Y)                     

out_csv <- file.path("<MINTCNA_ROOT>/Result/CRC/Benchmark/SCICoNE/input","Y_cells_by_bins.csv")
write.table(X, file=out_csv, sep=",", row.names=FALSE, col.names=FALSE, quote=FALSE)

meta <- data.frame(n_cells=ncells, n_bins=nbins)
write.table(meta, file=file.path("<MINTCNA_ROOT>/Result/CRC/Benchmark/SCICoNE/input","meta.txt"),
            sep="\\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

cat("Exported:", out_csv, "\\n")
cat("n_cells=", ncells, " n_bins=", nbins, "\\n")