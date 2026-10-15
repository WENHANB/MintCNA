library(data.table); library(GenomicRanges)

args <- commandArgs(trailingOnly = TRUE)
BASE <- if (length(args) >= 1) args[1] else "c3_N100_100kb_l1e4"
print(BASE)

CNA_DIR <- file.path("<MINTCNA_ROOT>/Simulation/cnasim", BASE)
load(file.path(CNA_DIR, "DNA_RC.RData"))  

nbins  <- nrow(Y)
ncells <- ncol(Y)
X <- t(Y)                     

out_csv <- file.path(CNA_DIR, "Calling/SCICoNE/input","Y_cells_by_bins.csv")
write.table(X, file=out_csv, sep=",", row.names=FALSE, col.names=FALSE, quote=FALSE)

meta <- data.frame(n_cells=ncells, n_bins=nbins)
write.table(meta, file=file.path(CNA_DIR, "Calling/SCICoNE/input","meta.txt"),
            sep="\\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

cat("Exported:", out_csv, "\\n")
cat("n_cells=", ncells, " n_bins=", nbins, "\\n")