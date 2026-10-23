library(data.table)

CNA_DIR <- "<MINTCNA_ROOT>/Result/CRC/Benchmark/rcCAE"

setwd("<MINTCNA_ROOT>/")
load("Data/CRC/DNAmeth/dnameth_chr_all_cells_100k_gc_mapp.RData")  

############## CNA result from rcCAE analysis ######################
CN <- fread(file.path(CNA_DIR, "out","copynumber.csv"))
P     <- fread(file.path(CNA_DIR, "out", "ploidy.csv"))

# copynumber: rows=cells, cols=bins
CN[, cell := paste0("cell", .I)]
setcolorder(CN, c("cell", setdiff(names(CN), "cell")))
bin_cols <- setdiff(names(CN), "cell")

# ploidy: typically 1 row with V1..Vn = per-cell ploidy

ploidy_vec <- as.numeric(P[1, ])
# baseline CN per cell from ploidy
base_ploidy <- as.integer(round(ploidy_vec))

# sanity check
n_cells <- nrow(CN)
if (length(base_ploidy) != n_cells) {
  warning(sprintf("ploidy length (%d) != #cells in copynumber (%d). Falling back to mode baseline.",
                  length(base_ploidy), n_cells))
  mode_int <- function(x) {
    tx <- table(x)
    as.integer(names(tx)[which.max(tx)])
  }
  base_ploidy <- CN[, mode_int(unlist(.SD)), .SDcols=bin_cols]
}

# compute del/nor/dup relative to baseline (round(ploidy))
cn_mat <- as.matrix(CN[, ..bin_cols])
delta  <- cn_mat - base_ploidy

# encode del/nor/dup as 1/2/3
cna_code <- matrix(2L, nrow=nrow(delta), ncol=ncol(delta))  # normal = 2
cna_code[delta <= -1] <- 1L
cna_code[delta >=  1] <- 3L

# output: bins as rows, cells as columns (cellXX)
cna_rccae <- (t(cna_code))
colnames(cna_rccae) <- colnames(Y)

save(cna_rccae, file = "<MINTCNA_ROOT>/Result/CRC/Benchmark/rcCAE.RData")

plot_S_heatmap(cna_rccae, ref_multicna_ae, code_map = code_map3, 
               file = "Result/CRC/Segmentation/CNV/rcCAE.pdf",
               width = 12, height = 20)
