setwd("<MINTCNA_ROOT>/")
source("<REPRODUCE_DIR>/04_benchmark_crc/SCCNV_Fun.R")
load("Data/CRC/DNAmeth/dnameth_chr_all_cells_100k_gc_mapp.RData")  
library(dplyr); library(GenomicRanges)

#nc_samples <- c("CRC01_NC_507","CRC01_NC_509","CRC01_NC_517",
#                "CRC01_NC_521","CRC01_NC_525","CRC01_NC_527","CRC01_NC_529")

#Y <- Y[,-which(colnames(Y)%in%nc_samples)]

chr   <- as.character(GenomicRanges::seqnames(gr))
start <- BiocGenerics::start(gr)
end   <- BiocGenerics::end(gr)

gc   <- S4Vectors::mcols(gr)$gc      # your GC column (looks like 0–100%)
mapp <- S4Vectors::mcols(gr)$mapp    # mappability (0–1)

# You do NOT currently have an "N fraction" column; set to 0 so nothing is filtered out by maxN
Nfrac <- rep(0, length(gr))

# SCCNV assumes:
# col1=chr, col2=mapp filter, col4=normalization divisor (they divide RC by col4), col6=gc, col7=Nfrac
dat1 <- data.frame(
  chr  = chr,        # col1
  mapp = mapp,       # col2
  start= start,      # col3 (not used by SCCNV, but fine)
  map_norm = pmax(mapp, 1e-6),  # col4; avoid division by 0
  end  = end,        # col5 (not used)
  gc   = gc,         # col6
  N    = Nfrac,      # col7
  stringsAsFactors = FALSE
)

# Add per-cell read counts as columns 8+
dat1 <- cbind(dat1, as.data.frame(Y, check.names = FALSE))

keep <- !is.na(dat1$gc) & !is.na(dat1$mapp) &
  dat1$gc >= 20 & dat1$gc <= 80 &
  dat1$mapp >= 0.90

dat1_f <- dat1[keep, , drop = FALSE]
cat("Before:", nrow(dat1), "After:", nrow(dat1_f), "\n")

binsize <- as.integer(stats::median(width(gr)))   # should be 100000
post_list <- SCCNV(dat1_f,
                   minmap = 0.90,   # consistent with your filter
                   maxN   = 0.1,
                   ploidy = 2,
                   windowsize = 20,
                   binsize = binsize)

# post_list: list of length = #cells
cn_mat <- sapply(post_list, function(pp) {
  as.integer(max.col(as.matrix(pp), ties.method = "first") - 1)  # CN0..CN20
})

# Make sure it is a matrix even if only 1 cell
cna_sccnv <- as.matrix(cn_mat)
colnames(cna_sccnv) <- names(post_list) %||% colnames(dat1_f)[8:ncol(dat1_f)]  # see note below


## Save the result
ref_sccnv <- as.data.frame(gr[keep]); colnames(ref_sccnv)[1] <- "chr";
cna_sccnv <- cna_sccnv + 1L; cna_sccnv[cna_sccnv >= 5] <- 5L 

save(cna_sccnv, ref_sccnv, file = "<MINTCNA_ROOT>/Result/CRC/Benchmark/SCCNV.RData")

source("<REPRODUCE_DIR>/lib/MultiCNA/CNAout.R")
code_map5 <- c(del = 1, loss = 2, nor = 3, gain = 4, amp = 5)
plot_S_heatmap(cna_sccnv, ref_sccnv, code_map = code_map5, 
               file = "Result/CRC/Segmentation/CNV/SCCNV.pdf",
               width = 12, height = 20)

