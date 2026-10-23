setwd("<MINTCNA_ROOT>/")
load("Data/CRC/DNAmeth/dnameth_chr_all_cells_100k_gc_mapp.RData")  

library(data.table); library(GenomicRanges)

keep <- !is.na(gr$gc) & !is.na(gr$mapp) & gr$gc >= 20 & gr$gc <= 80 & gr$mapp >= 0.90
#Y      <- Y[keep, c(cell_keep, nc_samples)]
Y      <- Y[keep, ]
gr <- gr[keep]

BIN_SIZE <- 100000L  # 100kb

############## Save data into rcCAE input format ######################

# chr: match rcCAE example (no "chr" prefix)
chr_vec <- as.numeric(seqnames(gr))
chr_vec <- sub("^chr", "", chr_vec, ignore.case = TRUE)

dt <- data.table(
  chr = as.integer(chr_vec),
  `gc-content` = as.numeric(gr$gc)*0.01,
  mappability  = as.numeric(gr$mapp)
)

Ydt <- as.data.table(Y)
dt <- cbind(dt, Ydt)

out_file <- "<MINTCNA_ROOT>/Result/CRC/Benchmark/rcCAE/rcCAE_input.txt"

# write comma-separated content into .txt, with the exact first line format
tmp <- tempfile()
fwrite(dt, tmp, sep = ",", col.names = TRUE)
writeLines(c(paste0("#bin size = ", BIN_SIZE), readLines(tmp)), out_file)
unlink(tmp)
cat("Wrote:", out_file, "\n")


