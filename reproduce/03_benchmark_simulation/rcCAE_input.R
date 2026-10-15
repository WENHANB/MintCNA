library(data.table); library(GenomicRanges)

args <- commandArgs(trailingOnly = TRUE)
BASE <- if (length(args) >= 1) args[1] else "c3_N100_100kb_l1e4"
print(BASE)

CNA_DIR <- file.path("<MINTCNA_ROOT>/Simulation/cnasim", BASE)
load(file.path(CNA_DIR, "DNA_RC.RData"))   # expects gr and Y

BIN_SIZE <- 100000L  # 100kb

############## Save data into rcCAE input format ######################

# chr: match rcCAE example (no "chr" prefix)
chr_vec <- as.numeric(seqnames(gr))
chr_vec <- sub("^chr", "", chr_vec, ignore.case = TRUE)

dt <- data.table(
  chr = as.integer(chr_vec),
  `gc-content` = as.numeric(gr$gc),
  mappability  = as.numeric(gr$mapp)
)

Ydt <- as.data.table(Y)
setnames(Ydt, paste0("cell", seq_len(ncol(Ydt))))
dt <- cbind(dt, Ydt)

out_dir  <- file.path(CNA_DIR, "Calling/rcCAE")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_file <- file.path(out_dir, "rcCAE_input.txt")

# write comma-separated content into .txt, with the exact first line format
tmp <- tempfile()
fwrite(dt, tmp, sep = ",", col.names = TRUE)
writeLines(c(paste0("#bin size = ", BIN_SIZE), readLines(tmp)), out_file)
unlink(tmp)
cat("Wrote:", out_file, "\n")


