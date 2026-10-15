setwd("<MINTCNA_ROOT>/")

suppressPackageStartupMessages({
  library(dplyr); library(GenomicRanges);library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
BASE <- if (length(args) >= 1) args[1] else "100kb_c2_N200_1"
print(BASE)
CNA_DIR <- file.path("<MINTCNA_ROOT>/Simulation/cnasim", BASE) 
base_dir <- file.path(CNA_DIR,"Calling/SCONCE")
load(file.path(CNA_DIR, "DNA_RC.RData"))  

############################## Read the result ##############################
k <- 5  # set 10 if you ran maxKploid=10

# master bin key in desired order
bins_df <- data.frame(
  chr   = as.character(seqnames(gr)),
  start = start(gr),
  end   = end(gr),
  stringsAsFactors = FALSE
)
bins_key <- paste0("chr",bins_df$chr, ":", bins_df$start+1, "-", bins_df$end)

# find output bed files
out_dir <- file.path(base_dir, "out")
bed_files <- list.files(out_dir, pattern = paste0("_k", k, "\\.bed$"), full.names = TRUE)
stopifnot(length(bed_files) > 0)

# cell IDs from filenames
cell_ids <- sub(paste0("_k", k, "\\.bed$"), "", basename(bed_files))

# allocate result matrix
CN <- matrix(NA_integer_, nrow = length(bins_key), ncol = length(bed_files),
             dimnames = list(bins_key, cell_ids))

# read each bed and align to master bins
for (i in seq_along(bed_files)) {
  dt <- fread(bed_files[i], header = FALSE)
  # expected: V1 chr, V2 start, V3 end, V4 copy_state
  key_i <- paste0(dt[[1]], ":", dt[[2]], "-", dt[[3]])
  idx <- match(key_i, bins_key)
  
  if (anyNA(idx)) {
    stop("Some bins in ", bed_files[i], " are not found in ref_gr bin set.")
  }
  
  CN[idx, i] <- as.integer(dt[[4]])
}

ref_sconce <- bins_df
cna_sconce <- CN + 1L; cna_sconce[cna_sconce >= 5] <- 5L 

save(cna_sconce, ref_sconce, file = file.path(CNA_DIR,"Calling/SCONCE.RData"))
