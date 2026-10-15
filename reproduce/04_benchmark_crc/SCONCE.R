setwd("<MINTCNA_ROOT>/")
load("Data/CRC/DNAmeth/dnameth_chr_all_cells_100k_gc_mapp.RData")  
load("Data/CRC/DNAmeth/cell_keep.RData")

suppressPackageStartupMessages({
  library(dplyr); library(GenomicRanges);library(data.table)
})

base_dir <- "<MINTCNA_ROOT>/Result/CRC/Benchmark/SCONCE"

nc_samples <- c("CRC01_NC_507","CRC01_NC_509","CRC01_NC_517",
                "CRC01_NC_521","CRC01_NC_525","CRC01_NC_527","CRC01_NC_529")


## --- 3) Bin QC (GC/mapp) + rebuild filtered GRanges ---
keep <- !is.na(gr$gc) & !is.na(gr$mapp) & gr$gc >= 20 & gr$gc <= 80 & gr$mapp >= 0.90
#Y      <- Y[keep, c(cell_keep, nc_samples)]
Y      <- Y[keep, ]
ref_gr <- gr[keep]

## Store as bed for SCONCE
## ---------- REQUIRED OBJECTS IN MEMORY ----------
## You must already have:
##   Y      : matrix (nBins x nCells), colnames(Y) are cell IDs
##   ref_gr : GRanges length nBins with bin coordinates
stopifnot(exists("Y"), exists("ref_gr"))

## ---------- CREATE FOLDERS ----------
dir.create(file.path(base_dir, "in"),  recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(base_dir, "ref"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(base_dir, "out"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(base_dir, "logs"),recursive = TRUE, showWarnings = FALSE)

## ---------- SANITY CHECKS ----------
stopifnot(!is.null(colnames(Y)))
stopifnot(nrow(Y) == length(ref_gr))

cell_ids <- colnames(Y)

missing_nc <- setdiff(nc_samples, cell_ids)
if (length(missing_nc) > 0) {
  stop("These nc_samples are not found in colnames(Y): ",
       paste(missing_nc, collapse = ", "))
}

labels <- ifelse(cell_ids %in% nc_samples, "normal", "tumor")
message("Label counts:")
print(table(labels))

## ---------- BUILD BINS FROM GRanges ----------
bins <- data.frame(
  chr   = as.character(GenomicRanges::seqnames(ref_gr)),
  start = BiocGenerics::start(ref_gr),
  end   = BiocGenerics::end(ref_gr)
)

## ---------- WRITE ONE BED PER CELL ----------
in_dir <- file.path(base_dir, "in")
bed_files <- character(ncol(Y))
names(bed_files) <- cell_ids

for (j in seq_len(ncol(Y))) {
  out_fn <- file.path(in_dir, paste0(cell_ids[j], ".bed"))
  df <- cbind(bins, readDepth = as.numeric(Y[, j]))
  write.table(df, out_fn, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE)
  bed_files[j] <- out_fn
}

## ---------- WRITE FILE LISTS ----------
diploid_beds <- bed_files[labels == "normal"]
tumor_beds   <- bed_files[labels == "tumor"]

writeLines(diploid_beds, file.path(base_dir, "ref", "diploidFileList.txt"))
writeLines(tumor_beds,   file.path(base_dir, "tumorFileList.txt"))

## Optional: save a small manifest
manifest <- data.frame(cell_id = cell_ids, label = labels, bed = bed_files, stringsAsFactors = FALSE)
write.csv(manifest, file.path(base_dir, "cell_manifest.csv"), row.names = FALSE)

message("Done. Wrote:")
message("  ", length(diploid_beds), " diploid BEDs  -> ", file.path(base_dir, "ref", "diploidFileList.txt"))
message("  ", length(tumor_beds),   " tumor BEDs    -> ", file.path(base_dir, "tumorFileList.txt"))
message("  BEDs are in             -> ", in_dir)

## Read the result
base_dir <- "./Result/CRC/Benchmark/SCONCE"
k <- 5  # set 10 if you ran maxKploid=10

# ref_gr must exist and define the bin order you want
stopifnot(exists("ref_gr"))

# master bin key in desired order
bins_df <- data.frame(
  chr   = as.character(seqnames(ref_gr)),
  start = start(ref_gr),
  end   = end(ref_gr),
  stringsAsFactors = FALSE
)
bins_key <- paste0(bins_df$chr, ":", bins_df$start, "-", bins_df$end)

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


ref_sconce <- as.data.frame(gr[keep]); colnames(ref_sconce)[1] <- "chr"
cna_sconce <- CN + 1L; cna_sconce[cna_sconce >= 5] <- 5L 

save(cna_sconce, ref_sconce, file = "<MINTCNA_ROOT>/Result/CRC/Benchmark/SCONCE.RData")

source("<MINTCNA_ROOT>/Code/MultiCNA/CNAout.R")
code_map5 <- c(del = 1, loss = 2, nor = 3, gain = 4, amp = 5)
plot_S_heatmap(cna_sconce, ref_sconce, code_map = code_map5, 
               file = "Result/CRC/Segmentation/CNV/SCONCE.pdf",
               width = 12, height = 20)
