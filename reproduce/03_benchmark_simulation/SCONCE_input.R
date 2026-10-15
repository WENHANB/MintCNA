setwd("<MINTCNA_ROOT>/")

suppressPackageStartupMessages({
  library(dplyr); library(GenomicRanges);library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
BASE <- if (length(args) >= 1) args[1] else "100kb_c2_N50_1"
print(BASE)
CNA_DIR <- file.path("<MINTCNA_ROOT>/Simulation/cnasim", BASE) 
base_dir <- file.path(CNA_DIR,"Calling/SCONCE")

load(file.path(CNA_DIR,"DNA_RC.RData"))
label <- fread(file.path(CNA_DIR,"cell_types.tsv"))

nc_samples <- label[which(label$normal == "normal"), root]

seqlevelsStyle(gr) <- "UCSC"
ranges(gr) <- IRanges(start = start(gr) + 1,
                      end   = end(gr))
gr$gc <- as.numeric(mcols(gr)$gc)*100 

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


