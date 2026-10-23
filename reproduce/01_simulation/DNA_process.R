library(data.table)

args <- commandArgs(trailingOnly = TRUE)
BASE <- if (length(args) >= 1) args[1] else "c3_N100_100kb_l1e4"
print(BASE)
CNA_DIR <- file.path("<MINTCNA_ROOT>/Simulation/cnasim", BASE)  

# ----------------------------
# 0) Read CNAsim readcounts.tsv (long format)
# ----------------------------
rc <- fread(file.path(CNA_DIR,"readcounts.tsv"))
label <- fread(file.path(CNA_DIR,"cell_types.tsv"))

# Standardize column names just in case
setnames(rc, tolower(names(rc)))
# expected: cell, chrom, start, end, acount, bcount

# ----------------------------
# 1) Create bin_id + total counts, then bin x cell matrix
# ----------------------------
rc[, bin_id := paste0(chrom, ":", start, "-", end)]
rc[, total := acount + bcount]

# aggregate safely (in case duplicates exist)
rc_sum <- rc[, .(total = sum(total)), by = .(bin_id, cell)]

# wide: rows=bins, cols=cells
RC_wide <- dcast(rc_sum, bin_id ~ cell, value.var = "total", fill = 0)

RCmat <- as.matrix(RC_wide[, -1])
rownames(RCmat) <- RC_wide$bin_id

# ----------------------------
# 2) Create ref dataframe (chr, start, end) aligned to RCmat rows
# ----------------------------
ref <- data.table(bin_id = rownames(RCmat))
ref[, c("chr", "pos") := tstrsplit(bin_id, ":", fixed = TRUE)]
ref[, c("start", "end") := tstrsplit(pos, "-", fixed = TRUE)]
ref[, `:=`(start = as.integer(start), end = as.integer(end))]
ref[, pos := NULL]

# ----------------------------
# 3) Reorder rows: chr1, chr2, ..., chr22, chrX, chrY; within chr by start/end
# ----------------------------
chr_levels <- c(paste0("chr", 1:22), "chrX", "chrY")
ref[, chr := factor(chr, levels = chr_levels)]
setorder(ref, chr, start, end)

RCmat <- RCmat[ref$bin_id, , drop = FALSE]

# keep only chr/start/end (same order as RCmat rows)
ref <- ref[, .(chr = as.character(chr), start, end)]

# ----------------------------
# 4) Reorder columns: cell1, cell2, ..., cell10, ...
# ----------------------------
cell_ids <- colnames(RCmat)
cell_num <- as.integer(sub("^cell", "", cell_ids))
RCmat <- RCmat[, order(cell_num), drop = FALSE]

# Outputs:
# - RCmat: bin x cell read count matrix (ordered)
# - ref  : data.table with chr/start/end aligned to RCmat rows
dim(RCmat)

######################################## Process the data ######################################## 
## CODEX2
library(CODEX2)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
ref <- GRanges(seqnames = gsub("chr", "", ref$chr),
               ranges = IRanges(start = ref$start,ref$end))
gc_codex2 <- getgc(ref,genome = BSgenome.Hsapiens.UCSC.hg38)
mapp_codex2 <- getmapp(ref,genome = BSgenome.Hsapiens.UCSC.hg38)

ref$gc_content <- gc_codex2*0.01
ref$mapp_score <- mapp_codex2

gc   <- mcols(ref)$gc_content ; mapp <- mcols(ref)$mapp_score

idx <- which(!is.na(gc) & mapp > 0.9 & gc > 0.2 & gc < 0.8)

reference_info <- data.frame(
  chr   = as.character(seqnames(ref)),
  start = start(ref),
  end   = end(ref),
  gc    = as.numeric(mcols(ref)$gc_content),
  mapp  = as.numeric(mcols(ref)$mapp_score),
  stringsAsFactors = FALSE)[idx,]
read_count_matrix <- RCmat[idx,]

Y <- read_count_matrix; gr <- GRanges(reference_info)
save(gr, Y, file = file.path(CNA_DIR,"DNA_RC.RData"))

### Two step normalization by GC and mappability
source("<REPRODUCE_DIR>/lib/MultiCNA/Normalization.functions.R")
read_count_matrix_norm_1 <- GC.MAP.normalization(Y = read_count_matrix,
                                                 ref_qc = reference_info)
### library size norm
read_count_matrix_norm <- libsize_normalize(read_count_matrix_norm_1$RC_norm, method="cpm")$X_norm
# read_count_matrix_norm <-Seurat::NormalizeData(read_count_matrix_norm_1$RC_norm, 
#                                                normalization.method = "RC",scale.factor = 1000000)
write.csv(as.matrix(read_count_matrix_norm), file = file.path(CNA_DIR,"RC_DNA.csv"))
# L2R calculation and process
nc_samples <- label[which(label$normal == "normal"), root]
#nc_mean <- rowMeans(as.matrix(read_count_matrix_norm[, nc_samples]), na.rm = TRUE)
nc_median <- apply(as.matrix(read_count_matrix_norm[, nc_samples]),1,median)
log2R_matrix <- apply(as.matrix(read_count_matrix_norm[, setdiff(colnames(read_count_matrix_norm), nc_samples)]), 
                      2, function(sample_counts) {log2(((sample_counts+0.01) / (nc_median+0.01)))})
rownames(log2R_matrix) <- rownames(read_count_matrix_norm)
colnames(log2R_matrix) <- setdiff(colnames(read_count_matrix_norm), nc_samples)

## Median shift normalization
log2R_matrix <- apply(log2R_matrix, 2, function(col) {col - median(col, na.rm = TRUE)})

## smooth() the l2r
log2R_matrix_smoothed <- modSaRa2::smooth(log2R_matrix,R = 10)
rownames(log2R_matrix_smoothed) <- rownames(log2R_matrix)
colnames(log2R_matrix_smoothed) <- colnames(log2R_matrix)

# apply(log2R_matrix_smoothed, 2, median)
# apply(log2R_matrix_smoothed, 2, mean)
# apply(log2R_matrix_smoothed, 2, sd)

## Order by genomic coordinate
ord <- with(reference_info, order(as.numeric(sub("chr", "", chr)), start))
log2R_matrix_smoothed <- log2R_matrix_smoothed[ord,]
reference_info <- reference_info[ord,]

## Save data on all chromosome
write.csv(log2R_matrix_smoothed, file = file.path(CNA_DIR,"l2r_smoothed_all_cells_100kb.csv"))
write.csv(reference_info, file = file.path(CNA_DIR,"ref_all_cells_100kb.csv"))

