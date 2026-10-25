#!/usr/bin/env Rscript
# AtaCNA (R package AtaCNV 0.1.0) on one multiome sample, package defaults: 1 Mb bins, GC correction,
# normalisation against the reference nuclei ("normal cells" mode), MBICseq segmentation.
# The package's fragment parser is replaced by an equivalent data.table binning (fragments assigned to the
# package's hg38 1 Mb bins by midpoint); samples with fewer than 100 reference nuclei borrow the reference
# nuclei of the same patient's other time point.
#   Rscript reproduce/analysis/multiome/methods/run_atacnv.R <care_id>
# Input : ATAC_RNA/Data/atac_filtered/<care_id>_fragments.tsv.gz, ATAC_RNA/Data/cell_lists/<care_id>_cells.tsv
# Output: ATAC_RNA/Result/benchmark/AtaCNV/<care_id>/atacnv_result.rds
# Run from the root of the MintCNA repository with MINTCNA_DATA set (see reproduce/README.md).
suppressPackageStartupMessages({ library(data.table); library(AtaCNV) })
ROOT <- Sys.getenv("MINTCNA_DATA"); if (!nzchar(ROOT)) stop("set MINTCNA_DATA (see reproduce/README.md)")
args <- commandArgs(trailingOnly = TRUE); if (!length(args)) stop("need <care_id>")
CARE <- args[1]
ROOT <- file.path(ROOT, "ATAC_RNA")
OUT  <- file.path(ROOT, "Result/benchmark/AtaCNV", CARE); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

data("bin_info_hg38", package = "AtaCNV")
bins <- as.data.table(bin_info_hg38); bins[, bin_id := .I]
bmap <- bins[, .(chr, start, bin_id)]; setkey(bmap, chr, start)

bin_fragments <- function(care, barcodes = NULL) {
  fr <- fread(cmd = paste("zcat", shQuote(file.path(ROOT, "Data/atac_filtered", paste0(care, "_fragments.tsv.gz")))),
              header = FALSE, select = 1:4, col.names = c("chr", "start", "end", "barcode"))
  if (!is.null(barcodes)) fr <- fr[barcode %in% barcodes]
  fr <- fr[chr %in% bins$chr]
  fr[, bin_start := as.numeric((((start + end) %/% 2L) %/% 1000000L) * 1000000L)]   # midpoint, 0-based 1 Mb bins
  fr <- bmap[fr, on = .(chr, start = bin_start), nomatch = 0]
  cnt <- fr[, .N, by = .(barcode, bin_id)]
  ids <- sort(unique(cnt$barcode))
  mat <- matrix(0L, nrow = length(ids), ncol = nrow(bins), dimnames = list(ids, NULL))
  mat[cbind(match(cnt$barcode, ids), cnt$bin_id)] <- cnt$N
  mat
}

# reference nuclei: this sample's, or the same patient's other time point when fewer than 100
cells <- fread(file.path(ROOT, "Data/cell_lists", paste0(CARE, "_cells.tsv")))
ref_src <- CARE
if (sum(cells$is_reference) < 100) {
  pat <- sub("T[0-9]+$", "", CARE)
  cand <- setdiff(grep(paste0("^", pat, "T"), sub("_cells.tsv", "", list.files(file.path(ROOT, "Data/cell_lists"), pattern = "_cells.tsv$")), value = TRUE), CARE)
  n_ref <- sapply(cand, function(cc) sum(fread(file.path(ROOT, "Data/cell_lists", paste0(cc, "_cells.tsv")))$is_reference))
  if (length(cand)) { ref_src <- cand[which.max(n_ref)]; message("reference nuclei borrowed from ", ref_src) }
}

mat <- bin_fragments(CARE)
if (ref_src == CARE) {
  is_ref <- rownames(mat) %in% cells[is_reference == TRUE, barcode]
} else {
  o_cells <- fread(file.path(ROOT, "Data/cell_lists", paste0(ref_src, "_cells.tsv")))
  omat <- bin_fragments(ref_src, o_cells[is_reference == TRUE, barcode])
  rownames(omat) <- paste0("REF_", rownames(omat))
  mat <- rbind(mat, omat)
  is_ref <- grepl("^REF_", rownames(mat))
}
message(sprintf("%s: %d nuclei x %d bins, %d reference nuclei (%s)", CARE, nrow(mat), ncol(mat), sum(is_ref), ref_src))

# MBICseq writes temporary files with fixed names into the working directory: isolate each run
scratch <- file.path(OUT, "mbicseq_scratch"); dir.create(scratch, showWarnings = FALSE, recursive = TRUE)
old_wd <- setwd(scratch); on.exit(setwd(old_wd), add = TRUE)

t0 <- Sys.time()
norm_re <- AtaCNV::normalize(mat, genome = "hg38", mode = "normal cells", normal_cells = is_ref, gc_correction = TRUE)
# a non-positive baseline entry (one cell in one bin, in two samples) makes MBICseq abort: replace it with
# the median of the positive values of that bin
bad <- which(norm_re$baseline <= 0 | !is.finite(norm_re$baseline), arr.ind = TRUE)
if (nrow(bad)) {
  for (j in unique(bad[, "col"])) { v <- norm_re$baseline[, j]; ok <- is.finite(v) & v > 0; norm_re$baseline[!ok, j] <- median(v[ok]) }
  message(sprintf("repaired %d non-positive baseline entries", nrow(bad)))
}
res <- AtaCNV::calculate_CNV(norm_count = norm_re$norm_count, baseline = norm_re$baseline)

saveRDS(list(care_id = CARE, ref_source = ref_src, n_ref = sum(is_ref), bins = bins, norm = norm_re, cnv = res),
        file.path(OUT, "atacnv_result.rds"))
setwd(old_wd); unlink(scratch, recursive = TRUE)
cat(sprintf("SUMMARY\t%s\t%d\t%d\t%s\t%.1f min\n", CARE, nrow(mat), sum(is_ref), ref_src,
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
