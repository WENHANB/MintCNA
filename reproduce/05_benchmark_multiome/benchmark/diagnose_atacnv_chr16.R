#!/usr/bin/env Rscript
# Why do P109T1 and P111T3 die inside AtaCNV, and what is the right fix?
#   Usage: Rscript diagnose_atacnv_chr16.R <care_id>
#
# Both fail identically:
#   Error, "read count" is negative (tumor_2543=99.0472, normal_2543 = -16.8362)
# thrown by the precompiled MBICseq binary, which then leaves no .bin_seg for R
# to read. The suspicion is that AtaCNV's LOESS GC correction drives some bins'
# NORMAL BASELINE negative on chr16, the most GC-rich autosome.
#
# This script stops one step short of segmentation: it bins the fragments, runs
# normalize(), and reports exactly which bins go non-positive and by how much.
# The normalize() output is saved so the fix can be applied without re-reading
# a multi-gigabyte fragments file.
#
# Two fixes are on the table and the numbers decide between them:
#   clamp  -- floor the baseline at a small positive value. Keeps the bin grid
#             identical across samples, but a clamped bin yields a huge copy
#             ratio and risks a spurious amplification call.
#   drop   -- exclude the offending bins. Honest (a bin whose GC correction went
#             negative is not trustworthy) but makes the grid sample-dependent.
# If only a handful of bins are affected, clamp-and-record is fine. If it is a
# contiguous stretch, dropping is the only defensible choice.

args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("need <care_id>")
CARE <- args[1]

suppressPackageStartupMessages({library(data.table); library(AtaCNV)})

ROOT <- "<MINTCNA_ROOT>/ATAC_RNA"
OUT  <- file.path(ROOT, "Result/benchmark/AtaCNV", CARE)
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

data("bin_info_hg38", package = "AtaCNV")
bins <- as.data.table(bin_info_hg38)
bins[, bin_id := .I]

cells <- fread(file.path(ROOT, "Data/cell_lists", paste0(CARE, "_cells.tsv")))

# ---- bin the fragments (identical to run_atacnv.R) -------------------------
frag_file <- file.path(ROOT, "Data/atac_filtered", paste0(CARE, "_fragments.tsv.gz"))
message("reading fragments ...")
fr <- fread(cmd = paste("zcat", shQuote(frag_file)), header = FALSE,
            select = c(1, 2, 3, 4), col.names = c("chr", "start", "end", "barcode"))
fr[, mid := (start + end) %/% 2L]
fr[, bin_start := as.numeric((mid %/% 1000000L) * 1000000L)]   # bin_info_hg38 is 0-based
fr <- fr[chr %in% bins$chr]

bmap <- bins[, .(chr, start, bin_id)]
fr <- bmap[fr, on = .(chr, start = bin_start), nomatch = 0]

cnt <- fr[, .N, by = .(barcode, bin_id)]
cell_ids <- sort(unique(cnt$barcode))
mat <- matrix(0L, nrow = length(cell_ids), ncol = nrow(bins), dimnames = list(cell_ids, NULL))
mat[cbind(match(cnt$barcode, cell_ids), cnt$bin_id)] <- cnt$N
is_ref <- rownames(mat) %in% cells[is_reference == TRUE, barcode]
message(sprintf("matrix: %d cells x %d bins, %d reference nuclei",
                nrow(mat), ncol(mat), sum(is_ref)))

# ---- normalize only, then stop and look ------------------------------------
scratch <- file.path(OUT, "diag_scratch"); dir.create(scratch, showWarnings = FALSE)
old_wd <- setwd(scratch); on.exit(setwd(old_wd), add = TRUE)

norm_re <- AtaCNV::normalize(mat, genome = "hg38", mode = "normal cells",
                             normal_cells = is_ref, gc_correction = TRUE)
setwd(old_wd)
saveRDS(list(care_id = CARE, norm = norm_re, bins = bins, is_ref = is_ref),
        file.path(OUT, "normalize_only.rds"))

b <- norm_re$baseline
cat("\n================ baseline diagnosis:", CARE, "================\n")
cat(sprintf("baseline dim: %s\n", paste(dim(b), collapse = " x ")))
cat(sprintf("range: %.3f to %.3f\n", min(b, na.rm = TRUE), max(b, na.rm = TRUE)))
cat(sprintf("non-positive entries: %d of %d (%.4f%%)\n",
            sum(b <= 0, na.rm = TRUE), length(b), 100 * mean(b <= 0, na.rm = TRUE)))

# AtaCNV drops bins during normalize, so map the surviving columns back to the grid
kept <- if (!is.null(colnames(b))) as.integer(colnames(b)) else seq_len(ncol(b))
bad_col <- which(apply(b, 2, function(v) any(v <= 0, na.rm = TRUE)))
cat(sprintf("\nbins containing at least one non-positive value: %d of %d\n",
            length(bad_col), ncol(b)))
if (length(bad_col)) {
  bb <- bins[kept[bad_col]]
  bb[, worst := apply(b[, bad_col, drop = FALSE], 2, min, na.rm = TRUE)]
  bb[, n_cells_bad := apply(b[, bad_col, drop = FALSE], 2, function(v) sum(v <= 0, na.rm = TRUE))]
  cat("\nby chromosome:\n"); print(bb[, .(n_bins = .N, worst = min(worst)), by = chr])
  cat("\nthe offending bins:\n")
  print(bb[, .(chr, start, end, gc = round(gc, 3), worst = round(worst, 2), n_cells_bad)], nrow = 60)
  cat(sprintf("\ncontiguous? bin_id span %d..%d over %d bins\n",
              min(kept[bad_col]), max(kept[bad_col]), length(bad_col)))
}
cat("\nsaved normalize output to", file.path(OUT, "normalize_only.rds"), "\n")
