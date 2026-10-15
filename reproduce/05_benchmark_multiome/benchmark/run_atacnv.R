#!/usr/bin/env Rscript
# AtaCNV on one sample.  Usage: Rscript run_atacnv.R <care_id>
#
# Two deviations from a naive "just call the package" run, both forced:
#
#  1. generate_input() is bypassed. It parses the fragments file in pure R, one
#     line at a time (readLines(con, n=1) inside while(TRUE)); on a 40-100M line
#     multiome fragments file that runs for days. The binning is reimplemented in
#     data.table here and produces the same cell x 1 Mb matrix the package expects.
#     Bin grid is taken from the package's own bin_info_hg38 so the columns line up.
#
#  2. Reference cells come from the authors' annotation (Myeloid + Endo/Mural +
#     Lymphocyte), which is exactly what normalize(mode="normal cells") wants.
#     Two samples fall below 100 reference nuclei and borrow from the same
#     patient's other timepoint -- recorded in the output so Methods can state it.
#
# Everything else is package default: 1 Mb bins (hard-coded K=1e6), gc_correction
# TRUE, mode "normal cells", MBICseq segmentation.

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
setkey(bins, chr, start, end)

# ---- reference cells ------------------------------------------------------
cells <- fread(file.path(ROOT, "Data/cell_lists", paste0(CARE, "_cells.tsv")))
ref_src <- CARE
if (sum(cells$is_reference) < 100) {
  # borrow the other timepoint of the same patient
  pat <- sub("T[0-9]+$", "", CARE)
  cand <- setdiff(grep(paste0("^", pat, "T"),
                       sub("_cells.tsv", "", list.files(file.path(ROOT,"Data/cell_lists"),
                                                        pattern="_cells.tsv$")), value = TRUE), CARE)
  best <- NULL; bestn <- 0
  for (cc in cand) {
    o <- fread(file.path(ROOT, "Data/cell_lists", paste0(cc, "_cells.tsv")))
    if (sum(o$is_reference) > bestn) { bestn <- sum(o$is_reference); best <- cc }
  }
  if (!is.null(best)) { ref_src <- best; message("borrowing reference from ", best) }
}

# ---- bin the fragments ----------------------------------------------------
frag_file <- file.path(ROOT, "Data/atac_filtered", paste0(CARE, "_fragments.tsv.gz"))
message("reading fragments ...")
fr <- fread(cmd = paste("zcat", shQuote(frag_file)), header = FALSE,
            select = c(1, 2, 3, 4), col.names = c("chr", "start", "end", "barcode"))
# assign by midpoint: unambiguous, and avoids double-counting a fragment that
# straddles a bin boundary (which findOverlaps(type="any") would do)
fr[, mid := (start + end) %/% 2L]
fr[, bin_start := as.numeric((mid %/% 1000000L) * 1000000L)]   # bin_info_hg38 is 0-based
fr <- fr[chr %in% bins$chr]
setkey(fr, chr, bin_start)

bmap <- bins[, .(chr, start, bin_id)]
setkey(bmap, chr, start)
fr <- bmap[fr, on = .(chr, start = bin_start), nomatch = 0]

cnt <- fr[, .N, by = .(barcode, bin_id)]
cell_ids <- sort(unique(cnt$barcode))
mat <- matrix(0L, nrow = length(cell_ids), ncol = nrow(bins),
              dimnames = list(cell_ids, NULL))
mat[cbind(match(cnt$barcode, cell_ids), cnt$bin_id)] <- cnt$N
message(sprintf("count matrix: %d cells x %d bins (1 Mb)", nrow(mat), ncol(mat)))

# ---- normalise + segment (package defaults) -------------------------------
if (ref_src == CARE) {
  is_ref <- rownames(mat) %in% cells[is_reference == TRUE, barcode]
} else {
  # external reference: bin the donor sample too and stack it on
  o_cells <- fread(file.path(ROOT, "Data/cell_lists", paste0(ref_src, "_cells.tsv")))
  o_frag  <- file.path(ROOT, "Data/atac_filtered", paste0(ref_src, "_fragments.tsv.gz"))
  ofr <- fread(cmd = paste("zcat", shQuote(o_frag)), header = FALSE, select = c(1,2,3,4),
               col.names = c("chr","start","end","barcode"))
  ofr <- ofr[barcode %in% o_cells[is_reference == TRUE, barcode] & chr %in% bins$chr]
  ofr[, mid := (start + end) %/% 2L][, bin_start := as.numeric((mid %/% 1000000L) * 1000000L)]
  ofr <- bmap[ofr, on = .(chr, start = bin_start), nomatch = 0]
  ocnt <- ofr[, .N, by = .(barcode, bin_id)]
  oid <- sort(unique(ocnt$barcode))
  omat <- matrix(0L, nrow = length(oid), ncol = nrow(bins), dimnames = list(paste0("REF_", oid), NULL))
  omat[cbind(match(ocnt$barcode, oid), ocnt$bin_id)] <- ocnt$N
  mat <- rbind(mat, omat)
  is_ref <- grepl("^REF_", rownames(mat))
}
message(sprintf("reference nuclei: %d (source %s)", sum(is_ref), ref_src))

# AtaCNV's segmentation shells out to the precompiled MBICseq binary, which reads
# and writes temp files in the PROCESS WORKING DIRECTORY under fixed names. Twenty
# array tasks sharing one cwd clobber each other's temp files -- the symptom is
# "Error in bin_read: unequal number of columns" and read.table connection
# failures on some chromosomes but not others. Isolate each task in its own dir.
scratch <- file.path(OUT, "mbicseq_scratch")
dir.create(scratch, showWarnings = FALSE, recursive = TRUE)
old_wd <- setwd(scratch); on.exit(setwd(old_wd), add = TRUE)

t0 <- Sys.time()
norm_re <- AtaCNV::normalize(mat, genome = "hg38", mode = "normal cells",
                             normal_cells = is_ref, gc_correction = TRUE)

# AtaCNV's LOESS GC correction can drive a baseline entry non-positive, and the
# precompiled MBICseq binary then aborts the whole sample with
#   Error, "read count" is negative (tumor_2543=99.05, normal_2543 = -16.84)
# leaving no .bin_seg for R to read. Measured on the two samples where it happened
# (P109T1, P111T3) this is ONE entry out of ~10 million -- one cell in one bin,
# the same bin in both -- not a contiguous region. Repair it with the median of
# the same bin across the cells where the baseline is positive, which is the
# least-damage option and keeps the bin grid identical across samples.
bad <- which(norm_re$baseline <= 0 | !is.finite(norm_re$baseline), arr.ind = TRUE)
if (nrow(bad)) {
  for (j in unique(bad[, "col"])) {
    v  <- norm_re$baseline[, j]
    ok <- is.finite(v) & v > 0
    norm_re$baseline[!ok, j] <- median(v[ok])
  }
  message(sprintf("repaired %d non-positive baseline entries in %d bin(s)",
                  nrow(bad), length(unique(bad[, "col"]))))
}

res <- AtaCNV::calculate_CNV(norm_count = norm_re$norm_count, baseline = norm_re$baseline)
el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))

saveRDS(list(care_id = CARE, ref_source = ref_src, n_ref = sum(is_ref),
             bins = bins, norm = norm_re, cnv = res),
        file.path(OUT, "atacnv_result.rds"))

# The scratch holds MBICseq's per-chromosome .bin / .bin_seg intermediates plus a
# CNV_result.rds that calculate_CNV writes as a side effect -- all of it is carried
# in the RDS above. Drop it on success; leave it in place on failure, where it is
# the only record of which chromosome MBICseq choked on.
setwd(old_wd); unlink(scratch, recursive = TRUE)

cat(sprintf("SUMMARY\t%s\t%d\t%d\t%s\t%.1f\n", CARE, nrow(mat), sum(is_ref), ref_src, el))
