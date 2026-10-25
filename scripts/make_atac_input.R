#!/usr/bin/env Rscript
# MintCNA - Build the MintCNA input from binned scATAC fragment counts
#
# Turns the long table written by scripts/atac_bin_fragments.awk into the three
# files the pipeline reads, in exactly the layout the DNA arm uses:
#
#   <out>/atac_counts.csv     bins x cells fragment counts; row names "chr:start-end"
#                             -> config input.dna_counts (with preprocessing.modality: atac)
#   <out>/bin_annotation.csv  chr, start, end, gc, mapp, frac_valid  (from the bundled track)
#                             -> config input.bin_annotation
#   <out>/cell_types.tsv      cell, normal ("normal" = reference nucleus), cell_type
#                             -> config input.cell_types
#
# Usage:
#   Rscript scripts/make_atac_input.R --counts <sample>_counts.tsv.gz \
#       --cells cells.tsv --out input/<sample> [--genome hg38] [--width 500000]
#
#   --counts  output of atac_bin_fragments.awk: barcode, chr, bin, count
#   --cells   TSV with columns  barcode, is_reference (TRUE/FALSE) [, cell_type]
#             is_reference marks the diploid reference nuclei (e.g. myeloid,
#             endothelial, lymphocyte -- non-malignant lineages unrelated to the tumour)
#   --width   bin width in bp; must equal the -v W used with the awk (default 500000)
#
# Reference nuclei are consumed by the denominator and dropped from the L2R
# output, exactly as in the DNA arm. Barcodes missing from --counts (e.g. those
# that failed the awk's MINFRAGS depth filter) are simply absent from the output.

suppressPackageStartupMessages({ library(optparse); library(data.table) })

option_list <- list(
  make_option("--counts",  type = "character", help = "long counts TSV(.gz) from atac_bin_fragments.awk"),
  make_option("--cells",   type = "character", help = "TSV: barcode, is_reference [, cell_type]"),
  make_option("--out",     type = "character", help = "output directory"),
  make_option("--genome",  type = "character", default = "hg38",  help = "bundled track: hg19 or hg38 [default %default]"),
  make_option("--width",   type = "double",    default = 500000, help = "bin width in bp [default %default]"),
  make_option("--ref_dir", type = "character", default = NULL,   help = "directory of bundled tracks [default <repo>/MintCNA/ref]"),
  make_option("--chrs",    type = "character", default = "chr1-22", help = "chromosomes: 'chr1-22' or comma list [default %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))
for (a in c("counts", "cells", "out")) if (is.null(opt[[a]])) stop("--", a, " is required.")

# Locate the repo root from this script's path, source the gc/mapp attachment.
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("--file=", args, value = TRUE)
repo <- if (length(file_arg)) dirname(dirname(normalizePath(sub("--file=", "", file_arg)))) else normalizePath(".")
source(file.path(repo, "MintCNA", "preprocess.R"))
ref_dir <- if (is.null(opt$ref_dir)) file.path(repo, "MintCNA", "ref") else opt$ref_dir

W    <- as.numeric(opt$width)
CHRS <- if (opt$chrs == "chr1-22") paste0("chr", 1:22) else strsplit(opt$chrs, ",")[[1]]

# ---- 1. Bin grid + GC / mappability / frac_valid from the bundled 5 kb track ----
track <- .load_gcmapp_track(opt$genome, ref_dir)
track <- track[track$chr %in% CHRS, ]
grid  <- unique(data.frame(chr = track$chr, bin = track$start %/% W, stringsAsFactors = FALSE))
grid$start <- grid$bin * W
grid$end   <- grid$start + W
grid <- grid[order(match(grid$chr, CHRS), grid$start), ]
ann  <- attach_gc_mapp(grid[, c("chr", "start", "end")], genome = opt$genome, ref_dir = ref_dir)
ann$id <- sprintf("%s:%d-%d", ann$chr, as.integer(ann$start), as.integer(ann$end))
key_b <- setNames(seq_len(nrow(ann)), paste(grid$chr, grid$bin))
message(sprintf("[make_atac_input] %d bins of %g bp on %d chromosomes (%s)",
                nrow(ann), W, length(CHRS), opt$genome))

# ---- 2. Cells: those present in the counts, split into reference / tumour ----
cnt   <- fread(opt$counts)
cells <- fread(opt$cells)
stopifnot(all(c("barcode", "chr", "bin", "count") %in% names(cnt)),
          all(c("barcode", "is_reference") %in% names(cells)))
cells$is_reference <- toupper(as.character(cells$is_reference)) %in% c("TRUE", "T", "1", "YES")
if (!"cell_type" %in% names(cells)) cells$cell_type <- ifelse(cells$is_reference, "reference", "tumor")

have <- unique(cnt$barcode)
use  <- intersect(cells$barcode, have)
ref  <- intersect(cells$barcode[cells$is_reference], have)
if (!length(use)) stop("No barcode in --cells appears in --counts.")
if (!length(ref)) warning("No reference nuclei found among the counted barcodes: ",
                          "the pipeline needs normal-control cells (is_reference = TRUE).")
message(sprintf("[make_atac_input] %d cells (%d reference, %d tumour); %d/%d barcodes in --cells had counts",
                length(use), length(ref), length(use) - length(ref), length(use), nrow(cells)))

# ---- 3. Pivot long -> bins x cells ----
Y  <- matrix(0L, nrow(ann), length(use), dimnames = list(ann$id, use))
bi <- key_b[paste(cnt$chr, cnt$bin)]
ci <- match(cnt$barcode, use)
ok <- !is.na(bi) & !is.na(ci)
Y[cbind(bi[ok], ci[ok])] <- as.integer(cnt$count[ok])
message(sprintf("[make_atac_input] %s fragments placed; %s dropped (off-grid chromosome or unlisted barcode)",
                format(sum(cnt$count[ok]), big.mark = ","), format(sum(cnt$count[!ok]), big.mark = ",")))

# ---- 4. Write ----
dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)
fwrite(data.table(bin = rownames(Y), as.data.table(Y)), file.path(opt$out, "atac_counts.csv"))
fwrite(as.data.table(ann[, c("chr", "start", "end", "gc", "mapp", "frac_valid")]),
       file.path(opt$out, "bin_annotation.csv"))
fwrite(data.table(cell = use,
                  normal = ifelse(use %in% ref, "normal", "tumor"),
                  cell_type = cells$cell_type[match(use, cells$barcode)]),
       file.path(opt$out, "cell_types.tsv"), sep = "\t")
message("[make_atac_input] wrote ", normalizePath(opt$out))
message("\nConfig snippet:\n",
        "input:\n",
        "  dna_counts:     \"", file.path(opt$out, "atac_counts.csv"), "\"\n",
        "  bin_annotation: \"", file.path(opt$out, "bin_annotation.csv"), "\"\n",
        "  cell_types:     \"", file.path(opt$out, "cell_types.tsv"), "\"\n",
        "  genome:         \"", opt$genome, "\"\n",
        "preprocessing:\n",
        "  modality: \"atac\"")
