#!/usr/bin/env Rscript
# Assemble the MintCNA input for one multiome sample: the ATAC bin-count matrix
# that stands in for DNA read counts, the bin annotation, and the cell table.
#   Usage: Rscript make_mintcna_input.R P111T3 [bin_width]

suppressPackageStartupMessages({library(data.table)})
ROOT <- "<MINTCNA_ROOT>/ATAC_RNA"; setwd(ROOT)
args <- commandArgs(trailingOnly = TRUE)
S <- if (length(args)) args[1] else "P111T3"
W <- if (length(args) > 1) as.numeric(args[2]) else 5e5
CHRS <- paste0("chr", 1:22)

OUT <- file.path("Data/mintcna_input", S); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

bins  <- fread("Data/ref/bins_500kb_hg38.tsv.gz")[chr %in% CHRS]
cells <- fread(sprintf("Data/cell_lists/%s_cells.tsv", S))
cnt   <- fread(sprintf("Data/mintcna_input/%s_counts_500kb.tsv.gz", S))

# Malignant nuclei are what gets called; reference nuclei are consumed by the
# denominator and dropped from the L2R output, exactly as in the DNA arm.
# Depth QC (> 3,000 fragments) is applied inside bin_fragments.awk, so a barcode
# that failed it is simply absent from the counts file. Intersect rather than
# re-deriving the threshold here, or those cells become all-zero columns.
have <- unique(cnt$barcode)
mal <- intersect(cells[cell_type == "Malignant", barcode], have)
ref <- intersect(cells[is_reference == TRUE,     barcode], have)
use <- c(mal, ref)
cat(sprintf("%s: depth QC keeps %d/%d malignant, %d/%d reference\n", S,
            length(mal), cells[cell_type == "Malignant", .N],
            length(ref), cells[is_reference == TRUE, .N]))

bins[, id := sprintf("%s:%d-%d", chr, start, end)]
setorderv(bins, c("chr", "start"))
bins[, ord := .I]
key_b <- setNames(bins$ord, paste(bins$chr, bins$bin))

Y <- matrix(0L, nrow(bins), length(use), dimnames = list(bins$id, use))
bi <- key_b[paste(cnt$chr, cnt$bin)]
ci <- match(cnt$barcode, use)
ok <- !is.na(bi) & !is.na(ci)
Y[cbind(bi[ok], ci[ok])] <- cnt$count[ok]
cat(sprintf("%s: %d bins x %d cells (%d malignant + %d reference); %s fragments placed\n",
            S, nrow(Y), ncol(Y), length(mal), length(ref),
            format(sum(cnt$count[ok]), big.mark = ",")))

ann <- bins[, .(chr, start, end, gc, mapp, frac_valid, frac_blacklist)]
ct  <- data.table(cell = use,
                  normal = ifelse(use %in% ref, "normal", "tumor"),
                  cell_type = cells$cell_type[match(use, cells$barcode)],
                  nFrags    = cells$nFrags[match(use, cells$barcode)])

fwrite(data.table(bin = rownames(Y), as.data.table(Y)),
       file.path(OUT, "atac_counts.csv"))
fwrite(ann, file.path(OUT, "bin_annotation.csv"))
fwrite(ct,  file.path(OUT, "cell_types.tsv"), sep = "\t")
cat("wrote", OUT, "\n")
