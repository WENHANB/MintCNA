#!/usr/bin/env Rscript
# Build the 500 kb bin annotation for the multiome samples: GC, mappability and
# ENCODE-blacklist overlap, aggregated up from the bundled MintCNA 5 kb hg38 track.
#
# GC and mappability are averaged over the 5 kb sub-bins that carry a value, and
# `frac_valid` records how much of the 500 kb window had one at all -- a bin
# sitting on a centromere or an assembly gap has GC/mapp defined on only a sliver
# of itself, and that sliver is not representative of the window.

suppressPackageStartupMessages(library(data.table))
ROOT <- "<MINTCNA_ROOT>"
args <- commandArgs(trailingOnly = TRUE)
W    <- if (length(args)) as.numeric(args[1]) else 5e5
CHRS <- paste0("chr", 1:22)
TAG  <- if (W >= 1e6) sprintf("%gmb", W / 1e6) else sprintf("%gkb", W / 1e3)

g <- fread("<REPRODUCE_DIR>/../ref/hg38_gcmapp_5kb.csv.gz")[chr %in% CHRS]
g[, bin := start %/% W]
b <- g[, .(gc         = mean(gc,   na.rm = TRUE),
           mapp       = mean(mapp, na.rm = TRUE),
           n_sub      = .N,
           frac_valid = mean(!is.na(gc))), by = .(chr, bin)]

# ENCODE blacklist v2: fraction of the 500 kb window inside a blacklisted region.
bl <- fread(file.path(ROOT, "ATAC_RNA/Data/ref/hg38-blacklist.v2.bed"),
            header = FALSE, col.names = c("chr", "start", "end", "reason"))[chr %in% CHRS]
bl_bp <- rbindlist(lapply(seq_len(nrow(bl)), function(i) {
  s <- bl$start[i]; e <- bl$end[i]
  k <- (s %/% W):(e %/% W)
  data.table(chr = bl$chr[i], bin = k,
             bp  = pmin(e, (k + 1) * W) - pmax(s, k * W))
}))[, .(bl_bp = sum(bp)), by = .(chr, bin)]

b <- merge(b, bl_bp, by = c("chr", "bin"), all.x = TRUE)
b[is.na(bl_bp), bl_bp := 0]
b[, `:=`(frac_blacklist = bl_bp / W,
         start = bin * W, end = (bin + 1) * W)]
b[, chr := factor(chr, levels = CHRS)]
setorder(b, chr, bin)
b[, chr := as.character(chr)]
b <- b[, .(chr, start, end, bin, gc, mapp, frac_valid, frac_blacklist)]

fwrite(b, file.path(ROOT, sprintf("ATAC_RNA/Data/ref/bins_%s_hg38.tsv.gz", TAG)), sep = "\t")

cat(sprintf("%d bins of %s on chr1-22\n\n", nrow(b), TAG))
cat("how many bins each filter would remove, applied on its own:\n")
cat(sprintf("  frac_valid  < 0.50 : %5d\n", b[frac_valid     < 0.50, .N]))
cat(sprintf("  mapp        < 0.90 : %5d\n", b[mapp           < 0.90, .N]))
cat(sprintf("  gc outside .2-.8   : %5d\n", b[gc < 0.2 | gc  > 0.80, .N]))
cat(sprintf("  blacklist   > 0.10 : %5d\n", b[frac_blacklist > 0.10, .N]))
cat(sprintf("  blacklist   > 0.50 : %5d\n", b[frac_blacklist > 0.50, .N]))
cat(sprintf("\nall four together keep %d bins\n",
            b[frac_valid >= .5 & mapp >= .9 & gc >= .2 & gc <= .8 & frac_blacklist <= .1, .N]))
cat("\nGC and mappability at 500 kb are averages of 100 sub-bins, so both are far\n",
    "narrower than at 5 kb -- the DNA thresholds bite differently here:\n", sep = "")
print(b[, .(gc = quantile(gc, c(.001,.01,.5,.99,.999), na.rm = TRUE),
            mapp = quantile(mapp, c(.001,.01,.5,.99,.999), na.rm = TRUE),
            q = c("0.1%","1%","50%","99%","99.9%"))])
