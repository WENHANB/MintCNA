#!/usr/bin/env Rscript
# Build ref/hg38_gcmapp_5kb.csv.gz, the hg38 twin of the bundled hg19 track.
#
# MintCNA's attach_gc_mapp() reads "<genome>_gcmapp_5kb.csv.gz" with base R and
# aggregates it onto whatever bins the user supplies, so the hg38 file only has
# to match the hg19 file's schema exactly:
#     chr,start,end,gc,mapp,n_valid
#     gc      fraction of valid bases that are G or C, 0-1, NA where n_valid = 0
#     mapp    50-mer mappability, 0-1
#     n_valid non-N bases in the 5 kb tile, 0-5000  (GC weight in aggregate_gcmapp)
#     start   0-BASED, end exclusive; autosomes only, as in hg19
#
# Source: QDNAseq.hg38 (asntech/QDNAseq.hg38), object hg38.5kbp.SR50. It carries
# GC, non-N fraction and GenMap 50-mer mappability precomputed on GRCh38, which
# is the same k-mer length as the wgEncodeCrgMapabilityAlign50mer track behind
# hg19 -- so GC/mappability correction stays comparable across the two builds.
# Building it from scratch would need GenMap plus a bedtools pass over the whole
# genome; there is no reason to when a maintained annotation exists.
#
# Three conversions are required and each one silently corrupts the track if
# missed: chromosome names gain a "chr" prefix, start goes from 1-based to
# 0-based, and bases/gc/mappability are PERCENTAGES in QDNAseq but FRACTIONS here.

suppressPackageStartupMessages({library(QDNAseq.hg38); library(Biobase)})

REF <- "/orange/feifeixiao/wenhan.bao/Multi-omics/MintCNA_git/ref"
OUT <- file.path(REF, "hg38_gcmapp_5kb.csv.gz")

data("hg38.5kbp.SR50", package = "QDNAseq.hg38")
b <- pData(hg38.5kbp.SR50)
message(sprintf("QDNAseq hg38 5kb bins: %d", nrow(b)))

AUTOSOMES <- as.character(1:22)
b <- b[b$chromosome %in% AUTOSOMES, ]

track <- data.frame(
  chr     = paste0("chr", b$chromosome),
  start   = as.integer(b$start - 1L),          # 1-based -> 0-based
  end     = as.integer(b$end),
  gc      = round(b$gc          / 100, 4),     # percent -> fraction
  mapp    = round(b$mappability / 100, 4),
  n_valid = as.integer(round(b$bases / 100 * 5000)),
  stringsAsFactors = FALSE
)

# Keep the hg19 convention: gc/mapp are blank where the tile has no valid bases.
track$gc[track$n_valid == 0]   <- NA
track$mapp[track$n_valid == 0] <- NA

# Order by chromosome as chr1..chr22, then position -- matching the hg19 file so
# any code that assumes contiguous per-chromosome blocks keeps working.
track <- track[order(match(track$chr, paste0("chr", 1:22)), track$start), ]

# Every chromosome ends in one short tile -- 22 of them, exactly as in the hg19
# file -- so widths are <= 5000, not == 5000, and n_valid is capped by the tile's
# own width rather than by 5000.
w <- track$end - track$start
stopifnot(!anyDuplicated(paste(track$chr, track$start)),
          all(w > 0 & w <= 5000L),
          sum(w != 5000L) == 22L,
          all(track$n_valid >= 0 & track$n_valid <= w),
          all(is.na(track$gc)   | (track$gc   >= 0 & track$gc   <= 1)),
          all(is.na(track$mapp) | (track$mapp >= 0 & track$mapp <= 1)))

con <- gzfile(OUT, "w"); write.csv(track, con, row.names = FALSE, na = ""); close(con)

hg19 <- read.csv(file.path(REF, "hg19_gcmapp_5kb.csv.gz"), stringsAsFactors = FALSE)
cat(sprintf("\nwrote %s\n", OUT))
cat(sprintf("  %-6s rows %7d   non-NA gc %7d (%.1f%%)   median gc %.4f   median mapp %.4f\n",
            "hg38", nrow(track), sum(!is.na(track$gc)), 100*mean(!is.na(track$gc)),
            median(track$gc, na.rm = TRUE), median(track$mapp, na.rm = TRUE)))
cat(sprintf("  %-6s rows %7d   non-NA gc %7d (%.1f%%)   median gc %.4f   median mapp %.4f\n",
            "hg19", nrow(hg19), sum(!is.na(hg19$gc)), 100*mean(!is.na(hg19$gc)),
            median(hg19$gc, na.rm = TRUE), median(hg19$mapp, na.rm = TRUE)))
cat("  columns identical to hg19: ", identical(names(track), names(hg19)), "\n")
