#!/usr/bin/env Rscript
# Turn the bulk-DNA GATK segments into a per-sample vector on the 100 kb hg38 grid.
#
# Two things have to happen and both are easy to get silently wrong:
#
#   1. LIFTOVER. The GATK segments are b37/hg19 -- the GLASS pipeline used
#      human_g1k_v37_decoy.fasta, and the table carries bare integer chromosome
#      names. The single-cell data is hg38 throughout. Verified empirically:
#      14 of 22 autosomes have maximum segment coordinates exceeding the hg38
#      chromosome length, and none exceed hg19. A coordinate join without
#      liftOver raises no error and silently misaligns every breakpoint.
#
#   2. LENGTH WEIGHTING. Segment lengths span 10 kb to 75 Mb, and 49% are shorter
#      than one 100 kb bin. Unweighted statistics are dominated by the short ones:
#      the unweighted median log2 across the cohort is -0.566, the length-weighted
#      median is +0.004. Only the latter is the true diploid centre.
#
# Output: Result/truth/<care_id>_truth_100kb.tsv.gz with columns
#   chr  start  end  log2  cnv_call  n_seg  frac_covered
# plus Result/truth/truth_summary.tsv and a liftOver audit.

suppressPackageStartupMessages({
  library(data.table); library(GenomicRanges); library(rtracklayer)
})

ROOT <- "<MINTCNA_ROOT>/ATAC_RNA"
OUT  <- file.path(ROOT, "Result/truth")
BIN  <- 1e5
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

SAMPLES <- c("P59T2","P59T3","P104T1","P104T2","P106T1","P106T2",
             "P100T1","P100T3","P103T1","P103T2","P105T1","P105T3",
             "P108T1","P108T2","P109T1","P109T2","P110T1","P110T2","P111T2","P111T3")

# ---- inputs ---------------------------------------------------------------
seg   <- fread(file.path(ROOT, "Data/synapse/variants_gatk_seg.csv"))
ali   <- fread(file.path(ROOT, "Data/synapse/biospecimen_aliquots.csv"))
smap  <- fread(file.path(ROOT, "Data/synapse/sample_identifier_mapping.csv"))
chain <- import.chain(file.path(ROOT, "Data/ref/hg19ToHg38.over.chain"))

# aliquot_barcode -> sample_barcode -> care_id
key <- merge(ali[, .(aliquot_barcode, sample_barcode, assay = aliquot_analysis_type)],
             smap[, .(sample_barcode, care_id)], by = "sample_barcode")
seg <- merge(seg, key, by = "aliquot_barcode")

# "[start,end)" -> integers.  Keep autosomes; 23/24 are X/Y in b37 naming and
# need their own baseline, so they are excluded rather than silently mixed in.
seg[, `:=`(seg_start = as.numeric(sub("\\[(\\d+),.*", "\\1", pos)),
           seg_end = as.numeric(sub(".*,(\\d+)\\)", "\\1", pos)))]
seg <- seg[suppressWarnings(!is.na(as.integer(chrom))) & as.integer(chrom) %in% 1:22]

# hg38 autosome lengths, for the analysis grid
hg38_len <- c(248956422,242193529,198295559,190214555,181538259,170805979,159345973,
              145138636,138394717,133797422,135086622,133275309,114364328,107043718,
              101991189,90338345,83257441,80373285,58617616,64444167,46709983,50818468)
names(hg38_len) <- paste0("chr", 1:22)

grid <- do.call(c, lapply(names(hg38_len), function(ch)
  GRanges(ch, IRanges(seq(1, hg38_len[[ch]], by = BIN), width = BIN))))
end(grid) <- pmin(end(grid), hg38_len[as.character(seqnames(grid))])
message(sprintf("analysis grid: %s bins of %.0f kb over 22 autosomes",
                format(length(grid), big.mark=","), BIN/1e3))

audit <- list(); summ <- list()

for (samp in SAMPLES) {
  # NB: loop and column names must not collide. data.table resolves a bare symbol
  # in i to the column first, so a loop variable named `s` next to a column named
  # `s` silently evaluates care_id == <the coordinate column> and matches nothing.
  d <- seg[care_id == samp]
  if (!nrow(d)) { message("no segments for ", samp); next }
  # Prefer WGS over WXS when both exist: genome-wide beats exon-restricted.
  want <- if ("WGS" %in% d$assay) "WGS" else "WXS"
  d <- d[assay == want]
  # if several aliquots of the same assay, take the one with most segments
  best <- d[, .N, by = aliquot_barcode][order(-N)][1, aliquot_barcode]
  d <- d[aliquot_barcode == best]

  gr19 <- GRanges(paste0("chr", d$chrom), IRanges(d$seg_start, d$seg_end - 1),
                  log2 = d$log2_copy_ratio, cnv_call = d$cnv_call)

  # ---- liftOver ----------------------------------------------------------
  lifted <- liftOver(gr19, chain)
  n_in   <- length(gr19)
  n_drop <- sum(lengths(lifted) == 0)
  n_split<- sum(lengths(lifted) > 1)
  gr38   <- unlist(lifted)
  # liftOver can emit tiny fragments and can reorder; drop sub-kb debris and
  # anything that landed on a different chromosome than it started on.
  gr38 <- gr38[width(gr38) >= 1000]
  gr38 <- sortSeqlevels(gr38); gr38 <- sort(gr38)
  gr38 <- gr38[seqnames(gr38) %in% names(hg38_len)]

  audit[[samp]] <- data.table(care_id = samp, assay = want, aliquot = best,
                           n_seg_hg19 = n_in, n_dropped = n_drop, n_split = n_split,
                           n_seg_hg38 = length(gr38),
                           mb_hg19 = sum(width(gr19))/1e6,
                           mb_hg38 = sum(width(gr38))/1e6)

  # ---- project onto the grid, weighted by overlap length -----------------
  hits <- findOverlaps(grid, gr38)
  w    <- width(pintersect(grid[queryHits(hits)], gr38[subjectHits(hits)]))
  v    <- gr38$log2[subjectHits(hits)]
  num  <- tapply(v * w, queryHits(hits), sum)
  den  <- tapply(w,     queryHits(hits), sum)
  nseg <- tapply(w,     queryHits(hits), length)

  log2v <- rep(NA_real_, length(grid)); cov <- rep(0, length(grid)); ns <- rep(0L, length(grid))
  idx <- as.integer(names(num))
  log2v[idx] <- num / den
  cov[idx]   <- pmin(den / width(grid)[idx], 1)
  ns[idx]    <- as.integer(nseg)

  # discrete call: majority by covered length
  callv <- rep(NA_character_, length(grid))
  cl <- data.table(bin = queryHits(hits), call = gr38$cnv_call[subjectHits(hits)], w = w)
  cl <- cl[, .(w = sum(w)), by = .(bin, call)][order(bin, -w)]
  cl <- cl[, .SD[1], by = bin]
  callv[cl$bin] <- cl$call

  # ---- recentre on the length-weighted diploid level ----------------------
  ok <- !is.na(log2v) & cov > 0.5
  ord <- order(log2v[ok])
  # as.numeric() is load-bearing: width() returns integer, and cumsum over 2.67 Gb
  # of autosome overflows .Machine$integer.max (2.147e9), silently yielding NA,
  # which then propagates through `log2v - centre` and blanks the whole column.
  cw <- cumsum(as.numeric(width(grid)[ok][ord]))
  centre <- log2v[ok][ord][which(cw >= tail(cw, 1) / 2)[1]]
  stopifnot(is.finite(centre))
  log2v <- log2v - centre

  out <- data.table(chr = as.character(seqnames(grid)), start = start(grid), end = end(grid),
                    log2 = round(log2v, 5), cnv_call = callv,
                    n_seg = ns, frac_covered = round(cov, 3))
  fwrite(out, file.path(OUT, paste0(samp, "_truth_100kb.tsv.gz")), sep = "\t")

  summ[[samp]] <- data.table(care_id = samp, assay = want,
                          bins_total = length(grid), bins_covered = sum(ok),
                          pct_covered = round(100*sum(ok)/length(grid), 1),
                          centre_applied = round(centre, 4),
                          frac_gain = round(mean(callv[ok] == "+", na.rm=TRUE), 3),
                          frac_loss = round(mean(callv[ok] == "-", na.rm=TRUE), 3))
  message(sprintf("  %-8s %s  %d->%d seg  %.1f%% bins covered  centre %+.3f",
                  samp, want, n_in, length(gr38), 100*sum(ok)/length(grid), centre))
}

A <- rbindlist(audit); S <- rbindlist(summ)
fwrite(A, file.path(OUT, "liftover_audit.tsv"), sep = "\t")
fwrite(S, file.path(OUT, "truth_summary.tsv"), sep = "\t")
cat("\n=== liftOver audit ===\n"); print(A)
cat("\n=== truth summary ===\n"); print(S)
cat(sprintf("\nsegments dropped by liftOver: %d / %d (%.2f%%)\n",
            sum(A$n_dropped), sum(A$n_seg_hg19), 100*sum(A$n_dropped)/sum(A$n_seg_hg19)))
