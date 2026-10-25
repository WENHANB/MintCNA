#!/usr/bin/env Rscript
# Bulk-DNA truth of the CARE samples on the 100 kb hg38 grid: the GATK copy-number segments (b37 coordinates)
# are lifted to hg38, projected onto the grid weighted by overlap length, and re-centred on the
# length-weighted median log2 ratio (the diploid level). WGS is preferred over WXS when both exist.
#   Rscript reproduce/analysis/multiome/evaluation/prepare_truth.R
# Input : ATAC_RNA/Data/synapse/{variants_gatk_seg.csv, biospecimen_aliquots.csv, sample_identifier_mapping.csv}
#         (download_synapse.py), ATAC_RNA/Data/ref/hg19ToHg38.over.chain
# Output: ATAC_RNA/Result/truth/<care_id>_truth_100kb.tsv.gz (chr, start, end, log2, cnv_call, n_seg, frac_covered),
#         liftover_audit.tsv, truth_summary.tsv
# Run from the root of the MintCNA repository with MINTCNA_DATA set (see reproduce/README.md).
suppressPackageStartupMessages({ library(data.table); library(GenomicRanges); library(rtracklayer) })
ROOT <- Sys.getenv("MINTCNA_DATA"); if (!nzchar(ROOT)) stop("set MINTCNA_DATA (see reproduce/README.md)")
ROOT <- file.path(ROOT, "ATAC_RNA")
OUT  <- file.path(ROOT, "Result/truth"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
BIN  <- 1e5
SAMPLES <- c("P59T2","P59T3","P104T1","P104T2","P106T1","P106T2","P100T1","P100T3","P103T1","P103T2",
             "P105T1","P105T3","P108T1","P108T2","P109T1","P109T2","P110T1","P110T2","P111T2","P111T3")

seg   <- fread(file.path(ROOT, "Data/synapse/variants_gatk_seg.csv"))
ali   <- fread(file.path(ROOT, "Data/synapse/biospecimen_aliquots.csv"))
smap  <- fread(file.path(ROOT, "Data/synapse/sample_identifier_mapping.csv"))
chain <- import.chain(file.path(ROOT, "Data/ref/hg19ToHg38.over.chain"))

key <- merge(ali[, .(aliquot_barcode, sample_barcode, assay = aliquot_analysis_type)], smap[, .(sample_barcode, care_id)], by = "sample_barcode")
seg <- merge(seg, key, by = "aliquot_barcode")
seg[, `:=`(seg_start = as.numeric(sub("\\[(\\d+),.*", "\\1", pos)), seg_end = as.numeric(sub(".*,(\\d+)\\)", "\\1", pos)))]
seg <- seg[suppressWarnings(!is.na(as.integer(chrom))) & as.integer(chrom) %in% 1:22]       # autosomes

hg38_len <- c(248956422,242193529,198295559,190214555,181538259,170805979,159345973,145138636,138394717,133797422,
              135086622,133275309,114364328,107043718,101991189,90338345,83257441,80373285,58617616,64444167,46709983,50818468)
names(hg38_len) <- paste0("chr", 1:22)
grid <- do.call(c, lapply(names(hg38_len), function(ch) GRanges(ch, IRanges(seq(1, hg38_len[[ch]], by = BIN), width = BIN))))
end(grid) <- pmin(end(grid), hg38_len[as.character(seqnames(grid))])

audit <- list(); summ <- list()
for (samp in SAMPLES) {
  d <- seg[care_id == samp]
  if (!nrow(d)) { message("no segments for ", samp); next }
  want <- if ("WGS" %in% d$assay) "WGS" else "WXS"
  d <- d[assay == want]
  best <- d[, .N, by = aliquot_barcode][order(-N)][1, aliquot_barcode]        # the aliquot with most segments
  d <- d[aliquot_barcode == best]
  gr19 <- GRanges(paste0("chr", d$chrom), IRanges(d$seg_start, d$seg_end - 1), log2 = d$log2_copy_ratio, cnv_call = d$cnv_call)

  lifted <- liftOver(gr19, chain)
  gr38 <- unlist(lifted)
  gr38 <- gr38[width(gr38) >= 1000]                                             # drop sub-kb liftOver debris
  gr38 <- sort(sortSeqlevels(gr38)); gr38 <- gr38[seqnames(gr38) %in% names(hg38_len)]
  audit[[samp]] <- data.table(care_id = samp, assay = want, aliquot = best, n_seg_hg19 = length(gr19),
                              n_dropped = sum(lengths(lifted) == 0), n_split = sum(lengths(lifted) > 1), n_seg_hg38 = length(gr38),
                              mb_hg19 = sum(width(gr19)) / 1e6, mb_hg38 = sum(width(gr38)) / 1e6)

  # overlap-weighted log2 ratio per grid bin; discrete call = the call covering most of the bin
  hits <- findOverlaps(grid, gr38)
  w <- width(pintersect(grid[queryHits(hits)], gr38[subjectHits(hits)]))
  num <- tapply(gr38$log2[subjectHits(hits)] * w, queryHits(hits), sum); den <- tapply(w, queryHits(hits), sum)
  idx <- as.integer(names(num))
  log2v <- rep(NA_real_, length(grid)); cov <- rep(0, length(grid)); ns <- rep(0L, length(grid))
  log2v[idx] <- num / den; cov[idx] <- pmin(den / width(grid)[idx], 1); ns[idx] <- as.integer(tapply(w, queryHits(hits), length))
  cl <- data.table(bin = queryHits(hits), call = gr38$cnv_call[subjectHits(hits)], w = w)[, .(w = sum(w)), by = .(bin, call)][order(bin, -w)]
  cl <- cl[, .SD[1], by = bin]
  callv <- rep(NA_character_, length(grid)); callv[cl$bin] <- cl$call

  # recentre on the length-weighted median of the covered bins (as.numeric: the cumulative width exceeds .Machine$integer.max)
  ok <- !is.na(log2v) & cov > 0.5
  ord <- order(log2v[ok])
  cw <- cumsum(as.numeric(width(grid)[ok][ord]))
  centre <- log2v[ok][ord][which(cw >= tail(cw, 1) / 2)[1]]
  stopifnot(is.finite(centre))
  log2v <- log2v - centre

  fwrite(data.table(chr = as.character(seqnames(grid)), start = start(grid), end = end(grid), log2 = round(log2v, 5),
                    cnv_call = callv, n_seg = ns, frac_covered = round(cov, 3)),
         file.path(OUT, paste0(samp, "_truth_100kb.tsv.gz")), sep = "\t")
  summ[[samp]] <- data.table(care_id = samp, assay = want, bins_total = length(grid), bins_covered = sum(ok),
                             pct_covered = round(100 * sum(ok) / length(grid), 1), centre_applied = round(centre, 4),
                             frac_gain = round(mean(callv[ok] == "+", na.rm = TRUE), 3), frac_loss = round(mean(callv[ok] == "-", na.rm = TRUE), 3))
  message(sprintf("%-8s %s  %d -> %d segments  %.1f%% bins covered  centre %+.3f", samp, want, length(gr19), length(gr38), 100 * sum(ok) / length(grid), centre))
}
fwrite(rbindlist(audit), file.path(OUT, "liftover_audit.tsv"), sep = "\t")
fwrite(rbindlist(summ),  file.path(OUT, "truth_summary.tsv"), sep = "\t")
