#!/usr/bin/env Rscript
# Post-QC of the CRC01 bulk WGS truth: altered events shorter than --min-bins are set to neutral.
# An event is a run of same-state bins, where runs separated by a neutral gap shorter than --gap-mb
# count as one event (CBS fragments real arm-level events). The default cut is 60 bins = 6 Mb; the
# filter only clears calls, never adds any.
#   Rscript reproduce/analysis/crc01/wgs_truth_postqc.R [--min-bins 60] [--gap-mb 2]
# Input : Result/CRC/Benchmark/DNAcopy.RData (wgs_truth.R)
# Output: Result/CRC/Benchmark/DNAcopy_qc.RData  (cna_dnacopy filtered, cna_dnacopy_raw, ref_dnacopy),
#         Result/CRC/Benchmark/wgs_postqc_removed.csv
# Run from the root of the MintCNA repository with MINTCNA_DATA set (see reproduce/README.md).
suppressPackageStartupMessages(library(data.table))
ROOT <- Sys.getenv("MINTCNA_DATA"); if (!nzchar(ROOT)) stop("set MINTCNA_DATA (see reproduce/README.md)")
BD   <- file.path(ROOT, "Result/CRC/Benchmark")
args <- commandArgs(trailingOnly = TRUE)
getn <- function(k, d) { i <- match(k, args); if (is.na(i)) d else as.numeric(args[i + 1]) }
MIN_BINS <- getn("--min-bins", 60)
GAP_MB   <- getn("--gap-mb", 2)
CHRS <- paste0("chr", 1:22)

e <- new.env(); load(file.path(BD, "DNAcopy.RData"), envir = e)
T0  <- e$cna_dnacopy
ref <- as.data.table(e$ref_dnacopy); setnames(ref, 1, "chr")
ref[, chr := factor(chr, levels = CHRS)]
stopifnot(nrow(ref) == nrow(T0), !any(is.na(ref$chr)), identical(order(ref$chr, ref$start), seq_len(nrow(ref))), all(T0 %in% 1:3))

rle_runs <- function(v) rbindlist(lapply(levels(ref$chr), function(cc) {
  i <- which(ref$chr == cc); r <- rle(v[i])
  en <- cumsum(r$lengths); st <- en - r$lengths + 1L
  data.table(chr = cc, state = r$values, n_bins = r$lengths, i_from = i[st], i_to = i[en],
             start = ref$start[i][st], end = ref$end[i][en])
}))

removed <- list(); T1 <- T0
for (s in colnames(T0)) {
  rn <- rle_runs(T0[, s])[state != 2L]
  setorder(rn, chr, start)
  rn[, ev := cumsum(c(1L, as.integer(head(state, -1) != tail(state, -1) |
                                     (tail(start, -1) - head(end, -1)) / 1e6 >= GAP_MB))), by = chr]
  ev <- rn[, .(state = state[1], n_bins = sum(n_bins), n_frag = .N, start = min(start), end = max(end)), by = .(chr, ev)]
  drop <- ev[n_bins < MIN_BINS]
  di <- merge(rn, drop[, .(chr, ev)], by = c("chr", "ev"))
  if (nrow(di)) T1[unlist(Map(`:`, di$i_from, di$i_to)), s] <- 2L
  removed[[s]] <- drop[, .(site = s, chr, state, n_bins, n_frag, start, end, span_mb = round((end - start) / 1e6, 2))]
  cat(sprintf("%-4s %3d events -> dropped %3d (%4d bins, %.2f%% of the genome)\n",
              s, nrow(ev), nrow(drop), sum(drop$n_bins), 100 * sum(drop$n_bins) / nrow(T0)))
}
removed <- rbindlist(removed)
cat(sprintf("altered bins: %d -> %d (%.1f%% removed)\n", sum(T0 != 2L), sum(T1 != 2L), 100 * (1 - sum(T1 != 2L) / sum(T0 != 2L))))
stopifnot(all(T1[T0 == 2L] == 2L), all(T1[T1 != 2L] == T0[T1 != 2L]))

cna_dnacopy_raw <- T0; cna_dnacopy <- T1; ref_dnacopy <- e$ref_dnacopy
attr(cna_dnacopy, "postQC") <- list(min_bins = MIN_BINS, gap_mb = GAP_MB, script = "reproduce/analysis/crc01/wgs_truth_postqc.R")
save(cna_dnacopy, cna_dnacopy_raw, ref_dnacopy, file = file.path(BD, "DNAcopy_qc.RData"))
fwrite(removed[order(-n_bins)], file.path(BD, "wgs_postqc_removed.csv"))
cat("wrote", file.path(BD, "DNAcopy_qc.RData"), "\n")
