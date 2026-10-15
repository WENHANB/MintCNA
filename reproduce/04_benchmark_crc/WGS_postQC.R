#!/usr/bin/env Rscript
# Post-QC of the bulk WGS truth: drop CBS segments too short to be believed.
#   Rscript WGS_postQC.R [--min-bins 60] [--gap-mb 2]
# Reads  Result/CRC/Benchmark/DNAcopy.RData  (+ wgs_fp_null.csv if present)
# Writes Result/CRC/Benchmark/DNAcopy_qc.RData
#        Result/CRC/Benchmark/wgs_postqc_removed.csv
# The original DNAcopy.RData is never modified.
#
# WHY 60 BINS (6 Mb)
#   WGS_postQC_calibrate.R runs the identical CBS on a residual that contains no
#   CNA by construction. Every segment it finds there is a false positive:
#   225 of them over 44 null runs, median 11 bins, 99th percentile 36, max 60.
#   The real truth has nothing at all between 4.8 Mb and 6.5 Mb, so ANY cut in
#   4.9-6.5 Mb produces the identical filtered truth -- the choice of 6 Mb is not
#   load-bearing. Below the gap sit 86 events holding 4.6% of the altered bins;
#   above it sit the arm-level and whole-chromosome events the biology is made of.
#
# WHY EVENTS, NOT SEGMENTS
#   CBS chops real arm-level events into pieces separated by short neutral gaps:
#   the 14q loss in ML1 arrives as six fragments of 67/26/131/40/25/47 bins
#   spanning 23.8-63.9 Mb. Judging each fragment on its own length would punch
#   holes through a real 40 Mb event. So same-state runs separated by a neutral
#   gap under --gap-mb are grouped into one event, and the length test is applied
#   to the event.
#
# DELETE ONLY
#   Grouping decides what to remove; it never writes a call. The neutral gaps
#   inside a surviving event stay neutral. This filter can only turn an altered
#   bin neutral, never the reverse, so it cannot invent support for any caller.

suppressPackageStartupMessages(library(data.table))
ROOT <- "<MINTCNA_ROOT>"
BD   <- file.path(ROOT, "Result/CRC/Benchmark")
args <- commandArgs(trailingOnly = TRUE)
getn <- function(k, d) { i <- match(k, args); if (is.na(i)) d else as.numeric(args[i + 1]) }
MIN_BINS <- getn("--min-bins", 60)
GAP_MB   <- getn("--gap-mb",    2)
CHRS     <- paste0("chr", 1:22)

e <- new.env(); load(file.path(BD, "DNAcopy.RData"), envir = e)
T0  <- e$cna_dnacopy
ref <- as.data.table(e$ref_dnacopy); setnames(ref, 1, "chr")
ref[, chr := factor(chr, levels = CHRS)]
stopifnot(nrow(ref) == nrow(T0), !any(is.na(ref$chr)))
stopifnot(identical(order(ref$chr, ref$start), seq_len(nrow(ref))))
stopifnot(all(T0 %in% 1:3))

nullf <- file.path(BD, "wgs_fp_null.csv")
if (file.exists(nullf)) {
  nl <- fread(nullf)
  cat(sprintf("calibration: %d null FP segments, p99 = %.0f bins, max = %.0f bins\n",
              nrow(nl), quantile(nl$num.mark, .99), max(nl$num.mark)))
  if (MIN_BINS <= quantile(nl$num.mark, .99))
    warning("--min-bins is at or below the 99th percentile of the null; FPs will survive")
} else cat("note: wgs_fp_null.csv not found -- using --min-bins without calibration evidence\n")

# ---- runs of identical state, then runs grouped into events across short gaps --
rle_runs <- function(v) {
  rbindlist(lapply(levels(ref$chr), function(cc) {
    i <- which(ref$chr == cc); r <- rle(v[i])
    en <- cumsum(r$lengths); st <- en - r$lengths + 1L
    data.table(chr = cc, state = r$values, n_bins = r$lengths,
               i_from = i[st], i_to = i[en], start = ref$start[i][st], end = ref$end[i][en])
  }))
}

removed <- list()
T1 <- T0
for (s in colnames(T0)) {
  rn <- rle_runs(T0[, s])[state != 2L]
  setorder(rn, chr, start)
  # one event = consecutive same-state runs whose neutral gap is under GAP_MB
  rn[, ev := cumsum(c(1L, as.integer(
        head(state, -1) != tail(state, -1) |
        (tail(start, -1) - head(end, -1)) / 1e6 >= GAP_MB))), by = chr]
  ev <- rn[, .(state = state[1], n_bins = sum(n_bins), n_frag = .N,
               start = min(start), end = max(end)), by = .(chr, ev)]
  drop <- ev[n_bins < MIN_BINS]
  di <- rn[ev %in% drop$ev & chr %in% drop$chr]
  di <- merge(rn, drop[, .(chr, ev)], by = c("chr", "ev"))
  if (nrow(di)) T1[unlist(Map(`:`, di$i_from, di$i_to)), s] <- 2L
  removed[[s]] <- drop[, .(site = s, chr, state, n_bins, n_frag,
                           start, end, span_mb = round((end - start) / 1e6, 2))]
  cat(sprintf("%-4s  %3d events -> dropped %3d (%4d bins, %.2f%% of the genome)\n",
              s, nrow(ev), nrow(drop), sum(drop$n_bins),
              100 * sum(drop$n_bins) / nrow(T0)))
}
removed <- rbindlist(removed)

cat(sprintf("\naltered bins: %d -> %d  (%.1f%% removed)\n",
            sum(T0 != 2L), sum(T1 != 2L), 100 * (1 - sum(T1 != 2L) / sum(T0 != 2L))))
stopifnot(all(T1[T0 == 2L] == 2L))                  # delete-only: no call invented
stopifnot(all(T1[T1 != 2L] == T0[T1 != 2L]))        # no call flipped, only cleared

cna_dnacopy_raw <- T0
cna_dnacopy     <- T1
ref_dnacopy     <- e$ref_dnacopy
attr(cna_dnacopy, "postQC") <- list(min_bins = MIN_BINS, gap_mb = GAP_MB,
                                    script = "<REPRODUCE_DIR>/04_benchmark_crc/WGS_postQC.R")
save(cna_dnacopy, cna_dnacopy_raw, ref_dnacopy, file = file.path(BD, "DNAcopy_qc.RData"))
fwrite(removed[order(-n_bins)], file.path(BD, "wgs_postqc_removed.csv"))
cat("wrote:", file.path(BD, "DNAcopy_qc.RData"), "\n")

cat("\nper-chromosome altered-bin fraction, before -> after\n")
b <- sapply(colnames(T0), function(s) sapply(split(T0[, s] != 2L, ref$chr), mean))
a <- sapply(colnames(T1), function(s) sapply(split(T1[, s] != 2L, ref$chr), mean))
chg <- which(rowSums(abs(a - b)) > 0.005)
print(round(cbind(before = rowMeans(b), after = rowMeans(a))[chg, ], 3))
