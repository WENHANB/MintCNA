# Drop-in fast replacement for RIDDLER's util/fast_window_matrix.R::window_matrix
#
# The original is unusable on multiome fragments: it does fread(data.table=FALSE)
# then frags[frags[,4] %in% barcodes,] on ~40M rows of a data.frame, then a
# further frags[frags[1] == chr,] scan per chromosome. A single sample did not
# finish the barcode filter in 8 hours. Same class of problem as AtaCNV's
# generate_input().
#
# This reimplements it in data.table. The window-assignment rule is copied line
# for line so the output is identical, not merely similar:
#
#   s_ids = findInterval(start, wnd_start)
#   e_ids = findInterval(end,   wnd_end) + 1
#   for fragments where s_ids != e_ids (they straddle a boundary):
#       offset = which.max(c(|start - wnd_end[s]|, |end - wnd_end[s]|)) - 1
#       s_ids  = s_ids + offset
#
# which.max returns the FIRST maximum, so a tie gives offset 0; the vectorised
# form therefore has to use a STRICT `>`, not `>=`.
#
# Output is byte-compatible with what riddler_fit() expects:
# write.table(cbind(window[,1:4], counts), sep = " ", header, no rownames.

window_matrix_dt <- function(bed_file, w_file, out_file, frag_min = 5000) {
  suppressPackageStartupMessages(library(data.table))
  chr_list <- c(paste0("chr", 1:22), "chrX", "chrY")

  message("reading ", bed_file)
  frags <- fread(bed_file, select = 1:4,
                 col.names = c("chr", "start", "end", "bc"), showProgress = FALSE)

  n_by_bc  <- frags[, .N, by = bc]
  barcodes <- sort(n_by_bc[N >= frag_min, bc])          # `>=`, as in the original
  frags    <- frags[bc %in% barcodes & chr %in% chr_list]
  message(sprintf("barcodes filtered: %d ; fragments to read: %d",
                  length(barcodes), nrow(frags)))

  window <- fread(w_file, header = FALSE, showProgress = FALSE)
  setnames(window, 1:4, c("wchr", "wstart", "wend", "wname"))
  window[, gid := .I]
  n <- nrow(window)

  assign_gid <- function(d, wnd) {
    s <- findInterval(d$start, wnd$wstart)
    e <- findInterval(d$end,   wnd$wend) + 1L
    straddle <- which(s != e)
    if (length(straddle)) {
      ws <- wnd$wend[s[straddle]]
      # strict `>` reproduces which.max's first-maximum tie rule
      s[straddle] <- s[straddle] +
        as.integer(abs(d$end[straddle] - ws) > abs(d$start[straddle] - ws))
    }
    ok <- s >= 1L & s <= nrow(wnd)
    list(gid = wnd$gid[s[ok]], bc = d$bc[ok])
  }

  parts <- lapply(chr_list, function(cc) {
    d <- frags[chr == cc]; if (!nrow(d)) return(NULL)
    w <- window[wchr == cc]; if (!nrow(w)) return(NULL)
    as.data.table(assign_gid(d, w))
  })
  hits <- rbindlist(parts)
  cnt  <- hits[, .N, by = .(gid, bc)]

  counts <- matrix(0L, nrow = n, ncol = length(barcodes),
                   dimnames = list(NULL, barcodes))
  counts[cbind(cnt$gid, match(cnt$bc, barcodes))] <- cnt$N
  message(sprintf("count matrix: %d windows x %d cells", n, length(barcodes)))

  # The original reads the window bed with header=FALSE, so its first four column
  # names come out V1..V4. riddler_fit() indexes them positionally, but match the
  # names anyway so the two outputs are byte-identical.
  out <- cbind(as.data.frame(window[, 1:4]), as.data.frame(counts))
  names(out)[1:4] <- paste0("V", 1:4)
  write.table(out, file = out_file, row.names = FALSE, col.names = TRUE, quote = FALSE)
  invisible(out_file)
}
