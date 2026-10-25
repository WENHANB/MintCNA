# data.table replacement for RIDDLER's util/fast_window_matrix.R::window_matrix(), which does not finish on a
# multiome fragments file. The window-assignment rule is reproduced exactly:
#   s_ids = findInterval(start, wnd_start); e_ids = findInterval(end, wnd_end) + 1
#   a fragment with s_ids != e_ids is assigned to the window nearer to its larger-distance end
#   (which.max over c(|start - wnd_end[s]|, |end - wnd_end[s]|); the first maximum wins ties)
# and the output has the layout riddler_fit() expects (window columns V1-V4, then one column per barcode).
window_matrix_dt <- function(bed_file, w_file, out_file, frag_min = 5000) {
  suppressPackageStartupMessages(library(data.table))
  chr_list <- c(paste0("chr", 1:22), "chrX", "chrY")

  frags <- fread(bed_file, select = 1:4, col.names = c("chr", "start", "end", "bc"), showProgress = FALSE)
  barcodes <- sort(frags[, .N, by = bc][N >= frag_min, bc])
  frags <- frags[bc %in% barcodes & chr %in% chr_list]
  message(sprintf("barcodes with >= %d fragments: %d; fragments: %d", frag_min, length(barcodes), nrow(frags)))

  window <- fread(w_file, header = FALSE, showProgress = FALSE)
  setnames(window, 1:4, c("wchr", "wstart", "wend", "wname"))
  window[, gid := .I]

  assign_gid <- function(d, wnd) {
    s <- findInterval(d$start, wnd$wstart)
    e <- findInterval(d$end, wnd$wend) + 1L
    straddle <- which(s != e)
    if (length(straddle)) {
      ws <- wnd$wend[s[straddle]]
      s[straddle] <- s[straddle] + as.integer(abs(d$end[straddle] - ws) > abs(d$start[straddle] - ws))
    }
    ok <- s >= 1L & s <= nrow(wnd)
    list(gid = wnd$gid[s[ok]], bc = d$bc[ok])
  }
  hits <- rbindlist(lapply(chr_list, function(cc) {
    d <- frags[chr == cc]; w <- window[wchr == cc]
    if (!nrow(d) || !nrow(w)) return(NULL)
    as.data.table(assign_gid(d, w))
  }))
  cnt <- hits[, .N, by = .(gid, bc)]
  counts <- matrix(0L, nrow = nrow(window), ncol = length(barcodes), dimnames = list(NULL, barcodes))
  counts[cbind(cnt$gid, match(cnt$bc, barcodes))] <- cnt$N
  message(sprintf("count matrix: %d windows x %d cells", nrow(window), length(barcodes)))

  out <- cbind(as.data.frame(window[, 1:4]), as.data.frame(counts))
  names(out)[1:4] <- paste0("V", 1:4)
  write.table(out, file = out_file, row.names = FALSE, col.names = TRUE, quote = FALSE)
  invisible(out_file)
}
