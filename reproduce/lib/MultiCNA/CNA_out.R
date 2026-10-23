reindex_cps <- function(cps_multi, chr_lengths, include_first_start = TRUE, include_ends = TRUE) {
  stopifnot(length(cps_multi) == length(chr_lengths))
  chr_lengths <- as.numeric(chr_lengths)
  
  # offsets: 0, L1, L1+L2, ...
  offsets <- c(0, cumsum(chr_lengths))[seq_along(chr_lengths)]
  
  out <- integer(0)
  if (include_first_start) out <- 1L  # only once
  
  for (i in seq_along(cps_multi)) {
    cp_i <- as.numeric(cps_multi[[i]])
    off  <- offsets[i]
    
    if (length(cp_i)) out <- c(out, off + cp_i)                 # CPs for chr i
    if (include_ends) out <- c(out, off + chr_lengths[i])       # end of chr i
  }
  
  sort(unique(out))
}


make_cna_table <- function(cnv.state, cnv.start, cnv.end, ref) {
  stopifnot(all(c("chr","start","end") %in% names(ref)))
  ref$chr <- as.character(ref$chr)
  n <- nrow(ref)
  
  # last row index for each chromosome (assumes ref is ordered by chr, then position)
  chr_ends <- tapply(seq_len(n), ref$chr, max)
  
  split_one <- function(st, s, e) {
    out <- list(); cur <- s
    while (TRUE) {
      chr_s <- ref$chr[cur]
      chr_end_idx <- chr_ends[[chr_s]]
      end_here <- min(e, chr_end_idx)
      out[[length(out)+1]] <- data.frame(
        state = st,
        chr   = chr_s,
        start = ref$start[cur],
        end   = ref$end[end_here],
        start_idx = cur,
        end_idx   = end_here,
        stringsAsFactors = FALSE
      )
      if (end_here == e) break
      cur <- end_here + 1
    }
    do.call(rbind, out)
  }
  
  pieces <- lapply(seq_along(cnv.state), function(i)
    split_one(cnv.state[i], cnv.start[i], cnv.end[i]))
  do.call(rbind, pieces)
}