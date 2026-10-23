multiInverse_chr <- function(x, ref, h = 10, hh = 2*h, precise = 1000, simT,
                             lambda = NULL, penalty = "l2", standardize.series = FALSE,
                             min_pool = 10 * precise) {
  stopifnot(ncol(x) == nrow(ref))
  n <- ncol(x); p <- nrow(x)
  
  # chromosome grouping (assumes ref rows align to x columns)
  chr <- as.character(ref$chr)
  chr_levels <- unique(chr)                # keep original order in ref
  idx_by_chr <- split(seq_len(n), chr)
  
  # preserve missingness pattern from real data
  miss <- is.na(x)
  
  # row-wise SD from real data (global across genome, as in your code)
  sd_row <- apply(x, 1, function(v) {
    v <- v[is.finite(v)]
    if (length(v) < 2) return(NA_real_)
    stats::sd(v)
  })
  
  empirical <- numeric(0)
  
  for (i in seq_len(simT)) {
    message("SimT: ", i)
    start_time <- Sys.time()
    
    # simulate null on full series then impose missingness
    x_null <- matrix(
      rnorm(p * n, mean = 0, sd = rep(sd_row, each = n)),
      nrow = p, ncol = n
    )
    x_null[miss] <- NA
    if (standardize.series) x_null <- rescale.variance(x_null)
    
    # run per chromosome, collect peak heights
    for (cc in chr_levels) {
      idx <- idx_by_chr[[cc]]
      n_chr <- length(idx)
      
      # localMax requires >= 2*hh+1 points
      if (n_chr < (2 * hh + 1)) next
      
      x_chr <- x_null[, idx, drop = FALSE]
      
      Dobj <- multiMissInspect_scores(
        x_chr, lambda = lambda, h = h,
        standardize.series = FALSE, penalty = penalty
      )$D
      
      s <- as.numeric(Dobj$score)
      
      pos_idx <- localMax(+s, span = hh)
      neg_idx <- localMax(-s, span = hh)
      
      id_all <- c(pos_idx, neg_idx)
      if (length(id_all)) {
        val <- abs(c(s[pos_idx], s[neg_idx]))
        o <- order(id_all); id_all <- id_all[o]; val <- val[o]
        
        # prune peaks closer than h (within this chromosome)
        if (length(id_all) > 1) {
          rm <- which(diff(id_all) < h)
          if (length(rm)) val <- val[-rm]
        }
        
        empirical <- c(empirical, val)
        
        if (length(empirical) >= min_pool) break
      }
    }
    
    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    message("Done, Seconds elapsed: ", round(elapsed, 3))
    
    if (length(empirical) >= min_pool) break
  }
  
  qgrid <- (0:precise) / precise
  qvec  <- as.numeric(stats::quantile(empirical, probs = qgrid, type = 8, na.rm = TRUE))
  data.frame(probs = qgrid, q = qvec)
}

localMax <- function(y, span) {
  if (length(y) < span * 2 + 1) return(NULL)
  n <- length(y)
  index <- integer(0)
  for (i in (span + 1):(n - span)) {
    if (y[i] == max(y[(i - span):(i + span)])) index <- c(index, i)
  }
  index
}
