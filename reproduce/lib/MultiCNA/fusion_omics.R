# Fuse scores from an omics-by-bins matrix and change-points 
# X   : matrix (n_omics x n_bins), rows = modalities (e.g., DNAm, RNA), cols = ordered bins (already centered)
# cps : integer vector of bin indices that are change-points (1..n_bins)
# method : "count" (n^alpha) or "v" (|v|^alpha, global per-omic)
# alpha  : exponent for raw weights
# v      : named numeric vector for method="v" (names must match rownames(X))
#
# Returns data.frame per segment with means, counts, weights, evidence, and fused score z
fusion_omics <- function(
    X, cps, 
    #ref,
    method = c("count", "v"),
    alpha = 1.0
) {
  stopifnot(is.matrix(X), ncol(X) >= 1, nrow(X) >= 1)
  
  if (is.null(rownames(X))) rownames(X) <- paste0("omic", seq_len(nrow(X)))
  
  n_bins <- ncol(X)
  cps <- sort(unique(as.integer(cps)))
  cps <- cps[cps >= 1 & cps <= n_bins]
  cps <- cps[cps != n_bins]  # drop terminal index if present
  
  starts <- c(1L, cps + 1L)
  ends   <- c(cps, n_bins)
  N <- length(starts)
  
  out <- data.frame(
    seg_id = seq_len(N),
    start_bin = starts,
    end_bin = ends,
    #start_genomic = ref[starts,]$start,
    #end_genomic = ref[ends,]$end,
    length_bins = ends - starts + 1L,
    stringsAsFactors = FALSE
  )
  
  # Per-omic means and counts (NA-safe)
  for (o in rownames(X)) {
    mns <- numeric(N); cts <- integer(N)
    for (i in seq_len(N)) {
      sli <- as.numeric(X[o, starts[i]:ends[i]])
      cts[i] <- sum(!is.na(sli))
      mns[i] <- if (cts[i] > 0) mean(sli, na.rm = TRUE) else NA_real_
    }
    out[[paste0("mean_", o)]] <- mns
    out[[paste0("n_", o)]]    <- cts
  }
  
  # Per-omic counts, means, sum of squares (NA-safe)
  n_mat   <- matrix(0L,  nrow = N, ncol = nrow(X), dimnames = list(NULL, rownames(X)))
  mean_mat<- matrix(NA_real_, nrow = N, ncol = nrow(X), dimnames = list(NULL, rownames(X)))
  ssq_mat <- matrix(0.0, nrow = N, ncol = nrow(X), dimnames = list(NULL, rownames(X)))
  
  for (o in rownames(X)) {
    nvec <- integer(N); mvec <- numeric(N); ssv <- numeric(N)
    for (i in seq_len(N)) {
      sli <- as.numeric(X[o, starts[i]:ends[i]])
      ok  <- !is.na(sli)
      nvec[i] <- sum(ok)
      mvec[i] <- if (nvec[i] > 0) mean(sli[ok]) else NA_real_
      ssv[i]  <- if (nvec[i] > 0) sum(sli[ok]^2) else 0.0
    }
    out[[paste0("n_", o)]]    <- nvec
    out[[paste0("mean_", o)]] <- mvec
    out[[paste0("sumsq_", o)]]<- ssv
    
    n_mat[,o]    <- nvec
    mean_mat[,o] <- mvec
    ssq_mat[,o]  <- ssv
  }
  
  # Raw weights \tilde w
  tilde <- matrix(0, nrow = N, ncol = nrow(X), dimnames = list(NULL, rownames(X)))
  if (method == "count") {
    for (o in rownames(X)) {
      nvec <- out[[paste0("n_", o)]]
      tilde[, o] <- ifelse(nvec > 0, nvec^alpha, 0)
    }
  } 
  if (method == "v") {
    for (i in 1:N) {
      x.cusum <- cusum(X[,starts[i]:ends[i]])
      p <- nrow(x.cusum)
      n <- ends[i] - starts[i]
      #n <- ncol(X)
      if (is.null(lambda)) lambda <- sqrt(n * log(p * log(n))) / 2
      vhat <- get.v.proj_l2(x.cusum, lambda)
      tilde[i,] <- abs(vhat)
    }
  }
  
  # Normalize weights per segment (sum to 1 over observed modalities)
  wsum <- rowSums(tilde)
  w <- tilde / wsum
  w[!is.finite(w)] <- NA_real_  # segments with no data in any omic
  
  # Store weights and evidence
  for (o in rownames(X)) out[[paste0("w_", o)]] <- w[, o]
  #out$evidence <- wsum
  
  # Fused score z_s
  z <- rep(NA_real_, N)
  for (i in seq_len(N)) {
    if (is.finite(wsum[i]) && wsum[i] > 0) {
      z_i <- 0
      for (o in rownames(X)) {
        m_oi <- out[[paste0("mean_", o)]][i]
        w_oi <- out[[paste0("w_", o)]][i]
        if (is.finite(w_oi) && !is.na(m_oi)) z_i <- z_i + w_oi * m_oi
      }
      z[i] <- z_i
    }
  }
  out$z <- z
  
  ## Weighted sufficient stats for EM on summaries
  w_filled <- w; w_filled[!is.finite(w_filled)] <- 0
  # sumsq_w_i = Σ_o w_{io}^2 * Σ_t x_{iot}^2
  sumsq_w <- rowSums((w_filled^2) * ssq_mat)

  out$sum.sq<- sumsq_w

  return(out)
}
