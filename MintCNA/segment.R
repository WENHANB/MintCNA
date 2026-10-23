# MintCNA - Segmentation: multivariate CUSUM scan, empirical null, segmentation driver
# Source-able R functions; sourced recursively by scripts/run_mintcna.R.


# ======================= scan.R =======================

# MintCNA - CUSUM-based Change-Point Scoring
# Core multivariate change-point detection via MissCUSUM + sparse rank-1 projection

library(softImpute)
library(InspectChangepoint)

# Source utilities (vector.norm, vector.normalise, vector.soft.thresh)
# These must be loaded before this file.

#' Change-point score computation
#'
#' Computes diagnostic scores by combining MissCUSUM with a sparse rank-1
#' projection. For each time point, a local CUSUM is calculated, projected to 1D
#' using either an l1 or l2 penalized vector, and the deviation at the center is
#' recorded.
#'
#' @param x Numeric matrix (p x n): modalities x bins, with possible NAs.
#' @param lambda Regularization parameter. Default: sqrt(n * log(p * log(n))) / 2, with n the number of bins in the series being scanned.
#' @param h Half-window size for local CUSUM and local maxima search. Default: 10.
#' @param standardize.series Logical; if TRUE, rescale variance of each row. Default: FALSE.
#' @param penalty Projection penalty: "l1" or "l2". Default: "l2".
#' @return List with:
#'   $D: data frame of all scores (position, score);
#'   $local_maxima: subset of D that are local maxima;
#'   $v_df: data frame of projection vectors per bin.
#' @param keep_v Return `v_df`, the per-position projection vectors. Building it
#'   costs ~40% of this function's runtime (one data.frame per scanned position,
#'   then a do.call(rbind)) and `scan_cells()` discards it immediately, so the
#'   pipeline sets FALSE. Default TRUE keeps the original return shape.
multiMissInspect_scores <- function(x, lambda = NULL, h = 10,
                                    standardize.series = FALSE, penalty = "l2",
                                    keep_v = TRUE) {
  p <- nrow(x)
  # n is the LENGTH OF THE SERIES, not the window. Setting n <- 2*h (as an
  # earlier version of this file did) makes the loop bound seq(h+1, n-h) run
  # backwards from h+1 to h, so the scan visits two positions instead of the
  # whole chromosome, and it changes lambda as well. Matches Code/MultiCNA.
  n <- ncol(x)
  if (is.null(lambda)) lambda <- sqrt(n * log(p * log(n))) / 2

  # Preallocate: growing a data.frame with rbind inside the loop reallocates and
  # copies on every iteration, which is O(n^2) in the number of bins and is run
  # once per chromosome per cell.
  ts    <- seq(h + 1, n - h)
  score <- numeric(length(ts))
  V_list <- if (keep_v) vector("list", length(ts)) else NULL

  for (j in seq_along(ts)) {
    t <- ts[j]
    x_window <- x[, (t - h):(t + h)]
    x.cusum <- cusum(x_window)

    vhat <- switch(penalty,
                   l1 = get.v.proj_l1(x.cusum, lambda),
                   l2 = get.v.proj_l2(x.cusum, lambda),
                   stop("penalty must be 'l1' or 'l2'"))

    score[j] <- as.vector(vhat %*% x.cusum)[h + 1]
    if (keep_v) V_list[[j]] <- data.frame(bin = t, t(vhat))
  }
  D <- data.frame(position = ts, score = score)

  v_df <- if (keep_v) do.call(rbind, V_list) else NULL

  is_local_max <- function(vec, idx, h) {
    left  <- max(1, idx - h)
    right <- min(length(vec), idx + h)
    all(abs(vec[idx]) >= abs(vec[left:right]))
  }

  local_maxima <- D[sapply(1:nrow(D), function(i) is_local_max(D$score, i, h)), ]

  list(D = D, local_maxima = local_maxima, v_df = v_df)
}

#' CUSUM transformation for a single time series with missing data
#'
#' @param x Numeric vector (may contain NAs).
#' @return CUSUM-transformed vector of length n-1.
cusum.vector <- function(x) {
  ob <- !is.na(x)
  z <- replace(x, !ob, 0)
  n <- length(x)
  leftsum  <- cumsum(z)
  rightsum <- leftsum[n] - leftsum
  L <- cumsum(ob)
  R <- L[n] - L
  x.cusum <- (leftsum[-n] / L[-n] - rightsum[-n] / R[-n]) * sqrt(L[-n] * R[-n] / L[n])
  x.cusum[is.nan(x.cusum)] <- 0
  x.cusum
}

#' Apply CUSUM to vector or matrix
#'
#' @param x Numeric vector or matrix (modalities x bins).
#' @return CUSUM-transformed matrix.
cusum <- function(x) {
  if (!is.matrix(x)) {
    return(cusum.vector(x))
  } else if (ncol(x) == 2) {
    return(as.matrix(apply(x, 1, cusum.vector)))
  } else {
    return(t(apply(x, 1, cusum.vector)))
  }
}

#' L1-penalized (Lasso) rank-1 projection
#'
#' @param x.cusum CUSUM matrix (p x n-1).
#' @param lambda Regularization parameter.
#' @param max_iter Maximum iterations.
#' @param tol Convergence tolerance.
#' @return Projection vector (length p).
get.v.proj_l1 <- function(x.cusum, lambda, max_iter = 1000, tol = 1e-10) {
  vhat <- svd(x.cusum)$u[, 1]
  if (sum(abs(x.cusum)) == 0) return(random.UnitVector(nrow(x.cusum)))
  for (iter in 1:max_iter) {
    vhat_old <- vhat
    what <- vector.normalise(t(x.cusum) %*% vhat)
    tmp <- x.cusum %*% what
    lambda_tmp <- min(lambda, max(abs(tmp)) - 1e-10)
    vhat <- vector.normalise(vector.soft.thresh(x.cusum %*% what, lambda_tmp))
    if (vector.norm(vhat_old - vhat) < tol) break
  }
  as.vector(vhat)
}

#' L2-penalized (Ridge) rank-1 projection
#'
#' @param x.cusum CUSUM matrix (p x n-1).
#' @param lambda Regularization parameter.
#' @param max_iter Maximum iterations.
#' @param tol Convergence tolerance.
#' @return Projection vector (length p).
get.v.proj_l2 <- function(x.cusum, lambda, max_iter = 1000, tol = 1e-10) {
  if (sum(abs(x.cusum)) == 0) return(random.UnitVector(nrow(x.cusum)))
  vhat <- svd(x.cusum)$u[, 1]
  for (iter in 1:max_iter) {
    vhat_old <- vhat
    what <- vector.normalise(t(x.cusum) %*% vhat)
    vhat <- (1 / (2 * lambda)) * (x.cusum %*% what)
    if (vector.norm(vhat_old - vhat) < tol) break
  }
  as.vector(vhat)
}

#' Generate a random unit vector
#'
#' @param n Dimension.
#' @return Random unit vector of length n.
random.UnitVector <- function(n) {
  v <- rnorm(n)
  v / vector.norm(v)
}

# MintCNA - Univariate Change-Point Scoring (DNA-only mode)
# CUSUM-based scoring for a single 1D signal (p=1)

#' Univariate change-point scoring for DNA-only mode
#'
#' Computes CUSUM-based scores for a single 1D signal. At each position t,
#' a local CUSUM is computed in a window of size 2h+1, and the absolute value
#' at the center is the score. No multivariate projection is needed.
#'
#' This is the DNA-only counterpart of multiMissInspect_scores().
#'
#' @param x Numeric vector of length n (one cell's signal along bins).
#' @param h Half-window size for local CUSUM and local maxima search. Default: 10.
#' @return List with:
#'   $D: data frame with columns position, score (all positions);
#'   $local_maxima: subset of D that are local maxima within +/- h.
univarInspect_scores <- function(x, h = 10) {
  if (is.matrix(x)) {
    if (nrow(x) == 1) {
      x <- as.numeric(x)
    } else {
      stop("univarInspect_scores expects a vector or 1-row matrix. Use multiMissInspect_scores for p > 1.")
    }
  }
  n <- length(x)
  if (n < 2 * h + 2) {
    warning("Signal too short for h = ", h, ". Returning empty results.")
    return(list(
      D = data.frame(position = integer(0), score = numeric(0)),
      local_maxima = data.frame(position = integer(0), score = numeric(0))
    ))
  }

  # Preallocate rather than rbind into a growing data.frame; see the note in
  # multiMissInspect_scores. Output is identical.
  ts    <- seq(h + 1, n - h)
  score <- numeric(length(ts))
  for (j in seq_along(ts)) {
    x_window <- x[(ts[j] - h):(ts[j] + h)]
    score[j] <- cusum.vector(x_window)[h + 1]
  }
  D <- data.frame(position = ts, score = score)

  # Local maxima of |score| within +/- h
  is_local_max <- function(vec, idx, h) {
    left  <- max(1, idx - h)
    right <- min(length(vec), idx + h)
    all(abs(vec[idx]) >= abs(vec[left:right]))
  }

  local_maxima <- D[sapply(1:nrow(D), function(i) is_local_max(D$score, i, h)), ]

  list(D = D, local_maxima = local_maxima)
}

# MintCNA - Math Utilities
# Vector operations used by the segmentation module

#' L_q norm of a vector
#'
#' @param v Numeric vector.
#' @param q Norm order (default 2). Use Inf for max norm, 0 for count of nonzero.
#' @param na.rm Remove NAs before computation.
#' @return Scalar norm value.
vector.norm <- function(v, q = 2, na.rm = FALSE) {
  if (na.rm) v <- na.omit(v)
  M <- max(abs(v))
  if (M == 0) return(0)
  v <- v / M
  if (q == Inf) {
    nm <- max(abs(v))
  } else if (q > 0) {
    nm <- (sum(abs(v)^q))^(1 / q)
  } else if (q == 0) {
    nm <- sum(v != 0)
  } else {
    return(NaN)
  }
  nm * M
}

#' Normalize a vector to unit L_q norm
#'
#' @param v Numeric vector.
#' @param q Norm order (default 2).
#' @param na.rm Remove NAs before computation.
#' @return Normalized vector.
vector.normalise <- function(v, q = 2, na.rm = FALSE) {
  v / vector.norm(v, q, na.rm)
}

#' Soft thresholding
#'
#' Entries are moved toward zero by lambda until they hit zero.
#'
#' @param x Numeric vector.
#' @param lambda Threshold value.
#' @return Soft-thresholded vector.
vector.soft.thresh <- function(x, lambda) {
  sign(x) * pmax(0, (abs(x) - lambda))
}

# ======================= null.R =======================

# MintCNA - Empirical Null Distribution for Change-Points
# Simulates null change-point scores per chromosome to build quantile lookup

#' Simulate null change-point distribution
#'
#' Generates Gaussian null data matching the variance and missingness pattern
#' of the real data, runs change-point scoring per chromosome, and collects
#' peak heights to build an empirical quantile table.
#'
#' @param x Numeric matrix (modalities x bins), with possible NAs.
#' @param ref Data frame with column 'chr' (one row per bin, aligned to x columns).
#' @param h Half-window size for CUSUM scan. Default: 10.
#' @param hh Span for local maxima detection. Default: 2*h.
#' @param precise Number of quantile grid points. Default: 1000.
#' @param simT Number of simulation iterations. Default: 10.
#' @param lambda Regularization parameter (NULL = auto).
#' @param penalty "l1" or "l2". Default: "l2".
#' @param standardize.series Logical. Default: FALSE.
#' @param min_pool Minimum number of null peaks to collect. Default: 10 * precise.
#' @return Data frame with columns 'probs' and 'q' (quantile values).
multiInverse_chr <- function(x, ref, h = 10, hh = 2 * h, precise = 1000, simT,
                             lambda = NULL, penalty = "l2",
                             standardize.series = FALSE,
                             min_pool = 10 * precise) {
  stopifnot(ncol(x) == nrow(ref))
  n <- ncol(x)
  p <- nrow(x)

  chr <- as.character(ref$chr)
  chr_levels <- unique(chr)
  idx_by_chr <- split(seq_len(n), chr)

  # Preserve missingness pattern
  miss <- is.na(x)

  # Row-wise SD from real data
  sd_row <- apply(x, 1, function(v) {
    v <- v[is.finite(v)]
    if (length(v) < 2) return(NA_real_)
    stats::sd(v)
  })

  empirical <- numeric(0)

  for (i in seq_len(simT)) {
    message("SimT: ", i)
    start_time <- Sys.time()

    # Simulate null data with matching variance and missingness
    # times = n, NOT each = n. matrix() fills COLUMN-major, so element (i,j)
    # takes the ((j-1)*p + i)-th value of the sd vector. `each = n` instead lays
    # the SDs out along the GENOME: with sd_row = (RNA 0.278, DNA 0.691) and
    # n = 5259, the first 2630 bins got 0.278 in BOTH rows and the rest 0.691 in
    # both, giving each row an effective SD of 0.527 that varies by position
    # rather than by modality. Measured effect on P111T3 at h = 15, simT = 10:
    # tau(alpha = 0.10) 0.585 -> 0.504, tau(alpha = 0.05) 0.722 -> 0.595, i.e.
    # the threshold was 14-18% too strict and change-points were being lost.
    x_null <- matrix(
      rnorm(p * n, mean = 0, sd = rep(sd_row, times = n)),
      nrow = p, ncol = n
    )
    x_null[miss] <- NA

    # Score per chromosome, collect peak heights
    for (cc in chr_levels) {
      idx   <- idx_by_chr[[cc]]
      n_chr <- length(idx)
      if (n_chr < (2 * hh + 1)) next

      x_chr <- x_null[, idx, drop = FALSE]
      Dobj  <- multiMissInspect_scores(
        x_chr, lambda = lambda, h = h,
        standardize.series = FALSE, penalty = penalty
      )$D

      s <- as.numeric(Dobj$score)
      pos_idx <- localMax(+s, span = hh)
      neg_idx <- localMax(-s, span = hh)

      id_all <- c(pos_idx, neg_idx)
      if (length(id_all)) {
        val <- abs(c(s[pos_idx], s[neg_idx]))
        o <- order(id_all)
        id_all <- id_all[o]
        val <- val[o]

        # Prune peaks closer than h
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

#' Find local maxima in a 1D vector
#'
#' @param y Numeric vector.
#' @param span Half-window: a point is a local max if it equals the max in y[(i-span):(i+span)].
#' @return Integer vector of indices that are local maxima.
localMax <- function(y, span) {
  if (length(y) < span * 2 + 1) return(NULL)
  n <- length(y)
  index <- integer(0)
  for (i in (span + 1):(n - span)) {
    if (y[i] == max(y[(i - span):(i + span)])) index <- c(index, i)
  }
  index
}

# MintCNA - Univariate Null Distribution (DNA-only mode)
# Simulates null change-point scores for a 1D signal per chromosome

#' Simulate null change-point distribution for univariate (DNA-only) mode
#'
#' Generates Gaussian null data matching the variance of the real signal,
#' runs univariate CUSUM scoring per chromosome, and collects peak heights
#' to build an empirical quantile table.
#'
#' This is the DNA-only counterpart of multiInverse_chr().
#'
#' @param x Numeric vector of length n (one cell's signal), or a 1-row matrix.
#' @param ref Data frame with column 'chr' (one row per bin, aligned to x).
#' @param h Half-window size for CUSUM scan. Default: 10.
#' @param hh Span for local maxima detection. Default: 2*h.
#' @param precise Number of quantile grid points. Default: 1000.
#' @param simT Number of simulation iterations. Default: 10.
#' @param min_pool Minimum number of null peaks to collect. Default: 10 * precise.
#' @return Data frame with columns 'probs' and 'q' (quantile values).
univarInverse_chr <- function(x, ref, h = 10, hh = 2 * h, precise = 1000, simT = 10,
                              min_pool = 10 * precise) {
  # Handle 1-row matrix input
  if (is.matrix(x)) {
    if (nrow(x) == 1) {
      x <- as.numeric(x)
    } else {
      stop("univarInverse_chr expects a vector or 1-row matrix. Use multiInverse_chr for p > 1.")
    }
  }

  stopifnot(length(x) == nrow(ref))
  n <- length(x)

  chr <- as.character(ref$chr)
  chr_levels <- unique(chr)
  idx_by_chr <- split(seq_len(n), chr)

  # Variance from real data (excluding NAs)
  x_finite <- x[is.finite(x)]
  sd_x <- if (length(x_finite) >= 2) sd(x_finite) else 1.0

  empirical <- numeric(0)

  for (i in seq_len(simT)) {
    message("SimT: ", i)
    start_time <- Sys.time()

    # Simulate null: Gaussian with same variance
    x_null <- rnorm(n, mean = 0, sd = sd_x)
    # Preserve missingness
    x_null[is.na(x)] <- NA

    for (cc in chr_levels) {
      idx   <- idx_by_chr[[cc]]
      n_chr <- length(idx)
      if (n_chr < (2 * hh + 1)) next

      x_chr <- x_null[idx]

      Dobj <- univarInspect_scores(x_chr, h = h)$D
      s <- as.numeric(Dobj$score)

      pos_idx <- localMax(+s, span = hh)
      neg_idx <- localMax(-s, span = hh)

      id_all <- c(pos_idx, neg_idx)
      if (length(id_all)) {
        val <- abs(c(s[pos_idx], s[neg_idx]))
        o <- order(id_all)
        id_all <- id_all[o]
        val <- val[o]

        # Prune peaks closer than h
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

# ======================= segment.R =======================

# MintCNA - Top-Level Segmentation Wrapper
# Dispatches between DNA-only and multi-omics segmentation pipelines
#
# The pipeline is split at the seam that matters: SCANNING is expensive and does
# not depend on `alpha`, `st` or the truncation bounds, while THRESHOLDING and
# CALLING are cheap and depend on nothing else. Keeping them in one function --
# as an earlier version of this file did -- means every threshold experiment
# re-runs the scan. `Simulation/code/MintCNA/MintCNA.R` already had this split,
# with a save() between the two loops; this restores it.
#
#   scan_cells()  -> CUSUM scores per cell per chromosome   [expensive, cache it]
#   call_cells()  -> threshold at tau(alpha) -> TGMM        [seconds]
#
# Dependencies: source all segmentation files before this one
# - utils.R, qn_transform.R, cusum_scan.R, univariate_scan.R
# - null_distribution.R, univariate_null.R
# - ../fusion/fuse_omics.R, ../calling/tgmm.R, ../calling/cnv_cluster.R
# - ../output/cna_output.R

#' Detect the segmentation mode from the shape of one cell's data
#'
#' @param x A cell's data: matrix (p x n_bins) for multi-omics, vector or 1-row
#'   matrix for DNA-only.
#' @return "multi" or "dna".
detect_mode <- function(x) if (is.matrix(x) && nrow(x) > 1) "multi" else "dna"

#' Coerce one cell's data to the p x n_bins matrix the scan expects
#' @param cell_data Vector, 1-row matrix, or p-row matrix.
#' @param mode "dna" or "multi".
as_cell_matrix <- function(cell_data, mode) {
  if (mode == "dna") {
    v <- if (is.matrix(cell_data)) {
      if (nrow(cell_data) != 1) stop("DNA-only mode expects a vector or 1-row matrix.")
      as.numeric(cell_data)
    } else if (is.numeric(cell_data) && is.null(dim(cell_data))) {
      cell_data
    } else stop("DNA-only mode expects a vector or 1-row matrix.")
    matrix(v, nrow = 1, dimnames = list("DNA", NULL))
  } else {
    if (!is.matrix(cell_data) || nrow(cell_data) < 2)
      stop("Multi-omics mode expects a matrix with >= 2 rows (modalities).")
    cell_data
  }
}

# ---------------------------------------------------------------------------
# STAGE 1 — scanning
# ---------------------------------------------------------------------------

#' CUSUM scan for every cell — the expensive, alpha-independent half
#'
#' Computes the shared empirical null once and then the per-chromosome
#' change-point scores for every cell. Nothing here depends on `alpha`, `st`, or
#' the TGMM truncation bounds, so the result is worth caching: every threshold
#' and calling experiment can start from it.
#'
#' @param cell_list Named list of cell data (already QN-transformed).
#' @param ref Data frame with chr, start, end — one row per bin.
#' @param mode "dna" or "multi". Auto-detected from the first cell if NULL.
#' @param h Half-window size for the CUSUM scan.
#' @param penalty Projection penalty for multi-omics: "l1" or "l2".
#' @param simT Simulation replicates for the empirical null. The null sets tau
#'   and therefore every change-point, and at simT = 1 its 95th percentile has
#'   ~16% Monte-Carlo CV; simT = 10 brings that to ~5% for a one-off cost.
#'   Default 1 matches Simulation/code/MintCNA/MintCNA.R.
#' @param seed RNG seed for the null simulation, so a scan is reproducible.
#' @param keep_D Keep the full score vector per chromosome as well as the local
#'   maxima. `D` is only needed for diagnostic plots and dominates the file size.
#' @param mc_cores Cells are independent once the null is fixed; > 1 runs them
#'   through parallel::mclapply.
#' @param file Optional path to save the result to (RDS).
#' @return List with $cps_track (cell -> chr -> scores), $res_null, $ref, $h,
#'   $penalty, $mode, $timing, and $cell_list (carried so calling needs nothing else).
scan_cells <- function(cell_list, ref, mode = NULL, h = 5, penalty = "l2",
                       simT = 1, seed = 1, keep_D = FALSE, mc_cores = 1,
                       file = NULL) {
  if (length(cell_list) == 0) stop("cell_list is empty.")
  if (is.null(mode)) mode <- detect_mode(cell_list[[1]])
  message(sprintf("[MintCNA] scan: %s mode, %d cells, h = %d.",
                  mode, length(cell_list), h))

  # The null depends on the variance structure, not on the signal, and QN has
  # already forced every cell to a common sd -- so one cell is representative.
  # The null is simulated, so without a seed two scans of the same data give
  # different thresholds -- 2.23 vs 1.69 on a toy example at simT = 1.
  if (!is.null(seed)) set.seed(seed)
  message(sprintf("[MintCNA] scan: empirical null (simT = %d, seed = %s)...",
                  simT, if (is.null(seed)) "none" else seed))
  ref_cell <- as_cell_matrix(cell_list[[length(cell_list)]], mode)
  res_null <- if (mode == "dna") {
    univarInverse_chr(as.numeric(ref_cell), ref, h = h, simT = simT)
  } else {
    multiInverse_chr(ref_cell, ref, h = h, simT = simT, penalty = penalty)
  }

  chr_levels <- unique(as.character(ref$chr))
  idx_by_chr <- split(seq_len(nrow(ref)), factor(ref$chr, levels = chr_levels))

  scan_one <- function(cell_name) {
    t0 <- Sys.time()
    cell_mat <- as_cell_matrix(cell_list[[cell_name]], mode)
    out <- list()
    # A chromosome shorter than the scan window yields no candidates. Store an
    # empty local_maxima rather than a sentinel, so the calling stage needs no
    # special case.
    empty <- list(local_maxima = data.frame(position = integer(0), score = numeric(0)))
    for (chr in chr_levels) {
      idx <- idx_by_chr[[chr]]
      if (length(idx) < 2 * h + 2) { out[[chr]] <- empty; next }
      sc <- if (mode == "dna") {
        univarInspect_scores(cell_mat[1, idx], h = h)
      } else {
        multiMissInspect_scores(cell_mat[, idx, drop = FALSE], h = h, penalty = penalty,
                                keep_v = FALSE)   # v_df is discarded below anyway
      }
      if (!keep_D) sc$D <- NULL
      out[[chr]] <- sc
    }
    list(scores = out, seconds = as.numeric(difftime(Sys.time(), t0, units = "secs")))
  }

  nm <- names(cell_list)
  cps_track <- if (mc_cores > 1) {
    parallel::mclapply(nm, scan_one, mc.cores = mc_cores)
  } else {
    lapply(nm, function(n) { message("[MintCNA] scan: ", n); scan_one(n) })
  }
  names(cps_track) <- nm

  # mclapply does not abort when a forked worker dies (OOM, segfault): it returns
  # a try-error in that slot. Silently writing such a checkpoint would drop cells
  # from every downstream result, so fail loudly here instead.
  bad <- !vapply(cps_track, function(x) is.list(x) && !is.null(x$seconds), logical(1))
  if (any(bad)) {
    msg <- vapply(cps_track[bad][seq_len(min(3, sum(bad)))],
                  function(x) paste(conditionMessage(attr(x, "condition")), collapse = " "),
                  character(1))
    stop(sprintf("%d of %d cells failed to scan (e.g. %s). First: %s",
                 sum(bad), length(nm), paste(nm[bad][seq_len(min(3, sum(bad)))], collapse = ", "),
                 paste(msg, collapse = " | ")))
  }
  timing <- data.frame(cell = nm,
                       seconds = vapply(cps_track, function(x) x$seconds, numeric(1)),
                       stringsAsFactors = FALSE)
  cps_track <- lapply(cps_track, function(x) x$scores)   # timing lives in `timing`
  names(cps_track) <- nm
  message(sprintf("[MintCNA] scan done: %.2f s per cell, %.1f min total.",
                  median(timing$seconds), sum(timing$seconds) / 60))

  out <- list(cps_track = cps_track, res_null = res_null, ref = ref, h = h,
              penalty = penalty, mode = mode, timing = timing,
              cell_list = cell_list, chr_levels = chr_levels,
              simT = simT, seed = seed)
  if (!is.null(file)) {
    saveRDS(out, file = file)
    message("[MintCNA] scan checkpoint written: ", file)
  }
  out
}

# ---------------------------------------------------------------------------
# STAGE 2 — thresholding and calling
# ---------------------------------------------------------------------------

#' Threshold a cached scan and call CNA states — the cheap half
#'
#' @param scan Output of scan_cells(), or a path to its RDS.
#' @param alpha Quantile threshold. tau is the (1 - alpha) quantile of the null,
#'   so SMALLER alpha means a stricter threshold and fewer change-points.
#' @param st Number of CNA states (3 or 5).
#' @param L Maximum EM iterations.
#' @param eps0,eps1 Neutral-band bounds for st = 3 (passed to updateEM).
#' @param lower,upper Explicit truncation bounds, overriding eps0/eps1.
#' @param mc_cores Cells are independent here too -- tau is a scalar fixed before
#'   the loop and each cell's TGMM sees only its own data -- so > 1 runs them
#'   through parallel::mclapply, as in scan_cells().
#' @param min_bins Drop segments whose (end_idx - start_idx) is not greater than
#'   this. Default 1, i.e. keep segments of at least 3 bins, matching
#'   Simulation/code/MintCNA/MintCNA.R.
#' @param file Optional path to save the result to (RData).
#' @return List with $cna_df_all, $cna_matrix, $converged_cells, $res_null,
#'   $per_cell, and the parameters used.
call_cells <- function(scan, alpha = 0.05, st = 3, L = 30,
                       eps0 = 0.30, eps1 = 0.30, lower = NULL, upper = NULL,
                       min_bins = 1, mc_cores = 1, file = NULL) {
  if (is.character(scan)) scan <- readRDS(scan)
  ref <- scan$ref; h <- scan$h; mode <- scan$mode
  chr_levels <- scan$chr_levels
  idx_by_chr <- split(seq_len(nrow(ref)), factor(ref$chr, levels = chr_levels))
  chr_lengths <- vapply(chr_levels, function(c) length(idx_by_chr[[c]]), integer(1))

  tau <- scan$res_null$q[scan$res_null$probs == (1 - alpha)]
  if (!length(tau) || !is.finite(tau))
    stop(sprintf("alpha = %s is not on the null's probability grid.", alpha))
  message(sprintf("[MintCNA] call: alpha = %s -> tau = %.3f; st = %d.", alpha, tau, st))

  call_one <- function(cell_name) {
    cell_mat <- as_cell_matrix(scan$cell_list[[cell_name]], mode)
    sc_by_chr <- scan$cps_track[[cell_name]]

    cps_by_chr <- list()
    for (chr in chr_levels) {
      lm <- sc_by_chr[[chr]]$local_maxima
      cps_by_chr[[chr]] <- sort(lm[abs(lm$score) > tau, "position"])
    }

    cps_global <- reindex_cps(cps_by_chr, chr_lengths)
    cps_global <- cps_global[cps_global > 1 & cps_global < ncol(cell_mat)]

    all_normal <- data.frame(state = "nor", chr = ref$chr[1], start = ref$start[1],
                             end = ref$end[nrow(ref)], start_idx = 1,
                             end_idx = nrow(ref), stringsAsFactors = FALSE)
    if (length(cps_global) == 0) {
      cna_table <- all_normal
      cnv_result <- NULL
    } else {
      cnv_result <- CNVcluster(Y = cell_mat, cp = cps_global, L = L, st = st,
                               eps0 = eps0, eps1 = eps1, lower = lower, upper = upper)
      # change-points were found but every segment was called neutral
      cna_table <- if (!length(cnv_result$cnv.state)) all_normal else
        make_cna_table(cnv_result$cnv.state, cnv_result$cnv.start,
                       cnv_result$cnv.end, ref)
    }
    # Segment-level statistics, exactly as Simulation/code/MintCNA/MintCNA.R
    # computes them: `mean` is over ALL modalities in the segment, and segments
    # spanning min_bins bins or fewer are dropped.
    cna_table$mean <- mapply(function(s, e) mean(cell_mat[, s:e], na.rm = TRUE),
                             cna_table$start_idx, cna_table$end_idx)
    cna_table$cell <- cell_name
    cna_table$bin_length <- cna_table$end_idx - cna_table$start_idx
    cna_table <- cna_table[cna_table$bin_length > min_bins, , drop = FALSE]

    list(cps = cps_by_chr, cps_global = cps_global, cnv_result = cnv_result,
         cna_table = cna_table,
         converged = !is.null(cnv_result) && cnv_result$iter < L)
  }

  nm <- names(scan$cps_track)
  per_cell <- if (mc_cores > 1) parallel::mclapply(nm, call_one, mc.cores = mc_cores)
              else lapply(nm, call_one)
  names(per_cell) <- nm
  bad <- !vapply(per_cell, function(x) is.list(x) && !is.null(x$cna_table), logical(1))
  if (any(bad)) stop(sprintf("%d of %d cells failed to call (first: %s)",
                             sum(bad), length(nm), nm[bad][1]))

  converged_cells <- nm[vapply(per_cell, function(x) isTRUE(x$converged), logical(1))]
  # one rbind, not 1,658 -- growing a data.frame in the loop is the same O(n^2)
  # antipattern that was fixed in scan.R
  cna_df_all <- data.table::rbindlist(lapply(per_cell, `[[`, "cna_table"))
  cna_df_all <- as.data.frame(cna_df_all[cna_df_all$cell %in% converged_cells, ])
  code_map <- if (st == 3) c(del = 1, nor = 2, dup = 3) else
                           c(del = 1, loss = 2, nor = 3, gain = 4, amp = 5)
  cna_matrix <- if (nrow(cna_df_all) > 0)
    build_cna_matrix(cna_df_all, ref, code_map = code_map) else NULL

  message(sprintf("[MintCNA] call done: %d/%d cells converged, %d segments.",
                  length(converged_cells), length(scan$cps_track), nrow(cna_df_all)))

  out <- list(cna_df_all = cna_df_all, cna_matrix = cna_matrix,
              converged_cells = converged_cells, res_null = scan$res_null,
              per_cell = per_cell, ref = ref,
              params = list(alpha = alpha, st = st, L = L, h = h,
                            eps0 = eps0, eps1 = eps1, min_bins = min_bins))
  if (!is.null(file)) {
    saveRDS(out, file = file)
    message("[MintCNA] call checkpoint written: ", file)
  }
  out
}

# ---------------------------------------------------------------------------
# Convenience wrapper — unchanged interface for callers that want one call
# ---------------------------------------------------------------------------

#' Run MintCNA segmentation and CNA calling on all cells
#'
#' Thin wrapper over scan_cells() + call_cells(), kept so existing callers do not
#' change. Prefer the two-stage form when experimenting with `alpha` or `st`.
#'
#' @inheritParams scan_cells
#' @inheritParams call_cells
#' @return List with $cna_df_all, $cna_matrix, $converged_cells, $res_null, $per_cell.
run_segmentation <- function(cell_list, ref, mode = NULL,
                             h = 5, penalty = "l2", alpha = 0.30,
                             simT = 1, seed = 1, L = 30, st = 3,
                             eps0 = 0.30, eps1 = 0.30, min_bins = 1,
                             mc_cores = 1, scan_file = NULL, call_file = NULL) {
  sc <- scan_cells(cell_list, ref, mode = mode, h = h, penalty = penalty,
                   simT = simT, seed = seed, mc_cores = mc_cores, file = scan_file)
  call_cells(sc, alpha = alpha, st = st, L = L, eps0 = eps0, eps1 = eps1,
             min_bins = min_bins, mc_cores = mc_cores, file = call_file)
}

