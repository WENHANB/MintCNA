# MintCNA - Calling: TGMM, cross-omics fusion, output tables, subclone clustering
# Source-able R functions; sourced recursively by scripts/run_mintcna.R.


# ======================= tgmm.R =======================

# MintCNA - Truncated Gaussian Mixture Model (TGMM) CNA Calling
# EM algorithm for assigning integer copy number states to segments

#' EM wrapper for CNA state assignment
#'
#' Runs iterative EM on segment-level fused means to assign copy number states.
#' Calls fusion_omics() for segment statistics and updateEM() for each iteration.
#'
#' @param x Matrix (modalities x bins) or vector of bin intensities.
#' @param cp Integer vector of change-point positions.
#' @param priors List with $p (mixing proportions), $mu (means), $sigma (SDs).
#' @param L Maximum EM iterations. Default: 100.
#' @param st Number of states (3 or 5).
#' @return List with:
#'   $p.final, $mu.final, $sigma.final: converged parameters;
#'   $cp.final: refined change-points after merging same-state neighbors;
#'   $index.final: bandwidth of change-points;
#'   $state.new: assigned state per segment;
#'   $iter: number of iterations to convergence.
#' @param eps0,eps1 Neutral-band bounds for st = 3. Neutral is (-eps0, eps1].
#' @param lower,upper Explicit truncation bounds of length `st`, overriding the
#'   defaults for either state count. They act on the QN-transformed scale.
gausianMixure <- function(x, cp, priors, L = 100, st,
                          eps0 = 0.30, eps1 = 0.30, lower = NULL, upper = NULL) {
  # Coerce vector to 1-row matrix (DNA-only mode)
  if (!is.matrix(x)) x <- matrix(x, nrow = 1, dimnames = list("DNA", NULL))
  T_len <- ncol(x)
  cp.new <- cp[which(cp != T_len)]
  index  <- names(cp.new)
  J      <- length(cp.new)
  start  <- c(1, cp.new + 1)
  end    <- c(cp.new, T_len)
  len    <- end - start + 1
  N      <- length(len)
  means    <- fusion_omics(x, cp)$z
  sum.x.sq <- fusion_omics(x, cp)$sum.sq
  state    <- rep(NA, N)

  # Truncation bounds. Caller-supplied bounds win; otherwise st = 3 falls through
  # to updateEM's eps0/eps1 defaults and st = 5 uses the fixed five-state grid.
  if (is.null(lower) || is.null(upper)) {
    if (st == 3) {
      lower <- upper <- NULL              # updateEM builds them from eps0/eps1
    } else if (st == 5) {
      lower <- c(-Inf, -0.5, -eps0, eps1, 0.5)
      upper <- c(-0.5, -eps0, eps1, 0.5, Inf)
    }
  }

  para.new <- para.old <- updateEM(priors$p, priors$mu, priors$sigma,
                                    means, sum.x.sq, N, len, st,
                                    lower = lower, upper = upper,
                                    eps0 = eps0, eps1 = eps1)
  for (iter in 1:L) {
    para.old <- para.new
    para.new <- updateEM(para.new$p.new, para.new$mu.new, para.new$sigma.new,
                          means, sum.x.sq, N, len, st,
                          lower = lower, upper = upper,
                          eps0 = eps0, eps1 = eps1)
    if (abs(sum(para.new$mu.new - para.old$mu.new)) < 1e-6) {
      message("Converged at iteration ", iter)
      break
    }
  }

  for (pt in 1:N) {
    state[pt] <- which.max(para.new$p[pt, ])
  }
  index.change <- which((state[2:(J + 1)] != state[1:J]) == TRUE)
  cp.final     <- c(cp.new[index.change])
  index.final  <- index[index.change]
  state.new    <- state[c(index.change, length(state))]

  list(p.final = para.new$p.new, mu.final = para.new$mu.new,
       sigma.final = para.new$sigma.new,
       cp.final = cp.final, index.final = index.final,
       state.new = state.new, iter = iter)
}

#' Single EM iteration (truncated Gaussian responsibilities + M-step)
#'
#' @param p.in Current mixing proportions (length st).
#' @param mu.in Current means (length st).
#' @param sigma.in Current variances (length st).
#' @param means Fused segment means z_s (length N).
#' @param sum.x.sq Precomputed sum(x^2) per segment (length N).
#' @param N Number of segments.
#' @param len Segment lengths in bins (length N).
#' @param st Number of states.
#' @param lower Truncation lower bounds (length st, or NULL for 3-state default).
#' @param upper Truncation upper bounds (length st, or NULL for 3-state default).
#' @param eps0 Neutral band half-width for 3-state default bounds. Default: 0.30.
#' @param eps1 Upper neutral bound for 3-state default. Default: 0.30.
#' @param eps_var Variance floor. Default: 1e-3.
#' @param eps_pi Mixing proportion floor. Default: 1e-8.
#' @param sigma_floor Minimum variance. Default: 1e-4.
#' @param pi_floor Minimum mixing proportion. Default: 1e-6.
#' @param reorder_by_mu Reorder components by mean. Default: FALSE.
#' @param enforce_mu2_near0 Clamp neutral state mean near 0. Default: TRUE.
#' @param mu2_max_abs Maximum absolute value for neutral mean. Default: 0.15.
#' @return List with $p.new, $mu.new, $sigma.new, $p (responsibilities),
#'   $lower, $upper.
updateEM <- function(p.in, mu.in, sigma.in, means, sum.x.sq, N, len, st,
                     lower = NULL, upper = NULL,
                     eps0 = 0.30, eps1 = 0.30,
                     eps_var = 1e-3, eps_pi = 1e-8,
                     sigma_floor = 1e-4, pi_floor = 1e-6,
                     reorder_by_mu = FALSE, enforce_mu2_near0 = TRUE,
                     mu2_max_abs = 0.15) {
  stopifnot(length(p.in) == st, length(mu.in) == st, length(sigma.in) == st,
            length(means) == N, length(sum.x.sq) == N, length(len) == N)

  # Default truncation bounds for 3-state
  if (is.null(lower) || is.null(upper)) {
    if (st == 3L) {
      lower <- c(-Inf, -eps0, eps1)
      upper <- c(-eps0, eps1, Inf)
    } else {
      stop("Provide 'lower' and 'upper' vectors of length st for st != 3.")
    }
  }
  stopifnot(length(lower) == st, length(upper) == st)

  # Keep parameters valid
  p.in     <- pmax(p.in, eps_pi); p.in <- p.in / sum(p.in)
  sigma.in <- pmax(sigma.in, eps_var)

  # --- E-step: truncated Gaussian responsibilities ---
  log_resp <- matrix(-Inf, N, st)
  for (k in 1:st) {
    sd_k  <- sqrt(sigma.in[k] / len)
    a     <- (upper[k] - mu.in[k]) / sd_k
    b     <- (lower[k] - mu.in[k]) / sd_k
    denom <- pmax(pnorm(a) - pnorm(b), 1e-12)

    ld <- dnorm(means, mean = mu.in[k], sd = sd_k, log = TRUE) - log(denom)
    inside <- (means >= lower[k]) & (means <= upper[k])
    ld[!inside] <- -Inf
    log_resp[, k] <- log(p.in[k]) + ld
  }

  # Row-stable normalization
  p <- matrix(0, N, st)
  for (i in 1:N) {
    m <- max(log_resp[i, ])
    den <- sum(exp(log_resp[i, ] - m))
    if (!is.finite(den) || den <= 0) {
      p[i, ] <- rep(1 / st, st)
    } else {
      p[i, ] <- exp(log_resp[i, ] - m) / den
    }
  }
  p[!is.finite(p)] <- 0
  rs <- rowSums(p); bad <- (!is.finite(rs)) | (rs <= 0)
  if (any(bad)) p[bad, ] <- 1 / st
  p <- pmax(p, 1e-12); p <- p / rowSums(p)

  # --- M-step ---
  N_eff <- sum(len)
  Nk    <- colSums(len * p)

  p.new <- mu.new <- sigma.new <- rep(NA_real_, st)
  for (k in 1:st) {
    den_len_w <- sum(len * p[, k])
    p.new[k] <- (Nk[k] + 1e-3) / (N_eff + 1e-3 * st)
    p.new[k] <- max(p.new[k], pi_floor)

    if (!is.finite(den_len_w) || den_len_w <= 1e-12) {
      mu.new[k]    <- as.numeric(quantile(means, (k - 0.5) / st, na.rm = TRUE))
      sigma.new[k] <- max(stats::var(means, na.rm = TRUE), eps_var)
      next
    }

    mu_k <- sum(means * len * p[, k]) / den_len_w
    neutral_k <- if (st == 5L) 3L else 2L
    if (enforce_mu2_near0 && k == neutral_k) {
      mu_k <- min(max(mu_k, -mu2_max_abs), mu2_max_abs)
    }
    mu.new[k] <- mu_k

    num  <- sum(p[, k] * (sum.x.sq - 2 * len * mu_k * means + len * mu_k^2))
    sig2 <- num / den_len_w
    if (!is.finite(sig2) || sig2 <= eps_var) sig2 <- eps_var
    sigma.new[k] <- max(sig2, sigma_floor)
  }

  p.new <- pmax(p.new, eps_pi); p.new <- p.new / sum(p.new)

  if (reorder_by_mu) {
    ord <- order(mu.new)
    p.new     <- p.new[ord]
    mu.new    <- mu.new[ord]
    sigma.new <- sigma.new[ord]
    p         <- p[, ord, drop = FALSE]
    lower     <- lower[ord]; upper <- upper[ord]
  }

  list(p.new = p.new, mu.new = mu.new, sigma.new = sigma.new,
       p = p, lower = lower, upper = upper)
}

# ======================= cnv_cluster.R =======================

# MintCNA - CNV Clustering Wrapper
# Top-level CNA calling: sets priors, runs EM, extracts states

#' CNV clustering via TGMM EM
#'
#' Computes fused segment means, sets up quantile-based priors, runs the
#' TGMM EM algorithm, and extracts CNA state labels.
#'
#' @param Y Matrix (modalities x bins) or vector of bin intensities.
#' @param cp Integer vector of change-point positions.
#' @param L Maximum EM iterations. Default: 100.
#' @param st Number of states (3 or 5).
#' @param eps0,eps1 Neutral-band bounds used when st = 3. The neutral state is
#'   (-eps0, eps1]; everything below is "del" and above is "dup". These act on
#'   the QN-transformed scale, whose SD is `sd_target`, not on the raw L2R.
#' @param lower,upper Explicit truncation bounds of length `st`, overriding
#'   eps0/eps1. Required for st = 5 if the defaults are not wanted.
#' @return List with:
#'   $newcp: refined change-points after merging same-state neighbors;
#'   $h: bandwidth indices;
#'   $cnv.state: CNA state labels per segment;
#'   $cnv.start, $cnv.end: segment boundaries;
#'   $iter: EM iterations to convergence.
CNVcluster <- function(Y, cp, L, st, eps0 = 0.30, eps1 = 0.30,
                       lower = NULL, upper = NULL) {
  m <- fusion_omics(Y, cp)$z

  if (st == 3) {
    priors <- list(
      p     = c(0.15, 0.7, 0.15),
      mu    = sort(c(quantile(m, .05, na.rm = TRUE),
                     quantile(m, .50, na.rm = TRUE),
                     quantile(m, .95, na.rm = TRUE))),
      sigma = rep(pmax(mad(m, na.rm = TRUE), 0.25), st)
    )
  }

  if (st == 5) {
    priors <- list(
      p     = c(.05, 0.1, 0.7, 0.1, .05),
      mu    = sort(c(quantile(m, .05, na.rm = TRUE),
                     quantile(m, .1, na.rm = TRUE),
                     quantile(m, .50, na.rm = TRUE),
                     quantile(m, .9, na.rm = TRUE),
                     quantile(m, .95, na.rm = TRUE))),
      sigma = rep(pmax(mad(m, na.rm = TRUE), 0.25), st)
    )
  }

  EM.f <- gausianMixure(Y, cp, priors, L, st,
                        eps0 = eps0, eps1 = eps1, lower = lower, upper = upper)
  newcp     <- EM.f$cp.final
  h         <- EM.f$index.final
  cnv.state <- getState(Y, EM = EM.f, st)

  list(newcp = newcp, h = h, cnv.state = cnv.state$cnv.state,
       cnv.start = cnv.state$cnv.start, cnv.end = cnv.state$cnv.end,
       iter = EM.f$iter)
}

#' Extract CNA state labels from EM results
#'
#' Maps integer state assignments to labels (del/nor/dup for 3-state,
#' del/loss/nor/gain/amp for 5-state) and determines segment boundaries.
#'
#' @param Y Matrix or vector (used for ncol to determine total length).
#' @param EM EM result list from gausianMixure().
#' @param st Number of states (3 or 5).
#' @return List with $cnv.state, $cnv.start, $cnv.end.
getState <- function(Y, EM, st) {
  state <- EM$state.new
  cp.f  <- EM$cp.final

  # A cell whose segments are ALL called neutral leaves end.index empty, and the
  # `is.na(cnv.end[length(cnv.end)])` guard below then indexes at position 0 and
  # errors on a zero-length condition. This is inherited from Code/MultiCNA and
  # never fired there because alpha = 0.30 always left something non-neutral; at
  # alpha = 0.05 on real data it does. Return empties and let the caller emit an
  # all-normal table.
  neutral <- if (st == 3) 2L else 3L
  if (!length(which(state != neutral)))
    return(list(cnv.state = character(0), cnv.start = integer(0), cnv.end = integer(0)))

  if (st == 3) {
    # Normal state = 2; non-normal segments are CNAs
    end.index <- which(state != 2)
    cnv.end   <- cp.f[end.index]
    if (is.na(cnv.end[length(cnv.end)])) cnv.end[length(cnv.end)] <- ncol(Y)
    start.index <- end.index - 1
    cnv.start   <- cp.f[start.index]
    if (start.index[1] == 0) cnv.start <- c(1, cnv.start)
    cnv.state <- state[end.index]
    cnv.state[cnv.state == 1] <- "del"
    cnv.state[cnv.state == 2] <- "nor"
    cnv.state[cnv.state == 3] <- "dup"
  }

  if (st == 5) {
    # Normal state = 3; non-normal segments are CNAs
    end.index <- which(state != 3)
    cnv.end   <- cp.f[end.index]
    if (is.na(cnv.end[length(cnv.end)])) cnv.end[length(cnv.end)] <- ncol(Y)
    start.index <- end.index - 1
    cnv.start   <- cp.f[start.index]
    if (start.index[1] == 0) cnv.start <- c(1, cnv.start)
    cnv.state <- state[end.index]
    cnv.state[cnv.state == 1] <- "del"
    cnv.state[cnv.state == 2] <- "loss"
    cnv.state[cnv.state == 3] <- "nor"
    cnv.state[cnv.state == 4] <- "gain"
    cnv.state[cnv.state == 5] <- "amp"
  }

  list(cnv.state = cnv.state, cnv.start = cnv.start, cnv.end = cnv.end)
}

# ======================= fuse_omics.R =======================

# MintCNA - Cross-Omics Fusion
# Coverage-proportional weighting of multi-omic segment means

#' Fuse multi-omic segment scores with coverage-proportional weights
#'
#' Given a multi-omic matrix and change-points, computes per-segment means,
#' coverage-proportional weights, and a weighted fused score z for each segment.
#'
#' @param X Numeric matrix (n_omics x n_bins), or a numeric vector (coerced to 1-row matrix).
#'   Rows = modalities, cols = bins.
#' @param cps Integer vector of bin indices that are change-points.
#' @param alpha Exponent for raw weights (n^alpha). Default: 1.0.
#' @return Data frame with one row per segment, containing:
#'   seg_id, start_bin, end_bin, length_bins,
#'   mean_<omic>, n_<omic>, sumsq_<omic>, w_<omic> for each modality,
#'   z (fused score), sum.sq (weighted sum of squares for EM).
fusion_omics <- function(X, cps, alpha = 1.0) {
  # Coerce vector to 1-row matrix (DNA-only mode)
  if (!is.matrix(X)) {
    X <- matrix(X, nrow = 1, dimnames = list("DNA", NULL))
  }
  stopifnot(ncol(X) >= 1, nrow(X) >= 1)

  if (is.null(rownames(X))) rownames(X) <- paste0("omic", seq_len(nrow(X)))

  n_bins <- ncol(X)
  cps <- sort(unique(as.integer(cps)))
  cps <- cps[cps >= 1 & cps <= n_bins]
  cps <- cps[cps != n_bins]

  starts <- c(1L, cps + 1L)
  ends   <- c(cps, n_bins)
  N <- length(starts)

  out <- data.frame(
    seg_id      = seq_len(N),
    start_bin   = starts,
    end_bin     = ends,
    length_bins = ends - starts + 1L,
    stringsAsFactors = FALSE
  )

  # Per-omic counts, means, sum of squares (NA-safe)
  n_mat    <- matrix(0L,       nrow = N, ncol = nrow(X), dimnames = list(NULL, rownames(X)))
  mean_mat <- matrix(NA_real_, nrow = N, ncol = nrow(X), dimnames = list(NULL, rownames(X)))
  ssq_mat  <- matrix(0.0,     nrow = N, ncol = nrow(X), dimnames = list(NULL, rownames(X)))

  for (o in rownames(X)) {
    nvec <- integer(N); mvec <- numeric(N); ssv <- numeric(N)
    for (i in seq_len(N)) {
      sli <- as.numeric(X[o, starts[i]:ends[i]])
      ok  <- !is.na(sli)
      nvec[i] <- sum(ok)
      mvec[i] <- if (nvec[i] > 0) mean(sli[ok]) else NA_real_
      ssv[i]  <- if (nvec[i] > 0) sum(sli[ok]^2) else 0.0
    }
    out[[paste0("n_", o)]]     <- nvec
    out[[paste0("mean_", o)]]  <- mvec
    out[[paste0("sumsq_", o)]] <- ssv

    n_mat[, o]    <- nvec
    mean_mat[, o] <- mvec
    ssq_mat[, o]  <- ssv
  }

  # Coverage-proportional weights: tilde_w = n^alpha
  tilde <- matrix(0, nrow = N, ncol = nrow(X), dimnames = list(NULL, rownames(X)))
  for (o in rownames(X)) {
    nvec <- out[[paste0("n_", o)]]
    tilde[, o] <- ifelse(nvec > 0, nvec^alpha, 0)
  }

  # Normalize weights per segment
  wsum <- rowSums(tilde)
  w <- tilde / wsum
  w[!is.finite(w)] <- NA_real_

  for (o in rownames(X)) out[[paste0("w_", o)]] <- w[, o]

  # Fused score z_s = sum_o(w_o * mean_o)
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

  # Weighted sufficient stats for EM
  w_filled <- w
  w_filled[!is.finite(w_filled)] <- 0
  sumsq_w <- rowSums((w_filled^2) * ssq_mat)
  out$sum.sq <- sumsq_w

  out
}

# ======================= cna_output.R =======================

# MintCNA - CNA Output Formatting
# Coordinate reindexing and CNA table construction

#' Reindex per-chromosome change-points to global genome coordinates
#'
#' Takes a list of per-chromosome change-point vectors and reindexes them
#' into a single global coordinate system.
#'
#' @param cps_multi List of integer vectors (one per chromosome), each containing
#'   within-chromosome change-point indices.
#' @param chr_lengths Numeric vector of chromosome lengths (in bins), same order as cps_multi.
#' @param include_first_start If TRUE, include position 1 as the first breakpoint. Default: TRUE.
#' @param include_ends If TRUE, include the end of each chromosome. Default: TRUE.
#' @return Sorted integer vector of global change-point positions.
reindex_cps <- function(cps_multi, chr_lengths, include_first_start = TRUE, include_ends = TRUE) {
  stopifnot(length(cps_multi) == length(chr_lengths))
  chr_lengths <- as.numeric(chr_lengths)

  offsets <- c(0, cumsum(chr_lengths))[seq_along(chr_lengths)]

  out <- integer(0)
  if (include_first_start) out <- 1L

  for (i in seq_along(cps_multi)) {
    cp_i <- as.numeric(cps_multi[[i]])
    off  <- offsets[i]
    if (length(cp_i)) out <- c(out, off + cp_i)
    if (include_ends) out <- c(out, off + chr_lengths[i])
  }

  sort(unique(out))
}

#' Build CNA table from state assignments
#'
#' Converts CNA state labels and segment boundaries into a genomic coordinate table.
#' Handles segments that span chromosome boundaries by splitting them.
#'
#' @param cnv.state Character vector of CNA state labels (e.g., "del", "nor", "dup").
#' @param cnv.start Integer vector of segment start indices (global bin coords).
#' @param cnv.end Integer vector of segment end indices.
#' @param ref Data frame with columns chr, start, end (one row per bin).
#' @return Data frame with columns: state, chr, start, end, start_idx, end_idx.
make_cna_table <- function(cnv.state, cnv.start, cnv.end, ref) {
  stopifnot(all(c("chr", "start", "end") %in% names(ref)))
  ref$chr <- as.character(ref$chr)
  n <- nrow(ref)

  chr_ends <- tapply(seq_len(n), ref$chr, max)

  split_one <- function(st, s, e) {
    out <- list(); cur <- s
    while (TRUE) {
      chr_s <- ref$chr[cur]
      chr_end_idx <- chr_ends[[chr_s]]
      end_here <- min(e, chr_end_idx)
      out[[length(out) + 1]] <- data.frame(
        state     = st,
        chr       = chr_s,
        start     = ref$start[cur],
        end       = ref$end[end_here],
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

# ======================= subclone.R =======================

# MintCNA - Subclone Clustering
# K-means clustering with stability selection and silhouette QC

library(mclust)
library(cluster)

#' K-means clustering with stability-selected K
#'
#' Performs K-means on a cell x feature matrix (e.g., latent embeddings from AG-DCAE)
#' with automatic K selection via subsampling stability (ARI) and a 1-SE rule.
#' Silhouette widths are computed as secondary QC.
#'
#' @param Z Numeric matrix (cells x features). Rownames = cell IDs (recommended).
#' @param K_min Minimum K to evaluate. Default: 2.
#' @param K_max Maximum K to evaluate. Default: min(10, floor(n/5)).
#' @param scale_Z Scale columns before clustering. Default: TRUE.
#' @param use_pca Reduce dimensions via PCA first. Default: TRUE.
#' @param pca_var Cumulative variance threshold for PCA. Default: 0.95.
#' @param pca_max Maximum PCA dimensions. Default: 20.
#' @param R Number of subsampling runs per K. Default: 200.
#' @param frac Subsample fraction. Default: 0.8.
#' @param nstart K-means random starts per run. Default: 50.
#' @param iter.max K-means maximum iterations. Default: 1000.
#' @param seed Random seed. Default: 1.
#' @return List with:
#'   $K_hat: selected K;
#'   $clusters: named integer vector (cell -> cluster);
#'   $stability: named numeric vector (K -> mean ARI);
#'   $silhouette: named numeric vector (K -> mean silhouette width);
#'   $km_fit: final kmeans object;
#'   $Z_used: processed feature matrix used for clustering.
cluster_latent_kmeans <- function(Z, K_min = 2, K_max = NULL,
                                  scale_Z = TRUE, use_pca = TRUE,
                                  pca_var = 0.95, pca_max = 20,
                                  R = 200, frac = 0.8,
                                  nstart = 50, iter.max = 1000, seed = 1) {
  set.seed(seed)
  Z <- as.matrix(Z)
  n <- nrow(Z); d <- ncol(Z)
  if (is.null(rownames(Z))) rownames(Z) <- paste0("cell", seq_len(n))

  # K range
  if (is.null(K_max)) {
    K_max <- min(10, floor(n / 5))
    K_max <- max(K_max, K_min)
    K_max <- min(K_max, n - 1)
  }
  Ks <- K_min:K_max
  if (length(Ks) < 1) stop("K range is empty. Need n >= 3.")

  # Preprocess
  Zp <- Z
  if (scale_Z) Zp <- scale(Zp)
  if (use_pca) {
    pca <- prcomp(Zp, center = FALSE, scale. = FALSE)
    var_prop <- (pca$sdev^2) / sum(pca$sdev^2)
    cum <- cumsum(var_prop)
    q <- which(cum >= pca_var)[1]
    q <- min(q, pca_max, n - 1, d)
    Zc <- pca$x[, seq_len(q), drop = FALSE]
  } else {
    Zc <- Zp
  }

  km_subset <- function(idx, K) {
    fit <- kmeans(Zc[idx, , drop = FALSE], centers = K, nstart = nstart, iter.max = iter.max)
    fit$cluster
  }

  # Stability selection
  stability <- setNames(rep(NA_real_, length(Ks)), Ks)
  for (K in Ks) {
    idx_list <- vector("list", R)
    cl_list  <- vector("list", R)
    for (r in seq_len(R)) {
      idx <- sort(sample.int(n, size = max(2, floor(frac * n)), replace = FALSE))
      idx_list[[r]] <- idx
      cl_list[[r]]  <- km_subset(idx, K)
    }
    npairs <- min(200, R * 2)
    aris <- numeric(0)
    for (t in seq_len(npairs)) {
      a <- sample.int(R, 1); b <- sample.int(R, 1)
      ia <- idx_list[[a]]; ib <- idx_list[[b]]
      common <- intersect(ia, ib)
      if (length(common) < max(5, floor(0.2 * length(ia)))) next
      pa <- match(common, ia); pb <- match(common, ib)
      aris <- c(aris, adjustedRandIndex(cl_list[[a]][pa], cl_list[[b]][pb]))
    }
    stability[as.character(K)] <- if (length(aris)) mean(aris) else NA_real_
  }

  # Pick K_hat: max stability with 1-SE rule
  bestK <- as.integer(names(which.max(stability)))
  bestS <- stability[as.character(bestK)]

  se_best <- NA_real_
  {
    idx_list <- vector("list", R); cl_list <- vector("list", R)
    for (r in seq_len(R)) {
      idx <- sort(sample.int(n, size = max(2, floor(frac * n)), replace = FALSE))
      idx_list[[r]] <- idx; cl_list[[r]] <- km_subset(idx, bestK)
    }
    npairs <- min(300, R * 2); aris <- numeric(0)
    for (t in seq_len(npairs)) {
      a <- sample.int(R, 1); b <- sample.int(R, 1)
      ia <- idx_list[[a]]; ib <- idx_list[[b]]
      common <- intersect(ia, ib)
      if (length(common) < max(5, floor(0.2 * length(ia)))) next
      pa <- match(common, ia); pb <- match(common, ib)
      aris <- c(aris, adjustedRandIndex(cl_list[[a]][pa], cl_list[[b]][pb]))
    }
    if (length(aris) > 5) se_best <- sd(aris) / sqrt(length(aris))
  }

  if (!is.na(se_best)) {
    cutoff <- bestS - se_best
    cand <- as.integer(names(stability)[which(stability >= cutoff)])
    K_hat <- min(cand, na.rm = TRUE)
  } else {
    K_hat <- bestK
  }

  # Final clustering
  km_fit <- kmeans(Zc, centers = K_hat, nstart = max(200, nstart), iter.max = iter.max)
  clusters <- km_fit$cluster
  names(clusters) <- rownames(Z)

  # Silhouette QC
  sil <- setNames(rep(NA_real_, length(Ks)), Ks)
  distZ <- dist(Zc)
  for (K in Ks) {
    fitK <- kmeans(Zc, centers = K, nstart = max(50, nstart), iter.max = iter.max)
    sil[K - K_min + 1] <- mean(silhouette(fitK$cluster, distZ)[, "sil_width"])
  }
  names(sil) <- as.character(Ks)

  list(K_hat = K_hat, clusters = clusters, stability = stability,
       silhouette = sil, km_fit = km_fit, Z_used = Zc)
}

# ======================= heatmap.R (build_cna_matrix only) =======================

#' Build CNA state matrix from CNA data frame
#'
#' Converts a CNA data frame (with cell, state, start_idx, end_idx) into a
#' bins x cells integer-coded matrix.
#'
#' @param cna_df_all Data frame with columns: cell, state, start_idx, end_idx.
#' @param ref Bin reference (data frame with chr, start, end).
#' @param code_map Named integer vector mapping state labels to codes.
#'   Default: c(del=1, nor=2, dup=3).
#' @param default_state Default state for bins not covered by any CNA. Default: "nor".
#' @return Integer matrix (bins x cells).
build_cna_matrix <- function(cna_df_all, ref,
                             code_map = c(del = 1, nor = 2, dup = 3),
                             default_state = "nor") {
  n_bins <- nrow(ref)
  cells  <- unique(cna_df_all$cell)
  M <- matrix(code_map[default_state], n_bins, length(cells),
              dimnames = list(sprintf("%s_%d_%d", ref$chr, ref$start, ref$end), cells))

  for (cell in cells) {
    segs <- cna_df_all[cna_df_all$cell == cell, ]
    for (i in seq_len(nrow(segs))) {
      M[segs$start_idx[i]:segs$end_idx[i], cell] <- code_map[segs$state[i]]
    }
  }
  M
}
