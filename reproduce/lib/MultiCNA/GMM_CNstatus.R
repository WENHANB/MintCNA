#' Clustering of CNVs using Expectation-Maximization algorithm
#'
#' This function clusters the identified change-points to make final CNV calling.The potential CNV segments between two neighbor candidate change-points are assigned to different copy number states according to the average intensities in the segment intervals. We use three clusters including duplication, normal state and deletion. Each cluster is presented by Gaussian distribution with unknown mean and variance. Expectation-Maximization (EM) algorithm is applied for the mixture of Gaussians to assign each segment to the most probable cluster/state. Two physically linked candidate CNV segments in the same group are merged to one unique CNV segment.
#' @param x the vector of the intensities of markers
#' @param cp the vector of the marker index of the identified change-points
#' @param priors given initial parameters for the EM algorithm
#' @param L repeat times in the EM algorithm. Defaults to 100
#' @param st number of assumed states in the EM algorithm
#' @return The return is the clustered CNV segments by presenting the start position and end position using SNP or CNV marker index, and the copy number states. It also returns a vector of final candidates of change-points.
#' @return p.final probability of falling into each state for each CNV segment after convergence
#' @return mu.final segment means of each state after convergence
#' @return cp.final list of change-points after EM algorithm
#' @return index.final the bandwidth of change-points
#' @return state.new assigned copy number state for each CNV
#' @export
gausianMixure <- function(x, cp, priors, L=100,st) {
  T      <- if (is.matrix(x)) ncol(x) else length(x)
  cp.new <- cp[which(cp!=T)]
  index  <- names(cp.new)
  J      <- length(cp.new)
  start  <- c(1, cp.new+1)
  end    <- c(cp.new, T)
  len    <- end - start + 1
  N      <- length(len)
  means  <- fusion_omics(x, cp, method = "count")$z
  sum.x.sq <- fusion_omics(x, cp, method = "count")$sum.sq
  state    <- rep(NA, N)
  
  # --- Setup TGMM parameters ---
  if (st == 3) {
    # nothing special needed — updateEM() handles truncation bounds internally
    lower <- upper <- NULL
  }
  
  if (st == 5) {
    # c0   <- mean(means, na.rm = TRUE)
    # beta <- max(0.25, 0.7 * mad(means, na.rm = TRUE))   # floor avoids ~0.01 drift
    # lower <- c(-Inf, c0 - 2*beta, c0 - beta, c0 + beta, c0 + 2*beta)
    # upper <- c(c0 - 2*beta, c0 - beta, c0 + beta, c0 + 2*beta, Inf)
    lower <- c(-Inf, -0.5, -0.30, 0.30, 0.5)
    upper <- c(-0.5, -0.30,0.30, 0.5, Inf)
  }
  
  para.new <- para.old <- updateEM(priors$p, priors$mu, priors$sigma, means, sum.x.sq, N, len,st,lower = lower, upper = upper)
  for (iter in 1:L){
    para.old <- para.new
    para.new <- updateEM(para.new$p.new, para.new$mu.new, para.new$sigma.new, means, sum.x.sq, N, len,st,lower = lower, upper = upper)
      if(abs(sum(para.new$mu.new - para.old$mu.new)) < 1e-6){
        message("Converged at iteration ", iter)
        break
      }
     #print(iter)
     #print(tail(para.new$p,10))
  }
  for (pt in 1:N) {
    state[pt] <- which.max(para.new$p[pt,])
  }
  index.change <- which((state[2:(J+1)] != state[1:J])== "TRUE")
  cp.final     <- c(cp.new[index.change])
  index.final  <- index[index.change]
  state.new    <- state[c(index.change,length(state))]
  return(list(p.final = para.new$p.new, mu.final = para.new$mu.new, sigma.final = para.new$sigma.new, 
              cp.final = cp.final, index.final = index.final, state.new = state.new, iter = iter))
}

########################################## TGMM Version 2############################################################################
# TGMM EM update (truncated responsibilities + "FastCall"-style M-step)
# - p.in, mu.in, sigma.in: current params (length st)
# - means: z_s per segment (length N)
# - sum.x.sq: precomputed sum over bins for each segment (∑ x^2), length N
# - N: number of segments; len: segment lengths (bins); st: #states
# - lower, upper: truncation bounds (length st)
# - eps0: small band around 0 for neutral when you don’t pass bounds (st=3 case)
updateEM <- function(
    p.in, mu.in, sigma.in, means, sum.x.sq, N, len, st,
    lower = NULL, upper = NULL,
    eps0 = 0.30,
    eps1 = 0.30,# used only if lower/upper are NULL and st==3
    eps_var = 1e-3, eps_pi = 1e-8,
    sigma_floor = 1e-4, pi_floor = 1e-6,
    reorder_by_mu = F, enforce_mu2_near0 = T, mu2_max_abs = 0.15
) {
  stopifnot(length(p.in)==st, length(mu.in)==st, length(sigma.in)==st,
            length(means)==N, length(sum.x.sq)==N, length(len)==N)
  
  # --- Truncation bounds (default 3-state CN: del, neu, dup) ---
  if (is.null(lower) || is.null(upper)) {
    if (st == 3L) {
      lower <- c(-Inf, -eps0,  eps1)
      upper <- c( -eps0, eps1,  Inf)
    } else {
      stop("Please provide 'lower' and 'upper' vectors of length st for st != 3.")
    }
  }
  stopifnot(length(lower)==st, length(upper)==st)

  # Keep parameters valid
  p.in     <- pmax(p.in, eps_pi); p.in <- p.in / sum(p.in)
  sigma.in <- pmax(sigma.in, eps_var)

  # ---------- E-step: truncated Gaussian responsibilities ----------
  # g_k(z) = [phi((z-mu_k)/sd_k) / sd_k] / [Phi((u_k-mu_k)/sd_k) - Phi((l_k-mu_k)/sd_k)]
  # If z is outside [l_k,u_k], density is 0.
  log_resp <- matrix(-Inf, N, st)
  for (k in 1:st) {
    sd_k <- sqrt(sigma.in[k] / len)                # segment-specific sd
    a    <- (upper[k] - mu.in[k]) / sd_k
    b    <- (lower[k] - mu.in[k]) / sd_k
    denom <- pmax(pnorm(a) - pnorm(b), 1e-12)      # truncation normalizer

    # base log-density for untruncated normal
    ld <- dnorm(means, mean = mu.in[k], sd = sd_k, log = TRUE) - log(denom)

    # zero outside the truncation interval
    inside <- (means >= lower[k]) & (means <= upper[k])
    ld[!inside] <- -Inf

    # add mixing weights in log-space
    log_resp[,k] <- log(p.in[k]) + ld
  }

  # Row-stable normalize to get responsibilities p[i,k]
  p <- matrix(0, N, st)
  for (i in 1:N) {
    m <- max(log_resp[i, ])
    den <- sum(exp(log_resp[i, ] - m))
    if (!is.finite(den) || den <= 0) {
      p[i, ] <- rep(1/st, st)
    } else {
      p[i, ] <- exp(log_resp[i, ] - m) / den
    }
  }
  # sanitize
  p[!is.finite(p)] <- 0
  rs <- rowSums(p); bad <- (!is.finite(rs)) | (rs <= 0)
  if (any(bad)) p[bad, ] <- 1/st
  p <- pmax(p, 1e-12); p <- p / rowSums(p)

  # ---------- M-step (approximate, treat as if untruncated) ----------
  N_eff <- sum(len)
  Nk    <- colSums(len * p)

  p.new     <- mu.new <- sigma.new <- rep(NA_real_, st)
  for (k in 1:st) {
    den_len_w <- sum(len * p[,k])

    # mixing proportions with tiny Dirichlet smoothing + floor
    p.new[k] <- (Nk[k] + 1e-3) / (N_eff + 1e-3 * st)
    p.new[k] <- max(p.new[k], pi_floor)

    if (!is.finite(den_len_w) || den_len_w <= 1e-12) {
      # re-init empty component near quantile
      mu.new[k]    <- as.numeric(quantile(means, (k - 0.5)/st, na.rm = TRUE))
      sigma.new[k] <- max(stats::var(means, na.rm = TRUE), eps_var)
      next
    }

    # weighted mean
    mu_k <- sum(means * len * p[,k]) / den_len_w

    # optional neutral constraint
    if (enforce_mu2_near0 && st >= 2 && k == 2) {
      mu_k <- min(max(mu_k, -mu2_max_abs), mu2_max_abs)
    }
    mu.new[k] <- mu_k

    # weighted variance in the "per-bin" parameterization
    num  <- sum(p[,k] * (sum.x.sq - 2 * len * mu_k * means + len * mu_k^2))
    den  <- den_len_w
    sig2 <- num / den
    if (!is.finite(sig2) || sig2 <= eps_var) sig2 <- eps_var
    sigma.new[k] <- max(sig2, sigma_floor)
  }

  # renormalize mixture weights
  p.new <- pmax(p.new, eps_pi); p.new <- p.new / sum(p.new)

  # optional ordering by means for interpretability
  if (reorder_by_mu) {
    ord <- order(mu.new)
    p.new     <- p.new[ord]
    mu.new    <- mu.new[ord]
    sigma.new <- sigma.new[ord]
    p         <- p[, ord, drop = FALSE]
    lower     <- lower[ord]; upper <- upper[ord]
  }

  list(
    p.new = p.new,mu.new = mu.new,sigma.new = sigma.new,p = p,lower = lower,upper = upper
  )
}

