# QN to N(0,1) per modality; optional per-chromosome
# X: r x n (rows = modalities). chr_ranges: list of integer vectors (columns per chr).
qn_modal <- function(X, chr_ranges = NULL, mean_target = 0, sd_target = 1) {
  X <- as.matrix(X)
  out <- X
  R <- nrow(X)
  if (length(sd_target) == 1) sd_target <- rep(sd_target, R)
  if (length(mean_target) == 1) mean_target <- rep(mean_target, R)
  
  if (is.null(chr_ranges)) {
    for (r in 1:R) {
      ok <- is.finite(X[r, ])
      m  <- sum(ok)
      if (m) {
        rk <- rank(X[r, ok], ties.method = "average")
        u  <- (rk - 0.5) / m
        z  <- qnorm(u)
        out[r, ok] <- mean_target[r] + sd_target[r] * z
      }
    }
  } else {
    for (r in 1:R) {
      for (rg in chr_ranges) {
        ok <- is.finite(X[r, rg])
        m  <- sum(ok)
        if (m) {
          rk <- rank(X[r, rg][ok], ties.method = "average")
          u  <- (rk - 0.5) / m
          z  <- qnorm(u)
          out[r, rg][ok] <- mean_target[r] + sd_target[r] * z
        }
      }
    }
  }
  out
}
