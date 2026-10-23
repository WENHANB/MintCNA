#' Change-points score computation
#'
#' Computes diagnostic scores by combining MissCUSUM with a sparse rank-1
#' projection. For each time point, a local CUSUM is calculated, projected to 1D
#' using either an l1 or l2 penalized vector, and the deviation at the center is
#' recorded. Returns the full score series, local maxima, and projection vectors.
#'
#' @param x A p × n numeric matrix (modalities × bins), with possible NAs.
#' @param lambda Regularization parameter for the projection step. Defaults to sqrt(n * log(p * log n)) / 2.
#' @param h Half-window size for local CUSUM computation and local maxima search. Default: 10
#' @param standardize.series Logical; if TRUE, rescale variance of each row. Default: FALSE.
#' @param penalty Projection penalty type, either "l1" or "l2". Default: "l2".
#'
#' @return A list with three elements:
#'   D: data frame of all scores (position, score);
#'   local_maxima: subset of D that are local maxima;
#'   v_df: data frame of projection vectors for each bin.
#'
#' @examples
#' set.seed(1)
#' x <- rbind(rnorm(200), rnorm(200))
#' res <- multiMissInspect_scores(x, h = 10, penalty = "l2")
#' head(res$D)
#'
#' @export
#' 
library(softImpute)
library(InspectChangepoint)

multiMissInspect_scores <- function(x, lambda = NULL, h = 10, standardize.series = FALSE, penalty = "l2") {
  p <- nrow(x)
  n <- ncol(x)
  if (is.null(lambda)) lambda <- sqrt(n * log(p * log(n))) / 2
  #if (standardize.series) x <- rescale.variance(x)
  
  D <- data.frame(position = integer(0), score = numeric(0))
  V_list <- list() 
  # Step: Compute local scores s_t
  for (t in seq(h + 1, n - h)) {
    x_window <- x[, (t - h):(t + h)]
    x.cusum <- cusum(x_window)
    
    # Choose projection method
    vhat <- switch(penalty,
                   l1 = get.v.proj_l1(x.cusum, lambda),
                   l2 = get.v.proj_l2(x.cusum, lambda),
                   stop("penalty must be 'l1' or 'l2'")
    )
    
    proj_scores <- as.vector(vhat %*% x.cusum)
    score <- proj_scores[h+1]
    D <- rbind(D, data.frame(position = t, score = score))
    
    V_list[[t]] <- data.frame(bin = t, t(vhat))
    #V_list[[t]] <- data.frame(bin = t, v1 = vhat[1], v2 = vhat[2])
  }
  
  v_df <- do.call(rbind, V_list)
  # Step: Identify local maxima of |s_t| within ±h
  is_local_max <- function(vec, idx, h) {
    left <- max(1, idx - h)
    right <- min(length(vec), idx + h)
    return(all(abs(vec[idx]) >= abs(vec[left:right])))
  }
  
  local_maxima <- D[sapply(1:nrow(D), function(i) is_local_max(D$score, i, h)), ]
  
  return(list(D = D, local_maxima = local_maxima, v_df = v_df))
}

## Computes a CUSUM-like transformation for a single time series, accounting for missing data 
cusum.vector <- function(x){
  ob <- !is.na(x)
  z <- replace(x, !ob, 0)
  n <- length(x)
  leftsum <- cumsum(z)
  rightsum <- leftsum[n] - leftsum
  L <- cumsum(ob)
  R <- L[n] - L
  x.cusum <- (leftsum[-n] / L[-n] - rightsum[-n] / R[-n]) * sqrt(L[-n] * R[-n] / L[n])
  x.cusum[is.nan(x.cusum)] <- 0
  x.cusum
}

## Applies cusum.vector() to vector or matrix
cusum <- function(x){
  if (!is.matrix(x)) {
    return(cusum.vector(x))
  } else if (ncol(x)==2){
    return(as.matrix(apply(x, 1, cusum.vector)))
  } else {
    return(t(apply(x, 1, cusum.vector)))
  }
}

## Applies rank-1 approximation with alternating updates using ℓ1-penalty (Lasso-style projection)
get.v.proj_l1 <- function(x.cusum, lambda, max_iter=1000, tol=1e-10){
  vhat <- svd(x.cusum)$u[, 1]
  if (sum(abs(x.cusum)) == 0) return(random.UnitVector(nrow(x.cusum)))
  for (iter in 1:max_iter){
    vhat_old <- vhat
    what <- vector.normalise(t(x.cusum) %*% vhat)
    tmp <- x.cusum %*% what
    lambda_tmp <- min(lambda, max(abs(tmp)) - 1e-10)
    vhat <- vector.normalise(vector.soft.thresh(x.cusum %*% what, lambda_tmp))
    if (vector.norm(vhat_old - vhat) < tol) break
  }
  return(as.vector(vhat))
}

## Applies rank-1 approximation with alternating updates using ℓ2-penalty (Ridge-style projection)
get.v.proj_l2 <- function(x.cusum, lambda, max_iter = 1000, tol = 1e-10) {
  if (sum(abs(x.cusum)) == 0) return(random.UnitVector(nrow(x.cusum)))
  vhat <- svd(x.cusum)$u[, 1]
  for (iter in 1:max_iter) {
    vhat_old <- vhat
    
    # Step 1: update w
    what <- vector.normalise(t(x.cusum) %*% vhat)
    
    # Step 2: update v with ridge-style shrinkage
    vhat <-(1 / (2 * lambda)) * (x.cusum %*% what)
    
    # Check convergence
    if (vector.norm(vhat_old - vhat) < tol) break
  }
  return(as.vector(vhat))
}


random.UnitVector <- function(n){
  v = rnorm(n)
  v/vector.norm(v)
}


## High-level function that Computes MissCUSUM, Computes vhat using the sparse rank-1 projection, 
### Projects to 1D and finds the time point with max deviation, Returns the single estimated change point. 
# !!Ignore this function
# locate.change <- function(x, lambda, standardize.series=FALSE, view.cusum=FALSE){
#   if (!is.matrix(x)) x <- t(x)
#   p <- nrow(x)
#   n <- ncol(x)
#   if (missing(lambda)) lambda <- sqrt(n*log(p*log(n)))/2
#   if (standardize.series) x <- rescale.variance(x)
#   
#   x.cusum <- cusum(x)
#   lambda <- min(lambda, max(apply(x.cusum, 1, vector.norm)) - 1e-10)
#   
#   vhat <- get.v.proj(x.cusum, lambda)
#   
#   x.cusum.proj <- colSums(as.vector(vhat) * x.cusum)
#   cp <- which.max(abs(x.cusum.proj) * (1 + rnorm(length(x.cusum.proj)) * 1e-10)) # random tie breaking
#   if (view.cusum) {
#     plot(x.cusum.proj, ylab = "projected cusum", pch = 20)
#   }
#   attr(cp, 'cusum') <- x.cusum.proj
#   attr(cp, 'vhat') <- vhat
#   return(cp)
# }
