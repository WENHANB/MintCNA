#' CNVcluster
#'
#' This function uses different priors settings to achieve convergence of the Expectation-Maximization algorithm, and then determine the best clustering results by applying the modified BIC.
#' @param Y the numeric vector of the intensities of markers
#' @param cp the numeric vector of the position index for the identified change-points
#' @param L repeat times in the EM algorithm, defaults to 100
#' @param method the fusion method for the omics, "count" is based on the number of observations for each omic or "v" is based on the weight learned from rank-1 approximation for each omic
#' @return The return is the clustered CNV segments by presenting the start position and end position, length of the CNV and the copy number states (duplication or deletion). It also returns a vector of final candidates of change-points.
#' @return  newcp the final list of change-points
#' @return  h the bandwidth used for the identification of change-points
#' @return  cnv.state copy number state for each CNV
#' @return  cnv.start start position of each CNV
#' @return  cnv.end end position of each CNV
#' @seealso \link{gausianMixure} for clustering of CNVs using Expectation-Maximization algorithm.
#' @export
CNVcluster <-function(Y, cp, L, method, st) {
  ##Gaucian mixture model based clustering
  m = fusion_omics(Y, cp, method)$z
  if (st==3) {
  # For three states including dup, del and normal
  # priors = get.priors(mean.vector, G = st)
  # priors <- list(
  #   p     = rep(1/st, st),
  #   mu    = c(-1,0.001,0.8),
  #   sigma = rep(0.5, st)
  # )
  priors <- list(
    p     = c(0.15, 0.7, 0.15),
    mu    = sort(c(quantile(m, .05, na.rm=TRUE),
                   quantile(m, .50, na.rm=TRUE),
                   quantile(m, .95, na.rm=TRUE))),
    sigma = rep(pmax(mad(m, na.rm=TRUE), 0.25), st)
  )}
  
  if (st==5) {
    # For three states including dup, del and normal
    # priors = get.priors(mean.vector, G = st)
    # priors <- list(
    #   p     = rep(1/st, st),
    #   mu    = c(-1,0.001,0.8),
    #   sigma = rep(0.5, st)
    # )
    priors <- list(
      p     = c(.05, 0.1, 0.7, 0.1, .05),
      mu    = sort(c(quantile(m, .05, na.rm=TRUE),
                     quantile(m, .1, na.rm=TRUE),
                     quantile(m, .50, na.rm=TRUE),
                     quantile(m, .9, na.rm=TRUE),
                     quantile(m, .95, na.rm=TRUE))),
      sigma = rep(pmax(mad(m, na.rm=TRUE), 0.25), st)
    )}
  # During EM after each M-step, reorder components by mu so del < neu < dup
  # and optionally clamp neutral mean near 0: mu[2] <- shrink(mu[2], target=0, kappa)
  EM.f = gausianMixure(Y, cp, priors, L, st)
  newcp = EM.f$cp.final
  h     = EM.f$index.final
  cnv.state <- getState(Y,EM = EM.f, st)
  return (list(newcp = newcp, h = h, cnv.state = cnv.state$cnv.state, 
               cnv.start = cnv.state$cnv.start, cnv.end = cnv.state$cnv.end, iter = EM.f$iter))
}

getState <-function (Y,EM = EM, st) {
    # state      = EM$state.new
    # cp.f       = EM$cp.final
    # end.index  = which(state!=2)
    # cnv.end     = cp.f[end.index]
    # if(end.index[1] == 1 | end.index[1] == 3) {
    #   end.index = end.index[-1]
    #   cnv.end   = cnv.end[-1]
    # } else if(end.index[length(end.index)] == 1 | end.index[length(end.index)] == 3) {
    #   end.index = end.index[-(length(end.index))]
    #   cnv.end   = cnv.end[-length(end.index)]
    # }
    # cnv.start = cp.f[end.index-1]
    # cnv.state = state[end.index]
  if (st==3) {
    state      = EM$state.new
    cp.f       = EM$cp.final
    end.index  = which(state!=2)
    cnv.end = cp.f[end.index]
    #cnv.end = cnv.end[!is.na(cnv.end)]
    if(is.na(cnv.end[length(cnv.end)])){cnv.end[length(cnv.end)] = ncol(Y)}
    start.index = end.index-1
    cnv.start     = cp.f[start.index]
    #cnv.start = cnv.start[!is.na(cnv.start)]
    if(start.index[1] == 0){cnv.start = c(1, cnv.start)}
    cnv.state = state[end.index]
    cnv.state[which(cnv.state == 1)] = "del"
    cnv.state[which(cnv.state == 2)] = "nor"
    cnv.state[which(cnv.state == 3)] = "dup"
  } 
  if (st==5) {
    state      = EM$state.new
    cp.f       = EM$cp.final
    end.index  = which(state!=3)
    cnv.end = cp.f[end.index]
    #cnv.end = cnv.end[!is.na(cnv.end)]
    if(is.na(cnv.end[length(cnv.end)])){cnv.end[length(cnv.end)] = ncol(Y)}
    start.index = end.index-1
    cnv.start     = cp.f[start.index]
    #cnv.start = cnv.start[!is.na(cnv.start)]
    if(start.index[1] == 0){cnv.start = c(1, cnv.start)}
    
    cnv.state = state[end.index]
    cnv.state[cnv.state == 1] <- "del"
    cnv.state[cnv.state == 2] <- "loss"
    cnv.state[cnv.state == 3] <- "nor"
    cnv.state[cnv.state == 4] <- "gain"
    cnv.state[cnv.state == 5] <- "amp"
  }

  return (list(cnv.state = cnv.state, cnv.start = cnv.start, cnv.end = cnv.end))
}

