#' calculate the prior information for the CNV states
#'
#' @param mean.vector the segmental means from all cells/samples
#' @param G number of CNV states, default 3 (del/nor/dup)
#' @return priors the proportion, mean, and standard deviation of the means for each state
#'
#'@export
#'
library(mclust)
get.priors<-function (mean.vector, G=3){

  yMclust <- Mclust(mean.vector,G=G,verbose = FALSE)
  init.clust<-yMclust$classification
  mus<-c()
  sds<-c()
  ps <-c()
  for (s1 in 1:G){
    mus<-c(mus,mean(mean.vector[init.clust==s1]))
    sds<-c(sds,sd(mean.vector[init.clust==s1]))
    ps<-c(ps,sum(init.clust==s1)/length(init.clust))
  }
  reorder<-order(mus)
  mu=mus[reorder]
  sd=sds[reorder]
  p=ps[reorder]
  if(anyNA(sd)){sd = rep(0.5, G)}
  if(anyNA(p)) {p = rep(1/G,G)}
  if(anyNA(mu)) {mu = c(-1, 0.001, 1)}
  priors <- list(p = p, mu = mu, sigma = sd)
  return(priors=priors)
  # return(list(priors=priors,init.clust=init.clust))
}
