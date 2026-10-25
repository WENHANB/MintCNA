# SCCNV (Wang et al.), re-implemented in R following the published algorithm.
# dat1: one row per bin; col 1 chr, col 2 mappability, col 4 normalisation divisor, col 6 GC (%),
# col 7 N fraction, cols 8+ read counts per cell. Returns, per cell, the posterior probability of
# copy number 0..20 for every bin (a data frame with columns CN0..CN20).
SCCNV <- function(dat1, minmap = 0.3, maxN = 0.1, ploidy = 2, windowsize = 10, binsize = 1000) {
  genomesize <- nrow(dat1) * binsize

  ## mappability filter and normalisation
  dat2 <- dat1[which(dat1[, 2] >= minmap & dat1[, 7] <= maxN), ]
  dat2[, 8:ncol(dat2)] <- dat2[, 8:ncol(dat2)] / dat2[, 4]

  ## GC normalisation by GC percentile
  print("Normalize by GC content")
  gc_percentile <- quantile(dat2[, 6], seq(0, 100) / 100)
  gc_pp <- rep(-1, nrow(dat2))
  for (i in 1:nrow(dat2)) {
    for (j in 1:(length(gc_percentile) - 1)) {
      if (dat2[i, 6] >= gc_percentile[j] & dat2[i, 6] <= gc_percentile[j + 1]) gc_pp[i] <- (j - 1)
    }
  }
  gcNorm <- function(RC, gc_pp = gc_pp) {
    tmp.ave.seq <- c()
    for (ii in 0:100) {
      tmp.ind <- which(gc_pp == ii)
      if (length(tmp.ind) == 0) tmp.ave.seq <- c(tmp.ave.seq, 1)
      if (length(tmp.ind) > 0)  tmp.ave.seq <- c(tmp.ave.seq, mean(RC[tmp.ind]))
    }
    # percentiles without bins take the mean of the nearest occupied ones
    com.seq <- which(tmp.ave.seq == 1)
    candi.seq <- which(tmp.ave.seq != 1)
    for (kk in com.seq) {
      tmp.ind2 <- which.min(abs(kk - candi.seq))
      tmp.ave.seq[kk] <- mean(tmp.ave.seq[c(tmp.ind2, tmp.ind2 + 1)])
    }
    overall.mean <- mean(RC)
    tmp.nr_n <- c()
    for (jj in 1:length(RC)) tmp.nr_n <- c(tmp.nr_n, overall.mean * RC[jj] / tmp.ave.seq[gc_pp[jj] + 1])
    tmp.nr_n
  }
  dat3 <- dat2
  for (i in 8:ncol(dat2)) dat3[, i] <- gcNorm(RC = dat2[, i], gc_pp = gc_pp)

  ## raw copy-number estimate and cross-sample normalisation
  print("Cross batch normalize")
  tmp.dat4 <- dat3[, 8:ncol(dat3)]
  cell_nr <- apply(tmp.dat4, 2, mean)
  dat4 <- dat3
  for (i in 1:ncol(tmp.dat4)) dat4[, i + 7] <- ploidy * tmp.dat4[, i] / cell_nr[i]

  dat5 <- dat4
  temp.dat5 <- dat4[, -c(1:7)]
  for (i in 1:nrow(temp.dat5)) {
    for (j in 1:ncol(temp.dat5)) {
      temp2 <- (sum(temp.dat5[i, ]) - temp.dat5[i, j]) / (ncol(temp.dat5) - 1)
      if (temp2 == 0) dat5[i, j + 7] <- 0 else dat5[i, j + 7] <- (temp.dat5[i, j] / temp2) * ploidy
    }
  }
  ## per cell, scale so that the median copy number is the ploidy
  nfactor <- apply(dat5[, -c(1:7)], 2, median)
  dat51 <- dat5
  for (i in 8:ncol(dat5)) dat51[, i] <- (ploidy / nfactor[i - 7]) * dat5[, i]
  dat5 <- dat51

  ## running-window smoothing within chromosome
  print("Futher smoothing")
  step <- windowsize / 2
  dat6 <- dat5
  for (i in 1:nrow(dat5)) {
    for (j in 8:ncol(dat5)) {
      p1 <- i - step; p2 <- i + step + 1
      if (p1 <= 0) p1 <- 1
      if (p2 > nrow(dat5)) p2 <- nrow(dat5)
      temp <- c()
      for (k in p1:p2) if (dat5[k, 1] == dat5[i, 1]) temp <- c(temp, dat5[k, j])
      dat6[i, j] <- mean(temp)
    }
  }

  ## copy-number calling: Gaussian posterior over CN 0..20 with a per-cell sigma
  print("CNV calling")
  nr <- dat6[, -c(1:7)]
  sigma <- c()
  for (i in 1:ncol(nr)) {
    sigma <- c(sigma, abs(quantile(nr[, i], 0.309) - ploidy) + abs(quantile(nr[, i], 0.691) - ploidy))
  }
  norm_pdf <- function(x, mm, sg) dnorm(x = (x - mm) / sg, mean = 0, sd = 1)

  post.prob.list <- list()
  for (j in 8:ncol(dat6)) {
    tmp.prob.data <- data.frame()
    for (i in 1:nrow(dat6)) {
      pdf <- norm_pdf(x = dat6[i, j], mm = seq(0, 20), sg = sigma[j - 7])
      tmp.prob.data <- rbind(tmp.prob.data, pdf / sum(pdf))
    }
    colnames(tmp.prob.data) <- paste0("CN", seq(0, 20))
    post.prob.list[[j - 7]] <- tmp.prob.data
  }
  post.prob.list
}
