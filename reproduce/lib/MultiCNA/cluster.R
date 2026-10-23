# -------------------------
# K-means clustering with stability-selected K (primary)
# + silhouette as secondary QC
# Input: Z (cells x latent_dim) matrix/data.frame with rownames = cell IDs (recommended)
# Output: list(K_hat, clusters, stability, silhouette, km_fit)
# -------------------------
library(mclust)
library(cluster)  
library(uwot) 

cluster_latent_kmeans <- function(
    Z,
    K_min = 2,
    K_max = NULL,
    scale_Z = TRUE,
    use_pca = TRUE,
    pca_var = 0.95,
    pca_max = 20,
    R = 200,                 # number of subsampling runs per K
    frac = 0.8,              # subsample fraction
    nstart = 50,             # kmeans random starts (per run)
    iter.max = 1000,
    seed = 1
) {
  set.seed(seed)
  
  Z <- as.matrix(Z)
  n <- nrow(Z)
  d <- ncol(Z)
  if (is.null(rownames(Z))) rownames(Z) <- paste0("cell", seq_len(n))
  
  # ---- K range ----
  if (is.null(K_max)) {
    K_max <- min(10, floor(n / 5))
    K_max <- max(K_max, K_min)          # handle small n
    K_max <- min(K_max, n - 1)
  }
  Ks <- K_min:K_max
  if (length(Ks) < 1) stop("K range is empty. Need n >= 3.")
  
  # ---- preprocess ----
  Zp <- Z
  if (scale_Z) Zp <- scale(Zp)
  
  if (use_pca) {
    # PCA on scaled Z
    pca <- prcomp(Zp, center = FALSE, scale. = FALSE)
    var <- (pca$sdev^2) / sum(pca$sdev^2)
    cum <- cumsum(var)
    q <- which(cum >= pca_var)[1]
    q <- min(q, pca_max, n - 1, d)
    Zc <- pca$x[, seq_len(q), drop = FALSE]
  } else {
    Zc <- Zp
  }
  
  # ---- helper: kmeans on subset ----
  km_subset <- function(idx, K) {
    fit <- kmeans(Zc[idx, , drop = FALSE], centers = K, nstart = nstart, iter.max = iter.max)
    fit$cluster
  }
  
  # ---- stability selection ----
  stability <- setNames(rep(NA_real_, length(Ks)), Ks)
  
  for (K in Ks) {
    # store partitions + indices
    idx_list <- vector("list", R)
    cl_list  <- vector("list", R)
    
    for (r in seq_len(R)) {
      idx <- sort(sample.int(n, size = max(2, floor(frac * n)), replace = FALSE))
      idx_list[[r]] <- idx
      cl_list[[r]]  <- km_subset(idx, K)
    }
    
    # compute mean ARI over random pairs (faster than all pairs)
    # You can increase npairs for more precision.
    npairs <- min(200, R * 2)
    aris <- numeric(0)
    
    for (t in seq_len(npairs)) {
      a <- sample.int(R, 1)
      b <- sample.int(R, 1)
      ia <- idx_list[[a]]
      ib <- idx_list[[b]]
      common <- intersect(ia, ib)
      if (length(common) < max(5, floor(0.2 * length(ia)))) next
      
      # map common cells to positions in each subset
      pa <- match(common, ia)
      pb <- match(common, ib)
      aris <- c(aris, adjustedRandIndex(cl_list[[a]][pa], cl_list[[b]][pb]))
    }
    
    stability[as.character(K)] <- if (length(aris)) mean(aris) else NA_real_
  }
  
  # pick K_hat by max stability; tie-break toward smaller K within 1-SE of best
  bestK <- as.integer(names(which.max(stability)))
  bestS <- stability[as.character(bestK)]
  
  # estimate SE at bestK using a fresh set of ARIs (quick bootstrap-ish)
  # (optional but helps "1-SE" rule)
  se_best <- NA_real_
  {
    # recompute aris for bestK to estimate se
    idx_list <- vector("list", R)
    cl_list  <- vector("list", R)
    for (r in seq_len(R)) {
      idx <- sort(sample.int(n, size = max(2, floor(frac * n)), replace = FALSE))
      idx_list[[r]] <- idx
      cl_list[[r]]  <- km_subset(idx, bestK)
    }
    npairs <- min(300, R * 2)
    aris <- numeric(0)
    for (t in seq_len(npairs)) {
      a <- sample.int(R, 1)
      b <- sample.int(R, 1)
      ia <- idx_list[[a]]
      ib <- idx_list[[b]]
      common <- intersect(ia, ib)
      if (length(common) < max(5, floor(0.2 * length(ia)))) next
      pa <- match(common, ia)
      pb <- match(common, ib)
      aris <- c(aris, adjustedRandIndex(cl_list[[a]][pa], cl_list[[b]][pb]))
    }
    if (length(aris) > 5) se_best <- sd(aris) / sqrt(length(aris))
  }
  
  if (!is.na(se_best)) {
    # smallest K within 1-SE of best stability
    cutoff <- bestS - se_best
    cand <- as.integer(names(stability)[which(stability >= cutoff)])
    K_hat <- min(cand, na.rm = TRUE)
  } else {
    K_hat <- bestK
  }
  
  # ---- final kmeans on all cells ----
  km_fit <- kmeans(Zc, centers = K_hat, nstart = max(200, nstart), iter.max = iter.max)
  clusters <- km_fit$cluster
  names(clusters) <- rownames(Z)
  
  # ---- silhouette QC (secondary) ----
  sil <- setNames(rep(NA_real_, length(Ks)), Ks)
  distZ <- dist(Zc)  # OK for n up to a few thousand; if huge n, skip or approximate
  for (K in Ks) {
    fitK <- kmeans(Zc, centers = K, nstart = max(50, nstart), iter.max = iter.max)
    sil[K - K_min + 1] <- mean(silhouette(fitK$cluster, distZ)[, "sil_width"])
  }
  names(sil) <- as.character(Ks)
  
  list(
    K_hat = K_hat,
    clusters = clusters,
    stability = stability,
    silhouette = sil,
    km_fit = km_fit,
    Z_used = Zc
  )
}