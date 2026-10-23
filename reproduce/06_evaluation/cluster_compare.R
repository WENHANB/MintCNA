args <- commandArgs(trailingOnly = TRUE)
BASE <- if (length(args) >= 1) args[1] else "c3_N100_100kb_l1e4"
print(BASE)
CNA_DIR <- file.path("<MINTCNA_ROOT>/Simulation/cnasim", BASE)  

library(data.table)
library(mclust)
library(cluster)  
library(uwot) 
library(ggplot2)

source("<REPRODUCE_DIR>/lib/MultiCNA/cluster.R")

l2r <- read.csv(file.path(CNA_DIR,"l2r_smoothed_all_cells_100kb.csv"), row.names = 1)

##################################################

library(data.table)
library(mclust)
library(Rtsne)

set.seed(1)

# must already exist in your environment:
# cluster_latent_kmeans()

# ------------------------------------------------------------
# 1) read truth
# ------------------------------------------------------------
truth <- fread(file.path(CNA_DIR, "cell_types.tsv"))
setnames(truth, c("root", "normal"), c("cell", "truth"))
truth <- truth[grepl("^cell\\d+$", cell)]

# ------------------------------------------------------------
# 2) read MintCNA latent
# ------------------------------------------------------------
latent <- read.csv(file.path(CNA_DIR, "latent_z_dna.csv"),
                   row.names = 1, check.names = FALSE)
latent <- as.matrix(latent)

# keep only cells present in latent
truth_latent <- truth[cell %in% rownames(latent)]
truth_latent <- truth_latent[match(rownames(latent), truth_latent$cell), ]

# safety
keep <- !is.na(truth_latent$truth)
latent_use <- latent[keep, , drop = FALSE]
truth_latent <- truth_latent[keep]

# for subclone recovery, exclude normal
tumor_idx <- truth_latent$truth != "normal"
latent_tumor <- latent_use[tumor_idx, , drop = FALSE]
truth_tumor_latent <- truth_latent[tumor_idx]

y_true_common <- factor(truth_tumor_latent$truth)
common_cells  <- truth_tumor_latent$cell

# ------------------------------------------------------------
# 3) MintCNA ARI
# ------------------------------------------------------------
res_mint <- cluster_latent_kmeans(latent_tumor, seed = 1, K_max = NULL)
yhat_mint <- factor(res_mint$clusters[common_cells])
ari_mint <- adjustedRandIndex(y_true_common, yhat_mint)

# ------------------------------------------------------------
# 4) read L2R and build PCA / tSNE on same cells
#    l2r is bins x cells, so transpose to cells x bins
# ------------------------------------------------------------
# if l2r already in memory, skip reading
# l2r <- read.csv(..., row.names = 1, check.names = FALSE)
l2r <- as.matrix(l2r)

# keep same cells and same order as MintCNA tumor cells
stopifnot(all(common_cells %in% colnames(l2r)))
X <- t(l2r[, common_cells, drop = FALSE])   # now cells x bins

# optional scaling by feature (bin)
X_scaled <- scale(X)

# ------------------------------------------------------------
# 5) PCA + kmeans
# ------------------------------------------------------------
pca_fit <- prcomp(X_scaled, center = FALSE, scale. = FALSE)

# keep PCs explaining 95% variance, max 20
pve <- pca_fit$sdev^2 / sum(pca_fit$sdev^2)
npc <- min(which(cumsum(pve) >= 0.95))
npc <- min(npc, 20)

Z_pca <- pca_fit$x[, 1:npc, drop = FALSE]
rownames(Z_pca) <- rownames(X)

res_pca <- cluster_latent_kmeans(Z_pca, seed = 1, K_max = NULL)
yhat_pca <- factor(res_pca$clusters[common_cells])
ari_pca <- adjustedRandIndex(y_true_common, yhat_pca)

# ------------------------------------------------------------
# 6) tSNE + kmeans
#    use PCA initialization dimension first for stability
# ------------------------------------------------------------
set.seed(1)
tsne_fit <- Rtsne(
  X_scaled,
  dims = 2,
  initial_dims = min(50, ncol(X_scaled)),
  perplexity = min(30, floor((nrow(X_scaled) - 1) / 3)),
  check_duplicates = FALSE,
  pca = TRUE,
  verbose = FALSE
)

Z_tsne <- tsne_fit$Y
rownames(Z_tsne) <- rownames(X)

res_tsne <- cluster_latent_kmeans(Z_tsne, seed = 1, K_max = NULL)
yhat_tsne <- factor(res_tsne$clusters[common_cells])
ari_tsne <- adjustedRandIndex(y_true_common, yhat_tsne)

# ------------------------------------------------------------
# 7) rcCAE native labels on same common cells
# ------------------------------------------------------------
read_rccae_label <- function(path) {
  x <- readLines(path, warn = FALSE)
  stopifnot(length(x) >= 2)
  cells <- trimws(unlist(strsplit(x[1], ",")))
  labs  <- as.integer(trimws(unlist(strsplit(x[2], ","))))
  data.table(cell = cells, pred = labs)
}

pred_rc <- read_rccae_label(file.path(CNA_DIR, "Calling/rcCAE/out/label.txt"))

dt_rc <- merge(
  data.table(cell = common_cells, truth = y_true_common),
  pred_rc,
  by = "cell",
  all.x = FALSE,
  all.y = FALSE
)

# keep same order as common_cells
dt_rc <- dt_rc[match(common_cells, dt_rc$cell), ]

ari_rc <- adjustedRandIndex(as.factor(dt_rc$truth), as.factor(dt_rc$pred))

# ------------------------------------------------------------
# 8) final 1 x 4 dataframe
# ------------------------------------------------------------
ari_res <- data.frame(
  MintCNA = ari_mint,
  PCA     = ari_pca,
  tSNE    = ari_tsne,
  rcCAE   = ari_rc,
  row.names = "ARI"
)

print(ari_res)

##################### All files ###################
library(data.table)
library(mclust)
library(Rtsne)

# ------------------------------------------------------------
# helper: read rcCAE labels
# ------------------------------------------------------------
read_rccae_label <- function(path) {
  x <- readLines(path, warn = FALSE)
  stopifnot(length(x) >= 2)
  cells <- trimws(unlist(strsplit(x[1], ",")))
  labs  <- as.integer(trimws(unlist(strsplit(x[2], ","))))
  data.table(cell = cells, pred = labs)
}

# ------------------------------------------------------------
# helper: parse scenario name
# no c => coverage = 0.1
# no n => noise = 0
# ------------------------------------------------------------
parse_scenario_info <- function(scn) {
  cov <- if (grepl("_c0.05", scn, fixed = TRUE)) {
    0.05
  } else if (grepl("_c0.2", scn, fixed = TRUE)) {
    0.2
  } else {
    0.1
  }
  
  noi <- if (grepl("_n0.05", scn, fixed = TRUE)) {
    0.05
  } else if (grepl("_n0.1", scn, fixed = TRUE)) {
    0.1
  } else {
    0
  }
  
  data.table(
    scenario = scn,
    coverage = cov,
    noise = noi
  )
}

# ------------------------------------------------------------
# main function: compute ARI for one scenario
# ------------------------------------------------------------
compute_ari_one <- function(CNA_DIR, seed = 1) {
  
  scenario_name <- basename(CNA_DIR)
  
  # -------- 1) read truth --------
  truth <- fread(file.path(CNA_DIR, "cell_types.tsv"))
  setnames(truth, c("root", "normal"), c("cell", "truth"))
  truth <- truth[grepl("^cell\\d+$", cell)]
  
  # -------- 2) MintCNA latent --------
  latent <- read.csv(file.path(CNA_DIR, "latent_z_dna.csv"),
                     row.names = 1, check.names = FALSE)
  latent <- as.matrix(latent)
  
  truth_latent <- truth[cell %in% rownames(latent)]
  truth_latent <- truth_latent[match(rownames(latent), truth_latent$cell), ]
  
  keep <- !is.na(truth_latent$truth)
  latent_use <- latent[keep, , drop = FALSE]
  truth_latent <- truth_latent[keep]
  
  # tumor only
  tumor_idx <- truth_latent$truth != "normal"
  latent_tumor <- latent_use[tumor_idx, , drop = FALSE]
  truth_tumor <- truth_latent[tumor_idx]
  
  common_cells <- truth_tumor$cell
  y_true <- factor(truth_tumor$truth)
  
  # -------- 3) MintCNA --------
  res_mint <- cluster_latent_kmeans(latent_tumor, seed = seed, K_max = NULL)
  yhat_mint <- factor(res_mint$clusters[common_cells])
  ari_mint <- adjustedRandIndex(y_true, yhat_mint)
  
  # -------- 4) read l2r --------
  # CHANGE THIS PART TO YOUR ACTUAL FILE
  # example if saved as csv with bins x cells:
  l2r <- read.csv(file.path(CNA_DIR, "l2r_smoothed_all_cells_100kb.csv"),
                  row.names = 1, check.names = FALSE)
  l2r <- as.matrix(l2r)
  
  stopifnot(all(common_cells %in% colnames(l2r)))
  
  # cells x bins
  X <- t(l2r[, common_cells, drop = FALSE])
  X_scaled <- scale(X)
  
  # -------- 5) PCA + kmeans --------
  pca_fit <- prcomp(X_scaled, center = FALSE, scale. = FALSE)
  pve <- pca_fit$sdev^2 / sum(pca_fit$sdev^2)
  npc <- min(which(cumsum(pve) >= 0.95))
  npc <- min(npc, 20)
  
  Z_pca <- pca_fit$x[, 1:npc, drop = FALSE]
  rownames(Z_pca) <- rownames(X)
  
  res_pca <- cluster_latent_kmeans(Z_pca, seed = seed, K_max = NULL)
  yhat_pca <- factor(res_pca$clusters[common_cells])
  ari_pca <- adjustedRandIndex(y_true, yhat_pca)
  
  # -------- 6) tSNE + kmeans --------
  set.seed(seed)
  perp <- min(30, floor((nrow(X_scaled) - 1) / 3))
  perp <- max(perp, 5)
  
  tsne_fit <- Rtsne(
    X_scaled,
    dims = 2,
    initial_dims = min(50, ncol(X_scaled)),
    perplexity = perp,
    check_duplicates = FALSE,
    pca = TRUE,
    verbose = FALSE
  )
  
  Z_tsne <- tsne_fit$Y
  rownames(Z_tsne) <- rownames(X)
  
  res_tsne <- cluster_latent_kmeans(Z_tsne, seed = seed, K_max = NULL)
  yhat_tsne <- factor(res_tsne$clusters[common_cells])
  ari_tsne <- adjustedRandIndex(y_true, yhat_tsne)
  
  # -------- 7) rcCAE native --------
  pred_rc <- read_rccae_label(file.path(CNA_DIR, "Calling/rcCAE/out/label.txt"))
  
  dt_rc <- merge(
    data.table(cell = common_cells, truth = y_true),
    pred_rc,
    by = "cell",
    all = FALSE
  )
  dt_rc <- dt_rc[match(common_cells, dt_rc$cell), ]
  
  ari_rc <- adjustedRandIndex(as.factor(dt_rc$truth), as.factor(dt_rc$pred))
  
  # -------- 8) output one row --------
  info <- parse_scenario_info(scenario_name)
  
  out <- cbind(
    info,
    data.table(
      MintCNA = ari_mint,
      PCA     = ari_pca,
      tSNE    = ari_tsne,
      rcCAE   = ari_rc
    )
  )
  
  return(out)
}

# root folder containing the 9 scenario folders
SIM_ROOT <- "<MINTCNA_ROOT>/Simulation/cnasim"

scenario_dirs <- file.path(
  SIM_ROOT,
  c(
    "c3_N100_100kb_l1e4",
    "c3_N100_100kb_l1e4_c0.2",
    "c3_N100_100kb_l1e4_c0.2_n0.1",
    "c3_N100_100kb_l1e4_c0.2_n0.05",
    "c3_N100_100kb_l1e4_c0.05",
    "c3_N100_100kb_l1e4_c0.05_n0.1",
    "c3_N100_100kb_l1e4_c0.05_n0.05",
    "c3_N100_100kb_l1e4_n0.1",
    "c3_N100_100kb_l1e4_n0.05"
  )
)

ari_table <- rbindlist(lapply(scenario_dirs, compute_ari_one), fill = TRUE)

# nicer ordering
setorder(ari_table, coverage, noise)

write.csv(ari_table,file = "<MINTCNA_ROOT>/Simulation/results/ari.csv")
