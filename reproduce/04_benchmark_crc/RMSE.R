load("<MINTCNA_ROOT>/Result/CRC/Benchmark/SCOPE.RData")
load("<MINTCNA_ROOT>/Result/CRC/Benchmark/MultiCNA_AE_197_d250.RData")
load("<MINTCNA_ROOT>/Result/CRC/Benchmark/HMMcopy.RData")
load("<MINTCNA_ROOT>/Result/CRC/Benchmark/DNAcopy.RData")
load("<MINTCNA_ROOT>/Result/CRC/Benchmark/SCCNV.RData")
load("<MINTCNA_ROOT>/Result/CRC/Benchmark/SCONCE.RData")
load("<MINTCNA_ROOT>/Result/CRC/Benchmark/rcCAE.RData")
load("<MINTCNA_ROOT>/Result/CRC/Benchmark/MultiCNA_AE_197_single_d250.RData")

colnames(ref_scope)[1] <- "chr";
cna_scope <- cna_scope[,colnames(cna_scope) %in% colnames(cna_multicna_ae)]
cna_scope <- cna_scope + 1L; cna_scope[cna_scope >= 5] <- 5L 
cna_hmmcopy <- cna_hmmcopy[,colnames(cna_hmmcopy) %in% colnames(cna_multicna_ae)]
cna_sconce <- cna_sconce[,colnames(cna_sconce) %in% colnames(cna_multicna_ae)]
cna_sccnv <- cna_sccnv[,colnames(cna_sccnv) %in% colnames(cna_multicna_ae)]
cna_rccae <- cna_rccae[,colnames(cna_rccae) %in% colnames(cna_multicna_ae)]

library(data.table)
library(viridis)
library(ggplot2)

## --- helpers from before (keep as-is) ---
coerce_123 <- function(mat) {
  v <- as.integer(mat)
  # If already 1/2/3 only, keep as-is
  if (all(v %in% c(1L,2L,3L, NA_integer_))) {
    return(matrix(v, nrow = nrow(mat), ncol = ncol(mat),
                  dimnames = dimnames(mat)))
  }
  # Otherwise map 1..5 -> 1/2/3
  out <- rep(NA_integer_, length(v))
  out[v %in% c(1L,2L)] <- 1L  # deletions
  out[v %in% 3L]        <- 2L  # neutral
  out[v %in% c(4L,5L)] <- 3L  # duplications
  matrix(out, nrow = nrow(mat), ncol = ncol(mat),
         dimnames = dimnames(mat))
}

cna_hmmcopy  = coerce_123(cna_hmmcopy)
cna_multicna_ae = coerce_123(cna_multicna_ae)
cna_scope    = coerce_123(cna_scope)
cna_sconce    = coerce_123(cna_sconce)
cna_sccnv    = coerce_123(cna_sccnv)
cna_rccae    = coerce_123(cna_rccae)

##### RMSE ### 

library(data.table)
library(ggplot2)

# -----------------------------
# 0) Put predictions into a list
# -----------------------------
methods <- list(
  "MintCNA"      = cna_multicna_ae,
  "MintCNA(DNA)"      = cna_multicna_ae_single,
  "SCOPE"        = cna_scope,
  "HMMcopy"      = cna_hmmcopy,
  "SCONCE"       = cna_sconce,
  "SCCNV"        = cna_sccnv,
  "rcCAE"        = cna_rccae,
  "SCICoNE"    = cna_scicone   # add if you have it
)

color_map <- c(
  "MintCNA"      = "#1B9E77",
  "MintCNA(DNA)" = "#D95F02",
  "SCOPE"        = "#7570B3",
  "HMMcopy"      = "#E7298A",
  "SCONCE"       = "#4DBBD5",  # changed (was #66A61E)
  "SCCNV"        = "#E6AB02",
  "rcCAE"        = "#A6761D",
  "SCICoNE"      = "#7A7A7A"
)

method_order <- c("MintCNA",
                  "MintCNA(DNA)",
                  #"MintCNA(AE)",
                  "SCOPE","HMMcopy","SCONCE","SCCNV","rcCAE","SCICoNE")

# (optional) ensure matrices
methods <- lapply(methods, as.matrix)
truth_mat <- as.matrix(cna_dnacopy)   # columns are site-level truth: PT, LN1, ML1, MP1

# -----------------------------
# 1) cell -> site -> truth column
# -----------------------------
cell_to_truthcol <- function(cell_names) {
  # extract site tag from cell id; expects "..._PT..." or "..._LN..." etc
  site <- sub("^.*_(PT|LN|ML|MP).*?$", "\\1", cell_names)
  site[!site %in% c("PT","LN","ML","MP")] <- NA_character_
  
  truth_col <- fifelse(site == "PT", "PT",
                       fifelse(site == "LN", "LN1",
                               fifelse(site == "ML", "ML1",
                                       fifelse(site == "MP", "MP1", NA_character_))))
  data.table(cell = cell_names, site = site, truth_col = truth_col)
}

# -----------------------------
# 2) per-cell RMSE (each cell uses its own site truth column)
# -----------------------------
rmse_per_cell_site_truth <- function(truth_mat, pred_mat) {
  pred_mat <- as.matrix(pred_mat)
  
  # align rows (bins) safely
  if (!is.null(rownames(truth_mat)) && !is.null(rownames(pred_mat))) {
    common_bins <- intersect(rownames(truth_mat), rownames(pred_mat))
    if (length(common_bins) == 0L) stop("No overlapping rownames (bins) between truth and prediction.")
    Tbin <- truth_mat[common_bins, , drop = FALSE]
    Pbin <- pred_mat[common_bins, , drop = FALSE]
  } else {
    stopifnot(nrow(truth_mat) == nrow(pred_mat))
    Tbin <- truth_mat
    Pbin <- pred_mat
  }
  
  stopifnot(!is.null(colnames(Tbin)), !is.null(colnames(Pbin)))
  
  map <- cell_to_truthcol(colnames(Pbin))
  map <- map[!is.na(truth_col) & truth_col %in% colnames(Tbin)]
  if (nrow(map) == 0L) stop("No cells in pred_mat could be mapped to a truth column (PT/LN/ML/MP).")
  
  rmse <- vapply(seq_len(nrow(map)), function(i) {
    cell_i <- map$cell[i]
    tc_i   <- map$truth_col[i]
    sqrt(mean((Tbin[, tc_i] - Pbin[, cell_i])^2, na.rm = TRUE))
  }, numeric(1))
  
  map[, rmse := as.numeric(rmse)]
  map
}

# -----------------------------
# 3) compute RMSE table for all methods
# -----------------------------
rmse_tab <- rbindlist(lapply(names(methods), function(m) {
  dt <- rmse_per_cell_site_truth(truth_mat, methods[[m]])
  dt[, method := m]
  dt
}))

# -----------------------------
# 4) boxplot (uses your color_map)
# -----------------------------
rmse_tab[, method := factor(method, levels = names(color_map))]

rmse_mean <- rmse_tab %>%
  group_by(method) %>%
  summarise(
    mean_rmse = median(rmse, na.rm = TRUE),
    .groups = "drop"
  )

rmse_mean

p_rmse <- ggplot(rmse_tab, aes(x = method, y = rmse, fill = method)) +
  geom_boxplot(width = 0.6, outlier.size = 0.35, linewidth = 0.3) +
  scale_fill_manual(values = color_map, guide = "none") +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 12),
    axis.ticks.y = element_line(linewidth = 0.4),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  coord_cartesian(ylim = c(0,1))

print(p_rmse)

ggsave("<MINTCNA_ROOT>/Result/CRC/Segmentation/Evaluation/wgs/RMSE.pdf", 
       plot = p_rmse, width = 6.5, height = 6, units = "in", dpi = 600)
