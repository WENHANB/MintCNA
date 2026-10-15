load("<MINTCNA_ROOT>/Result/CRC/Benchmark/SCOPE.RData")
load("<MINTCNA_ROOT>/Result/CRC/Benchmark/MultiCNA_AE_197_d250.RData")
load("<MINTCNA_ROOT>/Result/CRC/Benchmark/MultiCNA_AE_197_single_d250.RData")
load("<MINTCNA_ROOT>/Result/CRC/Benchmark/HMMcopy.RData")
load("<MINTCNA_ROOT>/Result/CRC/Benchmark/DNAcopy.RData")
load("<MINTCNA_ROOT>/Result/CRC/Benchmark/SCCNV.RData")
load("<MINTCNA_ROOT>/Result/CRC/Benchmark/SCONCE.RData")
load("<MINTCNA_ROOT>/Result/CRC/Benchmark/rcCAE.RData")
load("<MINTCNA_ROOT>/Result/CRC/Benchmark/SCICoNE.RData")

colnames(ref_scope)[1] <- "chr"
cna_scope <- cna_scope[, colnames(cna_scope) %in% colnames(cna_multicna_ae)]
cna_scope <- cna_scope + 1L
cna_scope[cna_scope >= 5] <- 5L
cna_hmmcopy <- cna_hmmcopy[, colnames(cna_hmmcopy) %in% colnames(cna_multicna_ae)]
cna_sconce <- cna_sconce[, colnames(cna_sconce) %in% colnames(cna_multicna_ae)]
cna_sccnv <- cna_sccnv[, colnames(cna_sccnv) %in% colnames(cna_multicna_ae)]
cna_rccae <- cna_rccae[, colnames(cna_rccae) %in% colnames(cna_multicna_ae)]
cna_scicone <- cna_scicone[,colnames(cna_scicone) %in% colnames(cna_multicna_ae)]

library(data.table)
library(ggplot2)

coerce_123 <- function(mat) {
  v <- as.integer(mat)
  if (all(v %in% c(1L, 2L, 3L, NA_integer_))) {
    return(matrix(v, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat)))
  }
  
  out <- rep(NA_integer_, length(v))
  out[v %in% c(1L, 2L)] <- 1L
  out[v %in% 3L] <- 2L
  out[v %in% c(4L, 5L)] <- 3L
  matrix(out, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
}

cna_hmmcopy <- coerce_123(cna_hmmcopy)
cna_multicna_ae <- coerce_123(cna_multicna_ae)
cna_multicna_ae_single <- coerce_123(cna_multicna_ae_single)
cna_scope <- coerce_123(cna_scope)
cna_sconce <- coerce_123(cna_sconce)
cna_sccnv <- coerce_123(cna_sccnv)
cna_rccae <- coerce_123(cna_rccae)
cna_scicone <- coerce_123(cna_scicone)

methods <- list(
  "MintCNA" = cna_multicna_ae,
  "MintCNA(DNA)" = cna_multicna_ae_single,
  "SCOPE" = cna_scope,
  "HMMcopy" = cna_hmmcopy,
  "SCONCE" = cna_sconce,
  "SCCNV" = cna_sccnv,
  "rcCAE" = cna_rccae,
  "SCICoNE" = cna_scicone
)

color_map <- c(
  "MintCNA" = "#1B9E77",
  "MintCNA(DNA)" = "#D95F02",
  "SCOPE" = "#7570B3",
  "HMMcopy" = "#E7298A",
  "SCONCE" = "#4DBBD5",
  "SCCNV" = "#E6AB02",
  "rcCAE" = "#A6761D",
  "SCICoNE" = "#7A7A7A"
)

methods <- lapply(methods, as.matrix)
truth_mat <- as.matrix(cna_dnacopy)

cell_to_truthcol <- function(cell_names) {
  site <- sub("^.*_(PT|LN|ML|MP).*?$", "\\1", cell_names)
  site[!site %in% c("PT", "LN", "ML", "MP")] <- NA_character_
  
  truth_col <- fifelse(
    site == "PT", "PT",
    fifelse(
      site == "LN", "LN1",
      fifelse(site == "ML", "ML1", fifelse(site == "MP", "MP1", NA_character_))
    )
  )
  
  data.table(cell = cell_names, site = site, truth_col = truth_col)
}

safe_pearson <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 2L) return(NA_real_)
  
  x <- x[keep]
  y <- y[keep]
  if (sd(x) == 0 || sd(y) == 0) return(NA_real_)
  
  cor(x, y, method = "pearson")
}

pearson_per_cell_site_truth <- function(truth_mat, pred_mat) {
  pred_mat <- as.matrix(pred_mat)
  
  if (!is.null(rownames(truth_mat)) && !is.null(rownames(pred_mat))) {
    common_bins <- intersect(rownames(truth_mat), rownames(pred_mat))
    if (length(common_bins) == 0L) {
      stop("No overlapping rownames (bins) between truth and prediction.")
    }
    truth_use <- truth_mat[common_bins, , drop = FALSE]
    pred_use <- pred_mat[common_bins, , drop = FALSE]
  } else {
    stopifnot(nrow(truth_mat) == nrow(pred_mat))
    truth_use <- truth_mat
    pred_use <- pred_mat
  }
  
  stopifnot(!is.null(colnames(truth_use)), !is.null(colnames(pred_use)))
  
  map <- cell_to_truthcol(colnames(pred_use))
  map <- map[!is.na(truth_col) & truth_col %in% colnames(truth_use)]
  if (nrow(map) == 0L) {
    stop("No cells in pred_mat could be mapped to a truth column (PT/LN/ML/MP).")
  }
  
  pearson <- vapply(seq_len(nrow(map)), function(i) {
    cell_i <- map$cell[i]
    truth_i <- map$truth_col[i]
    safe_pearson(
      as.numeric(truth_use[, truth_i]),
      as.numeric(pred_use[, cell_i])
    )
  }, numeric(1))
  
  map[, pearson := pearson]
  map
}

pearson_tab <- rbindlist(lapply(names(methods), function(m) {
  dt <- pearson_per_cell_site_truth(truth_mat, methods[[m]])
  dt[, method := m]
  dt
}))

pearson_tab[, method := factor(method, levels = names(color_map))]

pearson_mean <- pearson_tab %>%
  group_by(method) %>%
  summarise(
    mean_pearson = median(pearson, na.rm = TRUE),
    .groups = "drop"
  )

pearson_mean

p_pearson <- ggplot(pearson_tab, aes(x = method, y = pearson, fill = method)) +
  geom_boxplot(width = 0.6, outlier.size = 0.35, linewidth = 0.3, na.rm = TRUE) +
  scale_fill_manual(values = color_map, guide = "none") +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 12),
    axis.ticks.y = element_line(linewidth = 0.4),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  coord_cartesian(ylim = c(0, 1))

print(p_pearson)

ggsave(
  "<MINTCNA_ROOT>/Result/CRC/Segmentation/Evaluation/wgs/Pearson.pdf",
  plot = p_pearson,
  width = 6.5,
  height = 6,
  units = "in",
  dpi = 600
)
