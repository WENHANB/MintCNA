library(data.table)
library(ggplot2)

# per-cell RMSE across bins
rmse_per_cell <- function(truth_mat, pred_mat) {
  stopifnot(nrow(truth_mat) == nrow(pred_mat))
  stopifnot(!is.null(colnames(truth_mat)), !is.null(colnames(pred_mat)))
  
  keep_cells <- intersect(colnames(truth_mat), colnames(pred_mat))
  if (length(keep_cells) == 0L) stop("No overlapping cell names.")
  
  Tsub <- truth_mat[, keep_cells, drop = FALSE]
  Psub <- pred_mat[,  keep_cells, drop = FALSE]
  
  rmse_i <- sqrt(colMeans((Tsub - Psub)^2, na.rm = TRUE))
  data.table(cell = keep_cells, rmse = as.numeric(rmse_i), n_cells_used = length(keep_cells))
}

# compute for all methods
rmse_tab <- rbindlist(lapply(names(methods), function(m) {
  dt <- rmse_per_cell(truth_mat, methods[[m]])
  dt[, method := m]
  dt
}))

# map names/order
rmse_tab[, method_disp := name_map[method]]
rmse_tab[, method_disp := factor(method_disp, levels = method_order)]

# boxplot
p_rmse <- ggplot(rmse_tab, aes(x = method_disp, y = rmse, fill = method_disp)) +
  geom_boxplot(width = 0.6, outlier.size = 0.35, linewidth = 0.3) +
  scale_fill_manual(values = color_map, guide = "none") +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, size = 11),
    axis.text.y = element_text(size = 12),
    axis.ticks.y = element_line(linewidth = 0.4)
  ) +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  coord_cartesian(ylim = c(0, 1))   # better than ylim() (doesn't drop points)

print(p_rmse)
ggsave(file.path(CNA_DIR,"Result/RMSE.pdf"), plot = p_rmse, width = 6.5, height = 6, units = "in", dpi = 600)
