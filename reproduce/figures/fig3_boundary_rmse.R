#!/usr/bin/env Rscript
# Figure 3: the simulation benchmark as per-cell distributions. a = boundary bias (Mb), b = RMSE,
# each a 3 x 3 raincloud grid (rows noise 0 / 0.05 / 0.1, columns coverage 0.05 / 0.1 / 0.2) of
# eight methods: nudged half violin + box + jittered cells, a above b, one shared legend.
# Input : reproduce/figures/data/simulation_tables/{bb_percell,RMSE_percell}.csv
#         (reproduce/analysis/simulation/evaluate.R per setting, merged by merge_tables.R)
# Output: reproduce/figures/output/fig3/fig3_bb_rmse.pdf, a 600-dpi preview PNG, and the plotted
#         values as fig3a_boundary_bias_percell.csv and fig3b_rmse_percell.csv
# Run from the root of the MintCNA repository: Rscript reproduce/figures/fig3_boundary_rmse.R

suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(cowplot) })
FIG <- "reproduce/figures"
if (!file.exists(file.path(FIG, "shared/fig_palette.R"))) stop("run this script from the root of the MintCNA repository")
source(file.path(FIG, "shared/fig_palette.R"))
MERGED <- file.path(FIG, "data/simulation_tables")
OUT <- file.path(FIG, "output/fig3"); DATA <- OUT; QA <- OUT
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# Method order of the legend, the same as in Fig 2.
ORD <- c("MintCNA", "MintCNA(DNA)", "rcCAE", "SCOPE", "SCONCE", "SCCNV", "SCICoNE", "HMMcopy")
LAB <- method_labels(ORD)
COLS <- METHOD_HUES[ORD]
set.seed(1)                                                    # jitter is reproducible

prep <- function(d, ycol) {
  d <- d[method %in% ORD & !is.na(get(ycol))]
  d[, method := factor(method, levels = ORD)]
  d[, cov_f := factor(coverage, levels = c(0.05, 0.1, 0.2), labels = c("Coverage: 0.05", "Coverage: 0.1", "Coverage: 0.2"))]
  d[, noi_f := factor(noise, levels = c(0, 0.05, 0.1), labels = c("Noise: 0", "Noise: 0.05", "Noise: 0.1"))]
  d[, x_pt := as.numeric(method) - 0.18 + runif(.N, -0.06, 0.06)]
  d
}
bb <- prep(fread(file.path(MERGED, "bb_percell.csv")), "BB_cell_mb")
rm <- prep(fread(file.path(MERGED, "RMSE_percell.csv")), "rmse")
fwrite(bb[, .(coverage, noise, method, cell, BB_cell_mb)], file.path(DATA, "fig3a_boundary_bias_percell.csv"))
fwrite(rm[, .(coverage, noise, method, cell, rmse)],       file.path(DATA, "fig3b_rmse_percell.csv"))

rain <- function(d, ycol, ylim, ylab, log_y = FALSE) {
  ggplot(d, aes(x = method, y = .data[[ycol]])) +
    geom_violin(aes(fill = method), trim = TRUE, scale = "width", width = 0.85, alpha = 0.3,
                colour = NA, position = position_nudge(x = 0.18)) +
    geom_boxplot(aes(fill = method), width = 0.18, outlier.shape = NA, linewidth = 0.25,
                 colour = "grey30", alpha = 0.7) +
    geom_point(aes(x = x_pt, colour = method), size = 0.22, alpha = 0.45, stroke = 0) +
    scale_fill_manual(values = COLS, breaks = ORD, labels = LAB, name = NULL) +
    scale_colour_manual(values = COLS, breaks = ORD, labels = LAB, name = NULL) +
    (if (log_y) scale_y_log10(breaks = c(0.5, 1, 3, 10, 30), labels = c("0.5", "1", "3", "10", "30"))
     else scale_y_continuous(breaks = c(0.25, 0.5, 0.75, 1))) +
    coord_cartesian(ylim = ylim) +
    facet_grid(noi_f ~ cov_f) +
    labs(x = NULL, y = ylab) +
    theme_cowplot(BASE_PT, font_family = FIG_FAMILY) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = PT_SMALL, colour = "grey20"),
          axis.title.y = element_text(size = BASE_PT),
          axis.line = element_line(linewidth = 0.3, colour = "grey30"),
          axis.ticks = element_line(linewidth = 0.3, colour = "grey30"),
          panel.border = element_rect(colour = "grey80", fill = NA, linewidth = 0.3),
          panel.spacing = unit(3, "pt"), panel.grid = element_blank(),
          strip.background = element_rect(fill = "grey94", colour = NA),
          strip.text.x = element_text(size = BASE_PT - 0.5, face = "bold", margin = margin(t = 1.5, b = 1.5)),
          strip.text.y = element_text(size = BASE_PT - 0.5, face = "bold", angle = 270, margin = margin(l = 1.5, r = 1.5)),
          legend.position = "none", plot.margin = margin(3, 2, 1, 2))
}
# Panel a is drawn on a log10 axis clipped to 0.3-40 Mb, so that the large SCICoNE medians do not
# flatten the other methods; the cells outside the window are counted below and stated in the
# figure legend. Panel b starts at the data floor (rounded down to a multiple of 0.05) instead of 0.
BB_LIM <- c(0.3, 40)
out_n <- bb[BB_cell_mb < BB_LIM[1] | BB_cell_mb > BB_LIM[2], .N, by = method][order(-N)]
cat(sprintf("panel a log window [%.1f, %.0f] Mb: %d of %d cells outside -> %s\n", BB_LIM[1], BB_LIM[2],
            sum(out_n$N), nrow(bb), paste(sprintf("%s %d", out_n$method, out_n$N), collapse = ", ")))
RM_LO <- floor(min(rm$rmse) * 20) / 20
cat(sprintf("panel b y from %.2f (data minimum %.3f) to 1\n", RM_LO, min(rm$rmse)))
pa <- rain(bb, "BB_cell_mb", BB_LIM, "Boundary bias (Mb)", log_y = TRUE)
pb <- rain(rm, "rmse", c(RM_LO, 1), "Root mean square error")
leg <- get_plot_component(pa + theme(legend.position = "bottom", legend.justification = "center",
                                     legend.key.spacing.x = unit(3, "mm"),   # explicit key gap, as in fig2_radar.R
                                     legend.text = element_text(size = PT_SUP, margin = margin(l = 1.2, r = 1.5, unit = "mm")),
                                     legend.key.size = unit(3.2, "mm"), legend.margin = margin(t = 0)) +
                            guides(fill = guide_legend(nrow = 1, override.aes = list(alpha = 1)),
                                   colour = guide_legend(nrow = 1, override.aes = list(size = 2, alpha = 1))),
                          "guide-box-bottom", return_all = TRUE)
fig <- plot_grid(plot_grid(pa, pb, ncol = 1, align = "v", axis = "lr", labels = c("a", "b"), label_size = 8,
                           label_fontface = "bold", label_fontfamily = FIG_FAMILY, label_x = 0.002, label_y = 1,
                           hjust = 0, vjust = 1),
                 leg, ncol = 1, rel_heights = c(1, 0.05))
W_MM <- 180; H_MM <- 165; mm <- 1 / 25.4
cairo_pdf(file.path(OUT, "fig3_bb_rmse.pdf"), width = W_MM * mm, height = H_MM * mm, family = FIG_FAMILY)
print(fig); dev.off()
ragg::agg_png(file.path(QA, "fig3_bb_rmse_preview.png"), width = W_MM * mm, height = H_MM * mm, units = "in", res = 600)
print(fig); dev.off()
cat(sprintf("wrote fig3_bb_rmse.pdf  %d x %d mm\n", W_MM, H_MM))
