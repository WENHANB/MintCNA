#!/usr/bin/env Rscript
# Figure 2 -- simulation benchmark as a 3 x 3 radar grid: rows = noise 0 / 0.05 / 0.1,
# columns = coverage 0.05 / 0.1 / 0.2, five spokes P / R / F1 / DEL / DUP on the raw
# 0-1 scale. This is the author's original Fig 2
# (Simulation/code/Plot_multi_setting_radar_bb.R -> F1_plots/multi_setting/5_radar_3x3.pdf,
# 2026-07-24) rebuilt at Nature Communications size: 160 x 170 mm, Helvetica/Nimbus,
# 5.5-8 pt, lowercase letters. Same data (Simulation/results/tables/merged), same
# geometry (rings at 0.25/0.5/0.75/1, spokes, filled polygons), same two-green palette.
# RESTORED 2026-10-13: the merged "radar + boundary bias" Fig 2 of 2026-08-21 is kept in
# figure4paper/_archive/figure2-simulation-merged-2026-08-21/; the manuscript has Fig 2
# radar and Fig 3 boundary bias + RMSE as separate figures.
# Method order = the advisor's order used in Supp Figs 2-4 (two legend rows filled across).

.libPaths(c("<R_LIBS_USER>", .libPaths()))
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(cowplot) })
ROOT <- "<MINTCNA_ROOT>"
source(file.path(ROOT, "figure4paper/shared/fig_palette.R"))
MERGED <- file.path(ROOT, "Simulation/results/tables/merged")
HERE <- file.path(ROOT, "figure4paper/figure2-radar")
OUT <- file.path(HERE, "output"); DATA <- file.path(HERE, "data"); QA <- file.path(HERE, "qa")
for (d in c(OUT, DATA, QA)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

ORD <- c("MintCNA", "MintCNA(DNA)", "rcCAE", "SCOPE", "SCONCE", "SCCNV", "SCICoNE", "HMMcopy")
LAB <- method_labels(ORD)
COLS <- METHOD_HUES[ORD]

one <- fread(file.path(MERGED, "summary_one.csv"))[method %in% ORD]
del <- fread(file.path(MERGED, "summary_del.csv"))[method %in% ORD, .(coverage, noise, method, DEL = F1)]
dup <- fread(file.path(MERGED, "summary_dup.csv"))[method %in% ORD, .(coverage, noise, method, DUP = F1)]
src <- merge(merge(one[, .(coverage, noise, method, P = precision, R = recall, F1)], del,
                   by = c("coverage", "noise", "method")), dup, by = c("coverage", "noise", "method"))
fwrite(src[order(coverage, noise, method)], file.path(DATA, "fig2_radar_metrics.csv"))

SPOKES <- c("P", "R", "F1", "DEL", "DUP")
rad <- melt(src, id.vars = c("coverage", "noise", "method"), measure.vars = SPOKES,
            variable.name = "metric", value.name = "value")
rad[is.na(value), value := 0]                                  # as the original script
rad[, metric := factor(metric, levels = SPOKES)]
rad[, angle := (as.integer(metric) - 1) * 2 * pi / length(SPOKES)]   # P at 3 o'clock, anticlockwise
rad[, `:=`(x = value * cos(angle), y = value * sin(angle))]
rad[, method := factor(method, levels = ORD)]
rad[, cov_f := factor(coverage, levels = c(0.05, 0.1, 0.2), labels = c("Coverage: 0.05", "Coverage: 0.1", "Coverage: 0.2"))]
rad[, noi_f := factor(noise, levels = c(0, 0.05, 0.1), labels = c("Noise: 0", "Noise: 0.05", "Noise: 0.1"))]
setorder(rad, cov_f, noi_f, method, metric)
closed <- rad[, rbind(.SD, .SD[1]), by = .(cov_f, noi_f, method), .SDcols = c("x", "y")]

rings <- rbindlist(lapply(c(0.25, 0.5, 0.75, 1), function(r) {
  th <- seq(0, 2 * pi, length.out = 120); data.table(r = r, x = r * cos(th), y = r * sin(th)) }))
spokes <- data.table(metric = SPOKES, angle = (seq_along(SPOKES) - 1) * 2 * pi / length(SPOKES))
spokes[, `:=`(xend = cos(angle), yend = sin(angle), lx = 1.22 * cos(angle), ly = 1.22 * sin(angle))]

p <- ggplot() +
  geom_path(data = rings, aes(x, y, group = r), colour = "grey86", linewidth = 0.2) +
  geom_segment(data = spokes, aes(x = 0, y = 0, xend = xend, yend = yend), colour = "grey86", linewidth = 0.2) +
  geom_polygon(data = closed, aes(x, y, group = method, fill = method), alpha = 0.06) +
  geom_path(data = closed, aes(x, y, group = method, colour = method), linewidth = 0.35) +
  geom_point(data = rad, aes(x, y, colour = method), size = 0.7) +
  geom_text(data = spokes, aes(lx, ly, label = metric), size = (BASE_PT - 0.5) / .pt, fontface = "bold",
            colour = "grey30", family = FIG_FAMILY) +
  scale_colour_manual(values = COLS, breaks = ORD, labels = LAB, name = NULL) +
  scale_fill_manual(values = COLS, breaks = ORD, labels = LAB, name = NULL) +
  coord_fixed(xlim = c(-1.32, 1.32), ylim = c(-1.32, 1.32), expand = FALSE) +
  facet_grid(noi_f ~ cov_f) +
  theme_cowplot(BASE_PT, font_family = FIG_FAMILY) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "grey80", fill = NA, linewidth = 0.3),
        panel.spacing = unit(3, "pt"),
        strip.background = element_rect(fill = "grey94", colour = NA),
        strip.text.x = element_text(size = BASE_PT, face = "bold", margin = margin(t = 2, b = 2)),
        strip.text.y = element_text(size = BASE_PT, face = "bold", angle = 270, margin = margin(l = 2, r = 2)),
        legend.position = "bottom", legend.justification = "center", legend.byrow = TRUE,
        legend.key.spacing.x = unit(3, "mm"),               # explicit gap: cairo drew labels wider than the layout allowed (2026-10-13)
        legend.text = element_text(size = PT_SUP, margin = margin(l = 1.2, r = 1.5, unit = "mm")), legend.key.size = unit(3.2, "mm"),
        legend.margin = margin(t = 2), plot.margin = margin(2, 2, 1, 2)) +
  guides(colour = guide_legend(nrow = 2, override.aes = list(size = 1.8, linewidth = 0.6)),
         fill = guide_legend(nrow = 2))

W_MM <- 160; H_MM <- 170; mm <- 1 / 25.4
cairo_pdf(file.path(OUT, "fig2_radar.pdf"), width = W_MM * mm, height = H_MM * mm, family = FIG_FAMILY)
print(p); dev.off()
ragg::agg_png(file.path(QA, "fig2_radar_preview.png"), width = W_MM * mm, height = H_MM * mm, units = "in", res = 600)
print(p); dev.off()
cat(sprintf("wrote fig2_radar.pdf  %d x %d mm\n", W_MM, H_MM))
