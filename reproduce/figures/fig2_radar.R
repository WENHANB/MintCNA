#!/usr/bin/env Rscript
# Figure 2: the simulation benchmark as a 3 x 3 radar grid (rows noise 0 / 0.05 / 0.1, columns
# coverage 0.05 / 0.1 / 0.2) with five spokes, P / R / F1 / DEL / DUP, on the raw 0-1 scale.
# Input : reproduce/figures/data/simulation_tables/{summary_one,summary_del,summary_dup}.csv
#         (reproduce/analysis/simulation/evaluate.R per setting, merged by merge_tables.R)
# Output: reproduce/figures/output/fig2/fig2_radar.pdf, a 600-dpi preview PNG, and the plotted
#         values as fig2_radar_metrics.csv
# Run from the root of the MintCNA repository: Rscript reproduce/figures/fig2_radar.R

suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(cowplot) })
FIG <- "reproduce/figures"
if (!file.exists(file.path(FIG, "shared/fig_palette.R"))) stop("run this script from the root of the MintCNA repository")
source(file.path(FIG, "shared/fig_palette.R"))
MERGED <- file.path(FIG, "data/simulation_tables")
OUT <- file.path(FIG, "output/fig2"); DATA <- OUT; QA <- OUT
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# Method order of the legend (two rows, filled across); Fig 3 uses the same order.
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
rad[is.na(value), value := 0]                                  # an undefined metric is drawn at the centre
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
        legend.key.spacing.x = unit(3, "mm"),               # explicit key gap: cairo draws the labels wider than the layout allows for
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
