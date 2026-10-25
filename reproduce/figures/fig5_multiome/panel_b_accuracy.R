#!/usr/bin/env Rscript
# Figure 5a and 5b: accuracy of the compared methods on the 20 glioma multiome samples.
#   a  pseudobulk Pearson r against the matched bulk WGS profile, one violin per method (p_top)
#   b  precision-recall plane of the pooled bin-level calls, with constant-F1 contours (p_bot)
# Input (reproduce/figures/data/fig5_multiome/): pearson_vs_truth_own.tsv (reproduce/analysis/multiome/evaluation/pearson_vs_truth.R),
# f1_states.tsv (reproduce/analysis/multiome/evaluation/f1_states.R), sample_manifest.csv (subtype per sample).
# Output (reproduce/figures/output/fig5/): fig5a_pearson_source_data.csv, fig5b_f1_source_data.csv,
# fig5b_pooled_pr.csv and, when run on its own, fig5b_accuracy.pdf.
# Run from the root of the MintCNA repository: Rscript reproduce/figures/fig5_multiome/panel_b_accuracy.R
# (or through fig5_multiome/assemble_fig5.R).
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(cowplot)
})
FIG <- "reproduce/figures"   # run from the root of the MintCNA repository
if (!file.exists(file.path(FIG, "shared/fig_palette.R"))) stop("run this script from the root of the MintCNA repository")
source(file.path(FIG, "shared/fig_palette.R"))
DATA <- file.path(FIG, "data/fig5_multiome")
OUT  <- file.path(FIG, "output/fig5"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

tab <- fread(file.path(DATA, "pearson_vs_truth_own.tsv"))[!is.na(r)]        # reproduce/analysis/multiome/evaluation/pearson_vs_truth.R
man <- fread(file.path(DATA, "sample_manifest.csv")); setkey(man, care_id)
tab[, subtype := fifelse(grepl("^Oligo", man[care_id, subtype]), "Oligodendroglioma", "Astrocytoma")]
stopifnot(uniqueN(tab$care_id) == 20)
fwrite(tab[, .(care_id, subtype, method, pearson = r, n_bin, coverage)],
       file.path(OUT, "fig5a_pearson_source_data.csv"))

ORD <- tab[, .(m = median(r)), by = method][order(-m), method]
tab[, method := factor(method, levels = ORD)]
meds <- tab[, .(med = median(r)), by = method]

## ---- bin-level precision and recall --------------------------------------------
f1 <- fread(file.path(DATA, "f1_states.tsv"))                               # reproduce/analysis/multiome/evaluation/f1_states.R
f1[, method := factor(method, levels = ORD)]
f1 <- f1[!is.na(method)]
fwrite(f1, file.path(OUT, "fig5b_f1_source_data.csv"))
# Precision and recall are pooled over the 20 samples (bins summed), matching f1_states.R, rather
# than the median of per-sample values: epiAneufinder and inferCNV reach no modal call in most
# samples, so their per-sample F1 is undefined there and a median would describe how often they are
# undefined rather than how well they do.
pr <- f1[, {
  p <- if (sum(TP)+sum(FP) == 0) NA_real_ else sum(TP)/(sum(TP)+sum(FP))
  r <- if (sum(TP)+sum(FN) == 0) NA_real_ else sum(TP)/(sum(TP)+sum(FN))
  .(precision = p, recall = r,
    F1 = if (is.na(p)||is.na(r)||p+r == 0) NA_real_ else 2*p*r/(p+r))
}, by = method]
cat("\nprecision / recall / F1, POOLED over 20 samples:\n")
print(pr[order(-F1, na.last = TRUE)], row.names = FALSE)
fwrite(pr[order(-F1, na.last = TRUE)], file.path(OUT, "fig5b_pooled_pr.csv"))
stopifnot(nrow(pr[!is.na(precision)]) == length(unique(f1$method)))

## ---- theme and panels ------------------------------------------------------------
# Full-saturation method colours rather than the muted default: panel b identifies methods by colour
# alone and its points are about 2 mm, and panel a uses the same fills so the two share one colour key.
FILLS   <- method_fill(as.character(ORD), mode = "full")
lab_col <- label_ink(as.character(ORD))
W_MM <- 66; H_MM <- 93; MM <- 1/25.4

base <- theme_cowplot(BASE_PT, font_family = FIG_FAMILY) +
  theme(axis.text.y  = element_text(size = BASE_PT - 1, colour = "grey20"),
        axis.title   = element_text(size = BASE_PT),
        axis.line    = element_line(linewidth = 0.3, colour = "grey30"),
        axis.ticks   = element_line(linewidth = 0.3, colour = "grey30"),
        panel.grid   = element_blank(),
        plot.margin  = margin(2, 2, 1, 2),
        legend.position = "none")

## panel a: Pearson r. The y axis starts just below the data minimum (-0.28); only CONGAS+ and
## inferCNV reach below zero.
y_lo <- -0.30; y_hi <- 1.0; span <- y_hi - y_lo
i_ref <- which(ORD == "MintCNA")
top <- tab[, .(top = max(r, na.rm = TRUE)), by = method]
m2  <- merge(meds, top, by = "method")[, ytxt := pmin(top + 0.035*span, y_hi + 0.015*span)]

p_top <- ggplot(tab, aes(x = method, y = r, fill = method)) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey65", linetype = "22") +
  geom_violin(scale = "width", width = 0.86, colour = NA, alpha = 0.55, na.rm = TRUE) +
  geom_boxplot(width = 0.22, outlier.shape = NA, linewidth = 0.28,
               fill = "white", colour = "grey25", na.rm = TRUE) +
  # seeded jitter so the figure is reproducible run to run
  geom_point(position = position_jitter(width = 0.10, height = 0, seed = 1), size = 0.45, alpha = 0.5,
             colour = "grey30", na.rm = TRUE) +
  geom_text(data = m2, aes(x = method, y = ytxt, label = sprintf("%.2f", med)),
            inherit.aes = FALSE, size = PT_SMALL/.pt, family = FIG_FAMILY,
            fontface = "bold", colour = "grey15") +
  # No P values or significance stars are drawn: the medians above the violins carry the comparison.
  scale_fill_manual(values = FILLS, limits = ORD) +
  scale_x_discrete(limits = ORD, labels = method_labels(as.character(ORD))) +
  scale_y_continuous(breaks = seq(-0.25, 1, 0.25)) +
  labs(y = expression("Pearson"~italic(r)~"vs bulk WGS"), x = NULL) +
  coord_cartesian(ylim = c(y_lo, y_hi + 0.05*span), clip = "on") + base +
  # method names on panel a's own x axis, at 45 degrees; the key inside panel b repeats them
  theme(axis.text.x = element_text(size = PT_SUP, angle = 45, hjust = 1, vjust = 1,
                                   colour = lab_col))

## panel b: precision-recall plane with constant-F1 contours
fc <- rbindlist(lapply(c(0.2, 0.4, 0.6, 0.8), function(f) {
  x <- seq(f/2 + 1e-4, 1, length.out = 200)      # recall implied by precision at fixed F
  y <- f*x / (2*x - f)
  data.table(F = f, precision = x, recall = y)[recall <= 1 & recall >= 0]
}))
lab_fc <- fc[, .SD[which.min(abs(precision - 0.97))], by = F]   # contour labels at the right edge

# Diamonds for the two MintCNA configurations and circles for the other methods, with a white stroke
# so the two MintCNA marks, which nearly coincide, stay separate; MintCNA is drawn last so it sits on
# top. shape, fill and size are all mapped to method with identical limits and labels, so ggplot
# merges them into one key, which sits in the empty top-left corner instead of labels on the points.
p_bot <- ggplot() +
  geom_path(data = fc, aes(precision, recall, group = F),
            linewidth = 0.22, linetype = "22", colour = "grey72") +
  geom_text(data = lab_fc, aes(precision, recall, label = sprintf("F=%.1f", F)),
            size = PT_SMALL/.pt, colour = "grey55", hjust = 1, vjust = -0.35) +
  geom_point(data = pr[order(method == "MintCNA")], aes(precision, recall, fill = method, shape = method, size = method),
             stroke = 0.35, colour = "white") +
  scale_shape_manual(values = setNames(ifelse(as.character(ORD) %in% names(HERO), 23, 21), as.character(ORD)),
                     limits = ORD, name = NULL, labels = method_labels(as.character(ORD))) +
  scale_size_manual(values = setNames(ifelse(as.character(ORD) %in% names(HERO), 2.9, 2.4), as.character(ORD)),
                    limits = ORD, name = NULL, labels = method_labels(as.character(ORD))) +
  scale_fill_manual(values = FILLS, limits = ORD, name = NULL,
                    labels = method_labels(as.character(ORD)),
                    # one column: the labels are plotmath (MintCNA^ATAC) and two columns do not fit the square
                                         guide = guide_legend(ncol = 1,
                                         override.aes = list(stroke = 0.35))) +
  scale_x_continuous(limits = c(-0.03, 1.02), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-0.03, 1.02), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
  labs(x = "Precision", y = "Recall") + base +
  # The key sits inside the plot area, so a hairline box and a white fill make it one object that
  # is plainly not data.
  theme(legend.position = "inside", legend.position.inside = c(0.025, 0.975),
        legend.justification = c(0, 1), legend.direction = "vertical",
        legend.background = element_rect(fill = "white", colour = "grey72",
                                         linewidth = 0.25),
        legend.key = element_blank(),
        legend.margin = margin(1.2, 2.6, 1.2, 1.2), legend.spacing.y = unit(0.35, "mm"),
        legend.text = element_text(size = PT_SUP, colour = "grey20"),
        legend.key.size = unit(3.1, "mm")) +

  # Square: both axes are the same quantity on the same 0-1 scale, and the constant-F contours are
  # only correct curves at 1:1.
  theme(aspect.ratio = 1)

g_b <- plot_grid(p_top, p_bot, ncol = 1, rel_heights = c(1, 1.05))

if (!exists("FIG5_ASSEMBLING")) {
  cairo_pdf(file.path(OUT, "fig5b_accuracy.pdf"), width = W_MM*MM, height = H_MM*MM,
            family = FIG_FAMILY); print(g_b); dev.off()
  cat(sprintf("wrote fig5b_accuracy.pdf  (%d x %d mm)\n", W_MM, H_MM))
}
