#!/usr/bin/env Rscript
# Supplementary Fig. 5 -- what the RNA arm buys MintCNA on the simulation.
#   a  event-level F1 (DEL, DUP) against the RNA concordance percentage, baseline
#      simulation (0.1X, no noise); data from figS5a_concordance_compute.R
#   b-d precision / recall per gain state (CN3, CN4, CN5) for MintCNA, MintCNA(DNA)
#      and MintCNA(RNA) at 0.05X / 0.1X / 0.2X; data from
#      Simulation/cnasim/c3_N100_100kb_l1e4/RNA_benchmark/eval/s5_cn_pr_3cov.csv
#      (s5_cn_pr_3cov.R, 2026-07-16). dF1 in each panel = F1(MintCNA) - F1(MintCNA(DNA)),
#      recomputed from the stored TP/FP/FN so it is not a difference of rounded numbers.
#
# REBUILT 2026-10-09 FOR NATURE COMMUNICATIONS. The published S5 was four ggplots
# (Concor_F1.pdf + S5_B/C/D) pasted in an SVG: 160 x 250 mm, Helvetica + NotoSans-Bold
# + Liberation Sans, text 4.7-18 pt. Now 180 x 136 mm, Nimbus Sans only, 5.5-8 pt,
# lowercase letters. Layout: a as two stacked line panels on the left, b-d as a
# 3 x 3 precision-recall grid (rows = coverage, columns = CN state) on the right.
# The three MintCNA configurations use the hero family's greens (see fig_palette.R),
# not the orange/purple of the old build, so no competitor hue is reused here.
#
# NUMBERS: b-d unchanged from the published figure. Panel a shows THE PUBLISHED
# VALUES (author's decision, 2026-10-11): they are transcribed from the vector PDF
# Simulation/results/Concor_F1.pdf into data/figS5a_concordance_f1_published.csv,
# because the call files that produced the 0 % and 100 % points no longer exist on
# disk. The recomputed values (figS5a_concordance_compute.R; data/figS5a_concordance_f1.csv
# and the 2026-10-09 re-call) are kept for reference only. FIGS5A_SOURCE=computed
# switches panel a to them.

.libPaths(c("<R_LIBS_USER>", .libPaths()))
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(cowplot) })
ROOT <- "<MINTCNA_ROOT>"
source(file.path(ROOT, "figure4paper/shared/fig_palette.R"))
OUT <- Sys.getenv("FIGS5_OUT", file.path(ROOT, "figure4paper/figureS-supplementary"))
QA  <- Sys.getenv("FIGS5_QA",  file.path(ROOT, "figure4paper/figureS-supplementary/_work/qa"))
DAT <- file.path(ROOT, "figure4paper/figureS-supplementary/_work/data")
for (d in c(OUT, QA)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

# FIGS5_LAYOUT: "tall" (default since 2026-10-11, the published orientation) = a on top with
#   DEL | DUP side by side, then rows b / c / d (CN3 | CN4 | CN5), 130 x 175 mm; coverage is
#   NOT written on the figure -- the letters carry it and the legend spells it out.
#   "wide" = a as two stacked panels on the left, 3 x 3 grid with coverage strips on the right.
LAYOUT <- Sys.getenv("FIGS5_LAYOUT", "tall")
ORD <- c("MintCNA", "MintCNA(DNA)", "MintCNA(RNA)")
# Colours: the author's original S5 palette (RNA_benchmark/eval/s5_plot_BCD.R), kept at the
# author's request 2026-10-11 -- not the hero-family greens of fig_palette.R.
S5_PAL <- c("MintCNA" = "#1B9E77", "MintCNA(DNA)" = "#D95F02", "MintCNA(RNA)" = "#7570B3")

## ---- a: F1 vs concordance --------------------------------------------------------
fa <- if (Sys.getenv("FIGS5A_SOURCE", "published") == "computed") file.path(DAT, "figS5a_concordance_f1.csv") else
      file.path(DAT, "figS5a_concordance_f1_published.csv")
if (!file.exists(fa)) stop("missing: ", fa)
cat("panel a data:", basename(fa), "\n")
dA <- fread(fa)
dA[, state := factor(toupper(state), levels = c("DEL", "DUP"))]
yl <- c(floor(min(dA$F1) * 50) / 50, ceiling(max(dA$F1) * 50) / 50)
pA <- ggplot(dA, aes(percentage, F1, colour = state)) +
  geom_line(linewidth = 0.45) +
  geom_point(size = 1.3) +
  facet_wrap(~ state, ncol = if (LAYOUT == "tall") 2 else 1) +
  scale_colour_manual(values = c(DEL = LOSS, DUP = GAIN), guide = "none") +
  scale_x_continuous(breaks = c(0, 25, 50, 75, 100)) +
  scale_y_continuous(breaks = scales::breaks_pretty(n = 4)) +
  coord_cartesian(ylim = yl) +
  labs(x = "Concordance (%)", y = "F1 score") +
  theme_cowplot(BASE_PT, font_family = FIG_FAMILY) +
  theme(axis.text        = element_text(size = PT_SMALL, colour = "grey20"),
        axis.title       = element_text(size = BASE_PT),
        axis.line        = element_line(linewidth = 0.3, colour = "grey30"),
        axis.ticks       = element_line(linewidth = 0.3, colour = "grey30"),
        strip.background = element_rect(fill = "grey94", colour = NA),
        strip.text       = element_text(size = BASE_PT - 0.5, face = "bold", margin = margin(t = 1.5, b = 1.5)),
        panel.spacing.y  = unit(6, "pt"),
        panel.grid       = element_blank(),
        plot.margin      = margin(4, 4, 1, 2))

## ---- b-d: precision-recall per CN state ------------------------------------------
dB <- fread(file.path(ROOT, "Simulation/cnasim/c3_N100_100kb_l1e4/RNA_benchmark/eval/s5_cn_pr_3cov.csv"))
dB <- dB[cn %in% c("CN3", "CN4", "CN5")]
dB[, method := factor(c(`DNA+RNA` = "MintCNA", `DNA-only` = "MintCNA(DNA)", `RNA-only` = "MintCNA(RNA)")[method], levels = ORD)]
dB[, cov := factor(coverage, levels = c("0.05X", "0.1X", "0.2X"),
                   labels = if (LAYOUT == "tall") c("b", "c", "d") else c("0.05×", "0.1×", "0.2×"))]
# tall: the row strip is the panel letter, outside the panel (author, 2026-10-12); the
# coverage behind each letter is in the legend
dB[, cn  := factor(cn, levels = c("CN3", "CN4", "CN5"))]
dB[, `:=`(precision = TP / (TP + FP), recall = TP / (TP + FN))]
dB[, F1 := 2 * precision * recall / (precision + recall)]
delta <- dcast(dB, cov + cn ~ method, value.var = "F1")
delta[, d := MintCNA - `MintCNA(DNA)`]
delta[, lab := sprintf("ΔF1 = %+.3f", d)]
cat("dF1 (MintCNA - MintCNA(DNA)):\n"); print(dcast(delta, cov ~ cn, value.var = "d"), digits = 3)
fwrite(dB[, .(coverage, cn, method, TP, FP, FN, precision, recall, F1)], file.path(DAT, "figS5bcd_cn_state_pr.csv"))

iso <- rbindlist(lapply(c(0.2, 0.4, 0.6, 0.8), function(f) {
  P <- seq(f/2 + 1e-3, 1, length.out = 300)
  data.table(f = f, precision = P, recall = f*P/(2*P - f))
}))[recall <= 1]
tag <- data.table(cov = factor(levels(dB$cov), levels = levels(dB$cov)),
                  cn = factor("CN3", levels = levels(dB$cn)),
                  precision = 0.02, recall = 1.0, lab = c("b", "c", "d"))
pB <- ggplot(dB, aes(precision, recall)) +
  geom_line(data = iso, aes(precision, recall, group = f), inherit.aes = FALSE,
            colour = "grey82", linewidth = 0.22) +
  geom_point(aes(colour = method), shape = 17, size = 2.7) +   # solid triangles, no outline (as the published S5); 1.8 was too small
  geom_label(data = delta, aes(x = 0.98, y = 0.03, label = lab), inherit.aes = FALSE,
             hjust = 1, vjust = 0, size = PT_SMALL / .pt, family = FIG_FAMILY, colour = "grey20",
             fill = "white", linewidth = 0, label.padding = unit(0.5, "pt"), label.r = unit(0, "pt")) +
  (if (LAYOUT == "tall") NULL else
     geom_text(data = tag, aes(precision, recall, label = lab), inherit.aes = FALSE,
               hjust = 0, vjust = 1, size = 8/.pt, fontface = "bold", family = FIG_FAMILY)) +
  scale_colour_manual(values = S5_PAL, breaks = ORD, labels = method_labels(ORD), name = NULL) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
  ggh4x::facet_grid2(cov ~ cn, switch = "y", axes = "all", remove_labels = "none") +   # both axes on every panel
  labs(x = "Precision", y = "Recall") +
  theme_cowplot(BASE_PT, font_family = FIG_FAMILY) +
  theme(aspect.ratio      = 1,
        axis.text         = element_text(size = PT_SMALL, colour = "grey20"),
        axis.title        = element_text(size = BASE_PT),
        axis.line         = element_line(linewidth = 0.3, colour = "grey30"),
        axis.ticks        = element_line(linewidth = 0.3, colour = "grey30"),
        strip.background  = element_rect(fill = "grey94", colour = NA),
        strip.placement   = "outside",
        strip.text.x      = element_text(size = BASE_PT - 0.5, face = "bold", margin = margin(t = 1.5, b = 1.5)),
        strip.text.y.left = element_text(size = BASE_PT - 0.5, face = "bold", angle = 90, margin = margin(l = 1.5, r = 1.5)),
        panel.spacing.x   = unit(4, "pt"),
        panel.spacing.y   = unit(5, "pt"),
        panel.grid        = element_blank(),
        legend.position   = "bottom",
        legend.justification = "center",                    # theme_cowplot left-aligns it otherwise
        legend.key.spacing.x = unit(3, "mm"),               # explicit gap: cairo drew labels wider than the layout allowed (2026-10-13)
        legend.text       = element_text(size = PT_SUP, margin = margin(l = 1.2, r = 1.5, unit = "mm")),   # carries MintCNA^DNA / ^RNA
        legend.key.size   = unit(3.2, "mm"),
        legend.margin     = margin(t = 1),
        plot.margin       = margin(4, 3, 1, 2)) +
  guides(colour = guide_legend(nrow = 1, override.aes = list(size = 3)))
# tall layout: the left strip IS the panel letter -- 8 pt bold, top-left of each row, no box
# (must come AFTER the theme above)
if (LAYOUT == "tall") pB <- pB + theme(strip.text.y.left = element_text(size = 8, face = "bold", angle = 0,
                                                                        vjust = 1, hjust = 0, margin = margin(t = 0, r = 2)),
                                       strip.background.y = element_blank())

if (LAYOUT == "tall") {
  # "a" sits on the same vertical line as the b/c/d strips (which ggplot puts inside the
  # y-axis title column, ~5 mm in from the left edge)
  fig <- plot_grid(pA, pB, ncol = 1, rel_heights = c(0.25, 0.75), labels = c("a", ""), label_size = 8,
                   label_fontface = "bold", label_fontfamily = FIG_FAMILY, label_x = 5.3 / 130, label_y = 1,
                   hjust = 0, vjust = 1)
  W_MM <- 130; H_MM <- 175
} else {
  fig <- plot_grid(pA, pB, ncol = 2, rel_widths = c(0.30, 0.70), labels = c("a", ""), label_size = 8,
                   label_fontface = "bold", label_fontfamily = FIG_FAMILY, label_x = 0, label_y = 1,
                   hjust = 0, vjust = 1)
  W_MM <- 165; H_MM <- 125
}
mm <- 1 / 25.4
base <- "SuppFig05_RNA_concordance_CN"
cairo_pdf(file.path(OUT, paste0(base, ".pdf")), width = W_MM * mm, height = H_MM * mm, family = FIG_FAMILY)
print(fig); dev.off()
suppressWarnings(pdftools::pdf_convert(file.path(OUT, paste0(base, ".pdf")), format = "png", dpi = 200,
                                       filenames = file.path(QA, paste0(base, "_preview.png")), verbose = FALSE))
cat(sprintf("wrote %s.pdf  (%d x %d mm)\n", base, W_MM, H_MM))
