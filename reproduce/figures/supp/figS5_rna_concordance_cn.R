#!/usr/bin/env Rscript
# Supplementary Fig. 5 -- what the RNA arm adds to MintCNA on the simulation.
#   a    event-level F1 (DEL, DUP) against the RNA concordance percentage, baseline simulation
#        (0.1X, no noise). Panel a plots the values of the published figure, stored in
#        reproduce/figures/data/supp/figS5a_concordance_f1_published.csv.
#   b-d  precision and recall per gain state (CN3, CN4, CN5) for MintCNA, MintCNA(DNA) and MintCNA(RNA)
#        at 0.05X / 0.1X / 0.2X, from reproduce/figures/data/supp/s5_cn_pr_3cov.csv (produced by
#        reproduce/analysis/simulation/cn_state_precision_recall.R).
# Output: reproduce/figures/output/supp/SuppFig05_RNA_concordance_CN.pdf (130 x 175 mm), a PNG preview
#         and figS5bcd_cn_state_pr.csv.
# Run from the root of the MintCNA repository: Rscript reproduce/figures/supp/figS5_rna_concordance_cn.R

suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(cowplot) })
FIG <- "reproduce/figures"   # run from the root of the MintCNA repository
if (!file.exists(file.path(FIG, "shared/fig_palette.R"))) stop("run this script from the root of the MintCNA repository")
source(file.path(FIG, "shared/fig_palette.R"))
DAT <- file.path(FIG, "data/supp")
OUT <- file.path(FIG, "output/supp"); QA <- OUT; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# Arrangement: a on top with DEL | DUP side by side, then rows b / c / d (one coverage each) with
# CN3 | CN4 | CN5 as columns. Coverage is not written on the figure: the panel letters carry it
# and the legend spells it out.
ORD <- c("MintCNA", "MintCNA(DNA)", "MintCNA(RNA)")
# colours of the published figure; the three MintCNA configurations do not use the shared method palette
S5_PAL <- c("MintCNA" = "#1B9E77", "MintCNA(DNA)" = "#D95F02", "MintCNA(RNA)" = "#7570B3")

## ---- a: F1 vs concordance --------------------------------------------------------
fa <- file.path(DAT, "figS5a_concordance_f1_published.csv")
if (!file.exists(fa)) stop("missing: ", fa)
cat("panel a data:", basename(fa), "\n")
dA <- fread(fa)
dA[, state := factor(toupper(state), levels = c("DEL", "DUP"))]
yl <- c(floor(min(dA$F1) * 50) / 50, ceiling(max(dA$F1) * 50) / 50)
pA <- ggplot(dA, aes(percentage, F1, colour = state)) +
  geom_line(linewidth = 0.45) +
  geom_point(size = 1.3) +
  facet_wrap(~ state, ncol = 2) +
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
dB <- fread(file.path(DAT, "s5_cn_pr_3cov.csv"))   # reproduce/analysis/simulation/cn_state_precision_recall.R
dB <- dB[cn %in% c("CN3", "CN4", "CN5")]
dB[, method := factor(c(`DNA+RNA` = "MintCNA", `DNA-only` = "MintCNA(DNA)", `RNA-only` = "MintCNA(RNA)")[method], levels = ORD)]
dB[, cov := factor(coverage, levels = c("0.05X", "0.1X", "0.2X"),
                   labels = c("b", "c", "d"))]
# the row strip carries the panel letter; the coverage behind each letter is given in the legend
dB[, cn  := factor(cn, levels = c("CN3", "CN4", "CN5"))]
dB[, `:=`(precision = TP / (TP + FP), recall = TP / (TP + FN))]
dB[, F1 := 2 * precision * recall / (precision + recall)]
# dF1 = F1(MintCNA) - F1(MintCNA(DNA)), computed from the stored TP/FP/FN rather than from rounded F1 values
delta <- dcast(dB, cov + cn ~ method, value.var = "F1")
delta[, d := MintCNA - `MintCNA(DNA)`]
delta[, lab := sprintf("ΔF1 = %+.3f", d)]
cat("dF1 (MintCNA - MintCNA(DNA)):\n"); print(dcast(delta, cov ~ cn, value.var = "d"), digits = 3)
fwrite(dB[, .(coverage, cn, method, TP, FP, FN, precision, recall, F1)], file.path(OUT, "figS5bcd_cn_state_pr.csv"))

# iso-F1 contours at F1 = 0.2, 0.4, 0.6, 0.8
iso <- rbindlist(lapply(c(0.2, 0.4, 0.6, 0.8), function(f) {
  P <- seq(f/2 + 1e-3, 1, length.out = 300)
  data.table(f = f, precision = P, recall = f*P/(2*P - f))
}))[recall <= 1]
pB <- ggplot(dB, aes(precision, recall)) +
  geom_line(data = iso, aes(precision, recall, group = f), inherit.aes = FALSE,
            colour = "grey82", linewidth = 0.22) +
  geom_point(aes(colour = method), shape = 17, size = 2.7) +   # solid triangles, no outline
  geom_label(data = delta, aes(x = 0.98, y = 0.03, label = lab), inherit.aes = FALSE,
             hjust = 1, vjust = 0, size = PT_SMALL / .pt, family = FIG_FAMILY, colour = "grey20",
             fill = "white", linewidth = 0, label.padding = unit(0.5, "pt"), label.r = unit(0, "pt")) +
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
        legend.key.spacing.x = unit(3, "mm"),               # explicit key gap; cairo draws the labels wider than the default spacing allows
        legend.text       = element_text(size = PT_SUP, margin = margin(l = 1.2, r = 1.5, unit = "mm")),   # carries MintCNA^DNA / ^RNA
        legend.key.size   = unit(3.2, "mm"),
        legend.margin     = margin(t = 1),
        plot.margin       = margin(4, 3, 1, 2)) +
  guides(colour = guide_legend(nrow = 1, override.aes = list(size = 3)))
# the left strip is the panel letter: 8 pt bold, top-left of each row, no box (must come after the theme above)
pB <- pB + theme(strip.text.y.left = element_text(size = 8, face = "bold", angle = 0,
                                                  vjust = 1, hjust = 0, margin = margin(t = 0, r = 2)),
                 strip.background.y = element_blank())

# "a" sits on the same vertical line as the b/c/d strips, which ggplot puts inside the
# y-axis title column about 5 mm in from the left edge
fig <- plot_grid(pA, pB, ncol = 1, rel_heights = c(0.25, 0.75), labels = c("a", ""), label_size = 8,
                 label_fontface = "bold", label_fontfamily = FIG_FAMILY, label_x = 5.3 / 130, label_y = 1,
                 hjust = 0, vjust = 1)
W_MM <- 130; H_MM <- 175
mm <- 1 / 25.4
base <- "SuppFig05_RNA_concordance_CN"
cairo_pdf(file.path(OUT, paste0(base, ".pdf")), width = W_MM * mm, height = H_MM * mm, family = FIG_FAMILY)
print(fig); dev.off()
suppressWarnings(pdftools::pdf_convert(file.path(OUT, paste0(base, ".pdf")), format = "png", dpi = 200,
                                       filenames = file.path(QA, paste0(base, "_preview.png")), verbose = FALSE))
cat(sprintf("wrote %s.pdf  (%d x %d mm)\n", base, W_MM, H_MM))
