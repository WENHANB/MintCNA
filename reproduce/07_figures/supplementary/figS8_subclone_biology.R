#!/usr/bin/env Rscript
# Supplementary Fig. 8 (S7 until 2026-10-09) -- what the two subclones of P111T2 differ in, beyond
# the gene-level tests of main Figure 5E-F.
#
#   A  TF-motif accessibility volcano, 12 prespecified motifs (was main-figure E)
#   B  per-gene dosage: what copy number predicts vs what is observed, 21 lineage
#      genes (was main-figure G)
#   C  Fig 5F restricted to copy-neutral genes (|dCN| < 0.1): the accessibility-
#      expression correlation does not change (rho 0.513 vs 0.508), so the lineage
#      difference is not carried by gene dosage
#
# A and B are the panels panel_bio.R already builds (p_vol, p_dose); C is drawn
# here from the same per-gene table (dz) that B uses. The main-figure script is
# sourced under FIG5_ASSEMBLING so it builds its ggplot objects and skips its own
# standalone exports.

.libPaths(c("<R_LIBS_USER>", .libPaths()))
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(cowplot) })
ROOT <- "<MINTCNA_ROOT>"
FIG5_ASSEMBLING <- TRUE
source(file.path(ROOT, "figure4paper/figure5-multiome/code/panel_bio.R"), local = FALSE)
# AFTER the source: panel_bio.R defines its own OUT and DATA, and on the first
# build this figure was written into the main-figure folder because of it
OUT  <- file.path(ROOT, "figure4paper/figureS-supplementary/_work/output")
PDF_OUT <- file.path(ROOT, "figure4paper/figureS-supplementary")   # final PDFs live at the top level; everything else under _work/
QA   <- file.path(ROOT, "figure4paper/figureS-supplementary/_work/qa")
DATA <- file.path(ROOT, "figure4paper/figureS-supplementary/_work/data")

## ---- panel C: copy-neutral genes only ---------------------------------------------
NEU <- 0.1
neu <- dz[abs(d_cn) < NEU]
rho_neu <- cor(neu$d_atac, neu$d_rna, method = "spearman")
conc_neu <- 100 * neu[abs(d_rna) > 1, mean(sign(d_atac) == sign(d_rna))]
cat(sprintf("panel C: %d of %d genes copy-neutral (|dCN| < %.1f); Spearman rho = %.3f (all genes %.3f); %.1f%% concordant among |log2FC| > 1\n",
            nrow(neu), nrow(dz), NEU, rho_neu, rho_ar, conc_neu))
# 6,317 on the ATAC arm with strand-aware TSS windows (2026-09-14; 4,538 before the
# strand fix). The 5,399 / rho 0.505 quoted in earlier drafts came from
# the multi-arm dCN (Result/qc/cna_dosage_link_P111T2.tsv); B reads the ATAC arm on
# purpose (RNA must not help choose a gene's copy-number bin), so C does too, and
# the text follows the figure.
stopifnot(nrow(neu) == 6317)                          # strand-aware TSS windows, 2026-09-14
neu <- merge(neu, Lg[, .(gene, fam)], by = "gene", all.x = TRUE)
neu[is.na(fam), fam := "other"]
fwrite(neu[, .(gene, d_cn, d_atac, d_rna, fam)], file.path(DATA, "figS8c_copy_neutral_P111T2.csv"))

p_neu <- ggplot(neu, aes(d_atac, d_rna)) +
  geom_hline(yintercept = 0, linewidth = 0.22, colour = "grey75") +
  geom_vline(xintercept = 0, linewidth = 0.22, colour = "grey75") +
  geom_point(data = neu[fam == "other"], aes(colour = fam), size = 0.35, alpha = 0.45) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.3, linetype = "22", colour = "grey45") +
  geom_point(data = neu[fam != "other"], aes(colour = fam), size = 1.15) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.15, vjust = 1.5, size = PT_SUP/.pt, fontface = "bold",
           colour = "grey15", parse = TRUE, label = sprintf("rho == %.3f", rho_neu)) +
  # the gene count (6,317) is stated in the legend, not on the panel (author, 2026-10-11)
  scale_colour_manual(values = c(`OPC lineage` = STATE_COLORS[["OPC-like"]],
                                 `AC lineage`  = STATE_COLORS[["AC-like"]], other = "grey55"),
                      breaks = c("OPC lineage", "AC lineage", "other"),
                      labels = c("lineage gene, expected in C2 (OPC)", "lineage gene, expected in C1 (AC)", "all other genes"),
                      name = NULL,
                      guide = guide_legend(override.aes = list(size = c(1.6, 1.6, 1.0), alpha = 1))) +
  scale_x_continuous(breaks = seq(-6, 6, 2), expand = expansion(mult = 0.10)) +   # integer ticks, as Fig 5F
  scale_y_continuous(expand = expansion(mult = 0.05)) +
  labs(x = expression("Promoter accessibility  " * log[2] * " (C2/C1)"),
       y = expression("Expression  " * log[2] * " (C2/C1)")) +
  theme_cowplot(BASE_PT, font_family = FIG_FAMILY) +
  theme(axis.text = element_text(size = BASE_PT - 1, colour = "grey20"),
        axis.title = element_text(size = PT_SUP),
        axis.line = element_line(linewidth = 0.3, colour = "grey30"),
        axis.ticks = element_line(linewidth = 0.3, colour = "grey30"),
        panel.grid = element_blank(), aspect.ratio = 1,
        legend.position = "inside", legend.position.inside = c(0.98, 0.02),
        legend.justification = c(1, 0), legend.direction = "vertical",
        legend.background = element_blank(), legend.key = element_blank(),
        legend.text = element_text(size = BASE_PT - 1, colour = "grey20"),
        legend.key.size = unit(3, "mm"), legend.margin = margin(0, 0, 0, 0),
        plot.margin = margin(2, 2, 1, 2))

## ---- assemble ----------------------------------------------------------------------
W_MM <- 165; H_MM <- 70;  MM <- 1/25.4     # 165 = SI page text width (was 180)
pa <- p_vol  + theme(plot.margin = margin(7, 2, 1, 2))
pb <- p_dose + theme(plot.margin = margin(7, 2, 1, 2))
pc <- p_neu  + theme(plot.margin = margin(7, 2, 1, 2))
fig <- plot_grid(pa, pb, pc, ncol = 3, rel_widths = c(66, 50, 67), align = "h", axis = "b",
                 labels = c("a", "b", "c"), label_size = 8, label_fontfamily = FIG_FAMILY,
                 label_fontface = "bold", label_x = 0.005, label_y = 1.0, vjust = 1.1)
cairo_pdf(file.path(PDF_OUT, "SuppFig08_subclone_biology.pdf"), width = W_MM*MM, height = H_MM*MM, family = FIG_FAMILY)
print(fig); dev.off()
svglite::svglite(file.path(OUT, "figS8_subclone_biology.svg"), width = W_MM*MM, height = H_MM*MM)
print(fig); dev.off()
ragg::agg_png(file.path(QA, "figS8_subclone_biology_preview.png"), width = W_MM*MM, height = H_MM*MM, units = "in", res = 900)
print(fig); dev.off()
cat(sprintf("\nwrote figS8_subclone_biology.pdf/.svg  %d x %d mm\n", W_MM, H_MM))
