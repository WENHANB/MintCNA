#!/usr/bin/env Rscript
# Ad-hoc figure for the advisor (2026-09-13): main-figure panel F with the x axis
# swapped from promoter accessibility to COPY NUMBER -- per gene, the difference
# in copy-number state between the two subclones of P111T2 (C2 - C1, ATAC-only
# arm, so expression plays no part in a gene's copy-number assignment) against
# the expression shift log2(C2/C1). Same per-gene table as Supplementary Fig. S7B/C.
#
# Files:
#   extra_cnv_vs_expression_P111T2.pdf        copy number vs expression
#   extra_cnv_vs_accessibility_P111T2.pdf     copy number vs promoter accessibility
#   extra_cnv_vs_expression_side_by_side.pdf  A Fig 5F (acc vs expr), B CNV vs expr,
#                                             C CNV vs acc
.libPaths(c("<R_LIBS_USER>", .libPaths()))
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(cowplot) })
ROOT <- "<MINTCNA_ROOT>"
FIG5_ASSEMBLING <- TRUE
source(file.path(ROOT, "figure4paper/figure5-multiome/code/panel_bio.R"), local = FALSE)
OUT <- file.path(ROOT, "figure4paper/figure5-multiome/output")

cn <- merge(dz[, .(gene, d_cn, d_atac, d_rna)], Lg[, .(gene, fam)], by = "gene", all.x = TRUE)
cn[is.na(fam), fam := "other"]
rho_cn <- cor(cn$d_cn, cn$d_rna, method = "spearman")
r_cn   <- cor(cn$d_cn, cn$d_rna)
cat(sprintf("%d genes; Spearman rho(dCN, dExpr) = %.3f, Pearson r = %.3f, R^2 = %.3f; for comparison rho(dAcc, dExpr) = %.3f\n",
            nrow(cn), rho_cn, r_cn, r_cn^2, rho_ar))
# binned dose response: mean expression shift per copy-number-difference bin
BRK2 <- c(-Inf, -0.4, -0.2, 0.2, 0.4, Inf)
cn[, bin := cut(d_cn, BRK2)]
bm <- cn[, .(x = mean(d_cn), y = mean(d_rna), n = .N), by = bin][order(x)]
print(bm)

# One builder for both y axes. yvar = "d_rna" (expression) or "d_atac" (promoter
# accessibility). NOTE for the accessibility version: promoter accessibility and the
# copy-number call are both ATAC coverage, so that panel is an internal consistency
# check, not independent evidence -- the caption must say so.
LABS <- c("PDGFRA", "OLIG1", "DLL3", "CSPG4", "SOX10", "GFAP", "SLC1A3", "ID3")
cnv_panel <- function(yvar, ylab) {
  d <- copy(cn)[, y := get(yvar)]
  rho <- cor(d$d_cn, d$y, method = "spearman")
  bm  <- d[, .(x = mean(d_cn), y = mean(y), n = .N), by = bin][order(x)]
  cat(sprintf("  %s: Spearman rho(dCN, %s) = %.3f\n", yvar, yvar, rho)); print(bm)
  lab <- d[gene %in% LABS]
  XR <- range(d$d_cn) + c(-1, 1) * 0.10 * diff(range(d$d_cn))
  YR <- range(d$y) + c(-1, 1) * 0.05 * diff(range(d$y))
  LB <- place_labels(d$d_cn, d$y, d$fam != "other", match(lab$gene, d$gene), lab$gene,
                     XR, YR, c(47, 47), pt_mm = fifelse(d$fam != "other", 1.0, 0.4), text_pt = FS_LGT,
                     keepout = c(0, 13, 47 - 6, 47))
  LB <- merge(as.data.table(LB), d[, .(gene, fam)], by = "gene")
  ggplot(d, aes(d_cn, y)) +
    geom_hline(yintercept = 0, linewidth = 0.22, colour = "grey75") +
    geom_vline(xintercept = 0, linewidth = 0.22, colour = "grey75") +
    geom_point(data = d[fam == "other"], colour = "grey82", size = 0.35, alpha = 0.45) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.3, linetype = "22", colour = "grey45") +
    geom_line(data = bm, aes(x, y), colour = "grey20", linewidth = 0.35) +
    geom_point(data = bm, aes(x, y), colour = "grey20", shape = 18, size = 2.2) +
    geom_segment(data = LB[seg == TRUE], aes(x = ax, y = ay, xend = px, yend = py), linewidth = 0.18, colour = "grey60") +
    geom_point(data = d[fam != "other"], aes(colour = fam), size = 1.15) +
    geom_text(data = LB, aes(lx, ly, label = gene, colour = fam), size = FS_LGT/.pt, family = FIG_FAMILY, show.legend = FALSE) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.18, vjust = 1.5, size = PT_SUP/.pt, fontface = "bold",
             colour = "grey15", parse = TRUE, label = sprintf("rho == '%.2f'", rho)) +
    scale_colour_manual(values = c(`OPC lineage` = STATE_COLORS[["OPC-like"]],
                                   `AC lineage`  = STATE_COLORS[["AC-like"]]),
                        breaks = c("OPC lineage", "AC lineage"),
                        labels = c("expected in C2 (OPC lineage)", "expected in C1 (AC lineage)"),
                        name = NULL) +
    scale_x_continuous(expand = expansion(mult = 0.10)) +
    scale_y_continuous(expand = expansion(mult = 0.05)) +
    labs(x = "Copy-number difference  C2 \u2212 C1 (copies)", y = ylab) +
    theme_cowplot(BASE_PT, font_family = FIG_FAMILY) +
    theme(axis.text = element_text(size = BASE_PT - 1, colour = "grey20"),
          axis.title = element_text(size = PT_SUP),
          axis.line = element_line(linewidth = 0.3, colour = "grey30"),
          axis.ticks = element_line(linewidth = 0.3, colour = "grey30"),
          panel.grid = element_blank(), aspect.ratio = 1,
          legend.position = "inside", legend.position.inside = c(0.98, 0.02),
          legend.justification = c(1, 0), legend.background = element_blank(),
          legend.text = element_text(size = BASE_PT - 1, colour = "grey20"),
          legend.key.size = unit(3, "mm"), plot.margin = margin(2, 2, 1, 2))
}
p_cn  <- cnv_panel("d_rna",  expression("Expression  " * log[2] * " (C2/C1)"))
p_cna <- cnv_panel("d_atac", expression("Promoter accessibility  " * log[2] * " (C2/C1)"))

MM <- 1/25.4
cairo_pdf(file.path(OUT, "extra_cnv_vs_expression_P111T2.pdf"), width = 72*MM, height = 72*MM, family = FIG_FAMILY)
print(p_cn); dev.off()
ragg::agg_png(file.path(OUT, "extra_cnv_vs_expression_P111T2.png"), width = 72*MM, height = 72*MM, units = "in", res = 600)
print(p_cn); dev.off()

cairo_pdf(file.path(OUT, "extra_cnv_vs_accessibility_P111T2.pdf"), width = 72*MM, height = 72*MM, family = FIG_FAMILY)
print(p_cna); dev.off()
ragg::agg_png(file.path(OUT, "extra_cnv_vs_accessibility_P111T2.png"), width = 72*MM, height = 72*MM, units = "in", res = 600)
print(p_cna); dev.off()
side <- plot_grid(p_rep + theme(legend.position = "none", plot.margin = margin(7, 2, 1, 2)),
                  p_cn  + theme(legend.position = "none", plot.margin = margin(7, 2, 1, 2)),
                  p_cna + theme(plot.margin = margin(7, 2, 1, 2)),
                  ncol = 3, align = "h", axis = "tb", labels = c("a", "b", "c"), label_size = 8,
                  label_fontfamily = FIG_FAMILY, label_fontface = "bold", label_x = 0.005, label_y = 1.0, vjust = 1.1)
cairo_pdf(file.path(OUT, "extra_cnv_vs_expression_side_by_side.pdf"), width = 222*MM, height = 76*MM, family = FIG_FAMILY)
print(side); dev.off()
ragg::agg_png(file.path(OUT, "extra_cnv_vs_expression_side_by_side.png"), width = 222*MM, height = 76*MM, units = "in", res = 600)
print(side); dev.off()
fwrite(cn[, .(gene, d_cn, d_atac, d_rna, fam)], file.path(ROOT, "figure4paper/figure5-multiome/data/extra_cnv_vs_expression_P111T2.csv"))
cat("wrote extra_cnv_vs_expression_P111T2, extra_cnv_vs_accessibility_P111T2 and _side_by_side (A acc-expr, B CNV-expr, C CNV-acc)\n")
