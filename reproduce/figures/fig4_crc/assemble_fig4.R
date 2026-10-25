#!/usr/bin/env Rscript
# Figure 4 (CRC01 scTrio-seq2 benchmark): assembles panel a (single-cell copy-number heatmap),
# panel b (per-cell accuracy against the bulk WGS profile) and panel c (pseudobulk CNA tracks) into
# one 180 x 170 mm figure. Every panel is drawn at its final physical size inside one device, so
# the point sizes in the code are the point sizes on the page.
# Input : the panel scripts fig4_crc/panel_a_heatmap.R, panel_b_accuracy.R and panel_c_tracks.R,
#         which read reproduce/figures/data/fig4_crc/ (see fig4_crc/utils.R)
# Output: reproduce/figures/output/fig4/fig4_CRC.svg, .pdf and .tiff (600 dpi)
# Run from the root of the MintCNA repository: Rscript reproduce/figures/fig4_crc/assemble_fig4.R

FIG4_ASSEMBLING <- TRUE   # the panel scripts skip their standalone export when this exists
suppressPackageStartupMessages({ library(ggplot2); library(cowplot); library(grid)
                                 library(svglite); library(ragg) })
CODE <- "reproduce/figures/fig4_crc"

message("[1/3] panel a ..."); source(file.path(CODE, "panel_a_heatmap.R"), local = FALSE)
message("[2/3] panel b ..."); source(file.path(CODE, "panel_b_accuracy.R"), local = FALSE)
message("[3/3] panel c ..."); source(file.path(CODE, "panel_c_tracks.R"), local = FALSE)

g_a <- grid::grid.grabExpr(drawit())
g_b <- panel_b
g_c <- final

# 180 mm is the double-column width. The first row (a, b) is 93 mm high with panel a taking
# 117/183 of the width; the second row (c) fills the rest.
W_MM <- 180; H_MM <- 170; mm <- 1/25.4
wa <- 117/183; hr1 <- 93/H_MM

# The composition is built inside the output device: cowplot measures text as it assembles, and
# on R's default null pdf() device that measurement fails for these fonts. Each format is
# therefore drawn afresh in its own device rather than converted from another.
build <- function() {
  ggdraw() +
    draw_grob(g_a, x = 0.005, y = 1 - hr1, width = wa - 0.01,      height = hr1 - 0.004) +
    draw_plot(g_b, x = wa,    y = 1 - hr1, width = 1 - wa - 0.004, height = hr1 - 0.004) +
    draw_plot(g_c, x = 0.004, y = 0.002,  width = 0.994,           height = 1 - hr1 - 0.004) +
    draw_label("a", x = 0.004,     y = 0.997,          fontface = "bold", size = 8,
               fontfamily = FIG_FAMILY, hjust = 0, vjust = 1) +
    draw_label("b", x = wa + 0.004, y = 0.997,         fontface = "bold", size = 8,
               fontfamily = FIG_FAMILY, hjust = 0, vjust = 1) +
    draw_label("c", x = 0.004,     y = 1 - hr1 + 0.004, fontface = "bold", size = 8,
               fontfamily = FIG_FAMILY, hjust = 0, vjust = 1)
}

OUTF <- file.path(OUT, "fig4_CRC")
svglite::svglite(paste0(OUTF, ".svg"), width = W_MM*mm, height = H_MM*mm)
print(build()); dev.off()

cairo_pdf(paste0(OUTF, ".pdf"), width = W_MM*mm, height = H_MM*mm, family = FIG_FAMILY)
print(build()); dev.off()

ragg::agg_tiff(paste0(OUTF, ".tiff"), width = W_MM*mm, height = H_MM*mm,
               units = "in", res = 600, compression = "lzw")
print(build()); dev.off()
cat(sprintf("\nwrote %s .svg/.pdf/.tiff   %d x %d mm\n", OUTF, W_MM, H_MM))
