#!/usr/bin/env Rscript
# Figure 3 -- assembly.  One script, one vector PDF, no panel rescaled.
#
# The original figure was assembled from panels authored at three different
# physical sizes (a 180x150 mm, b 254x102 mm, c 330x127 mm) and then fitted to a
# single page. Whatever scale each panel received, its "13 pt" or "14 pt" text
# arrived on the page as something else, and the three panels disagreed about
# what a point was. Here every panel is drawn at its final size inside one
# device, so the 7 pt in the code is the 7 pt on the page -- inside the 5-7 pt
# Nature Communications in-figure range, at the 183 mm double-column width.
#
# Panel letters are uppercase bold A/B/C, matching the rest of this manuscript's
# main figures. Nature Communications' own style note says lowercase, but internal
# consistency across Figs 1-4 outweighs a soft house rule, and production
# relabels display items anyway.

.libPaths(c("<R_LIBS_USER>", .libPaths()))
FIG4_ASSEMBLING <- TRUE
suppressPackageStartupMessages({ library(ggplot2); library(cowplot); library(grid)
                                 library(svglite); library(ragg) })
CODE <- "<REPRODUCE_DIR>/07_figures/figure4-CRC"

message("[1/3] panel a ..."); source(file.path(CODE, "panel_a_heatmap.R"), local = FALSE)
message("[2/3] panel b ..."); source(file.path(CODE, "panel_b_accuracy.R"), local = FALSE)
message("[3/3] panel c ..."); source(file.path(CODE, "panel_c_tracks.R"), local = FALSE)

g_a <- grid::grid.grabExpr(drawit())
g_b <- panel_b
g_c <- final

# 183 x 170 mm. Width is the Nature Portfolio double-column measure; the 170 mm
# height ceiling is the flagship Nature production spec ("do not exceed 170 mm
# height, leaving room for the legend"). Nature Communications states only the
# widths, but 170 costs nothing here and keeps the figure legal for either.
W_MM <- 180; H_MM <- 170; mm <- 1/25.4          # 93 (row 1) + 75 (row 2) + letters
wa <- 117/183; hr1 <- 93/H_MM

# The composition is BUILT INSIDE the target device. cowplot::draw_plot measures
# text as it assembles, and on R's default null pdf() device that measurement
# fails for these fonts ("failed to find or load PDF CID font"). Opening cairo
# first means the measuring device and the drawing device are the same one.
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
# PDF only. A QA raster can be written by setting FIG4_QA_PNG to a path outside
# the deliverable folder; nothing but the PDF is produced by default.
# Full export bundle -- editable SVG, vector PDF, 600-dpi TIFF -- matching Fig 4.
# Rebuilt inside each device rather than converted, because cowplot measures text
# as it assembles and the measuring device must be the drawing device.
svglite::svglite(paste0(OUTF, ".svg"), width = W_MM*mm, height = H_MM*mm)
print(build()); dev.off()

cairo_pdf(paste0(OUTF, ".pdf"), width = W_MM*mm, height = H_MM*mm, family = FIG_FAMILY)
print(build()); dev.off()

ragg::agg_tiff(paste0(OUTF, ".tiff"), width = W_MM*mm, height = H_MM*mm,
               units = "in", res = 600, compression = "lzw")
print(build()); dev.off()
qa <- Sys.getenv("FIG4_QA_PNG", "")
if (nzchar(qa)) {
  ragg::agg_png(qa, width = W_MM*mm, height = H_MM*mm, units = "in", res = 900)
  print(build()); dev.off(); cat("QA raster:", qa, "\n")
}
cat(sprintf("\nwrote %s .svg/.pdf/.tiff   %d x %d mm   colour mode: %s\n", OUTF, W_MM, H_MM, B_COLOR_MODE))
