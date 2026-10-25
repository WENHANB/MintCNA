#!/usr/bin/env Rscript
# Figure 5 (glioma single-nucleus multiome benchmark): assembles the panels into one page.
# Left column: a Pearson r vs bulk WGS (p_top), c per-nucleus CNA heatmap of P111T2 (drawit),
# d RNA UMAP pair (p_umap_clone, p_umap_state). Right column: b precision-recall plane (p_bot),
# e subclone differential expression (p_de), f promoter accessibility vs expression (p_rep).
# The panel scripts are sourced with FIG5_ASSEMBLING defined, so they build their plot objects
# and skip their own exports; the page grid comes from grid_fig5.R. Panels c-f are sample P111T2.
# Output: reproduce/figures/output/fig5/fig5_multiome.svg, .pdf and .tiff (180 x 163.5 mm).
# Run from the root of the MintCNA repository: Rscript reproduce/figures/fig5_multiome/assemble_fig5.R
FIG5_ASSEMBLING <- TRUE
suppressPackageStartupMessages({ library(ggplot2); library(cowplot); library(grid)
                                 library(svglite); library(ragg) })
CODE <- "reproduce/figures/fig5_multiome"   # run from the root of the MintCNA repository

source(file.path(CODE, "grid_fig5.R"))     # page geometry in mm
mm <- 1/25.4

# Every panel is drawn at its final size inside one device, so the point sizes in the panel
# scripts are the sizes on the page. The heatmap's internals are in absolute millimetres, so
# panel_a_heatmap.R takes its size from the grid (LEFT_W x C_HM_H) when sourced here and the
# drawing is captured at exactly that size instead of being rescaled.
message("[1/5] panel f and the Supplementary Fig. 8 panels (panel_bio.R) ..."); source(file.path(CODE, "panel_bio.R"),      local = FALSE)
message("[2/5] panels a,b (accuracy) ..."); source(file.path(CODE, "panel_b_accuracy.R"), local = FALSE)
message("[3/5] panel e (RNA DE) ...");        source(file.path(CODE, "panel_de.R"),         local = FALSE)
message("[4/5] panel c (heatmap) ...")
source(file.path(CODE, "panel_a_heatmap.R"), local = FALSE)
g_c <- grid::grid.grabExpr(drawit(), width = LEFT_W*mm, height = C_HM_H*mm)

message("[5/5] panel d (UMAP pair) ...")
source(file.path(CODE, "panel_umap.R"), local = FALSE)
# the two UMAPs share a row, aligned so their square panels match exactly
g_umap <- plot_grid(plotlist = align_plots(p_umap_clone, p_umap_state,
                                           align = "h", axis = "tb"),
                    ncol = 2, rel_widths = c(1, 1))

# Right column. align_plots pads the axis strips of the three plots to a common size but not the
# axis-title rows, and the x titles of e and f carry plotmath while b's does not, so the plot areas
# would still differ by a fraction of a millimetre. .eq() sets every fixed (non-null) gtable row and
# column to the largest value across the plots, so the three gtables have identical non-panel
# geometry and therefore identical panel rectangles when drawn into equal boxes.
.eq <- function(gs, what) {
  ish <- function(u) as.character(grid::unitType(u)) == "null"
  if (what == "w") {
    vs <- lapply(gs, function(g) grid::convertWidth(g$widths, "mm", valueOnly = TRUE))
  } else {
    vs <- lapply(gs, function(g) grid::convertHeight(g$heights, "mm", valueOnly = TRUE))
  }
  # pmax over vectors of different lengths would recycle silently, so the gtables must agree in size
  stopifnot(length(unique(lengths(vs))) == 1L)
  nul <- ish(if (what == "w") gs[[1]]$widths else gs[[1]]$heights)
  mx  <- do.call(pmax, vs)
  for (i in seq_along(gs)) {
    if (what == "w") gs[[i]]$widths [!nul] <- grid::unit(mx[!nul], "mm")
    else             gs[[i]]$heights[!nul] <- grid::unit(mx[!nul], "mm")
  }
  gs
}
col_r <- align_plots(p_bot, p_de, p_rep, align = "v", axis = "tblr")
col_r <- .eq(col_r, "w")
col_r <- .eq(col_r, "h")
col_grid <- plot_grid(plotlist = col_r, ncol = 1,
                      rel_heights = c(B_H, E_H, F_H))

fx <- function(x) x/PAGE_W; fy <- function(y) y/PAGE_H   # mm -> npc on this page

build <- function() {
  # every position and size comes from grid_fig5.R; nothing is nudged by hand
  p <- ggdraw() +
    draw_plot(p_top,      x = fx(LEFT_X),  y = fy(A_Y),   width = fx(LEFT_W),  height = fy(A_H)) +
    draw_plot(col_grid,   x = fx(RIGHT_X), y = 0,       width = fx(RIGHT_W), height = 1) +
    draw_grob(g_c,        x = fx(LEFT_X),  y = fy(C_HM_Y), width = fx(LEFT_W), height = fy(C_HM_H)) +
    draw_plot(g_umap,     x = fx(LEFT_X),  y = fy(C_UM_Y), width = fx(LEFT_W), height = fy(C_UM_H))
    # Panel letters run down the left column (a, c, d), then down the right (b, e, f); each sits
    # LTR_DROP mm below the top of its own region.
    ltr <- function(l, x, ytop_mm)
      draw_label(l, x = fx(x), y = fy(ytop_mm - LTR_DROP), fontface = "bold",
                 size = 8, fontfamily = FIG_FAMILY, hjust = 0, vjust = 1)   # 8 pt bold panel letters
    p + ltr("a", LEFT_X,  PAGE_H) + ltr("b", RIGHT_X, PAGE_H) +
        ltr("c", LEFT_X,  C_H)    + ltr("d", LEFT_X,  C_UM_H) +
        ltr("e", RIGHT_X, E_Y + E_H) + ltr("f", RIGHT_X, F_Y + F_H)
}

OUTF <- file.path(OUT, "fig5_multiome")

# Editable SVG, vector PDF and 600-dpi TIFF. The figure is rebuilt inside each device rather than
# drawn once and converted, because cowplot measures text as it assembles and the measuring device
# must be the drawing device.
svglite::svglite(paste0(OUTF, ".svg"), width = PAGE_W*mm, height = PAGE_H*mm)
print(build()); dev.off()

cairo_pdf(paste0(OUTF, ".pdf"), width = PAGE_W*mm, height = PAGE_H*mm, family = FIG_FAMILY)
print(build()); dev.off()

ragg::agg_tiff(paste0(OUTF, ".tiff"), width = PAGE_W*mm, height = PAGE_H*mm,
               units = "in", res = 600, compression = "lzw")
print(build()); dev.off()
cat(sprintf("\nwrote %s .svg/.pdf/.tiff   %g x %g mm\n", OUTF, PAGE_W, PAGE_H))   # %g: PAGE_H is not an integer
