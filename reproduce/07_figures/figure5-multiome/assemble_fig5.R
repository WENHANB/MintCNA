#!/usr/bin/env Rscript
# Figure 4 -- assembly.  One script, one vector PDF, no panel rescaled.
#
# LAYOUT, and why it is not Figure 3's.
# Figure 3 is [wide heatmap | one metric] over [full-width track row]: the data
# picture leads, the numbers follow. Figure 4 has a second thing to say -- the
# subclones mean something biologically -- so it is built the other way up:
#
#     row 1 (benchmark, 20 samples)   A  Pearson vs bulk WGS  |  B  precision-recall
#     row 2 (case study, P111T2)      C  per-nucleus CNA heat  |  E  TF motif volcano
#                                     D  clone / state UMAPs   |  F  ATAC vs RNA
#                                                              |  G  CNA dosage
#
# One vertical seam at x = 123 runs the whole page: A/C/D left, B/E/F/G right. The
# reading order is the argument: it works across the cohort (A, B) -> here is what
# it sees in one tumour, and its subclones line up with the RNA cell states (C) ->
# and reproduce that partition in the transcriptome's own embedding (D) -> the
# subclones differ at lineage TF motifs (E) -> that difference reaches the
# transcriptome (F) -> and copy number moves both modalities but explains only
# 1.5% of it (G). Panels C-G are all sample P111T2.
#
# Uppercase bold letters, matching Figures 1-3 of this manuscript. Nature
# Communications' style note says lowercase; internal consistency wins and
# production relabels display items anyway.
#
# Every panel is drawn at its FINAL size inside one device, so the 7 pt in the
# code is the 7 pt on the page. The heatmap is the one panel whose internals are
# in absolute millimetres rather than relative units, so it is re-authored at the
# region size via FIG4A_W/FIG4A_H instead of being stretched into place.
.libPaths(c("<R_LIBS_USER>", .libPaths()))
FIG5_ASSEMBLING <- TRUE
suppressPackageStartupMessages({ library(ggplot2); library(cowplot); library(grid)
                                 library(svglite); library(ragg) })
CODE <- "<MINTCNA_ROOT>/figure4paper/figure5-multiome/code"

# The page grid lives in grid_fig5.R and is shared with qa_panel_geometry.R, so
# the assembly and the geometry check can never disagree about the layout.
source(file.path(CODE, "grid_fig5.R"))
mm <- 1/25.4

# The heatmap is the one panel whose internals are in absolute millimetres, so it
# is re-authored at its grid slot rather than stretched into it.

message("[1/3] panels d,e (motif) ..."); source(file.path(CODE, "panel_bio.R"),      local = FALSE)
message("[2/3] panels a,b (accuracy) ..."); source(file.path(CODE, "panel_b_accuracy.R"), local = FALSE)
message("[2b]  panel e (RNA DE) ...");        source(file.path(CODE, "panel_de.R"),         local = FALSE)
message("[3/4] panel c (heatmap) ...")
Sys.setenv(FIG5_SAMPLE = "P111T2", FIG4A_W = LEFT_W, FIG4A_H = C_HM_H, FIG4C_LEGEND = "TRUE")
source(file.path(CODE, "panel_a_heatmap.R"), local = FALSE)
g_c <- grid::grid.grabExpr(drawit(), width = LEFT_W*mm, height = C_HM_H*mm)

message("[4/4] panel d (umap pair) ...")
Sys.setenv(FIG5_UMAP_SAMPLE = "P111T2")
source(file.path(CODE, "panel_umap.R"), local = FALSE)
# the two UMAPs share a row, aligned so their square panels match exactly
g_umap <- plot_grid(plotlist = align_plots(p_umap_clone, p_umap_state,
                                           align = "h", axis = "tb"),
                    ncol = 2, rel_widths = c(1, 1))

# Snapping the DRAW BOXES to a grid is not enough: B, D and E carry y titles of
# very different lengths ("Recall", "-log10 FDR", "P111T2 OPC - AC"), so their
# plotting areas still started up to 9 mm apart inside boxes that shared an edge.
# align_plots pads every y-axis strip to the widest one, so the three panels in
# the right column share a plot-area left edge and a right edge, not just an
# outer box. Returned as gtables, hence draw_grob below.
# axis = "tblr", not "lr": padding only left/right left the four plot areas
# 0.71 mm apart in height, because the x-axis titles are not all one line high.
# align_plots pads the AXIS strips to a common size, but not the axis TITLE rows,
# and E and F carry plotmath subscripts in their x titles (log[2], Delta) that B
# and D do not. That left the four plot areas 0.71 mm apart in height, measured on
# the exported page. equalize_fixed() finishes the job: every gtable row and
# column that is not the panel's null unit is set to the largest such value across
# the four, so all four gtables have identical non-panel geometry and therefore
# identical panel rectangles when drawn into identical boxes.
# WIDTHS and HEIGHTS are equalised over DIFFERENT SETS, and that is deliberate.
# Widths: all five plots, so every plot area in the right column shares a left and
# a right edge. Heights: only B, E and F, so their non-panel geometry matches and
# the gaps between them are uniform. G's two rows must NOT join the height group --
# they are a stacked pair that has to sit closer together than the column gap, and
# the upper one has no x axis at all, so forcing it to carry the maximum x-axis
# height would open a 7.7 mm hole in the middle of one panel.
.eq <- function(gs, what) {
  ish <- function(u) as.character(grid::unitType(u)) == "null"
  if (what == "w") {
    vs <- lapply(gs, function(g) grid::convertWidth(g$widths, "mm", valueOnly = TRUE))
  } else {
    vs <- lapply(gs, function(g) grid::convertHeight(g$heights, "mm", valueOnly = TRUE))
  }
  # pmax over vectors of different lengths recycles silently and lands fixed sizes
  # on panel rows; a faceted plot in this list collapsed the panel to a sliver.
  stopifnot(length(unique(lengths(vs))) == 1L)
  nul <- ish(if (what == "w") gs[[1]]$widths else gs[[1]]$heights)
  mx  <- do.call(pmax, vs)
  for (i in seq_along(gs)) {
    if (what == "w") gs[[i]]$widths [!nul] <- grid::unit(mx[!nul], "mm")
    else             gs[[i]]$heights[!nul] <- grid::unit(mx[!nul], "mm")
  }
  gs
}
# The TF-motif volcano (p_vol) left the main figure on 2026-09-05 and is
# Supplementary. Its slot holds the cohort panel: without it every biology panel
# was one tumour chosen for showing the effect, and the figure could not carry the
# claim its caption makes. See panel_cohort.R.
col_r <- align_plots(p_bot, p_de, p_rep, align = "v", axis = "tblr")
col_r <- .eq(col_r, "w")
col_r <- .eq(col_r, "h")
col_grid <- plot_grid(plotlist = col_r, ncol = 1,
                      rel_heights = c(B_H, E_H, F_H))

fx <- function(x) x/PAGE_W; fy <- function(y) y/PAGE_H   # mm -> npc on this page

build <- function() {
  # every x/y/width/height comes from grid_fig5.R and NOTHING is nudged. The
  # earlier "- 0.006" / "- 0.012" npc corrections were eyeball fudges from an
  # earlier layout, and -0.006 npc is 1.02 mm on this page -- exactly the
  # top-edge mismatch that qa_svg_geometry.py measured between A and B.
  p <- ggdraw() +
    draw_plot(p_top,      x = fx(LEFT_X),  y = fy(A_Y),   width = fx(LEFT_W),  height = fy(A_H)) +
    draw_plot(col_grid,   x = fx(RIGHT_X), y = 0,       width = fx(RIGHT_W), height = 1) +
    draw_grob(g_c,        x = fx(LEFT_X),  y = fy(C_HM_Y), width = fx(LEFT_W), height = fy(C_HM_H)) +
    draw_plot(g_umap,     x = fx(LEFT_X),  y = fy(C_UM_Y), width = fx(LEFT_W), height = fy(C_UM_H))
    # Letters run down the LEFT column (A, C, D), then down the right (B, E, F, G).
    # Each sits LTR_DROP mm below the top of its own region -- in MILLIMETRES, so
    # the inset does not silently change when the page height does.
    ltr <- function(l, x, ytop_mm)
      draw_label(l, x = fx(x), y = fy(ytop_mm - LTR_DROP), fontface = "bold",
                 size = 8, fontfamily = FIG_FAMILY, hjust = 0, vjust = 1)   # 8 pt: Nature figure guide (2026-10-09)
    p + ltr("a", LEFT_X,  PAGE_H) + ltr("b", RIGHT_X, PAGE_H) +
        ltr("c", LEFT_X,  C_H)    + ltr("d", LEFT_X,  C_UM_H) +
        ltr("e", RIGHT_X, E_Y + E_H) + ltr("f", RIGHT_X, F_Y + F_H)
}

OUTF <- file.path(OUT, "fig5_multiome")

# Full export bundle: editable SVG, vector PDF, 600-dpi TIFF. The figure is
# rebuilt inside each device rather than drawn once and converted, because
# cowplot measures text as it assembles and the measuring device must be the
# drawing device.
svglite::svglite(paste0(OUTF, ".svg"), width = PAGE_W*mm, height = PAGE_H*mm)
print(build()); dev.off()

cairo_pdf(paste0(OUTF, ".pdf"), width = PAGE_W*mm, height = PAGE_H*mm, family = FIG_FAMILY)
print(build()); dev.off()

ragg::agg_tiff(paste0(OUTF, ".tiff"), width = PAGE_W*mm, height = PAGE_H*mm,
               units = "in", res = 600, compression = "lzw")
print(build()); dev.off()
qa <- Sys.getenv("FIG5_QA_PNG", "")
if (nzchar(qa)) {
  ragg::agg_png(qa, width = PAGE_W*mm, height = PAGE_H*mm, units = "in", res = 700)
  print(build()); dev.off(); cat("QA raster:", qa, "\n")
}
# %g not %d: PAGE_H is 163.5 and %d errors on a non-integer -- after the files
# were written, so the figure looked fine while the exit code broke the SI chain.
cat(sprintf("\nwrote %s .svg/.pdf/.tiff   %g x %g mm\n", OUTF, PAGE_W, PAGE_H))
