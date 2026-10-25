#!/usr/bin/env Rscript
# Page grid of Figure 5, in mm, sourced by fig5_multiome/assemble_fig5.R; reads no data and writes nothing.
# The grid is defined on plot areas rather than on panel boxes: a panel with rotated tick labels spends
# more of its box on the axis strip than one with horizontal labels, so equal boxes would put the data
# rectangles at different heights on the page.
#   left column 120 mm: a Pearson vs bulk WGS, c per-nucleus CNA heatmap, d UMAP pair
#   right column 60 mm: b precision-recall, e differential expression, f ATAC vs RNA (all three square)
# Run from the root of the MintCNA repository through fig5_multiome/assemble_fig5.R.

PAGE_W <- 180      # double-column page width
PAGE_H <- 163.5    # three square right-column panels, 3 x SQ_H

MARG    <- 1                                        # outer margin
GUT     <- 2                                        # gutter between columns
# The right column is narrow because b and f are square and take their size from the smaller box
# dimension; a wider column would only add empty width beside them.
LEFT_W  <- 120
RIGHT_W <- PAGE_W - MARG - LEFT_W - GUT - MARG      # 60
LEFT_X  <- MARG
RIGHT_X <- MARG + LEFT_W + GUT                      # 123

# A square panel is capped by the plot width available in the column (46.90 mm); its box height is
# that cap, less a 0.1 mm margin, plus the 7.70 mm x-axis strip.
SQ_H <- 54.5      # square box: width cap 46.90 - 0.1 + 7.70 x strip
B_H <- SQ_H       # precision-recall   (square)
E_H <- SQ_H       # differential expression   (square: both axes are unit-free scores)
F_H <- SQ_H       # ATAC vs RNA        (square)
stopifnot(abs(B_H + E_H + F_H - PAGE_H) < 1e-9)

# a's box is taller than b's by exactly the height its 45-degree rotated method labels add over b's
# horizontal ones, so the two plot areas share the same top and bottom edges.
A_EXTRA <- 7.285
A_H <- B_H + A_EXTRA
C_H <- PAGE_H - A_H
# The left column below a is the heatmap over the UMAP pair. The two UMAPs are square, each with
# half the column, so they are height-limited by C_UM_H.
C_HM_H  <- 60                     # heatmap, carries the shared CN + state legend
C_UM_H  <- C_H - C_HM_H           # the UMAP pair

# Panel letters are inset from the top-left of their region by this much, in mm rather than npc,
# so the inset does not change with the page size.
LTR_DROP <- 0.6

# y origins (bottom-left, as cowplot wants them)
A_Y    <- PAGE_H - A_H
C_HM_Y <- C_UM_H
C_UM_Y <- 0
F_Y <- 0
E_Y <- F_H
B_Y <- F_H + E_H              # == PAGE_H - B_H
stopifnot(abs(B_Y - (PAGE_H - B_H)) < 1e-9)
