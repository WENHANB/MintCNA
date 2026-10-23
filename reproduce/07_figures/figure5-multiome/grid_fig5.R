#!/usr/bin/env Rscript
# Figure 4's page grid, in one place, sourced by BOTH the assembly and the
# geometry QA script so the two can never disagree about what the layout is.
#
# The grid is defined on PLOT AREAS, not on boxes. Equal boxes were the mistake
# that made this figure look uneven through several rounds: a panel with
# 45-degree rotated tick labels spends 14.3 mm of its box on those labels while a
# panel with horizontal labels spends 6.8 mm, so two panels sharing a box height
# end up with plot areas 7.5 mm apart and sitting at different heights on the
# page. The eye compares the data rectangles, not the boxes.
#
#   left column 120 mm        right column 60 mm
#   A  Pearson vs bulk WGS    B  precision-recall   (SQUARE)
#   C  per-nucleus CNA heat   E  TF motif volcano
#   D  clone / state UMAPs    F  ATAC vs RNA        (SQUARE)
#                             G  CNA dosage
#
# A's box is 7.285 mm taller than B's, which is exactly the extra height A's
# rotated method labels consume, so the two plot areas land on the same top and
# the same bottom. The cost is that the A/C seam does not line up with the B/E
# seam; a seam offset is far less visible than two data rectangles of different
# height side by side, and qa_svg_geometry.py asserts the plot areas, not boxes.

PAGE_W <- 180      # Nature Communications double-column measure (183 until 2026-10-11)
PAGE_H <- 163.5    # three square right-column panels: 3 x SQ_H. Fig 3 is 170.

MARG    <- 1                                        # outer margin
GUT     <- 2                                        # gutter between columns
# The right column is NARROW, and that is deliberate. B and F are square, so each
# one's size is its SMALLER dimension. At a 73 mm column their squares were
# pinned by the box HEIGHT while 26 mm of column width sat empty beside them.
# Widening the column would only have grown that empty band.
LEFT_W  <- 120
RIGHT_W <- PAGE_W - MARG - LEFT_W - GUT - MARG      # 60
LEFT_X  <- MARG
RIGHT_X <- MARG + LEFT_W + GUT                      # 123

# ---- WHY THE PAGE IS 215 mm, AND WHY NOT MORE -------------------------------
# Measured on the 170 mm page: B's box was 50 mm, of which 7.70 mm went to the
# x-axis strip, leaving a 42.30 mm plot HEIGHT; the width available to it in this
# column is 46.90 mm (what the free-aspect panels E and G actually occupy). A
# square takes the smaller of the two, so B and F were HEIGHT-limited at 42.30 and
# 4.6 mm of their width sat unused. Growing the page therefore does grow them --
# but only up to 46.90 mm, where the width takes over as the binding constraint.
#
# That cap also fixes the MAXIMUM page height, because A and B share a row and
# qa_svg_geometry.py asserts their plot heights are equal:
#
#     A plot height = A_H - A_OVERHEAD = (B_H + 7.285) - 14.985 = B_H - 7.70
#     B plot height = min(B_H - 7.70, 46.90)          <- the square's cap
#
# The two agree only while B_H - 7.70 <= 46.90, i.e. B_H <= 54.60. Past that, A
# keeps growing and B stops, and the row breaks. B_H = 54.5 sits 0.1 mm inside
# that limit rather than on the knife edge.
#
# So A, B and F are DONE at any page height >= 195, and the remaining height goes
# where it still buys something: the UMAP pair, and E and G. The UMAPs set the
# height. Each has half the column, but the LEFT one gives up ~5.5 mm to its
# "UMAP 2" title, so the binding width is 54.5 mm, not 60. They are square and
# height-limited at C_UM_H - 5.435, and pushing them all the way to 54.5 would
# leave a 2.8 mm gutter between two UMAPs -- they would read as one panel. At
# 47.8 mm the gutter is 10.4 mm and the pair spans the column the way the heatmap
# above it does, which is what 200 mm buys over 195 mm.
#   A  109.05 x 42.30 -> 109.05 x 46.80      B, F  42.29 sq -> 46.80 sq
#   C  120 x 72       -> 120 x 98            D     35.28 sq -> ~48.8 sq
#   E     46.90 x 27.30 -> 46.90 x 38.30     G     27.30 -> 52.30 over two rows
# 52.15, down from 54.5. G's y ticks are now "-0.50"-wide, and the assembly
# equalises non-panel WIDTHS across the whole right column, so the column's plot
# width fell 46.90 -> 44.55 mm. B and F are square, so their cap fell with it, and
# B stopped matching A's plot height in row 1. The square's box height has to
# follow its width cap, not the other way round: 44.55 cap - 0.1 margin + 7.70 mm
# of x-axis strip = 52.15.
# 49.6. Third retune, same mechanism each time: G's y-axis labels set the width of
# the strip that the assembly shares across the whole right column, so whenever G
# changes what it puts on its y axis the column's PLOT WIDTH moves, and the square
# panels -- capped by that width -- stop matching A's plot height in row 1. G now
# carries 21 gene NAMES rather than numbers, the strip grew again, and the column
# width fell 44.55 -> 41.99. SQ_H = 41.99 cap - 0.1 margin + 7.70 x-axis strip.
SQ_H <- 54.5      # square box: width cap 46.90 - 0.1 + 7.70 x strip (see above)
# 2026-09-06: the figure is A-F. The cohort scatter and the per-gene dosage chart
# are Supplementary; they were defensive panels (sample selection, "is it dosage")
# and the standard find-clones -> compare-clones -> interpret narrative has no slot
# for defence in the main figure. Right column is B, E (motif volcano), F.
B_H <- SQ_H       # precision-recall   (square)
E_H <- SQ_H       # TF motif volcano   (square: both axes are unit-free scores)
F_H <- SQ_H       # ATAC vs RNA        (square)
stopifnot(abs(B_H + E_H + F_H - PAGE_H) < 1e-9)

# A's extra height is exactly what its 45-degree rotated method labels cost over
# B's horizontal ones. TUNED AGAINST A MEASUREMENT, not guessed: run
# qa_svg_geometry.py on the exported SVG, which reads the axis-line coordinates
# out of the vector file and reports the A-vs-B plot-area mismatch in mm; move
# A_EXTRA by that amount and re-run until it reports OK.
A_EXTRA <- 7.285
A_H <- B_H + A_EXTRA
C_H <- PAGE_H - A_H
# The left column below A is the heatmap over the UMAP pair. The UMAPs are the
# accuracy evidence for the AG-DCAE partition -- the same embedding coloured by
# our clone and by the study's own cell state -- and they are square. Each has
# half the column, so ~54 mm of width is available and they are HEIGHT-limited:
# their size is C_UM_H less ~5.4 mm of title and axis title.
C_HM_H  <- 60                     # heatmap, carries the shared CN + state legend
C_UM_H  <- C_H - C_HM_H           # the UMAP pair

# Panel letters are inset from the top-left of their region by this much, in mm.
# NOT in npc: 0.004 npc was 0.68 mm on a 170 mm page and 0.78 mm on this one, so
# an npc inset silently changes size whenever the page does.
LTR_DROP <- 0.6

# y origins (bottom-left, as cowplot wants them)
A_Y    <- PAGE_H - A_H
C_HM_Y <- C_UM_H
C_UM_Y <- 0
F_Y <- 0
E_Y <- F_H
B_Y <- F_H + E_H              # == PAGE_H - B_H
stopifnot(abs(B_Y - (PAGE_H - B_H)) < 1e-9)
