#!/usr/bin/env Rscript
# Supplementary Fig. 1 -- AG-DCAE denoising example: two simulated cells, three
# rows each (raw input, CAE + MSE baseline, AG-DCAE), chromosomes 1-6 plus a zoom
# on one focal event.
#
# REBUILT 2026-10-09 FOR NATURE COMMUNICATIONS. The published S1 was two runs of
# Simulation/cnasim/c3_N100_100kb_l1e4_n0.1/FigureS1/scripts/plot.py (matplotlib,
# 13 x 7 in per cell) pasted together in an SVG: 290 x 330 mm, all axis text
# rasterised into the bitmaps, the only vector text the 18 pt panel letters in
# Arial-Bold + Liberation Sans. NC wants 180 mm wide, <= 170 mm tall, 5-8 pt text
# in ONE embedded typeface, so the figure is regenerated here in R at final size
# with the same Helvetica/Nimbus family as every other figure. The scatter clouds
# (11,513 bins x 6 rows) are the only rasterised objects (ggrastr, 600 dpi); every
# line and glyph is vector.
#
# Numbers are reproduced from plot.py, not re-derived: same chr1-6 window, same
# zoom rule (first run of the target level 25-100 bins long on chr1, +-30 bins),
# same RMSE and the same transition-width metric (W = 8, central 50% band between
# the locally recovered levels; not shown for the raw row). Cells: panel a is
# cell68 (focal deletion, CN = 1) and panel b is cell40 (focal duplication,
# CN = 3) -- identified by matching the RMSE values printed on the published
# figure (0.309/0.231/0.226 and 0.319/0.279/0.264) against every cell.

.libPaths(c("<R_LIBS_USER>", .libPaths()))
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(cowplot); library(ggrastr)
})
ROOT <- "<MINTCNA_ROOT>"
source("<REPRODUCE_DIR>/07_figures/shared/fig_palette.R")
SIM <- file.path(ROOT, "Simulation/cnasim/c3_N100_100kb_l1e4_n0.1")
S1D <- file.path(SIM, "FigureS1/data")
OUT <- Sys.getenv("FIGS1_OUT", file.path(ROOT, "figure4paper/figureS-supplementary"))
QA  <- Sys.getenv("FIGS1_QA",  file.path(ROOT, "figure4paper/figureS-supplementary/_work/qa"))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE); dir.create(QA, showWarnings = FALSE, recursive = TRUE)

PANELS <- list(a = list(cell = "cell68", type = "deletion",    lvl = -1.0,  cn = 1),
               b = list(cell = "cell40", type = "duplication", lvl = 0.585, cn = 3))
ROWS      <- c("Raw input", "CAE + MSE", "AG-DCAE")
ROW_COL   <- c("Raw input" = "#8A9199", "CAE + MSE" = "#D08C2E", "AG-DCAE" = unname(HERO["MintCNA"]))
TRUTH_RED <- "#B2182B"; ZOOM_FILL <- "#FFD84D"
YL <- c(-1.6, 1.6); ZOOM_PAD <- 30L; W_TR <- 8L

## ---- data -------------------------------------------------------------------
read_mat <- function(f) { d <- fread(f); setnames(d, 1, "bin"); d }
raw <- read_mat(file.path(SIM, "l2r_smoothed_all_cells_100kb.csv"))
tru <- read_mat(file.path(S1D, "truth_log2.csv"))
cae <- read_mat(file.path(S1D, "caemse_denoised.csv"))
ag  <- read_mat(file.path(S1D, "agdcae_denoised.csv"))
stopifnot(nrow(tru) == nrow(raw), nrow(cae) == nrow(raw), nrow(ag) == nrow(raw))
for (d in list(tru, cae, ag)) if (is.character(d$bin)) stopifnot(identical(d$bin, raw$bin))
chrom <- sub(":.*", "", raw$bin)
keep  <- chrom %in% paste0("chr", 1:6)
chrom_v <- chrom[keep]; n <- sum(keep)
chr_start <- which(chrom_v != shift(chrom_v, fill = chrom_v[1])) - 1L      # 0-based, as plot.py
chr_lab <- data.table(chr = unique(chrom_v))
chr_lab[, mid := sapply(chr, function(c) { i <- which(chrom_v == c); (min(i) + max(i)) / 2 - 1 })]

find_focal_run <- function(m, min_len = 25L, max_len = 100L) {
  r <- rle(m); e <- cumsum(r$lengths); s <- e - r$lengths + 1L
  i <- which(r$values & r$lengths >= min_len & r$lengths <= max_len)
  if (!length(i)) return(NULL)
  c(s[i[1]], e[i[1]])                                                        # 1-based, inclusive
}
rmse <- function(a, b) sqrt(mean((a - b)^2))
# plot.py transition_width(): per true breakpoint, the number of bins inside
# +-W whose value lies in the central 50% band between the locally recovered
# (median) levels on either side. Lower = sharper edge.
transition_width <- function(s, truth_z, W = W_TR) {
  bk <- which(abs(diff(truth_z)) > 0.05)                                     # left element of each step
  if (!length(bk)) return(NA_real_)
  widths <- c()
  for (k in bk) {
    l_lo <- max(1L, k - W); l_hi <- k - 1L; r_lo <- k + 1L; r_hi <- min(length(s), k + W)
    if (l_hi - l_lo + 1L < 3L || r_hi - r_lo + 1L < 3L) next
    L <- median(s[l_lo:l_hi]); R <- median(s[r_lo:r_hi])
    mid <- (L + R) / 2; half <- 0.25 * abs(R - L)
    if (half < 1e-6) next
    widths <- c(widths, sum(abs(s[l_lo:r_hi] - mid) < half))
  }
  if (length(widths)) mean(widths) else NA_real_
}

build_cell <- function(spec) {
  v <- lapply(list(raw = raw, truth = tru, cae = cae, ag = ag), function(d) d[[spec$cell]][keep])
  is_ev <- abs(v$truth - spec$lvl) < 0.05
  run <- find_focal_run(is_ev & chrom_v == "chr1"); if (is.null(run)) run <- find_focal_run(is_ev)
  zi <- max(1L, run[1] - ZOOM_PAD):min(n, run[2] + ZOOM_PAD)
  z0 <- min(zi) - 1; z1 <- max(zi)                                           # 0-based half-open, as plot.py
  sig <- list("Raw input" = v$raw, "CAE + MSE" = v$cae, "AG-DCAE" = v$ag)
  long <- rbindlist(lapply(ROWS, function(r) data.table(row = r, x = seq_len(n) - 1, value = sig[[r]], truth = v$truth)))
  long[, row := factor(row, levels = ROWS)]
  met <- rbindlist(lapply(ROWS, function(r) data.table(
    row = factor(r, levels = ROWS),
    rmse_g = rmse(sig[[r]], v$truth), rmse_z = rmse(sig[[r]][zi], v$truth[zi]),
    tw = if (r == "Raw input") NA_real_ else transition_width(sig[[r]][zi], v$truth[zi]))))
  met[, lab_g := sprintf("RMSE = %.3f", rmse_g)]
  met[, lab_z := ifelse(is.na(tw), sprintf("RMSE = %.3f", rmse_z),
                        sprintf("RMSE = %.3f\ntransition = %.1f bins", rmse_z, tw))]
  cat(sprintf("%s %s: zoom bins [%d, %d) on %s, event %d bins\n", spec$cell, spec$type, z0, z1,
              chrom_v[run[1]], run[2] - run[1] + 1L))
  print(met[, .(row, rmse_g = round(rmse_g, 3), rmse_z = round(rmse_z, 3), tw)])
  list(long = long, zoom = long[x >= z0 & x < z1], met = met, z0 = z0, z1 = z1,
       # two lines: at 165 mm the zoom column is ~37 mm, and the one-line title was clipped
       ztitle = sprintf("Zoom on %s\nfocal %s, CN = %d", chrom_v[run[1]], spec$type, spec$cn),
       cell = spec$cell)
}

## ---- plotting ---------------------------------------------------------------
base_theme <- function() {
  theme_cowplot(BASE_PT, font_family = FIG_FAMILY) +
    theme(axis.text        = element_text(size = PT_SMALL, colour = "grey20"),
          axis.title       = element_text(size = BASE_PT),
          axis.line        = element_line(linewidth = 0.3, colour = "grey30"),
          axis.ticks       = element_line(linewidth = 0.3, colour = "grey30"),
          axis.ticks.length= unit(1.5, "pt"),
          strip.background = element_rect(fill = "grey94", colour = NA),
          strip.text       = element_text(size = BASE_PT - 0.5, face = "bold"),
          strip.placement  = "outside",
          panel.spacing.y  = unit(4, "pt"),
          panel.grid       = element_blank(),
          legend.position  = "none",
          plot.margin      = margin(1, 2, 1, 1))
}
ybreaks <- scale_y_continuous(breaks = c(-1, 0, 1), labels = c("−1", "0", "1"))
ann_box <- function(d, x, y, lab, hj = 0) {
  geom_label(data = d, aes(x = {{x}}, y = {{y}}, label = {{lab}}), inherit.aes = FALSE,
             hjust = hj, vjust = 1, size = PT_SMALL / .pt, family = FIG_FAMILY, lineheight = 0.95,
             fill = "white", alpha = 0.85, linewidth = 0, label.padding = unit(0.6, "pt"),
             label.r = unit(0, "pt"))
}

genome_plot <- function(cd) {
  ggplot(cd$long, aes(x, value)) +
    annotate("rect", xmin = cd$z0, xmax = cd$z1, ymin = -Inf, ymax = Inf, fill = ZOOM_FILL, alpha = 0.45) +
    geom_vline(xintercept = chr_start, colour = "grey72", linewidth = 0.25, linetype = "12") +
    geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.2) +
    rasterise(geom_point(aes(colour = row), size = 0.28, alpha = 0.45, shape = 16, stroke = 0),
              dpi = 1200, dev = "ragg") +   # 600 until 2026-10-13; the author asked for higher resolution
    geom_step(aes(y = truth), direction = "hv", colour = TRUTH_RED, linewidth = 0.3) +
    ann_box(cd$met, 90, 1.5, lab_g) +
    facet_grid(row ~ ., switch = "y") +
    coord_cartesian(xlim = c(0, n), ylim = YL, expand = FALSE) +
    scale_colour_manual(values = ROW_COL) + ybreaks +
    # chromosome names ride on a secondary axis: ggplot draws it above the top
    # row only, outside the panel, so the panel can stay clipped (a truth spike
    # of CN = 7 and the noisiest raw bins exceed the +-1.6 window)
    scale_x_continuous(breaks = seq(0, 10000, 2000),
                       sec.axis = dup_axis(breaks = chr_lab$mid, labels = chr_lab$chr, name = NULL)) +
    labs(x = "Genomic bin (100 kb), chromosomes 1–6", y = expression(log[2]~ratio)) +
    base_theme() +
    theme(axis.title.y = element_text(size = PT_SUP),                       # carries a subscript
          axis.text.x.top  = element_text(size = PT_SMALL, colour = "grey35", margin = margin(b = 1)),
          axis.ticks.x.top = element_blank(), axis.line.x.top = element_blank(),
          strip.text.y.left = element_text(angle = 90, size = BASE_PT - 0.5, face = "bold",
                                           margin = margin(l = 1.5, r = 1.5)))
}
zoom_plot <- function(cd) {
  ggplot(cd$zoom, aes(x, value)) +
    geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.2) +
    geom_point(aes(colour = row), size = 0.75, shape = 16, alpha = 0.9, stroke = 0) +
    geom_step(aes(y = truth), direction = "hv", colour = TRUTH_RED, linewidth = 0.4) +
    ann_box(cd$met, cd$z0 + 0.02 * (cd$z1 - cd$z0), 1.5, lab_z) +
    facet_grid(row ~ ztitle, labeller = labeller(ztitle = function(x) cd$ztitle)) +
    coord_cartesian(xlim = c(cd$z0, cd$z1), ylim = YL, expand = FALSE) +
    scale_colour_manual(values = ROW_COL) + ybreaks +
    scale_x_continuous(breaks = scales::breaks_pretty(n = 3)) +
    labs(x = "Genomic bin (100 kb)", y = NULL) +
    base_theme() +
    theme(strip.text.y = element_blank(), strip.background.y = element_blank(),
          strip.text.x = element_text(size = BASE_PT - 0.5, face = "bold", lineheight = 0.95,
                                      margin = margin(t = 1.2, b = 1.2)))
}
panel <- function(cd) {
  cd$zoom[, ztitle := "z"]
  patchworkGrob((genome_plot(cd) | zoom_plot(cd)) + plot_layout(widths = c(3, 1.05)))
}

cells <- lapply(PANELS, build_cell)
key <- ggdraw() +
  draw_line(x = c(0.300, 0.330), y = c(0.5, 0.5), colour = TRUTH_RED, linewidth = 0.4) +
  draw_label(expression("Ground truth, " * log[2] * "(CN/2)"), x = 0.338, y = 0.5, hjust = 0,
             size = PT_SUP, fontfamily = FIG_FAMILY) +
  draw_grob(grid::rectGrob(gp = grid::gpar(fill = ZOOM_FILL, alpha = 0.45, col = NA)),
            x = 0.555, y = 0.22, width = 0.030, height = 0.56) +
  draw_label("Zoom region", x = 0.593, y = 0.5, hjust = 0, size = BASE_PT, fontfamily = FIG_FAMILY)

fig <- plot_grid(panel(cells$a), panel(cells$b), key, ncol = 1, rel_heights = c(1, 1, 0.075),
                 labels = c("a", "b", ""), label_size = 8, label_fontface = "bold",
                 label_fontfamily = FIG_FAMILY, label_x = 0.002, label_y = 0.995, hjust = 0, vjust = 1)

W_MM <- 165; H_MM <- 150; mm <- 1 / 25.4   # 165 = text width of the Letter/1-inch SI page; Word would scale 180 to 92 %
base <- "SuppFig01_denoising_example"
cairo_pdf(file.path(OUT, paste0(base, ".pdf")), width = W_MM * mm, height = H_MM * mm, family = FIG_FAMILY)
print(fig); dev.off()
suppressWarnings(pdftools::pdf_convert(file.path(OUT, paste0(base, ".pdf")), format = "png", dpi = 200,
                      filenames = file.path(QA, paste0(base, "_preview.png")), verbose = FALSE))
cat(sprintf("wrote %s.pdf  (%d x %d mm)\n", base, W_MM, H_MM))
