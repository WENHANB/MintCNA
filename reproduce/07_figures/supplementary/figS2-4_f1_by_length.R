#!/usr/bin/env Rscript
# Supplementary Figs S2 / S3 / S4 -- event-level precision-recall stratified by CNA
# length, one figure per coverage (S2 = 0.1X, S3 = 0.05X, S4 = 0.2X).
#
# REBUILT 2026-10-09 FOR NATURE COMMUNICATIONS. The previous build
# (Simulation/results/S2_S4_F1_length/assemble_S2_S4.R) laid the 18 panels out as
# 3 columns x 6 rows on a 170 x 305 mm canvas, and the PDF-level lowercase-letter
# patch then re-wrapped it to 210 x 357 mm and embedded Verdana-Bold next to
# Helvetica and Liberation Sans. Against NC's rules -- 180 mm double column, a
# 210 x 276 mm page, ONE typeface (Arial/Helvetica) embedded, 5-8 pt text -- that
# fails on size and on typeface, and scaling it down to fit would put the text at
# 3-4 pt. So the layout is rebuilt rather than rescaled:
#   columns = length stratum (outer strip) x DEL/DUP (inner strip)   -> 6 across
#   rows    = noise level a / b / c (0, 0.05, 0.1)                     -> 3 down
# 180 mm wide, ~27 mm per square panel, text 5.5-7.4 pt, letters 8 pt bold
# lowercase, everything in the shared Helvetica/Nimbus family and embedded.
# The numbers are untouched: the per-setting P/R/F1 come from the cached
# summary_mss CSVs written by compute_one.R (TP = >=20% overlap, greedy 1:1,
# pooled over cells, by direction and by length stratum).

.libPaths(c("<R_LIBS_USER>", .libPaths()))
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(ggh4x); library(cowplot)
})
ROOT  <- "<MINTCNA_ROOT>"
source(file.path(ROOT, "figure4paper/shared/fig_palette.R"))
CACHE <- file.path(ROOT, "Simulation/results/S2_S4_F1_length/cache")
OUT   <- Sys.getenv("FIGS24_OUT", file.path(ROOT, "figure4paper/figureS-supplementary"))
QA    <- Sys.getenv("FIGS24_QA",  file.path(ROOT, "figure4paper/figureS-supplementary/_work/qa"))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
dir.create(QA,  showWarnings = FALSE, recursive = TRUE)

# FIGS24_LAYOUT: "tall" (default since 2026-10-11, the author's original orientation) =
#                3 columns (stratum) x 6 rows (noise block a/b/c x DEL/DUP), portrait,
#                130 x 200 mm, 30 mm squares; the noise levels are NOT written on the
#                figure -- the letters a/b/c carry them and the legend spells them out.
#                "wide" = 6 columns (stratum x DEL/DUP) x 3 rows (noise), 165 x 103 mm.
LAYOUT <- Sys.getenv("FIGS24_LAYOUT", "tall")
# legend order and fill: advisor's request 2026-10-12 -- two rows filled ACROSS,
#   MintCNA  MintCNA^DNA  rcCAE  SCOPE
#   SCONCE   SCCNV        SCICoNE HMMcopy
ORD <- c("MintCNA", "MintCNA(DNA)", "rcCAE", "SCOPE", "SCONCE", "SCCNV", "SCICoNE", "HMMcopy")
PT_SZ <- 1.9   # point size: 1.35 was invisible, 2.3 too big (author, 2026-10-12)
PM  <- c("MintCNA" = "MintCNA", "MintCNA(DNA)" = "MintCNA^DNA", "rcCAE" = "rcCAE",
         "SCOPE" = "SCOPE", "SCONCE" = "SCONCE", "HMMcopy" = "HMMcopy",
         "SCCNV" = "SCCNV", "SCICoNE" = "SCICoNE")
lab_by_name <- function(x) parse(text = PM[as.character(x)])     # keyed by name, never by position
SHAPES <- setNames(rep(21L, length(ORD)), ORD); SHAPES[c("MintCNA", "MintCNA(DNA)")] <- 24L

cov_map <- list(
  S2 = list(cov = "0.1",  bases = c("c3_N100_100kb_l1e4",       "c3_N100_100kb_l1e4_n0.05",       "c3_N100_100kb_l1e4_n0.1")),
  S3 = list(cov = "0.05", bases = c("c3_N100_100kb_l1e4_c0.05", "c3_N100_100kb_l1e4_c0.05_n0.05", "c3_N100_100kb_l1e4_c0.05_n0.1")),
  S4 = list(cov = "0.2",  bases = c("c3_N100_100kb_l1e4_c0.2",  "c3_N100_100kb_l1e4_c0.2_n0.05",  "c3_N100_100kb_l1e4_c0.2_n0.1")))

read_block <- function(base) {
  f <- file.path(CACHE, paste0("summary_mss_", base, ".csv"))
  if (!file.exists(f)) stop("missing cache: ", f)
  d <- fread(f)
  d[, method  := factor(method,  levels = ORD)]
  d[, state   := factor(state,   levels = c("del", "dup"), labels = c("DEL", "DUP"))]
  d[, stratum := factor(stratum, levels = c("<3Mb", "3-10Mb", ">10Mb"),
                        labels = c("Short (<3 Mb)", "Medium (3–10 Mb)", "Long (>10 Mb)"))]
  d[!is.na(precision) & !is.na(recall)]
}

# iso-F1 contours: F1 = 2PR/(P+R)  ->  R = F1*P / (2P - F1)
iso <- rbindlist(lapply(c(0.2, 0.4, 0.6, 0.8), function(f) {
  P <- seq(f/2 + 1e-3, 1, length.out = 300)
  data.table(f = f, precision = P, recall = f*P/(2*P - f))
}))[recall <= 1]
iso_lab <- iso[, .SD[which.min(abs(precision - 0.97))], by = f]

mk_fig <- function(spec) {
  dt <- rbindlist(lapply(seq_along(spec$bases), function(i)
    read_block(spec$bases[i])[, block := c("a", "b", "c")[i]]))
  dt[, block := factor(block, levels = c("a", "b", "c"))]
  # tall layout: one row per block x direction; the strip shows only DEL / DUP
  dt[, rowfac := factor(paste(block, state, sep = "."), levels = as.vector(t(outer(c("a", "b", "c"), c("DEL", "DUP"), paste, sep = "."))))]
  # panel letter: top-left of each row's first panel. 8 pt bold lowercase (NC).
  tag <- data.table(block = factor(c("a", "b", "c"), levels = c("a", "b", "c")),
                    state = factor("DEL", levels = c("DEL", "DUP")),
                    rowfac = factor(paste(c("a", "b", "c"), "DEL", sep = "."), levels = levels(dt$rowfac)),
                    stratum = factor("Short (<3 Mb)", levels = levels(dt$stratum)),
                    precision = 0.02, recall = 1.0, lab = c("a", "b", "c"))
  ggplot(dt, aes(precision, recall)) +
    geom_line(data = iso, aes(precision, recall, group = f), inherit.aes = FALSE,
              colour = "grey82", linewidth = 0.22) +
    geom_point(aes(colour = method, fill = method, shape = method), size = PT_SZ, stroke = 0.3) +
    geom_text(data = tag, aes(precision, recall, label = lab), inherit.aes = FALSE,
              hjust = 0, vjust = 1, size = 8/.pt, fontface = "bold", family = FIG_FAMILY) +
    scale_colour_manual(values = METHOD_HUES[ORD], breaks = ORD, labels = lab_by_name, name = NULL) +
    scale_fill_manual  (values = METHOD_HUES[ORD], breaks = ORD, labels = lab_by_name, name = NULL) +
    scale_shape_manual (values = SHAPES,           breaks = ORD, labels = lab_by_name, name = NULL) +
    scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
    ggh4x::facet_nested(block ~ stratum + state, axes = "all", remove_labels = "none") +   # both axes on every panel
    # No contour labels and no key line (author, 2026-10-11): the iso-F1 levels are
    # stated in the legend.
    labs(x = "Precision", y = "Recall") +
    theme_cowplot(BASE_PT, font_family = FIG_FAMILY) +
    theme(aspect.ratio      = 1,
          axis.text         = element_text(size = PT_SMALL, colour = "grey20"),
          axis.title        = element_text(size = BASE_PT),
          axis.line         = element_line(linewidth = 0.3, colour = "grey30"),
          axis.ticks        = element_line(linewidth = 0.3, colour = "grey30"),
          strip.background  = element_rect(fill = "grey94", colour = NA),
          strip.text.x      = element_text(size = BASE_PT - 0.5, face = "bold",
                                           margin = margin(t = 1.5, b = 1.5)),
          strip.text.y.left = element_text(size = BASE_PT - 0.5, face = "bold", angle = 90,
                                           margin = margin(l = 1.5, r = 1.5)),
          strip.placement   = "outside",
          panel.spacing.x   = unit(3.5, "pt"),
          panel.spacing.y   = unit(5, "pt"),
          panel.grid        = element_blank(),
          legend.position   = "bottom",
          legend.justification = "center",                    # theme_cowplot left-aligns it otherwise
          legend.byrow      = TRUE,                           # fill the two rows across, not down
          legend.text       = element_text(size = PT_SUP),   # carries MintCNA^DNA
          legend.key.size   = unit(3.2, "mm"),
          legend.margin     = margin(t = 1),
          plot.margin       = margin(2, 3, 1, 2)) +
    guides(colour = guide_legend(nrow = 1, override.aes = list(size = PT_SZ + 0.3)),
           fill = guide_legend(nrow = 1), shape = guide_legend(nrow = 1))
}

## ---- tall layout: three stacked blocks (a / b / c = noise 0 / 0.05 / 0.1) ---------------
# Each block is its own 2 x 3 facet plot (DEL / DUP rows, length-stratum columns) with
# ordinary DEL/DUP strips; the block letter is drawn by plot_grid at the top-left corner of
# the block, outside the panels and outside the strips (author, 2026-10-12: no nested
# letter strip -- it made one large grey box). Column strips only on the first block, the
# y title on the middle one, the x title and the legend on the last one.
mk_block <- function(dt, blk, col_strips = FALSE, x_title = FALSE, y_title = FALSE, legend = FALSE) {
  ggplot(dt[block == blk], aes(precision, recall)) +
    geom_line(data = iso, aes(precision, recall, group = f), inherit.aes = FALSE,
              colour = "grey82", linewidth = 0.22) +
    geom_point(aes(colour = method, fill = method, shape = method), size = PT_SZ, stroke = 0.3) +
    scale_colour_manual(values = METHOD_HUES[ORD], breaks = ORD, labels = lab_by_name, name = NULL) +
    scale_fill_manual  (values = METHOD_HUES[ORD], breaks = ORD, labels = lab_by_name, name = NULL) +
    scale_shape_manual (values = SHAPES,           breaks = ORD, labels = lab_by_name, name = NULL) +
    scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
    ggh4x::facet_grid2(state ~ stratum, switch = "y", axes = "all", remove_labels = "none") +
    # every block carries a y title so the three blocks have identical widths and their
    # panels line up; only the middle one is visible
    labs(x = if (x_title) "Precision" else NULL, y = "Recall") +
    theme_cowplot(BASE_PT, font_family = FIG_FAMILY) +
    theme(aspect.ratio      = 1,
          axis.text         = element_text(size = PT_SMALL, colour = "grey20"),
          axis.title        = element_text(size = BASE_PT),
          axis.title.y      = element_text(size = BASE_PT, colour = if (y_title) "black" else NA),
          axis.line         = element_line(linewidth = 0.3, colour = "grey30"),
          axis.ticks        = element_line(linewidth = 0.3, colour = "grey30"),
          strip.background  = element_rect(fill = "grey94", colour = NA),
          strip.text.x      = if (col_strips) element_text(size = BASE_PT - 0.5, face = "bold",
                                                           margin = margin(t = 1.5, b = 1.5)) else element_blank(),
          strip.background.x= if (col_strips) element_rect(fill = "grey94", colour = NA) else element_blank(),
          strip.text.y.left = element_text(size = BASE_PT - 0.5, face = "bold", angle = 90,
                                           margin = margin(l = 1.5, r = 1.5)),
          strip.placement   = "outside",
          panel.spacing.x   = unit(3.5, "pt"),
          panel.spacing.y   = unit(5, "pt"),
          panel.grid        = element_blank(),
          legend.position   = if (legend) "bottom" else "none",
          legend.justification = "center",
          legend.byrow      = TRUE,                           # two rows filled across (advisor)
          legend.key.spacing.x = unit(3, "mm"),               # explicit gap: cairo drew labels wider than the layout allowed (2026-10-13)
          legend.text       = element_text(size = PT_SUP, margin = margin(l = 1.2, r = 1.5, unit = "mm")),
          legend.key.size   = unit(3.2, "mm"),
          legend.margin     = margin(t = 1),
          # room above the DEL row for the block letter (blocks b and c have no column strips)
          plot.margin       = margin(if (col_strips) 2 else 11, 3, 1, 2)) +
    guides(colour = guide_legend(nrow = 2, override.aes = list(size = PT_SZ + 0.3)),
           fill = guide_legend(nrow = 2), shape = guide_legend(nrow = 2))
}
mk_fig_tall <- function(spec) {
  dt <- rbindlist(lapply(seq_along(spec$bases), function(i)
    read_block(spec$bases[i])[, block := c("a", "b", "c")[i]]))
  dt[, block := factor(block, levels = c("a", "b", "c"))]
  pa <- mk_block(dt, "a", col_strips = TRUE)
  pb <- mk_block(dt, "b", y_title = TRUE)
  pc <- mk_block(dt, "c", x_title = TRUE, legend = TRUE)
  plot_grid(pa, pb, pc, ncol = 1, align = "v", axis = "lr", rel_heights = c(1, 1.02, 1.22),
            labels = c("a", "b", "c"), label_size = 8, label_fontface = "bold",
            label_fontfamily = FIG_FAMILY, label_x = 8 / 120, label_y = 1, hjust = 0, vjust = 1)   # 8 mm in, just left of the DEL strip (author: 稍微右移)
}

if (LAYOUT == "tall") { W_MM <- 120; H_MM <- 205 } else { W_MM <- 165; H_MM <- 103 }; mm <- 1/25.4   # 165 = SI page text width
for (nm in names(cov_map)) {
  fig <- if (LAYOUT == "tall") mk_fig_tall(cov_map[[nm]]) else mk_fig(cov_map[[nm]])
  # file names follow NC's "Supplementary Fig. N": zero-padded so the folder lists in order
  base <- sprintf("SuppFig%02d_F1_by_length_c%s", as.integer(sub("S", "", nm)), cov_map[[nm]]$cov)
  cairo_pdf(file.path(OUT, paste0(base, ".pdf")), width = W_MM*mm, height = H_MM*mm,
            family = FIG_FAMILY)
  print(fig); dev.off()
  ragg::agg_png(file.path(QA, paste0(base, "_preview.png")), width = W_MM*mm,
                height = H_MM*mm, units = "in", res = 600)
  print(fig); dev.off()
  cat(sprintf("wrote %s.pdf  (%d x %d mm)\n", base, W_MM, H_MM))
}
