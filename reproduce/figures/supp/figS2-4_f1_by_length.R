#!/usr/bin/env Rscript
# Supplementary Figs. 2-4: event-level precision and recall by CNA length (short < 3 Mb,
# medium 3-10 Mb, long > 10 Mb) and direction (DEL / DUP), one figure per coverage
# (S2 = 0.1X, S3 = 0.05X, S4 = 0.2X), three blocks a / b / c for noise 0 / 0.05 / 0.1.
# An event counts as recovered when a called segment of the same direction overlaps >= 20 %
# of it (greedy one-to-one matching, pooled over cells).
# Input : reproduce/figures/data/simulation_tables/by_setting/summary_mss_<setting>.csv
#         (reproduce/analysis/simulation/evaluate.R)
# Output: reproduce/figures/output/supp/SuppFig02-04_F1_by_length_c<coverage>.pdf (+ PNG previews)
# Run from the root of the MintCNA repository: Rscript reproduce/figures/supp/figS2-4_f1_by_length.R

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(ggh4x); library(cowplot)
})
FIG <- "reproduce/figures"   # run from the root of the MintCNA repository
if (!file.exists(file.path(FIG, "shared/fig_palette.R"))) stop("run this script from the root of the MintCNA repository")
source(file.path(FIG, "shared/fig_palette.R"))
CACHE <- file.path(FIG, "data/simulation_tables/by_setting")   # reproduce/analysis/simulation/evaluate.R
OUT <- file.path(FIG, "output/supp"); QA <- OUT; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# legend order, two rows filled across: MintCNA, MintCNA^DNA, rcCAE, SCOPE / SCONCE, SCCNV, SCICoNE, HMMcopy
ORD <- c("MintCNA", "MintCNA(DNA)", "rcCAE", "SCOPE", "SCONCE", "SCCNV", "SCICoNE", "HMMcopy")
PT_SZ <- 1.9
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

## ---- three stacked blocks (a / b / c = noise 0 / 0.05 / 0.1) ---------------------------
# Each block is its own 2 x 3 facet plot (DEL / DUP rows, length-stratum columns); the block letter
# is drawn by plot_grid at the top-left corner of the block, outside the panels and the strips.
# Column strips only on the first block, the y title on the middle one, the x title and the legend
# on the last one. The noise levels are not written on the figure: the letters carry them and the
# legend spells them out.
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
          legend.justification = "center",                    # theme_cowplot left-aligns it otherwise
          legend.byrow      = TRUE,                           # fill the two rows across, not down
          legend.key.spacing.x = unit(3, "mm"),               # explicit gap: cairo renders the labels wider than ggplot measured them
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
            label_fontfamily = FIG_FAMILY, label_x = 8 / 120, label_y = 1, hjust = 0, vjust = 1)   # 8 mm in, just left of the DEL strip
}

W_MM <- 120; H_MM <- 205; mm <- 1/25.4
for (nm in names(cov_map)) {
  fig <- mk_fig_tall(cov_map[[nm]])
  # zero-padded so the folder lists in figure order
  base <- sprintf("SuppFig%02d_F1_by_length_c%s", as.integer(sub("S", "", nm)), cov_map[[nm]]$cov)
  cairo_pdf(file.path(OUT, paste0(base, ".pdf")), width = W_MM*mm, height = H_MM*mm,
            family = FIG_FAMILY)
  print(fig); dev.off()
  ragg::agg_png(file.path(QA, paste0(base, "_preview.png")), width = W_MM*mm,
                height = H_MM*mm, units = "in", res = 600)
  print(fig); dev.off()
  cat(sprintf("wrote %s.pdf  (%d x %d mm)\n", base, W_MM, H_MM))
}
