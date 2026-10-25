#!/usr/bin/env Rscript
# Figure 4b: per-cell accuracy of every caller on CRC01 against the matched bulk WGS profile,
# Pearson r (top) and RMSE (bottom) per cell as a violin and box per method.
# Input (reproduce/figures/data/fig4_crc/): the CRC01 call set of every method and the post-QC bulk
#   WGS truth (DNAcopy_qc.RData), loaded through load_all_methods() in fig4_crc/utils.R.
# Output: reproduce/figures/output/fig4/fig4b_accuracy.pdf (66 x 93 mm).
# Run from the root of the MintCNA repository:
#   Rscript reproduce/figures/fig4_crc/panel_b_accuracy.R (or through fig4_crc/assemble_fig4.R)

suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(data.table); library(cowplot)
})
FIG <- "reproduce/figures"   # run from the root of the MintCNA repository
if (!file.exists(file.path(FIG, "shared/fig_palette.R"))) stop("run this script from the root of the MintCNA repository")
source(file.path(FIG, "fig4_crc/utils.R"))
source(file.path(FIG, "shared/fig_palette.R"))
OUT <- file.path(FIG, "output/fig4"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

## ------------------------------------------------------------------ metrics ---
d <- load_all_methods(); methods <- d$methods; truth <- d$truth_mat

# The truth matrix carries no rownames, so a caller is aligned to it by bin position. For every
# caller whose rownames carry coordinates, check_bin_order() verifies that they match the reference
# grid, so a data swap cannot silently break the pairing.
check_bin_order <- function(mats, ref) {
  ref <- as.data.frame(ref)
  parse_rn <- function(rn) {
    if (all(grepl("^chr[^_]*_\\d+_\\d+$", rn)))
      data.frame(chr = sub("_.*", "", rn),
                 start = as.numeric(sub("^[^_]+_(\\d+)_.*", "\\1", rn)))
    else if (all(grepl("^chr[^:]*:\\d+-\\d+$", rn)))
      data.frame(chr = sub(":.*", "", rn),
                 start = as.numeric(sub(".*:(\\d+)-.*", "\\1", rn)))
    else NULL
  }
  for (nm in names(mats)) {
    rn <- rownames(mats[[nm]]); if (is.null(rn)) next
    p <- parse_rn(rn); if (is.null(p)) next
    stopifnot(identical(p$chr, as.character(ref$chr)),
              all(p$start == ref$start | p$start == ref$start + 1))
  }
  stopifnot(all(vapply(mats, nrow, integer(1)) == nrow(ref)))
  invisible(TRUE)
}
check_bin_order(methods, d$ref)
cat("bin-order guard: all callers with coordinates match the reference grid\n")

align_bins <- function(T, P) {
  P <- as.matrix(P); cm <- intersect(rownames(T), rownames(P))
  if (length(cm) > 0) list(T = T[cm, , drop = FALSE], P = P[cm, , drop = FALSE]) else list(T = T, P = P)
}
per_cell <- function(T, P) {
  a <- align_bins(T, P); T <- a$T; P <- a$P
  m <- cell_to_truthcol(colnames(P)); m <- m[!is.na(truth_col) & truth_col %in% colnames(T)]
  m[, rmse := vapply(seq_len(.N), function(i)
        sqrt(mean((T[, truth_col[i]] - P[, cell[i]])^2, na.rm = TRUE)), numeric(1))]
  m[, pearson := vapply(seq_len(.N), function(i) {
        x <- as.numeric(T[, truth_col[i]]); y <- as.numeric(P[, cell[i]])
        k <- is.finite(x) & is.finite(y)
        if (sum(k) < 2L || sd(x[k]) == 0 || sd(y[k]) == 0) NA_real_ else cor(x[k], y[k])
      }, numeric(1))]
  m
}
tab <- rbindlist(lapply(names(methods), function(m) per_cell(truth, methods[[m]])[, method := m]))

ORD <- rank_methods(tab, "pearson")
tab[, method := factor(as.character(method), levels = ORD)]
sum_p <- tab[, .(med = median(pearson, na.rm = TRUE)), by = method][match(ORD, method)]
sum_r <- tab[, .(med = median(rmse,    na.rm = TRUE)), by = method][match(ORD, method)]

## cells with negative Pearson r ---------------------------------------------------
cat(sprintf("cells with Pearson < 0: %d of %d\n",
            tab[pearson < 0, .N], tab[, .N]))
print(tab[, .(min = round(min(pearson, na.rm = TRUE), 3),
              n_neg = sum(pearson < 0, na.rm = TRUE)), by = method][n_neg > 0])

## --------------------------------------------------------------------- theme ---
FILLS   <- method_fill(ORD)
lab_col <- label_ink(ORD)

base <- theme_cowplot(BASE_PT, font_family = FIG_FAMILY) +
  theme(
    axis.title.x       = element_blank(),
    # PT_SUP rather than BASE_PT - 1: this row carries MintCNA^DNA, and plotmath renders the
    # superscript at about 0.7x the parent size.
    axis.text.x        = element_text(size = PT_SUP, angle = 45, hjust = 1,
                                      vjust = 1, colour = lab_col),
    axis.text.y        = element_text(size = BASE_PT - 1, colour = "grey20"),
    axis.title.y       = element_text(size = BASE_PT),
    axis.line          = element_line(linewidth = 0.3, colour = "grey30"),
    axis.ticks         = element_line(linewidth = 0.3, colour = "grey30"),
    # No gridlines; the zero reference line of the Pearson panel is added separately because it
    # marks a meaningful value.
    panel.grid         = element_blank(),
    plot.margin        = margin(2, 2, 1, 2),
    legend.position    = "none"
  )

vb <- function(dt, yv, ylim, ylab, meds, show_x, clipped = NULL) {
  y_lo <- ylim[1]; y_hi <- ylim[2]; span <- y_hi - y_lo
  # Median labels sit directly above their own violin.
  top <- dt[, .(top = max(get(yv), na.rm = TRUE)), by = method]
  meds2 <- merge(meds, top, by = "method")
  # clip = "on" below cuts the low tails at the floor, so the median labels have to fit inside
  # head_room above the data.
  head_room <- 0.17
  meds2[, ytxt := pmin(top + 0.035 * span, y_hi + 0.02 * span)]
  p <- ggplot(dt, aes(x = method, y = .data[[yv]], fill = method)) +
    geom_violin(scale = "width", width = 0.86, colour = NA, alpha = 0.55, na.rm = TRUE) +
    geom_boxplot(width = 0.22, outlier.shape = NA, linewidth = 0.28,
                 fill = "white", colour = "grey25", na.rm = TRUE) +
    geom_text(data = meds2, aes(x = method, y = ytxt, label = sprintf("%.2f", med)),
              inherit.aes = FALSE, size = PT_SMALL / .pt,
              family = FIG_FAMILY, fontface = "bold", colour = "grey15", vjust = 0) +
    scale_fill_manual(values = FILLS, limits = ORD) +
    scale_x_discrete(limits = ORD, labels = method_labels(ORD)) +
    # Breaks come from the data range, not the padded range, so the headroom never adds a tick
    # above the metric's maximum.
    scale_y_continuous(breaks = scales::breaks_pretty(5)(c(y_lo, y_hi))) +
    labs(y = ylab) + coord_cartesian(ylim = c(y_lo, y_hi + head_room * span), clip = "on") + base
  if (y_lo < 0) p <- p + geom_hline(yintercept = 0, colour = "grey45",
                                    linewidth = 0.25, linetype = "22")
  if (!is.null(clipped) && nrow(clipped)) {
    # Cells cut by the floor are marked with a filled triangle on the axis; shape 25 is drawn by
    # the device, whereas a text glyph such as U+25BC would pull in a second font. Their count and
    # minimum are given in the figure legend.
    p <- p +
      geom_point(data = clipped, aes(x = method, y = y_lo), inherit.aes = FALSE,
                 shape = 25, size = 1.1, fill = "grey30", colour = "grey30") +
      annotate("blank", x = 0.55, y = y_lo,
               size = PT_SMALL / .pt)
  }
  # No P values are drawn: the cells come from a single tumour, so the replicate unit for a
  # between-method claim is the patient (n = 1). The medians above each violin carry the comparison.
  if (!show_x) p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  p
}

# The Pearson axis has a fixed floor; anything it cuts is counted here and marked in the panel, so
# the clip is disclosed inside the figure.
PEARSON_FLOOR <- -0.1
y_pear <- c(PEARSON_FLOOR, 1)
below <- tab[pearson < PEARSON_FLOOR, .(n_below = .N, min_r = min(pearson)), by = method][n_below > 0]
cat(sprintf("cells below the Pearson floor (%.2f): %d of %d\n",
            PEARSON_FLOOR, sum(below$n_below), nrow(tab)))
print(below, row.names = FALSE)
y_rmse <- c(0, ceiling(max(tab$rmse, na.rm = TRUE) * 10) / 10)
cat(sprintf("Pearson axis %.1f..%.1f   RMSE axis %.1f..%.1f\n", y_pear[1], y_pear[2], y_rmse[1], y_rmse[2]))

p1 <- vb(tab, "pearson", y_pear, expression("Pearson"~italic(r)~"vs bulk WGS"),
         sum_p, show_x = FALSE, clipped = below)
p2 <- vb(tab, "rmse",    y_rmse, "RMSE vs bulk WGS", sum_r, show_x = TRUE)

panel_b <- p1 / p2 + plot_layout(heights = c(1, 1))

## ---------------------------------------------------------------- final size ---
W_MM <- 66; H_MM <- 93                     # right column of the top row of Fig. 4
mm <- 1 / 25.4
if (!exists("FIG4_ASSEMBLING")) {
ggsave(file.path(OUT, "fig4b_accuracy.pdf"), panel_b,
       width = W_MM * mm, height = H_MM * mm, device = cairo_pdf, family = FIG_FAMILY)
cat(sprintf("wrote fig4b_accuracy.pdf  (%d x %d mm, base %g pt)\n", W_MM, H_MM, BASE_PT))
}
