#!/usr/bin/env Rscript
# Figure 3c -- pseudobulk CNA tracks vs matched bulk WGS.  REVISION of 13_pseudobulk_tracks_grid.R
#
# The one structural change: the bulk WGS profile is no longer only the top-left
# tile of an 8-tile grid. It is repeated as a grey silhouette BEHIND every
# caller's track. The panel exists to answer "does this caller put its mass where
# the bulk does", and in the original that comparison is made across tiles -- the
# reader holds the truth in working memory and scans seven panels against it.
# Drawn behind each caller, the same comparison becomes local: any red or blue
# spilling outside the grey is a false call, any grey left uncovered is a miss.
# The first tile still shows the bulk alone, so the grey is DEFINED before it is
# used; that is what buys the right to repeat it silently in the other seven.
#
# Second change: the encoding key moved into the panel. The original set
# legend.position = "none" and left the colours to the caption; a reader who
# starts at the figure then has no way to know grey is the reference.
#
# Kept from the original, deliberately: 500 kb display bins with every statistic
# computed on the raw 100 kb bins; alternating chromosome bands instead of 22
# crowded ticks; the y axis as "fraction of tumor altered" rather than a copy
# number the collapsed 1/2/3 codes no longer carry; and the Fig 3b panel order.

.libPaths(c("<R_LIBS_USER>", .libPaths()))
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(cowplot)
})
ROOT <- "<MINTCNA_ROOT>"
source("<REPRODUCE_DIR>/07_figures/figure4-CRC/utils.R")
source("<REPRODUCE_DIR>/07_figures/shared/fig_palette.R")
OUT <- file.path(ROOT, "figure4paper/figure4-CRC/output")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

args <- commandArgs(trailingOnly = TRUE)
# No printed r in the strips (pass --with-r to restore). Per-cell r duplicated
# panel B and was removed earlier; the pseudobulk r on its own reads as a ranking
# -- HMMcopy 0.97 next to MintCNA 0.98 -- which is the opposite of what this panel
# is for. It is here to be read as SHAPE: does the caller put its mass where the
# bulk does. The numbers stay in data/fig4c_track_metrics.csv.
SHOW_R <- "--with-r" %in% args

DISP <- 5e5
SITE_TO_TRUTHCOL <- c(PT = "PT", LN = "LN1", ML = "ML1", MP = "MP1")

## ------------------------------------------------------------------- load ---
d       <- load_all_methods()
methods <- d$methods[setdiff(names(d$methods), "MintCNA(DNA)")]
truth   <- d$truth_mat
ref     <- as.data.table(d$ref)
ref[, chr := factor(chr, levels = CHR_ORDER)]
stopifnot(nrow(truth) == 25262L, nrow(ref) == 25262L, !any(is.na(ref$chr)))
stopifnot(identical(order(ref$chr, ref$start), seq_len(nrow(ref))))

cells  <- colnames(methods$MintCNA)
site   <- sub("^CRC01_([A-Z]+)[0-9]*_.*", "\\1", cells)
w_site <- table(factor(site, levels = SITE_ORDER)) / length(cells)

freq_method <- function(m) list(gain = rowMeans(m == 3L), loss = rowMeans(m == 1L))
freq_truth <- function(T) {
  g <- l <- numeric(nrow(T))
  for (s in SITE_ORDER) { tc <- SITE_TO_TRUTHCOL[[s]]
    g <- g + w_site[[s]] * (T[, tc] == 3L); l <- l + w_site[[s]] * (T[, tc] == 1L) }
  list(gain = g, loss = l)
}
tr_truth <- freq_truth(truth)
tracks   <- lapply(methods, freq_method)
net      <- lapply(c(list(`bulk WGS` = tr_truth), tracks), function(x) x$gain - x$loss)

## -------------------------------------------------------------- statistics ---
percell <- rbindlist(lapply(names(methods), function(k) {
  m <- methods[[k]]
  r <- suppressWarnings(vapply(cells, function(cc)
    cor(truth[, SITE_TO_TRUTHCOL[[site[match(cc, cells)]]]], m[, cc]), numeric(1)))
  data.table(method = k, pearson = r[is.finite(r)])
}))
stats <- rbindlist(lapply(names(methods), function(k) data.table(
  method = k,
  pseudobulk_r     = cor(net[[k]], net[["bulk WGS"]], use = "complete.obs"),
  pseudobulk_rmse  = sqrt(mean((net[[k]] - net[["bulk WGS"]])^2)),
  percell_median_r = percell[method == k, median(pearson)])))
setorder(stats, -pseudobulk_r)
print(stats)
# The discordance this panel exists to expose: HMMcopy is 2nd by pseudobulk r and
# 6th by per-cell r. Its error is per-cell noise that averages out of a population
# profile. That is the panel's one inference that Fig 3b cannot make, and in the
# original it is invisible because no number is printed anywhere on the figure.
cat("\nrank by pseudobulk r :", paste(stats[order(-pseudobulk_r), method], collapse = " > "), "\n")
cat("rank by per-cell   r :", paste(stats[order(-percell_median_r), method], collapse = " > "), "\n")

## ------------------------------------------------------- genome coordinates ---
ref[, gb := start %/% DISP]
len <- ref[, .(m = max(end)), by = chr][order(chr)]
len[, `:=`(off = cumsum(c(0, head(m, -1))), band = seq_len(.N) %% 2 == 1)]

to_long <- function(tr, nm) data.table(method = nm, chr = ref$chr, gb = ref$gb,
                                       gain = tr$gain, loss = tr$loss)
long <- rbindlist(c(list(to_long(tr_truth, "bulk WGS")),
                    lapply(names(tracks), function(k) to_long(tracks[[k]], k))))
agg  <- long[, .(gain = mean(gain), loss = mean(loss)), by = .(method, chr, gb)]
agg  <- merge(agg, len[, .(chr, off)], by = "chr")
agg[, x := (off + gb * DISP + DISP / 2) / 1e6]

ord <- c("bulk WGS", rank_methods(percell, "pearson"))
LAB <- c(`bulk WGS` = "Bulk WGS", setNames(names(methods), names(methods)))
if (SHOW_R) {
  # Only the PSEUDOBULK r. The median per-cell r used to sit here too, but that is
  # exactly the number panel B already plots as a distribution, and design-theory.md
  # section 10's anti-redundancy check asks that panel c "adds a dimension absent
  # from a and b". The pseudobulk r is that dimension: it is computed on the
  # population-level track, not per cell.
  #
  # Consequence for the text: the discordance this panel exposes -- HMMcopy reaches
  # r_pb 0.97 while its per-cell median is 0.67, because its error is cell-level
  # noise that averages away in the pseudobulk -- is no longer readable by comparing
  # two numbers inside one strip. It must be stated in the legend or Results instead;
  # the per-method values for both are in data/fig4c_track_metrics.csv.
  rr <- setNames(sprintf("~~italic(r)[pb]==%.2f", stats$pseudobulk_r), stats$method)
}

# One red/blue layer for all eight facets, the bulk WGS track included. The bulk
# used to be drawn grey and repeated as a silhouette behind every caller; that
# bought a within-panel comparison at the cost of a second encoding the reader
# had to learn. Now every track uses the same gain/loss colours and the reader
# compares shapes across facets against the bulk facet at top-left.
caller <- melt(agg, id.vars = c("method", "x"),
               measure.vars = c("gain", "loss"), variable.name = "type", value.name = "f")
caller[, y := fifelse(type == "loss", -f, f)]
caller[, type := factor(fifelse(type == "gain", "Gain", "Loss"), levels = c("Gain", "Loss"))]

panel_lv <- ord
mk_lab <- function(k) {
  nm <- sprintf('bold("%s")', LAB[[k]])
  if (SHOW_R && k != "bulk WGS") paste0(nm, rr[[k]]) else nm
}
caller[, panel := factor(method, levels = panel_lv)]

bands <- len[band == TRUE, .(xmin = off / 1e6, xmax = (off + m) / 1e6)]
mid   <- len[, .(chr, x = (off + m / 2) / 1e6, w_mb = m / 1e6)][, txt := sub("chr", "", chr)]

## label thinning, from the FINAL rendered width -------------------------------
W_MM <- 183; H_MM <- 75; mm <- 1 / 25.4
panel_in <- (W_MM * mm - 0.55) / 4
lbl_in   <- 1.6 * PT_SMALL / 72
gap_mb   <- max(mid$x + mid$w_mb / 2) / panel_in * lbl_in
last <- -Inf
for (i in seq_len(nrow(mid))) if (mid$x[i] - last < gap_mb) mid$txt[i] <- "" else last <- mid$x[i]
cat("chromosome ticks shown:", paste(mid$txt[mid$txt != ""], collapse = " "), "\n")

strip_labs <- setNames(vapply(panel_lv, mk_lab, character(1)), panel_lv)

p <- ggplot() +
  geom_rect(data = bands, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = "grey55", alpha = 0.07, inherit.aes = FALSE) +
  geom_col(data = caller, aes(x, y, fill = type), width = DISP / 1e6,
           position = "identity", alpha = 0.85) +
  geom_hline(yintercept = 0, colour = "grey30", linewidth = 0.2) +
  facet_wrap(~ panel, nrow = 2, labeller = labeller(panel = as_labeller(strip_labs, label_parsed))) +
  scale_fill_manual(values = c(Gain = GAIN, Loss = LOSS), guide = "none") +
  scale_x_continuous(breaks = mid$x, labels = mid$txt, expand = c(0.004, 0.004)) +
  scale_y_continuous(breaks = c(-1, -0.5, 0, 0.5, 1), labels = c("1.0", "0.5", "0", "0.5", "1.0")) +
  coord_cartesian(ylim = c(-1, 1)) +
  labs(x = "Chromosome", y = "Fraction of tumor altered") +
  theme_cowplot(BASE_PT, font_family = FIG_FAMILY) +
  theme(
    strip.background   = element_blank(),
    # PT_SUP: the strips carry an r[pb] subscript, reduced to ~0.7x
    # The strip carries an r[pb] subscript, so the line box is taller than
    # a plain 7.4 pt label and it sits directly under the row above. A top margin
    # plus a wider inter-row gap keeps the descenders clear of that panel's floor.
    strip.text         = element_text(size = PT_SUP, hjust = 0,
                                      margin = margin(t = 0.5, b = 1.5)),
    panel.spacing.x    = unit(4, "pt"), panel.spacing.y = unit(7, "pt"),
    axis.text.x        = element_text(size = PT_SMALL, colour = "grey25"),
    axis.text.y        = element_text(size = PT_SMALL, colour = "grey25"),
    axis.title         = element_text(size = BASE_PT),
    axis.ticks         = element_line(colour = "grey70", linewidth = 0.2),
    axis.line          = element_line(colour = "grey70", linewidth = 0.25),
    # No gridlines. design-theory.md sections 2 and 6 both say "no grid; sparse
    # y-ticks". Eight facets x three rules was 24 extra lines competing with the
    # bars, and these tracks are read as shape against the zero line, not off a
    # value; the zero line itself stays.
    panel.grid         = element_blank(),
    plot.margin        = margin(2, 2, 1, 2)
  )

## direct labels instead of a legend ------------------------------------------
# The encoding is named inside the first facet: Gain/Loss are written in their
# own colours on the two halves of the axis. Direct labels beat a legend here because the categories are
# spatially fixed -- up is gain, down is loss, always. This also avoids overlaying
# a second grob on the facet grid, which is what the original legend-in-caption
# choice was working around.
xmax <- max(agg$x)
# Coloured swatch + black text, not coloured text: the Nature figure guide asks
# for "coloured boxes combined with black text" and says coloured text may be
# recoloured to black in production (2026-10-09). The swatch is a filled square
# (shape 15) placed 0.9 mm into the facet; the word starts 2.4 mm in.
mm_mb <- xmax / (panel_in * 25.4)   # Mb per mm of facet width
keydt <- data.table(
  panel = factor("bulk WGS", levels = panel_lv),
  x_sw  = xmax * 0.015 + 0.9 * mm_mb,   # chr1-3: empty in the bulk track, unlike the chr17-21 tail
  x_txt = xmax * 0.015 + 2.4 * mm_mb,
  y     = c(0.80, -0.80),
  lab   = c("Gain", "Loss"),
  col   = c(GAIN, LOSS))
p <- p +
  geom_point(data = keydt, aes(x_sw, y, colour = col), shape = 15, size = 1.9,
             inherit.aes = FALSE, show.legend = FALSE) +
  geom_text(data = keydt, aes(x_txt, y, label = lab), colour = "grey15",
            hjust = 0, size = PT_SMALL / .pt, family = FIG_FAMILY,
            fontface = "bold", inherit.aes = FALSE, show.legend = FALSE) +
  scale_colour_identity()

final <- p

if (!exists("FIG4_ASSEMBLING")) {
ggsave(file.path(OUT, "fig4c_tracks.pdf"), final,
       width = W_MM * mm, height = H_MM * mm, device = cairo_pdf, family = FIG_FAMILY)
fwrite(stats, file.path(ROOT, "figure4paper/figure4-CRC/data/fig4c_track_metrics.csv"))
cat(sprintf("wrote fig4c_tracks.pdf  (%d x %d mm, base %g pt)\n", W_MM, H_MM, BASE_PT))
}
