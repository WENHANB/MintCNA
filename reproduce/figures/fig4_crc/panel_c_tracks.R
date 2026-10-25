#!/usr/bin/env Rscript
# Figure 4c: pseudobulk CNA tracks of every caller on CRC01 against the matched bulk WGS profile.
# Each track is the fraction of cells with a gain (up) or a loss (down) along the genome; the bulk
# WGS profile, weighted by the number of cells per site, is the first facet.
# Input (reproduce/figures/data/fig4_crc/): the CRC01 call set of every method and the post-QC bulk
#   WGS truth (DNAcopy_qc.RData), loaded through load_all_methods() in fig4_crc/utils.R.
# Output: reproduce/figures/output/fig4/fig4c_tracks.pdf (183 x 75 mm) and fig4c_track_metrics.csv
#   (pseudobulk r and RMSE and the median per-cell r of every method).
# Run from the root of the MintCNA repository:
#   Rscript reproduce/figures/fig4_crc/panel_c_tracks.R (or through fig4_crc/assemble_fig4.R)

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(cowplot)
})
FIG <- "reproduce/figures"   # run from the root of the MintCNA repository
if (!file.exists(file.path(FIG, "shared/fig_palette.R"))) stop("run this script from the root of the MintCNA repository")
source(file.path(FIG, "fig4_crc/utils.R"))
source(file.path(FIG, "shared/fig_palette.R"))
OUT <- file.path(FIG, "output/fig4"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

DISP <- 5e5   # 500 kb display bins; every statistic below is computed on the raw 100 kb bins
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

# Calls are collapsed to loss / neutral / gain (1/2/3), so a track is the fraction of cells altered
# rather than a copy number. The bulk WGS track weights the four site profiles by their cell counts.
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
# The discordance this panel adds to panel b: HMMcopy ranks 2nd by pseudobulk r but 6th by
# per-cell r, because its error is per-cell noise that averages out of a population profile.
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

# One gain/loss colour pair for all eight facets, the bulk WGS track included; the reader compares
# the shape of each caller's track with the bulk facet at top left.
caller <- melt(agg, id.vars = c("method", "x"),
               measure.vars = c("gain", "loss"), variable.name = "type", value.name = "f")
caller[, y := fifelse(type == "loss", -f, f)]
caller[, type := factor(fifelse(type == "gain", "Gain", "Loss"), levels = c("Gain", "Loss"))]

panel_lv <- ord
mk_lab <- function(k) sprintf('bold("%s")', LAB[[k]])
caller[, panel := factor(method, levels = panel_lv)]

bands <- len[band == TRUE, .(xmin = off / 1e6, xmax = (off + m) / 1e6)]
mid   <- len[, .(chr, x = (off + m / 2) / 1e6, w_mb = m / 1e6)][, txt := sub("chr", "", chr)]

## chromosome label thinning from the rendered width ----------------------------
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
    # Strip labels sit directly under the row above; the top margin and the wider inter-row gap
    # keep them clear of that row's axis.
    strip.text         = element_text(size = PT_SUP, hjust = 0,
                                      margin = margin(t = 0.5, b = 1.5)),
    panel.spacing.x    = unit(4, "pt"), panel.spacing.y = unit(7, "pt"),
    axis.text.x        = element_text(size = PT_SMALL, colour = "grey25"),
    axis.text.y        = element_text(size = PT_SMALL, colour = "grey25"),
    axis.title         = element_text(size = BASE_PT),
    axis.ticks         = element_line(colour = "grey70", linewidth = 0.2),
    axis.line          = element_line(colour = "grey70", linewidth = 0.25),
    # No gridlines: the tracks are read as shape against the zero line, which is kept.
    panel.grid         = element_blank(),
    plot.margin        = margin(2, 2, 1, 2)
  )

## direct labels instead of a legend ------------------------------------------
# Gain and Loss are named inside the first facet, each on its own half of the axis: up is always
# gain and down is always loss, so a legend would add nothing.
xmax <- max(agg$x)
# A coloured swatch (shape 15) with black text rather than coloured text; the swatch sits 0.9 mm
# into the facet and the word starts 2.4 mm in.
mm_mb <- xmax / (panel_in * 25.4)   # Mb per mm of facet width
keydt <- data.table(
  panel = factor("bulk WGS", levels = panel_lv),
  x_sw  = xmax * 0.015 + 0.9 * mm_mb,   # over chr1-3, where the bulk track is empty
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
fwrite(stats, file.path(OUT, "fig4c_track_metrics.csv"))     # source data, written in every run

if (!exists("FIG4_ASSEMBLING")) {
ggsave(file.path(OUT, "fig4c_tracks.pdf"), final,
       width = W_MM * mm, height = H_MM * mm, device = cairo_pdf, family = FIG_FAMILY)
cat(sprintf("wrote fig4c_tracks.pdf  (%d x %d mm, base %g pt)\n", W_MM, H_MM, BASE_PT))
}
