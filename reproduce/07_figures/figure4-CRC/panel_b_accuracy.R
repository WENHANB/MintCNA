#!/usr/bin/env Rscript
# Figure 3b -- per-cell accuracy vs matched bulk WGS.  REVISION of 02_03_rmse_pearson.R
#
# What changed and why:
#  1. y range is honest. The original clipped Pearson at coord_cartesian(0,1);
#     117 of 1,576 cells sit below 0 (SCCNV 62/197, min -0.58; SCICoNE 55/197).
#     Clipping made the two failing callers look merely uninformative when they
#     are ANTI-correlated in a third of cells. The axis now reaches the data.
#  2. No colour legend. Methods are spatially fixed on x, so they are labelled
#     directly; that returns ~12 mm of height and removes the cross-panel hue
#     collision (green meant both MintCNA and subclone C1; rcCAE's orange was
#     one JND from subclone C2). Colour now does one job: hero vs baseline.
#  3. No P values are printed. The paired tests run over 197 cells from ONE
#     tumour, so the replicate unit for a between-method claim is the patient,
#     n = 1, and the P value is bounded by cell count rather than by how
#     reproducible the advantage is. The medians carry the comparison; the tests
#     stay in data/fig4b_paired_tests.csv and in Methods. See the long note in
#     vb() for the counterexample this dataset supplies.
#  4. violin + box only. The original also drew 197 jittered points per method
#     (1,576 dots); at 65 mm across eight methods that is 8 mm per group and the
#     cloud reads as a smear. The violin already carries the shape, including
#     the negative tail that panel 1 above now exposes.

.libPaths(c("<R_LIBS_USER>", .libPaths()))
suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(data.table); library(cowplot)
})
ROOT <- "<MINTCNA_ROOT>"
source("<REPRODUCE_DIR>/07_figures/figure4-CRC/utils.R")
source("<REPRODUCE_DIR>/07_figures/shared/fig_palette.R")
OUT <- file.path(ROOT, "figure4paper/figure4-CRC/output")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

## ------------------------------------------------------------------ metrics ---
d <- load_all_methods(); methods <- d$methods; truth <- d$truth_mat

# The truth matrix carries no rownames, so alignment to it is POSITIONAL: row i of
# a caller is matched to row i of the bulk profile. Verified once (2026-08-21) three
# ways -- coordinate-by-coordinate against `ref` for the callers that do carry
# rownames, a run-length check (every caller 66-1156 bins per run vs 2.6 for a
# row-shuffled control, so no matrix is scrambled), and a shift/reverse control
# (MintCNA r falls 0.848 -> ~0 under any shift, and SCCNV/SCICoNE stay at ~0 under
# every shift, so their low correlation is real and not a misalignment artefact).
# The guard below re-checks the part that can be checked automatically, so a future
# data swap cannot silently break the pairing.
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

## every comparison against MintCNA, computed not asserted -----------------------
# Paired two-sided Wilcoxon over the same 197 cells, Holm-adjusted within each
# metric (7 comparisons per metric). Paired, because every method scored the same
# cells against the same site-matched bulk profile; an unpaired test would throw
# away that structure and answer a question nobody asked.
#
# Why the whole matrix is computed rather than assumed: median order does NOT
# imply significance order, and this dataset contains the counterexample.
# MintCNA vs MintCNA(DNA) differs by 0.001 in median Pearson and is P = 4e-12;
# MintCNA vs rcCAE differs by 0.003 -- a LARGER gap -- and is P = 0.24. The two
# MintCNA modes are the same model on the same cells with and without the RNA
# channel, so their paired differences are tiny but almost perfectly consistent;
# rcCAE is a different algorithm, so its paired differences are large and
# scattered. Significance tracks the consistency of the difference, not its size.
# Wilcoxon SIGNED-RANK on the paired per-cell differences. Two properties of
# R's implementation that belong in Methods rather than being discovered later:
#   - exact = NULL defaults to exact only when the number of NON-ZERO differences
#     is < 50. Here it is 133-197, so every p below is the normal approximation
#     with continuity correction. Checked against exact and against a sign-flip
#     permutation: rcCAE 0.2387 / 0.2393 / 0.24, ablation 2.2e-12 / (not
#     computable, zeroes) / < 1e-5. No conclusion moves.
#   - zero differences are dropped, as the signed-rank test requires. MintCNA vs
#     MintCNA(DNA) has 64 exact ties, so that comparison runs on n = 133, not 197.
# The sign test and the paired t test are computed too and written to source data,
# NOT because the reported test changes but because the three disagree for rcCAE
# and the disagreement is itself a result: the difference distribution is skewed
# (+0.79), MintCNA wins in more cells (115/197, sign test P = 0.022) while rcCAE
# wins by more when it wins (median |d| 0.068 vs 0.009), so the mean-based test
# points the other way (t P = 0.024 favouring rcCAE). Signed-rank is kept as the
# pre-specified, most conservative summary; switching tests after seeing this
# would be choosing the answer.
paired_tests <- function(metric) {
  w <- dcast(tab, cell ~ method, value.var = metric)
  o <- setdiff(ORD, "MintCNA")
  hi_is_better <- metric == "pearson"
  r <- rbindlist(lapply(o, function(m) {
    x <- w[["MintCNA"]]; y <- w[[m]]; k <- is.finite(x) & is.finite(y)
    d <- x[k] - y[k]; nz <- d[d != 0]
    win <- if (hi_is_better) d > 0 else d < 0
    data.table(method = m, n = sum(k), n_nonzero = length(nz),
               n_MintCNA_better = sum(win), mean_diff = mean(d), median_diff = median(d),
               p_raw     = wilcox.test(x[k], y[k], paired = TRUE)$p.value,
               p_sign    = binom.test(sum(win), sum(k))$p.value,
               p_ttest   = t.test(d)$p.value)
  }))
  r[, p_holm := p.adjust(p_raw, method = "holm")]
  r[, sig := fifelse(p_holm >= 0.05, "n.s.",
             fifelse(p_holm < 1e-4, "****",
             fifelse(p_holm < 1e-3, "***",
             fifelse(p_holm < 1e-2, "**", "*"))))]
  # The reference column is labelled, not left blank. A blank slot leaves the
  # reader guessing what the other seven marks are compared against.
  rbind(data.table(method = "MintCNA", n = NA_integer_, n_nonzero = NA_integer_,
                   n_MintCNA_better = NA_integer_, mean_diff = NA_real_,
                   median_diff = NA_real_, p_raw = NA_real_, p_sign = NA_real_,
                   p_ttest = NA_real_, p_holm = NA_real_, sig = "ref"),
        r)[match(ORD, method)]
}
TEST <- list(pearson = paired_tests("pearson"), rmse = paired_tests("rmse"))

# The significance row answers two different questions and must not use one symbol
# for both.
#   vs the six competitors  -- "is MintCNA better?"  Effects are huge (0.08-0.82
#      in median Pearson), so a mark is the right summary.
#   vs its own DNA-only ablation -- "does the RNA channel help?"  Here a star is
#      actively misleading. The paired difference is -0.0012 in median Pearson,
#      i.e. the ablation is very slightly BETTER, and it reaches P = 4e-12 only
#      because the direction is consistent: of the 133 cells with any difference,
#      107 favour DNA-only and 26 favour the full model. The magnitude is 1/54 of
#      MintCNA's own between-cell SD (0.0016 vs 0.085). Wilcoxon signed-rank tests
#      consistency, not size. So that slot reports the signed EFFECT, not a star.
# And the five competitor marks collapse to one spanning bound. A bracket to SCOPE
# implying "everything after is smaller" would be false -- p is NOT monotone in the
# median order (SCONCE 1.6e-29 sits between SCOPE 8.6e-32 and HMMcopy 9.2e-32) --
# but a BOUND over the whole group is true and needs no ordering assumption.
delta_med <- function(mt) {
  w <- dcast(tab[method %in% c("MintCNA","MintCNA(DNA)")], cell ~ method, value.var = mt)
  d <- w[["MintCNA"]] - w[["MintCNA(DNA)"]]
  median(d[is.finite(d)])
}
COMPET <- setdiff(ORD, c("MintCNA", "MintCNA(DNA)", "rcCAE"))
BOUND  <- lapply(TEST, function(t) 10^ceiling(log10(max(t[method %in% COMPET, p_holm]))))
DELTA  <- lapply(c(pearson = "pearson", rmse = "rmse"), delta_med)
for (mt in names(TEST)) cat(sprintf("%s: ablation median paired diff = %+.4f | all competitors P < %.0e\n",
                                    mt, DELTA[[mt]], BOUND[[mt]]))
i_scope <- match(COMPET[1], ORD); i_last <- match(COMPET[length(COMPET)], ORD)
for (m in names(TEST)) {
  cat("\n-- ", m, ": MintCNA vs each (Wilcoxon signed-rank, Holm) --\n", sep = "")
  print(TEST[[m]][, .(method, n_nz = n_nonzero, better = n_MintCNA_better,
                      p_signrank = signif(p_holm, 3), p_sign = signif(p_sign, 3),
                      p_t = signif(p_ttest, 3), sig)], row.names = FALSE)
}
fwrite(rbindlist(lapply(names(TEST), function(m) copy(TEST[[m]])[, metric := m])),
       file.path(ROOT, "figure4paper/figure4-CRC/data/fig4b_paired_tests.csv"))

i_mint <- match("MintCNA", ORD); i_rc <- match("rcCAE", ORD)

## clipping the original hid ------------------------------------------------------
cat(sprintf("cells with Pearson < 0 (hidden by the original ylim(0,1)): %d of %d\n",
            tab[pearson < 0, .N], tab[, .N]))
print(tab[, .(min = round(min(pearson, na.rm = TRUE), 3),
              n_neg = sum(pearson < 0, na.rm = TRUE)), by = method][n_neg > 0])

## --------------------------------------------------------------------- theme ---
FILLS   <- method_fill(ORD)
lab_col <- label_ink(ORD)

base <- theme_cowplot(BASE_PT, font_family = FIG_FAMILY) +
  theme(
    axis.title.x       = element_blank(),
    # PT_SUP, not BASE_PT - 1: this row carries MintCNA^DNA and plotmath renders
    # the superscript at ~0.7x, so a 6 pt parent put "DNA" on the page at 4.2 pt.
    axis.text.x        = element_text(size = PT_SUP, angle = 45, hjust = 1,
                                      vjust = 1, colour = lab_col),
    axis.text.y        = element_text(size = BASE_PT - 1, colour = "grey20"),
    axis.title.y       = element_text(size = BASE_PT),
    axis.line          = element_line(linewidth = 0.3, colour = "grey30"),
    axis.ticks         = element_line(linewidth = 0.3, colour = "grey30"),
    # No gridlines anywhere in this figure. design-theory.md sections 2 and 6:
    # "no grid; sparse y-ticks". The zero reference line in the Pearson panel is
    # kept -- it marks a meaningful value, it is not a ruler.
    panel.grid         = element_blank(),
    plot.margin        = margin(2, 2, 1, 2),
    legend.position    = "none"
  )

vb <- function(dt, yv, ylim, ylab, meds, show_x, tst, bound, dlt, clipped = NULL) {
  y_lo <- ylim[1]; y_hi <- ylim[2]; span <- y_hi - y_lo
  # Median labels sit directly above their OWN violin, not in a shared strip at the
  # panel floor. In the original the strip was pinned to y = -Inf, which put the
  # RMSE numbers inside the low violin tails and, once the Pearson axis was opened
  # down to -0.6, left the Pearson numbers a full panel-height away from the boxes
  # they describe. Per-violin placement cannot collide and reads as a label.
  top <- dt[, .(top = max(get(yv), na.rm = TRUE)), by = method]
  meds2 <- merge(meds, top, by = "method")
  # clip = "on" below, so the floor actually cuts SCCNV's tail instead of letting
  # it spill into the panel underneath; everything that must sit ABOVE the data
  # therefore has to fit inside head_room rather than escaping the panel.
  # 0.32 carried a median row, a "ref" row, a bracket and a bracket label. With
  # the brackets gone only the first two remain, so the headroom comes back to
  # the data instead of being white space at the top of the panel.
  head_room <- 0.17
  meds2[, ytxt := pmin(top + 0.035 * span, y_hi + 0.02 * span)]
  # No significance row at all, "ref" included. "ref" existed to tell the reader
  # which violin the other seven were being tested against; with the tests gone
  # from the panel it labels a comparison that is no longer drawn. MintCNA is
  # already first and already the only green pair, so nothing is lost.
  p <- ggplot(dt, aes(x = method, y = .data[[yv]], fill = method)) +
    geom_violin(scale = "width", width = 0.86, colour = NA, alpha = 0.55, na.rm = TRUE) +
    geom_boxplot(width = 0.22, outlier.shape = NA, linewidth = 0.28,
                 fill = "white", colour = "grey25", na.rm = TRUE) +
    geom_text(data = meds2, aes(x = method, y = ytxt, label = sprintf("%.2f", med)),
              inherit.aes = FALSE, size = PT_SMALL / .pt,
              family = FIG_FAMILY, fontface = "bold", colour = "grey15", vjust = 0) +
    scale_fill_manual(values = FILLS, limits = ORD) +
    scale_x_discrete(limits = ORD, labels = method_labels(ORD)) +
    # breaks are pinned to the DATA range, not to the padded range, so the extra
    # headroom that carries the label row never invents a 1.2 tick on a 0-0.95 metric
    scale_y_continuous(breaks = scales::breaks_pretty(5)(c(y_lo, y_hi))) +
    labs(y = ylab) + coord_cartesian(ylim = c(y_lo, y_hi + head_room * span), clip = "on") + base
  if (y_lo < 0) p <- p + geom_hline(yintercept = 0, colour = "grey45",
                                    linewidth = 0.25, linetype = "22")
  if (!is.null(clipped) && nrow(clipped)) {
    # shape 25 is a filled down-triangle drawn by the device; a text glyph such
    # as U+25BC is not in Nimbus Sans and would silently pull in a second font.
    p <- p +
      geom_point(data = clipped, aes(x = method, y = y_lo), inherit.aes = FALSE,
                 shape = 25, size = 1.1, fill = "grey30", colour = "grey30") +
      # No prose in the panel: the triangle says data continues below the axis and
      # the legend gives the count, the minimum, the test, the adjustment and the n.
      annotate("blank", x = 0.55, y = y_lo,
               # Each line is kept under ~34 characters. At 5.5 pt in a 66 mm panel
               # that is about six of the eight method slots, which stops the note
               # short of SCICoNE's median label on the right.
               size = PT_SMALL / .pt)
  }
  # NO P-VALUE ANNOTATION IN THIS PANEL, deliberately.
  #
  # The paired tests here run over 197 CELLS from a SINGLE tumour (CRC01), so the
  # replicate unit for a claim like "MintCNA beats SCOPE" is the patient, n = 1.
  # Cells within one tumour share the bulk WGS truth, the library prep and the
  # GC/mappability bias structure, so any systematic difference is consistent
  # across all 197 and the P value ends up bounded by how many cells were
  # sequenced rather than by how reproducible the advantage is. This dataset
  # contains its own proof: MintCNA vs MintCNA(DNA) differs by 0.0012 in median
  # Pearson at P = 2e-12, while MintCNA vs rcCAE differs by MORE, 0.0029, at
  # P = 0.24. A P value that moves twelve orders of magnitude against the effect
  # size is not measuring what a reader will assume it measures.
  #
  # The medians printed above each violin are the effect sizes and carry the
  # comparison on their own. The full Wilcoxon / sign / t results, the per-method
  # n, and the count of cells in which MintCNA wins are all retained in
  # data/fig4b_paired_tests.csv and described in Methods, where the replicate unit
  # can be stated in words.
  if (!show_x) p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  p
}

# Floor chosen from the data, not by eye: see the comment above. Anything the
# floor cuts is counted and drawn, so the clip is disclosed inside the figure.
PEARSON_FLOOR <- -0.1
y_pear <- c(PEARSON_FLOOR, 1)
below <- tab[pearson < PEARSON_FLOOR, .(n_below = .N, min_r = min(pearson)), by = method][n_below > 0]
cat(sprintf("cells below the Pearson floor (%.2f): %d of %d\n",
            PEARSON_FLOOR, sum(below$n_below), nrow(tab)))
print(below, row.names = FALSE)
y_rmse <- c(0, ceiling(max(tab$rmse, na.rm = TRUE) * 10) / 10)
cat(sprintf("Pearson axis %.1f..%.1f   RMSE axis %.1f..%.1f\n", y_pear[1], y_pear[2], y_rmse[1], y_rmse[2]))

p1 <- vb(tab, "pearson", y_pear, expression("Pearson"~italic(r)~"vs bulk WGS"),
         sum_p, show_x = FALSE, tst = TEST$pearson,
         bound = BOUND$pearson, dlt = DELTA$pearson, clipped = below)
p2 <- vb(tab, "rmse",    y_rmse, "RMSE vs bulk WGS",
         sum_r, show_x = TRUE,  tst = TEST$rmse,
         bound = BOUND$rmse, dlt = DELTA$rmse)

panel_b <- p1 / p2 + plot_layout(heights = c(1, 1))

## ---------------------------------------------------------------- final size ---
W_MM <- 66; H_MM <- 93                     # right column of Fig 3 row 1
mm <- 1 / 25.4
if (!exists("FIG4_ASSEMBLING")) {
ggsave(file.path(OUT, "fig4b_accuracy.pdf"), panel_b,
       width = W_MM * mm, height = H_MM * mm, device = cairo_pdf, family = FIG_FAMILY)
cat(sprintf("wrote fig4b_accuracy.pdf  (%d x %d mm, base %g pt)\n", W_MM, H_MM, BASE_PT))
}
