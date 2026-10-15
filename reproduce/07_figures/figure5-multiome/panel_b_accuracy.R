#!/usr/bin/env Rscript
# Figure 4b -- pseudobulk Pearson against matched bulk WGS, 20 samples, 8 methods.
#
# Grammar is Figure 3b's: one violin per method at fixed x with direct labels, so
# colour carries nothing position does not already carry and the competitors are
# muted to one grey family; the hero pair keeps the green reserved for our method
# across Figures 2, 3 and 4. Medians sit above their own violin. One significance
# row, referenced to MintCNA, with a bracket only where the comparison is a tie --
# "n.s." alone does not say n.s. against what.
#
# ROW 2 IS A PRECISION-RECALL SCATTER, matching AtaCNV's own Figure 3I/3J:
# methods placed by precision and recall with dashed constant-F contours. A
# single F1 number hides the trade-off that separates these methods -- AtaCNV
# reaches precision 0.92 at recall 0.43 while MintCNA is 0.68 / 0.60, and both
# land near F 0.6. On the plane they are in different places; as one number they
# are not.
#
# PENDING: AtaCNV is currently discretised with the 0.75/1.25 rule, but that rule
# is the one AtaCNV's paper applied to Copy-scAT, inferCNV and CopyKAT -- to its
# OWN output it applies a Bayesian MAP state caller with a copy-neutral rescale.
# Scoring it by someone else's threshold is the same class of error as running
# epiAneufinder at tuned parameters, and its precision/recall will move once
# estimate_cnv_state_cluster() is used instead.
#
# TWO ROWS, two questions. Pearson asks whether the consensus profile has the
# right SHAPE; it is affine-invariant and weights each bin by its squared
# deviation, so a low-amplitude false positive is nearly free. F1 counts BINS, so
# it charges for every one. The two rows do not rank the methods the same way,
# and that disagreement is the point: AtaCNV reaches precision 0.92 at recall
# 0.43 -- it calls almost nothing wrong and misses more than half.
#
# F1 follows AtaCNV's own published procedure, on a truth that needs no threshold
# of ours: GATK's cnv_call. Per cell, "elements larger than 1.25 were defined as
# gains, those less than 0.75 were defined as losses"; then "the common result
# across all tumor cells by taking the mode for each bin"; then bins are counted.
# A bin whose direction disagrees with the truth counts as both FP and FN.
#
# inferCNV is N/A on the F1 row, not zero. That threshold is defined on copy
# RATIOS; inferCNV emits a relative expression score, and its own documentation
# does not call it a copy ratio, so the rule has no meaning there. Reporting the
# 0.000 that the category error produces would be reporting our own mistake as
# the method's result. epiAneufinder IS zero and is drawn: the rule applies to it
# exactly (its per-cell values are -1 / 0 / +0.585, the same three the state
# methods produce), and it reaches zero because at its published default it calls
# 4.0% of bins per cell with at most 38% of nuclei ever agreeing on one bin, so
# no bin reaches a modal call.
#
# NOT faceted by subtype. The oligodendroglioma / astrocytoma split is real and
# large (median r 0.830 against 0.629) but it is a property of OUR method's
# weakness, not of the between-method comparison this panel makes; splitting the
# panel would halve every violin's width to show a second claim. It goes in the
# text and in panel d, which is inherently subtype-resolved.
#
# MintCNA(ATAC) is the ablation, not a competitor. It is drawn because the honest
# result is that it ties: on this modality pair RNA adds no independent
# information, and a figure that omitted the arm would be hiding that.
.libPaths(c("<R_LIBS_USER>", .libPaths()))
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(cowplot)
})
ROOT <- "<MINTCNA_ROOT>"
source(file.path(ROOT, "figure4paper/shared/fig_palette.R"))
OUT  <- file.path(ROOT, "figure4paper/figure5-multiome/output")
DATA <- file.path(ROOT, "figure4paper/figure5-multiome/data")
for (d in c(OUT, DATA)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

tab <- fread(file.path(ROOT, "ATAC_RNA/Result/qc/pearson_vs_truth_own.tsv"))[!is.na(r)]
man <- fread(file.path(ROOT, "ATAC_RNA/Data/meta/sample_manifest.csv")); setkey(man, care_id)
tab[, subtype := fifelse(grepl("^Oligo", man[care_id, subtype]), "Oligodendroglioma", "Astrocytoma")]
stopifnot(uniqueN(tab$care_id) == 20)
fwrite(tab[, .(care_id, subtype, method, pearson = r, n_bin, coverage)],
       file.path(DATA, "fig5b_source_data.csv"))

ORD <- tab[, .(m = median(r)), by = method][order(-m), method]
tab[, method := factor(method, levels = ORD)]
meds <- tab[, .(med = median(r)), by = method]

# Paired Wilcoxon against MintCNA on the same 20 samples, Holm-adjusted. Paired,
# because every method saw the same samples and the between-sample spread (0.37
# to 0.94 for MintCNA alone) is far larger than the between-method differences.
W <- dcast(tab, care_id ~ method, value.var = "r")
ref_m <- "MintCNA"
tst <- rbindlist(lapply(setdiff(ORD, ref_m), function(m) {
  x <- W[[ref_m]]; y <- W[[m]]; k <- is.finite(x) & is.finite(y)
  data.table(method = m, n = sum(k),
             p_raw = wilcox.test(x[k], y[k], paired = TRUE)$p.value,
             delta = median(x[k] - y[k]))
}))
tst[, p_holm := p.adjust(p_raw, "holm")]
tst[, sig := fifelse(p_holm < 1e-3, "***", fifelse(p_holm < 1e-2, "**",
             fifelse(p_holm < 0.05, "*", "n.s.")))]
tst <- rbind(data.table(method = ref_m, n = nrow(W), p_raw = NA_real_, delta = 0,
                        p_holm = NA_real_, sig = "ref"), tst)
tst[, method := factor(method, levels = ORD)]
cat("paired Wilcoxon vs MintCNA (Holm):\n")
print(tst[order(method), .(method, n, delta = round(delta,3),
                           p_holm = signif(p_holm,3), sig)], row.names = FALSE)
fwrite(tst, file.path(DATA, "fig5b_statistics.csv"))

# Equivalence, not just absence of a difference. The panel claims MintCNA is not
# worse than MintCNA(ATAC) or AtaCNV; a non-significant Wilcoxon cannot support
# that, so TOST is computed against a +/-0.05 Pearson margin and written out.
# MARGIN is the largest gap we would call scientifically unimportant; it is stated
# rather than derived, and the minimum detectable difference at this n is reported
# beside it so the reader can see the design's resolution.
MARGIN <- 0.05
eq <- rbindlist(lapply(setdiff(names(W), c("care_id", ref_m)), function(m) {
  x <- W[[ref_m]]; y <- W[[m]]; k <- is.finite(x) & is.finite(y); d <- x[k] - y[k]
  wt <- wilcox.test(x[k], y[k], paired = TRUE, conf.int = TRUE)
  data.table(method = m, n = sum(k), hodges_lehmann = unname(wt$estimate),
             HL_lo = wt$conf.int[1], HL_hi = wt$conf.int[2], margin = MARGIN,
             p_TOST = max(t.test(d + MARGIN, alternative = "greater")$p.value,
                          t.test(d - MARGIN, alternative = "less")$p.value))
}))
eq[, equivalent := p_TOST < 0.05]
eq[, mdd_80pct := (qt(.975, n-1) + qt(.80, n-1)) * sd(W[[ref_m]] - W[["AtaCNV"]], na.rm = TRUE) / sqrt(n)]
cat("\nequivalence vs MintCNA (TOST, margin", MARGIN, "Pearson r):\n")
print(eq[order(p_TOST), .(method, hodges_lehmann = round(hodges_lehmann,4),
                          HL_lo = round(HL_lo,4), HL_hi = round(HL_hi,4),
                          p_TOST = signif(p_TOST,3), equivalent)], row.names = FALSE)
fwrite(eq, file.path(DATA, "fig5b_equivalence.csv"))

## --------------------------------------------------------------------- theme ---
f1 <- fread(file.path(ROOT, "ATAC_RNA/Result/qc/f1_states.tsv"))
f1[, method := factor(method, levels = ORD)]
f1 <- f1[!is.na(method)]
fwrite(f1, file.path(DATA, "fig5b_f1_source_data.csv"))
# POOLED (micro) over the 20 samples, matching f1_states.R -- NOT the median of
# per-sample values. epiAneufinder and inferCNV are undefined in 17 of 20 samples
# because they reach no modal call at all, so their median is a statement about
# how often they are undefined rather than about how well they do; pooling counts
# bins, which is the quantity the panel claims to show. Taking the median here
# also silently disagreed with every F1 number quoted in the text.
pr <- f1[, {
  p <- if (sum(TP)+sum(FP) == 0) NA_real_ else sum(TP)/(sum(TP)+sum(FP))
  r <- if (sum(TP)+sum(FN) == 0) NA_real_ else sum(TP)/(sum(TP)+sum(FN))
  .(precision = p, recall = r,
    F1 = if (is.na(p)||is.na(r)||p+r == 0) NA_real_ else 2*p*r/(p+r))
}, by = method]
cat("\nprecision / recall / F1, POOLED over 20 samples:\n")
print(pr[order(-F1, na.last = TRUE)], row.names = FALSE)
fwrite(pr[order(-F1, na.last = TRUE)], file.path(DATA, "fig5b_pooled_pr.csv"))
stopifnot(nrow(pr[!is.na(precision)]) == length(unique(f1$method)))

# FULL saturation, not the "muted" default. Muting desaturates by 0.55 and
# lightens by 0.30, which is fine under a large violin but leaves panel B's 2 mm
# points impossible to tell apart -- and B now identifies methods by colour alone.
# A uses the same fills so the two panels share one colour key.
FILLS   <- method_fill(as.character(ORD), mode = "full")
lab_col <- label_ink(as.character(ORD))
W_MM <- 66; H_MM <- 93; MM <- 1/25.4

base <- theme_cowplot(BASE_PT, font_family = FIG_FAMILY) +
  theme(axis.text.y  = element_text(size = BASE_PT - 1, colour = "grey20"),
        axis.title   = element_text(size = BASE_PT),
        axis.line    = element_line(linewidth = 0.3, colour = "grey30"),
        axis.ticks   = element_line(linewidth = 0.3, colour = "grey30"),
        panel.grid   = element_blank(),
        plot.margin  = margin(2, 2, 1, 2),
        legend.position = "none")

## row 1 -- Pearson. y starts at the data minimum (-0.278), not at a round -0.35:
## only CONGAS+ and inferCNV reach below zero and the extra band was a quarter of
## the panel height carrying two tails.
y_lo <- -0.30; y_hi <- 1.0; span <- y_hi - y_lo
i_ref <- which(ORD == "MintCNA")
top <- tab[, .(top = max(r, na.rm = TRUE)), by = method]
m2  <- merge(meds, top, by = "method")[, ytxt := pmin(top + 0.035*span, y_hi + 0.015*span)]

p_top <- ggplot(tab, aes(x = method, y = r, fill = method)) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey65", linetype = "22") +
  geom_violin(scale = "width", width = 0.86, colour = NA, alpha = 0.55, na.rm = TRUE) +
  geom_boxplot(width = 0.22, outlier.shape = NA, linewidth = 0.28,
               fill = "white", colour = "grey25", na.rm = TRUE) +
  geom_jitter(width = 0.10, height = 0, size = 0.45, alpha = 0.5,
              colour = "grey30", na.rm = TRUE) +
  geom_text(data = m2, aes(x = method, y = ytxt, label = sprintf("%.2f", med)),
            inherit.aes = FALSE, size = PT_SMALL/.pt, family = FIG_FAMILY,
            fontface = "bold", colour = "grey15") +
  # NO P VALUES IN THIS PANEL, matching Fig 3b.
  #
  # Unlike Fig 3b the test here is sound -- the unit is the SAMPLE and there are
  # 20 independent tumours -- but the claim the panel actually makes about the
  # first three violins is a claim of NO DIFFERENCE, and a large P value cannot
  # establish that. "paired P > 0.99" says only that n = 20 did not detect a
  # difference; read as equivalence it is absence of evidence dressed as
  # evidence of absence. The honest version is an equivalence test, which is
  # computed below and written to source data rather than printed here:
  #   vs MintCNA(ATAC)  Hodges-Lehmann +0.003, 95% CI [-0.010, +0.020], TOST P = 6e-07
  #   vs AtaCNV         Hodges-Lehmann -0.006, 95% CI [-0.034, +0.022], TOST P = 0.009
  # both equivalent within a +/-0.05 Pearson margin, against a minimum detectable
  # difference of 0.045 at n = 20 and 80% power.
  #
  # The per-method stars go with it: the medians printed above each violin are the
  # effect sizes, and 0.80 against 0.11 does not need an asterisk to be read.
  # "ref" goes too -- it named the comparison the tests made, and with no tests
  # drawn it labels something that is not in the panel. Matches Fig 3b.
  scale_fill_manual(values = FILLS, limits = ORD) +
  scale_x_discrete(limits = ORD, labels = method_labels(as.character(ORD))) +
  scale_y_continuous(breaks = seq(-0.25, 1, 0.25)) +
  # italic r, matching Fig 3b: the statistic symbol is italic by convention and
  # the two figures carry the same axis.
  labs(y = expression("Pearson"~italic(r)~"vs bulk WGS"), x = NULL) +
  coord_cartesian(ylim = c(y_lo, y_hi + 0.05*span), clip = "on") + base +
  # Method names ON A'S OWN X AXIS. A shared key docked under A but serving B
  # 120 mm away is the one arrangement none of the six benchmarked papers uses --
  # every shared key in that literature is adjacent to every panel it serves. With
  # A self-decoding, the small key inside B becomes redundant reassurance rather
  # than the sole decoder, which is the CopyKAT pattern: method encoded twice, by
  # x position here and by fill there.
  # 45 degrees, not 90: epiAneufinder and RIDDLER both angle at 45 for narrow
  # panels; only Copy-scAT goes to 90 and it is not worth copying.
  theme(axis.text.x = element_text(size = PT_SUP, angle = 45, hjust = 1, vjust = 1,
                                   colour = lab_col))

## row 2 -- precision-recall plane with constant-F contours (AtaCNV Fig 3I/J)
fc <- rbindlist(lapply(c(0.2, 0.4, 0.6, 0.8), function(f) {
  x <- seq(f/2 + 1e-4, 1, length.out = 200)      # recall implied by precision at fixed F
  y <- f*x / (2*x - f)
  data.table(F = f, precision = x, recall = y)[recall <= 1 & recall >= 0]
}))
# Contour labels return to the RIGHT edge. They were moved to the top when the
# panel carried direct method labels that collided with them; the labels are gone
# now, and the top-left is wanted for the method key.
lab_fc <- fc[, .SD[which.min(abs(precision - 0.97))], by = F]

# NO aspect.ratio here. Both axes run 0-1 so a square panel is the natural shape,
# but aspect.ratio forces the panel square INSIDE its slot: at a 63.0 x 49.2 mm
# slot the drawn panel became 49.2 x 49.2 with ~7 mm of white on each side, so B
# read as visibly narrower than D and E directly below it in the same column.
# Column alignment wins -- the axes still both span 0-1 with matching ticks, and
# the constant-F contours are guides, not data. Same call as panel E, which gave
# up coord_fixed() for the same reason.
p_bot <- ggplot() +
  geom_path(data = fc, aes(precision, recall, group = F),
            linewidth = 0.22, linetype = "22", colour = "grey72") +
  # PT_SMALL (5.5 pt), not 4.6 pt. Every rendered glyph has to clear the 5 pt
  # floor; these contour labels were the only text in the figure below it.
  geom_text(data = lab_fc, aes(precision, recall, label = sprintf("F=%.1f", F)),
            size = PT_SMALL/.pt, colour = "grey55", hjust = 1, vjust = -0.35) +
  # 2.4, up from 2.1: the panel grew from 42.3 to 46.8 mm and the marks did not.
  # MintCNA and MintCNA^ATAC sit 0.014 recall apart -- 0.65 mm here -- so they
  # overlap at any readable size. That is the honest picture (panel A's TOST says
  # they are equivalent); the white stroke is what keeps them as two marks.
  # Diamonds for our two marks, circles for the competitors (advisor request,
  # 2026-09-21): colour already separates them, the shape makes the pair findable
  # without the key. A diamond of the same nominal size draws smaller than a circle,
  # so the heroes get 2.9 to the circles' 2.4. shape / fill / size are all mapped
  # to `method` with identical limits and labels, so ggplot merges them into ONE key.
  # MintCNA drawn LAST so the paper's own method is the mark on top where the two
  # diamonds overlap (they are 0.014 recall apart).
  geom_point(data = pr[order(method == "MintCNA")], aes(precision, recall, fill = method, shape = method, size = method),
             stroke = 0.35, colour = "white") +
  scale_shape_manual(values = setNames(ifelse(as.character(ORD) %in% names(HERO), 23, 21), as.character(ORD)),
                     limits = ORD, name = NULL, labels = method_labels(as.character(ORD))) +
  scale_size_manual(values = setNames(ifelse(as.character(ORD) %in% names(HERO), 2.9, 2.4), as.character(ORD)),
                    limits = ORD, name = NULL, labels = method_labels(as.character(ORD))) +
  # A KEY, not eight labels on the points. At a 42 mm square the repelled labels
  # piled onto each other and onto the data; AtaCNV Fig 3G-J solves the same
  # problem the same way. Two columns fit the top-left, which is empty because
  # every method sits at recall < 0.6 and the high-precision ones are on the right.
  scale_fill_manual(values = FILLS, limits = ORD, name = NULL,
                    labels = method_labels(as.character(ORD)),
                    # ONE column, not two. The labels are plotmath (MintCNA^ATAC), so the text has
                                         # to sit at PT_SUP or the superscript renders at
                                         # 3.85 pt; at that size two columns need 46 mm and
                                         # the square is 42.3 mm. One column fits the empty
                                         # top-left, which spans precision < 0.54 above
                                         # recall 0.34 -- no method is in there.
                                         # the key was SMALLER than the mark it
                                         # stands for (1.6 vs 2.1); a legend key
                                         # must never be harder to see than the data
                                         guide = guide_legend(ncol = 1,
                                         override.aes = list(stroke = 0.35))) +
  scale_x_continuous(limits = c(-0.03, 1.02), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-0.03, 1.02), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
  labs(x = "Precision", y = "Recall") + base +
  # The key is BOXED. It sits inside the plot area, so without a rule it reads as
  # eight loose labels floating in the data region; a hairline box and a white fill
  # make it one object that is plainly not data. Inset from the corner rather than
  # flush against it, so the box does not touch the axis lines. Row spacing and key
  # size both go up -- at 2.6 mm and 0.1 mm the eight entries ran together.
  theme(legend.position = "inside", legend.position.inside = c(0.025, 0.975),
        legend.justification = c(0, 1), legend.direction = "vertical",
        legend.background = element_rect(fill = "white", colour = "grey72",
                                         linewidth = 0.25),
        legend.key = element_blank(),
        legend.margin = margin(1.2, 2.6, 1.2, 1.2), legend.spacing.y = unit(0.35, "mm"),
        legend.text = element_text(size = PT_SUP, colour = "grey20"),
        legend.key.size = unit(3.1, "mm")) +

  # SQUARE. Both axes are the same quantity on the same 0-1 scale, so a
  # rectangular panel stretches one of them and invents a visual difference the
  # data does not contain. It costs width -- a square plot in this column is
  # 34.8 mm wide against D's and F's 60.9 mm -- but a distorted metric plane is
  # a worse defect than an uneven column. The constant-F contours are only
  # correct curves at 1:1 as well.
  theme(aspect.ratio = 1)

g_b <- plot_grid(p_top, p_bot, ncol = 1, rel_heights = c(1, 1.05))

if (!exists("FIG5_ASSEMBLING")) {
  cairo_pdf(file.path(OUT, "fig5b_accuracy.pdf"), width = W_MM*MM, height = H_MM*MM,
            family = FIG_FAMILY); print(g_b); dev.off()
  cat(sprintf("wrote fig5b_accuracy.pdf  (%d x %d mm)\n", W_MM, H_MM))
}
