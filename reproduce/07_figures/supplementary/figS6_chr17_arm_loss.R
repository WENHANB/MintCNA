#!/usr/bin/env Rscript
# Supplementary Fig. 6 -- share of CRC01 cells carrying a chromosome-17 arm-level
# loss (17p, 17q) at each sampling site, MintCNA calls, 197 cells.
#
# REBUILT 2026-10-09 FOR NATURE COMMUNICATIONS. The published S6 was two
# 5.5 x 4.5 in ggplots from Code/CRC/Evaluation/chr17.R (2026-03-20) pasted on a
# 250 x 160 mm page with 14-16 pt text, Helvetica + Verdana-Bold + Liberation Sans.
# Same data, same rule, redrawn at 88 x 48 mm (single column) in the shared
# Helvetica/Nimbus family, 5.5-8 pt, lowercase letters, site colours taken from
# the Fig. 3a site strip so PT/LN/ML/MP mean the same colour in both figures.
# n per site is shown under the site label (LN has 7 cells). Wilson 95% intervals were
# drawn 2026-10-09 and removed at the author's request 2026-10-12 (they stay in the CSV).
#
# Rule (chr17.R): hg19 17p = chr17 bins with start <= 25,248,165; 17q = bins with
# end >= 25,248,166. A cell carries an arm loss when > 50 % of that arm's bins are
# called loss (state 1). Numbers reproduce the published 62.7/85.7/98.2/98.7 (17p)
# and 76.3/100/98.2/98.7 (17q).

.libPaths(c("<R_LIBS_USER>", .libPaths()))
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(cowplot) })
ROOT <- "<MINTCNA_ROOT>"
source("<REPRODUCE_DIR>/07_figures/shared/fig_palette.R")
OUT <- Sys.getenv("FIGS6_OUT", file.path(ROOT, "figure4paper/figureS-supplementary"))
QA  <- Sys.getenv("FIGS6_QA",  file.path(ROOT, "figure4paper/figureS-supplementary/_work/qa"))
DAT <- file.path(ROOT, "figure4paper/figureS-supplementary/_work/data")
for (d in c(OUT, QA, DAT)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

SFX <- Sys.getenv("CRC_SUFFIX", "_d250")   # default "_d250" = 2026-10-14 re-run (AG-DCAE latent_dim 250, 197 cells); "" = the pre-2026-10-14 files
load(file.path(ROOT, paste0("Result/CRC/Benchmark/MultiCNA_AE_197", SFX, ".RData")))   # cna_multicna_ae, ref_multicna_ae
EXCL <- strsplit(Sys.getenv("CRC_EXCLUDE", "CRC01_ML3_085"), ",")[[1]]; EXCL <- EXCL[nzchar(EXCL)]   # per-cell QC exclusion, see Result/CRC/figures/scripts/utils.R
if (length(EXCL)) cna_multicna_ae <- cna_multicna_ae[, !colnames(cna_multicna_ae) %in% EXCL, drop = FALSE]
P_END <- 25248165L
SITES <- c("PT", "LN", "ML", "MP")
site <- regmatches(colnames(cna_multicna_ae), regexpr("LN|PT|ML|MP", colnames(cna_multicna_ae)))
stopifnot(length(site) == ncol(cna_multicna_ae))
arm_idx <- list(`17p` = ref_multicna_ae$chr == "chr17" & ref_multicna_ae$start <= P_END,
                `17q` = ref_multicna_ae$chr == "chr17" & ref_multicna_ae$end   >= P_END + 1L)
wilson <- function(k, n, z = qnorm(0.975)) {
  p <- k / n; den <- 1 + z^2 / n; ctr <- (p + z^2 / (2 * n)) / den
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / den
  data.table(lo = pmax(0, ctr - half), hi = pmin(1, ctr + half))
}
tab <- rbindlist(lapply(names(arm_idx), function(a) {
  loss <- colMeans(cna_multicna_ae[arm_idx[[a]], , drop = FALSE] == 1, na.rm = TRUE) > 0.5
  d <- data.table(arm = a, site = factor(site, levels = SITES), loss = loss)[, .(k = sum(loss), n = .N), by = .(arm, site)]
  cbind(d, wilson(d$k, d$n))
}))
tab[, pct := 100 * k / n]
setorder(tab, arm, site)
tab[, site_lab := sprintf("%s\n(n = %d)", site, n)]
print(tab[, .(arm, site, k, n, pct = round(pct, 1), lo = round(100 * lo, 1), hi = round(100 * hi, 1))], row.names = FALSE)
fwrite(tab, file.path(DAT, "figS6_chr17_arm_loss_by_site.csv"))

mk <- function(a) {
  d <- tab[arm == a]
  ggplot(d, aes(site, pct, fill = site)) +
    geom_col(width = 0.68) +
    geom_text(aes(y = pct, label = sprintf("%.1f%%", pct)), vjust = -0.5, size = PT_SMALL / .pt,
              family = FIG_FAMILY, colour = "grey10") +
    scale_fill_manual(values = SITE_COLORS2, guide = "none") +
    scale_x_discrete(labels = setNames(d$site_lab, as.character(d$site))) +
    scale_y_continuous(limits = c(0, 112), breaks = seq(0, 100, 25), expand = c(0, 0)) +
    labs(x = "Tumor site", y = sprintf("Cells with %s loss (%%)", a)) +   # US spelling, as the manuscript
    theme_cowplot(BASE_PT, font_family = FIG_FAMILY) +
    theme(axis.text  = element_text(size = PT_SMALL, colour = "grey20"),
          axis.title = element_text(size = BASE_PT),
          axis.line  = element_line(linewidth = 0.3, colour = "grey30"),
          axis.ticks = element_line(linewidth = 0.3, colour = "grey30"),
          axis.ticks.length = unit(1.5, "pt"),
          plot.margin = margin(4, 3, 1, 2))
}
fig <- plot_grid(mk("17p"), mk("17q"), ncol = 2, labels = c("a", "b"), label_size = 8,
                 label_fontface = "bold", label_fontfamily = FIG_FAMILY, label_x = 0, label_y = 1,
                 hjust = 0, vjust = 1)
W_MM <- 88; H_MM <- 48; mm <- 1 / 25.4
base <- "SuppFig06_CRC_chr17_by_site"
cairo_pdf(file.path(OUT, paste0(base, ".pdf")), width = W_MM * mm, height = H_MM * mm, family = FIG_FAMILY)
print(fig); dev.off()
suppressWarnings(pdftools::pdf_convert(file.path(OUT, paste0(base, ".pdf")), format = "png", dpi = 300,
                                       filenames = file.path(QA, paste0(base, "_preview.png")), verbose = FALSE))
cat(sprintf("wrote %s.pdf  (%d x %d mm)\n", base, W_MM, H_MM))
