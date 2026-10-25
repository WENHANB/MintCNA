#!/usr/bin/env Rscript
# Supplementary Fig. 6 -- share of CRC01 cells carrying a chromosome 17 arm-level loss
# (a: 17p, b: 17q) at each sampling site, from the MintCNA calls.
# Input:  reproduce/figures/data/fig4_crc/MultiCNA_AE_197_d250.RData (cna_multicna_ae: copy-number
#         state per bin and cell, 1 = loss; ref_multicna_ae: hg19 bin coordinates)
# Output: reproduce/figures/output/supp/SuppFig06_CRC_chr17_by_site.pdf (88 x 48 mm), a PNG preview
#         and figS6_chr17_arm_loss_by_site.csv (counts per arm and site, with Wilson 95% intervals
#         that are written but not drawn)
# Run from the root of the MintCNA repository: Rscript reproduce/figures/supp/figS6_chr17_arm_loss.R

suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(cowplot) })
FIG <- "reproduce/figures"   # run from the root of the MintCNA repository
if (!file.exists(file.path(FIG, "shared/fig_palette.R"))) stop("run this script from the root of the MintCNA repository")
source(file.path(FIG, "shared/fig_palette.R"))
DAT <- file.path(FIG, "data/fig4_crc")
OUT <- file.path(FIG, "output/supp"); QA <- OUT; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

load(file.path(DAT, "MultiCNA_AE_197_d250.RData"))   # cna_multicna_ae, ref_multicna_ae
EXCL <- "CRC01_ML3_085"   # per-cell QC exclusion, as in fig4_crc/utils.R
if (length(EXCL)) cna_multicna_ae <- cna_multicna_ae[, !colnames(cna_multicna_ae) %in% EXCL, drop = FALSE]
# hg19 17p = chr17 bins with start <= 25,248,165, 17q = bins with end >= 25,248,166; a cell carries
# an arm loss when more than 50% of that arm's bins are called loss (state 1)
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
fwrite(tab, file.path(OUT, "figS6_chr17_arm_loss_by_site.csv"))

mk <- function(a) {
  d <- tab[arm == a]
  ggplot(d, aes(site, pct, fill = site)) +
    geom_col(width = 0.68) +
    geom_text(aes(y = pct, label = sprintf("%.1f%%", pct)), vjust = -0.5, size = PT_SMALL / .pt,
              family = FIG_FAMILY, colour = "grey10") +
    scale_fill_manual(values = SITE_COLORS2, guide = "none") +   # site colours shared with the CRC01 site strip of Fig. 4
    scale_x_discrete(labels = setNames(d$site_lab, as.character(d$site))) +
    scale_y_continuous(limits = c(0, 112), breaks = seq(0, 100, 25), expand = c(0, 0)) +
    labs(x = "Tumor site", y = sprintf("Cells with %s loss (%%)", a)) +   # US spelling, as in the manuscript
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
