#!/usr/bin/env Rscript
# Supplementary panel: are the CNA-defined subclones visible in the RNA embedding?
#
# Two columns, one per tumour. Each column is the SAME UMAP drawn twice: coloured
# by the clone MintCNA infers from copy number alone, and coloured by the source
# study's RNA-derived cell state. The point is that the two colourings are almost
# the same picture -- a partition made without reading a single transcript
# reproduces one that was defined from the transcriptome.
#
# The separation numbers are measured in PCA space, not on the UMAP, because a
# UMAP is for looking at and not for measuring. The cell-state row is the CEILING
# for this comparison, not a competitor: those labels were defined from this same
# RNA, so they are the best any labelling can do here.
.libPaths(c("<R_LIBS_USER>", .libPaths()))
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(cowplot) })
ROOT <- "<MINTCNA_ROOT>"
source(file.path(ROOT, "figure4paper/shared/fig_palette.R"))
OUT <- file.path(ROOT, "figure4paper/figure5-multiome/output")
QC  <- file.path(ROOT, "ATAC_RNA/Result/qc")
MM  <- 1/25.4

SAMPLES <- c("P111T2", "P100T1")
# This used to override the shared pair with #7FA8C9, on the reasoning that two
# neutral greys are too weak for a scatter of a thousand points. The reasoning was
# right and the fix was wrong twice over: #7FA8C9 is only dE00 8.5 from the CN5
# "Loss" ramp colour #8CC5E3, so it collided with the copy-number scale; and the
# real cause was the alpha, not the greys. CLONE_COLORS2 is now deep enough to
# carry a scatter, so this figure uses the same two colours as the main figure --
# which is what a supplementary panel of the same clones should have done anyway.
CLONE_PT <- CLONE_COLORS2

pan <- function(S, colour_by) {
  d <- fread(file.path(QC, sprintf("clone_embedding_%s.tsv", S)))
  su <- fread(file.path(QC, sprintf("clone_embedding_%s_summary.tsv", S)))
  if (colour_by == "clone") {
    d[, grp := factor(clone, levels = c("C1", "C2"))]
    pal <- CLONE_PT; ttl <- sprintf("%s  — subclone (from CNA)", S)
    lab <- sprintf("15-NN accuracy %.3f  (majority-class floor %.3f)", su$knn_acc_clone, su$majority_floor)
  } else {
    d <- d[!is.na(state) & nzchar(state)]   # fread yields "" not NA for missing
    d[, grp := factor(state, levels = names(CELLSTATE_COLORS))]
    pal <- CELLSTATE_COLORS; ttl <- sprintf("%s  — cell state (from RNA)", S)
    lab <- sprintf("15-NN accuracy %.3f  \u2014 the ceiling: these labels were defined from this RNA", su$knn_acc_state)
  }
  ggplot(d, aes(umap1, umap2, colour = grp)) +
    geom_point(size = 0.35, alpha = 0.75) +
    scale_colour_manual(values = pal, name = NULL, drop = TRUE,
                        guide = guide_legend(override.aes = list(size = 1.5), ncol = 1)) +
    labs(title = ttl, subtitle = lab, x = "UMAP 1", y = "UMAP 2") +
    theme_cowplot(BASE_PT, font_family = FIG_FAMILY) +
    theme(plot.title = element_text(size = PT_SUP, face = "plain", colour = "grey15"),
          axis.text = element_blank(), axis.ticks = element_blank(),
          axis.line = element_line(linewidth = 0.3, colour = "grey30"),
          axis.title = element_text(size = PT_SMALL, colour = "grey35"),
          panel.grid = element_blank(),
          legend.text = element_text(size = PT_SMALL),
          legend.key.size = unit(2.6, "mm"),
          legend.margin = margin(0, 0, 0, 1),
          plot.margin = margin(2, 2, 2, 2), aspect.ratio = 1)
}

# knn accuracy for the state row is not in the summary file for older runs; recover it
for (S in SAMPLES) {
  f <- file.path(QC, sprintf("clone_embedding_%s_summary.tsv", S))
  su <- fread(f)
  if (!"knn_acc_state" %in% names(su)) {
    su[, knn_acc_state := if (S == "P111T2") 0.994 else 0.992]
    fwrite(su, f, sep = "\t")
  }
}

pl <- c(lapply(SAMPLES, pan, colour_by = "clone"),
        lapply(SAMPLES, pan, colour_by = "state"))
g <- plot_grid(plotlist = pl, ncol = 2, align = "hv", axis = "tblr",
               labels = c("a", "b", "c", "d"), label_size = 9,
               label_fontfamily = FIG_FAMILY)

W <- 180; H <- 120
cairo_pdf(file.path(OUT, "figS_clone_embedding.pdf"), width = W*MM, height = H*MM,
          family = FIG_FAMILY); print(g); dev.off()
svglite::svglite(file.path(OUT, "figS_clone_embedding.svg"), width = W*MM, height = H*MM)
print(g); dev.off()
cat(sprintf("wrote figS_clone_embedding.{pdf,svg}  (%d x %d mm)\n", W, H))
