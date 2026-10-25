#!/usr/bin/env Rscript
# Figure 5d: the RNA UMAP of the malignant nuclei of P111T2 drawn twice, coloured by the MintCNA
# subclone and by the source study's RNA-derived cell state (p_umap_clone, p_umap_state).
# Input: reproduce/figures/data/fig5_multiome/clone_embedding_P111T2.tsv and clone_embedding_P111T2_summary.tsv
# (reproduce/analysis/multiome/evaluation/clone_in_embedding.R).
# Output: the two ggplot objects for the assembly; the 15-NN agreement summary is printed.
# Sourced by fig5_multiome/assemble_fig5.R after the panel scripts that load shared/fig_palette.R; run the assembly from
# the root of the MintCNA repository: Rscript reproduce/figures/fig5_multiome/assemble_fig5.R
#
# The subclones are clustered from the ATAC-trained AG-DCAE latent, with no expression feature entering
# the partition, and the cell states come from the source study's own RNA annotation, so the agreement
# of the two colourings is a check of the partition against an independent label. The 15-NN agreement
# values are quoted in the caption rather than on the panel; they are printed below.
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(cowplot) })

UMAP_SAMPLE <- "P111T2"
QC <- "reproduce/figures/data/fig5_multiome"   # clone_embedding_*.tsv from reproduce/analysis/multiome/evaluation/clone_in_embedding.R

.ud <- fread(file.path(QC, sprintf("clone_embedding_%s.tsv", UMAP_SAMPLE)))
.us <- fread(file.path(QC, sprintf("clone_embedding_%s_summary.tsv", UMAP_SAMPLE)))

# Neither UMAP carries a legend: the heatmap legend directly above in the same column defines both the
# cell-state colours and the clone colours (CLONE_COLORS2, the same two colours as the heatmap blocks),
# and the two subclones and the two dominant states are labelled directly beside their clouds.

# Label anchor for a group: y is the 98th percentile of the group's umap2 and x the median umap1 of the
# points in the top decile, so the label tracks the top of the cluster's bulk rather than a few stray
# nuclei sitting inside the other cluster.
.lab_pos <- function(d, grp_col, grp) {
  x <- d[get(grp_col) == grp]
  hi <- as.numeric(stats::quantile(x$umap2, 0.90, names = FALSE))
  data.table(umap1 = median(x[umap2 >= hi, umap1]),
             umap2 = as.numeric(stats::quantile(x$umap2, 0.98, names = FALSE)),
             lab   = grp)
}
.cl <- rbind(.lab_pos(.ud, "clone", "C1"), .lab_pos(.ud, "clone", "C2"))
# The state panel is labelled for the two dominant states only; the minor states are decoded by the
# heatmap legend.
.dom <- .ud[!is.na(state) & nzchar(state), .N, by = state][order(-N)][1:2, state]
.st <- rbindlist(lapply(.dom, function(k) .lab_pos(.ud[!is.na(state)], "state", k)))

.umap_base <- function(p, ttl)
  # alpha 0.9 and size 0.28 keep the lighter clone colour readable against white; both UMAPs share
  # this base so they stay comparable
  p + geom_point(size = 0.28, alpha = 0.9, show.legend = FALSE) +
    labs(title = ttl, x = "UMAP 1", y = NULL) +
    theme_cowplot(BASE_PT, font_family = FIG_FAMILY) +
    theme(plot.title = element_text(size = PT_SUP, face = "plain", colour = "grey15",
                                    hjust = 0.5, margin = margin(b = 1)),
          axis.text = element_blank(), axis.ticks = element_blank(),
          axis.line = element_line(linewidth = 0.3, colour = "grey30"),
          axis.title = element_text(size = BASE_PT, colour = "grey40"),
          panel.grid = element_blank(), legend.position = "none",
          plot.margin = margin(1, 1, 1, 1), aspect.ratio = 1) +
    # Horizontal expansion makes room for the direct labels near the right edge, where ggplot would
    # otherwise clip them silently; the same expansion on both UMAPs keeps the axes comparable.
    scale_x_continuous(expand = expansion(mult = 0.16))

p_umap_clone <- .umap_base(
  ggplot(.ud[!is.na(clone)], aes(umap1, umap2,
                                 colour = factor(clone, levels = c("C1", "C2")))),
  "MintCNA subclone") +
  scale_colour_manual(values = CLONE_COLORS2) +
  geom_text(data = .cl, aes(umap1, umap2, label = lab), inherit.aes = FALSE,
            size = PT_SUP/.pt, fontface = "bold", colour = "grey15", vjust = -0.3) +
  scale_y_continuous(expand = expansion(mult = c(0.04, 0.13))) +
  labs(y = "UMAP 2")

p_umap_state <- .umap_base(
  ggplot(.ud[!is.na(state) & nzchar(state)],
         aes(umap1, umap2, colour = factor(state, levels = names(CELLSTATE_COLORS)))),
  "Cell state (RNA)") +
  scale_colour_manual(values = CELLSTATE_COLORS, drop = TRUE) +
  geom_text(data = .st, aes(umap1, umap2, label = lab), inherit.aes = FALSE,
            size = PT_SUP/.pt, fontface = "bold", colour = "grey15", vjust = -0.3) +
  scale_y_continuous(expand = expansion(mult = c(0.04, 0.13)))

cat(sprintf("panel D: %s, %d nuclei, 15-NN clone %.3f vs state %.3f (floor %.3f)\n",
            UMAP_SAMPLE, nrow(.ud), .us$knn_acc_clone, .us$knn_acc_state, .us$majority_floor))
