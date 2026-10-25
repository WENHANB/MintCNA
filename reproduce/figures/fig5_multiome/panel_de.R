#!/usr/bin/env Rscript
# Figure 5e: gene-level differential expression between the two MintCNA subclones of P111T2, drawn as
# pseudobulk log2(C2/C1) expression against -log10 BH-FDR of a per-gene Wilcoxon test across nuclei.
# Input (reproduce/figures/data/fig5_multiome/): fig5e_rna_de_P111T2.tsv (reproduce/analysis/multiome/evaluation/rna_de_clones.R),
# atac_rna_clone_P111T2.tsv (reproduce/analysis/multiome/evaluation/atac_to_rna_clone.R).
# Output: p_de for the assembly, and reproduce/figures/output/fig5/fig5_rna_de.pdf when FIG5_ASSEMBLING is not set.
# Sourced by fig5_multiome/assemble_fig5.R after panel_bio.R and panel_b_accuracy.R, which define DATA, S and
# OUT and load shared/fig_palette.R; run the assembly from the root of the MintCNA repository:
# Rscript reproduce/figures/fig5_multiome/assemble_fig5.R
#
# The x axis is the same quantity as panel f's y axis, the pseudobulk log2(C2/C1) from atac_rna_clone_P111T2.tsv,
# so a gene can be carried from e to f on one scale. All expressed genes are tested; the 22 lineage genes
# are only highlighted, and the legend states the direction expected for each family.
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(cowplot); library(ggrepel) })

de <- fread(file.path(DATA, "fig5e_rna_de_P111T2.tsv"))[, .(gene, fdr, rank)]
ar <- fread(sprintf("%s/atac_rna_clone_%s.tsv", DATA, S))[, .(gene, d_rna)]
dd <- merge(ar, de, by = "gene")
LIN22 <- c(OLIG1="OPC",OLIG2="OPC",SOX10="OPC",PDGFRA="OPC",CSPG4="OPC",`NKX2-2`="OPC",SOX8="OPC",ASCL1="OPC",DLL3="OPC",BCAN="OPC",
           GFAP="AC",AQP4="AC",ALDH1L1="AC",S100B="AC",SLC1A2="AC",SLC1A3="AC",NFIA="AC",NFIB="AC",NFIX="AC",SOX9="AC",CLU="AC",ID3="AC")
dd[, fam := fifelse(gene %in% names(LIN22)[LIN22 == "OPC"], "expected in C2 (OPC lineage)",
             fifelse(gene %in% names(LIN22)[LIN22 == "AC"],  "expected in C1 (AC lineage)", "other"))]
dd[, fam := factor(fam, levels = c("expected in C2 (OPC lineage)", "expected in C1 (AC lineage)", "other"))]
dd[, lp := -log10(pmax(fdr, 1e-300))]
cat(sprintf("panel E (DE): %d genes with both a test and a pseudobulk ratio; FDR<0.05 %d; lineage genes significant %d of %d\n",
            nrow(dd), dd[fdr < .05, .N], dd[fam != "other" & fdr < .05, .N], dd[fam != "other", .N]))
# Labels: the three strongest genes per lineage family plus the two strongest overall, placed by
# place_labels() (shared/fig_palette.R). Every point in the panel is an obstacle, a name goes to the
# nearest clear spot with a leader line when it cannot sit beside its point, and the legend corner
# is a keep-out box.
lab_de <- rbind(dd[fam != "other"][, rk := frank(fdr), by = fam][rk <= 3],
                dd[fam == "other"][order(fdr)][1:2][, rk := 0L])[order(fdr)]
E_XR <- range(dd$d_rna) + c(-1, 1) * 0.08 * diff(range(dd$d_rna))
E_YR <- c(0, max(dd$lp)) + c(-0.02, 0.14) * max(dd$lp)
E_MM <- c(47, 47)                                  # plot area, mm (60 mm column, square)
LB <- place_labels(dd$d_rna, dd$lp, dd$fam != "other", match(lab_de$gene, dd$gene), lab_de$gene,
                   E_XR, E_YR, E_MM, pt_mm = fifelse(dd$fam != "other", 1.0, 0.45), text_pt = FS_LGT,
                   keepout = c(0, 30, E_MM[2] - 9, E_MM[2]))   # legend, top-left
LB <- merge(as.data.table(LB), dd[, .(gene, fam)], by = "gene")
FAMC <- c(`expected in C2 (OPC lineage)` = STATE_COLORS[["OPC-like"]],
          `expected in C1 (AC lineage)`  = STATE_COLORS[["AC-like"]],
          other = "grey55")
p_de <- ggplot(dd, aes(d_rna, lp)) +
  geom_hline(yintercept = -log10(0.05), linewidth = 0.22, linetype = "22", colour = "grey70") +
  geom_vline(xintercept = 0, linewidth = 0.22, colour = "grey70") +
  geom_point(data = dd[fam == "other"], aes(colour = fam), size = 0.45, alpha = 0.45) +
  geom_segment(data = LB[seg == TRUE], aes(x = ax, y = ay, xend = px, yend = py),
               linewidth = 0.18, colour = "grey60") +
  geom_point(data = dd[fam != "other"], aes(colour = fam), size = 1.15) +
  # gene names in near-black; the point beside each name carries the lineage colour
  geom_text(data = LB, aes(lx, ly, label = gene), colour = "grey15", size = FS_LGT/.pt,
            family = FIG_FAMILY, show.legend = FALSE) +
  # "other" is in the scale and the legend through the same mapping as the lineage genes
  scale_colour_manual(values = FAMC, breaks = names(FAMC),
                      labels = c("expected in C2 (OPC lineage)", "expected in C1 (AC lineage)", "all other genes"),
                      name = NULL,
                      guide = guide_legend(override.aes = list(size = c(1.6, 1.6, 1.0), alpha = 1))) +
  scale_x_continuous(expand = expansion(mult = 0.08)) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.14))) +
  labs(x = expression("Expression  "*log[2]*" (C2/C1)"), y = expression(-log[10]*" FDR")) +
  theme_cowplot(BASE_PT, font_family = FIG_FAMILY) +
  theme(aspect.ratio = 1,
        axis.text = element_text(size = BASE_PT - 1, colour = "grey20"),
        axis.title = element_text(size = PT_SUP),            # both titles carry plotmath
        axis.line = element_line(linewidth = 0.3, colour = "grey30"),
        axis.ticks = element_line(linewidth = 0.3, colour = "grey30"),
        panel.grid = element_blank(),
        # legend in the top-left corner, where the plane is sparse
        legend.position = "inside", legend.position.inside = c(0.01, 0.995),
        legend.justification = c(0, 1), legend.direction = "vertical",
        legend.background = element_blank(), legend.key = element_blank(),
        legend.margin = margin(0, 0, 0, 0), legend.spacing.y = unit(0.1, "mm"),
        legend.text = element_text(size = BASE_PT - 1, colour = "grey20"),
        legend.key.size = unit(3.0, "mm"), plot.margin = margin(2, 2, 1, 2))
if (!exists("FIG5_ASSEMBLING")) {
  cairo_pdf(file.path(OUT, "fig5_rna_de.pdf"), width = 70/25.4, height = 70/25.4, family = FIG_FAMILY)
  print(p_de); dev.off(); cat("wrote fig5_rna_de.pdf\n")
}
