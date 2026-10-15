#!/usr/bin/env Rscript
# Panel E -- gene-level differential expression between the two subclones.
#
# Replaces the TF-motif volcano (2026-09-06), which is now Supplementary. The
# motif test was the only chromatin-level "regulatory" comparison available for a
# cohort with no peak calls, but it is the weaker instrument: of 12 prespecified
# motifs 6 reach FDR < 0.05 and on the oligodendroglial side only SOX10 does. The
# field's standard "compare the subclones" panel is gene-level DE (CopyKAT, Numbat),
# and here it is overwhelming: 4,289 of 10,756 genes at FDR < 0.05, 1,147 below
# 1e-10, with the canonical OPC genes (LHFPL3, PCDH15, OPCML, SEMA5A) at the top.
#
# x is the SAME quantity as panel F's y axis -- pseudobulk log2(C2/C1) from
# Result/qc/atac_rna_clone_*.tsv -- so a reader can carry a gene from E to F on one
# scale. y is -log10 BH-FDR from a per-gene Wilcoxon across nuclei
# (Code/eval/rna_de_clones.R). All expressed genes are tested; the 22 lineage
# genes are only HIGHLIGHTED, and the legend states the prediction each carries.
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(cowplot); library(ggrepel) })

de <- fread(file.path(DATA, "fig5e_rna_de_P111T2.tsv"))[, .(gene, fdr, rank)]
ar <- fread(sprintf("Result/qc/atac_rna_clone_%s.tsv", S))[, .(gene, d_rna)]
dd <- merge(ar, de, by = "gene")
LIN22 <- c(OLIG1="OPC",OLIG2="OPC",SOX10="OPC",PDGFRA="OPC",CSPG4="OPC",`NKX2-2`="OPC",SOX8="OPC",ASCL1="OPC",DLL3="OPC",BCAN="OPC",
           GFAP="AC",AQP4="AC",ALDH1L1="AC",S100B="AC",SLC1A2="AC",SLC1A3="AC",NFIA="AC",NFIB="AC",NFIX="AC",SOX9="AC",CLU="AC",ID3="AC")
dd[, fam := fifelse(gene %in% names(LIN22)[LIN22 == "OPC"], "expected in C2 (OPC lineage)",
             fifelse(gene %in% names(LIN22)[LIN22 == "AC"],  "expected in C1 (AC lineage)", "other"))]
dd[, fam := factor(fam, levels = c("expected in C2 (OPC lineage)", "expected in C1 (AC lineage)", "other"))]
dd[, lp := -log10(pmax(fdr, 1e-300))]
cat(sprintf("panel E (DE): %d genes with both a test and a pseudobulk ratio; FDR<0.05 %d; lineage genes significant %d of %d\n",
            nrow(dd), dd[fdr < .05, .N], dd[fam != "other" & fdr < .05, .N], dd[fam != "other", .N]))
# labels: 3 strongest per lineage family plus the 3 strongest genes overall. Fewer
# than before because the labels grew from 5.5 to 6.5 pt. EVERY coloured point is
# passed to the repel layer -- the unlabelled ones with an empty label -- because
# ggrepel only avoids points in its own data; without this, "DLL3" printed across
# two unlabelled OPC genes.
# Three per lineage plus the two strongest genes overall (LHFPL3, XKR4), placed by
# place_labels() (shared/fig_palette.R): every point in the panel, grey included, is
# an obstacle, and a name goes to the cheapest clear spot near its point with a
# leader when it cannot sit beside it. The legend's corner is a keep-out box.
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
  # Black names, not lineage-coloured: the Nature figure guide asks authors to avoid
  # coloured text, and the point beside each name already carries the lineage colour
  # (2026-10-09). The "other" names stay grey through the same fixed colour.
  geom_text(data = LB, aes(lx, ly, label = gene), colour = "grey15", size = FS_LGT/.pt,
            family = FIG_FAMILY, show.legend = FALSE) +
  # "other" IS in the scale and IS in the legend, drawn through the same mapping
  # as the coloured genes -- a literal-colour cloud is what produced empty legend
  # keys twice before.
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
        # top-LEFT: the significant-and-negative (C1) side is sparse above y ~ 90,
        # while the C2 side carries the extreme points
        legend.position = "inside", legend.position.inside = c(0.01, 0.995),
        legend.justification = c(0, 1), legend.direction = "vertical",
        legend.background = element_blank(), legend.key = element_blank(),
        legend.margin = margin(0, 0, 0, 0), legend.spacing.y = unit(0.1, "mm"),
        # legend text is 6.0 pt (the house legend size). At 6.5 the three-line key grew
        # wide enough to meet the top-right gene labels.
        legend.text = element_text(size = BASE_PT - 1, colour = "grey20"),
        legend.key.size = unit(3.0, "mm"), plot.margin = margin(2, 2, 1, 2))
if (!exists("FIG5_ASSEMBLING")) {
  cairo_pdf(file.path(OUT, "fig5_rna_de.pdf"), width = 70/25.4, height = 70/25.4, family = FIG_FAMILY)
  print(p_de); dev.off(); cat("wrote fig5_rna_de.pdf\n")
}
