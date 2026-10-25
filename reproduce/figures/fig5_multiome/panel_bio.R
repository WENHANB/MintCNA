#!/usr/bin/env Rscript
# Figure 5 panel f and Supplementary Fig. 8 panels a and b, on the two MintCNA subclones of
# P111T2: p_rep (Fig. 5f) promoter accessibility against expression between the clones, p_vol
# (Supplementary Fig. 8a) the TF-motif volcano, p_dose (Supplementary Fig. 8b) the per-gene
# copy-number dosage.
# Input : reproduce/figures/data/fig5_multiome/ (P111T2 MintCNA results, cell list and cell states;
#         tf_motif_P111T2.tsv, atac_rna_clone_P111T2.tsv, cna_dosage_link_P111T2_atac.tsv from
#         reproduce/analysis/multiome/evaluation/)
# Output: reproduce/figures/output/fig5/ (fig5c_*.csv, fig5f_*.csv source data) and output/supp/
#         (figS8a_*.csv, figS8b_*.csv); one PDF per panel when run on its own; assemble_fig5.R and
#         supp/figS8_subclone_biology.R source it with FIG5_ASSEMBLING defined
# Run from the root of the MintCNA repository: Rscript reproduce/figures/fig5_multiome/panel_bio.R
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(cowplot); library(ggrepel)
})
FIG <- "reproduce/figures"   # run from the root of the MintCNA repository
if (!file.exists(file.path(FIG, "shared/fig_palette.R"))) stop("run this script from the root of the MintCNA repository")
source(file.path(FIG, "shared/fig_palette.R"))
DATA <- file.path(FIG, "data/fig5_multiome")
OUT  <- file.path(FIG, "output/fig5"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
OUT_S8 <- file.path(FIG, "output/supp"); dir.create(OUT_S8, showWarnings = FALSE, recursive = TRUE)
S <- "P111T2"; TAG <- "h15a005e020"

# P111T2 is used because its subclones are recoverable from the discrete copy-number calls alone:
# k-means on the 1/2/3 state matrix reproduces the AG-DCAE latent clustering (ARI 0.831). The
# subclones come from ATAC coverage only, the cell states from the study's RNA annotation and the
# motif accessibility from ATAC at promoter resolution, so no step sees the next step's data.
e <- new.env(); load(sprintf("%s/%s_results_%s.RData", DATA, S, TAG), e)
Sm <- base::get("cna_state_matrix", e); sub <- base::get("subclones", e)
keep <- fread(sprintf("%s/%s_cells.tsv", DATA, S))[cell_type == "Malignant", barcode]
i <- colnames(Sm) %in% keep; bc <- colnames(Sm)[i]; sub <- sub[i]
tb <- sort(table(sub), decreasing = TRUE); o <- as.integer(names(tb)[1:2])
cl <- ifelse(sub == o[1], "C1", ifelse(sub == o[2], "C2", NA_character_))

## ---- clone x cell state -----------------------------------------------------
# The clone-by-state counts are source data; the correspondence itself is shown by the
# cell-state strip of the heatmap (Fig. 5c), so no bar chart is drawn here.
md <- fread(sprintf("%s/snrna_malignant_meta_%s.tsv", DATA, S))[care_id == S]
md[, bare := sub("^.*-([ACGT]+-[0-9]+)$", "\\1", CellID)]
st <- md$cell_state[match(bc, md$bare)]
d <- data.table(clone = cl, state = st)[!is.na(clone) & !is.na(state)]
fwrite(d[, .N, by = .(clone, state)], file.path(OUT, sprintf("fig5c_clone_by_state_%s.csv", S)))
d[, st3 := fifelse(state %in% c("OPC-like","AC-like"), state, "other")]
d[, st3 := factor(st3, levels = c("OPC-like","AC-like","other"))]
tab <- d[, .N, by = .(clone, st3)][, pct := 100*N/sum(N), by = clone]
cat("clone x state:\n"); print(dcast(tab, clone ~ st3, value.var = "pct", fill = 0))

## ---- Supplementary Fig. 8a: TF motif volcano --------------------------------
# The highlighted motifs are the prespecified lineage TF sets of
# reproduce/analysis/multiome/evaluation/tf_motif_subclone.R, including STAT3, which came out null.
LINEAGE <- c(OPC = "OLIG1|OLIG2|SOX10|SOX8|NKX2-2",
             AC  = "NFIA|NFIB|NFIC|NFIX|SOX9|STAT3|RFX4")
tf <- fread(sprintf("%s/tf_motif_%s.tsv", DATA, S))          # reproduce/analysis/multiome/evaluation/tf_motif_subclone.R
# Axes are named for the clones (C1, C2), not for the cell states. The sign is oriented so that
# positive means more accessible in the OPC-like clone, which in this sample makes x = C2 - C1.
opc_is <- tab[st3 == "OPC-like"][which.max(pct), clone]
tf[, x := if (opc_is == "C1") diff else -diff]
tf[, fam := fifelse(grepl(LINEAGE["OPC"], TF), "OPC lineage",
             fifelse(grepl(LINEAGE["AC"], TF), "AC lineage", "other"))]
tf[, lp := -log10(pmax(fdr, 1e-300))]
# Labelled: every prespecified motif at FDR < 0.05, whichever direction it moved, the three
# strongest motifs that were not named in advance, and ASCL1 (the canonical OPC bHLH factor;
# it was not prespecified, so it keeps the grey of the other unnamed motifs).
show <- rbind(tf[fam != "other" & fdr < 0.05][order(-abs(x))],
              head(tf[fam == "other"][order(-lp)], 3),
              tf[fam == "other" & TF == "ASCL1"])
show <- show[!duplicated(motif)]
# explicit nudges keep the NFIX(var.2) and NFIA labels (either side of zero, just above the FDR line) apart
show[, `:=`(nx = 0, ny = 0)]
show[TF == "NFIX(var.2)", `:=`(nx = -0.26, ny = 0)]
show[TF == "NFIA",        `:=`(nx =  0.15, ny = -0.9)]
p_vol <- ggplot(tf, aes(x, lp)) +
  geom_hline(yintercept = -log10(0.05), linewidth = 0.22, linetype = "22", colour = "grey70") +
  geom_vline(xintercept = 0, linewidth = 0.22, colour = "grey70") +
  # the background cloud is not mapped to colour; the labelled unnamed motifs are drawn through
  # the colour scale so that the "other motifs" legend entry has a glyph
  geom_point(data = tf[fam == "other"], colour = "grey82", size = 0.45, alpha = 0.55,
             show.legend = FALSE) +
  geom_point(data = tf[fam != "other"], aes(colour = fam), size = 1.15) +
  geom_point(data = show[fam == "other"], aes(colour = fam), size = 1.15) +
  # one repel layer for every label, so prespecified and unnamed labels avoid each other
  geom_text_repel(data = show, aes(label = TF, colour = fam), size = PT_SMALL/.pt,
                  position = position_nudge_repel(x = show$nx, y = show$ny),
                  segment.size = 0.16, segment.colour = "grey65", min.segment.length = 0,
                  box.padding = 0.38, point.padding = 0.12, max.overlaps = Inf,
                  force = 5, force_pull = 0.5, max.iter = 50000, seed = 11, show.legend = FALSE) +
  # colour encodes provenance (named in tf_motif_subclone.R before the test, or not), so the
  # grey class has its own legend entry
  scale_colour_manual(values = c(`OPC lineage` = STATE_COLORS[["OPC-like"]],
                                 `AC lineage`  = STATE_COLORS[["AC-like"]],
                                 `other` = "grey55"),
                      breaks = c("OPC lineage", "AC lineage", "other"),
                      labels = c("expected in C2 (OPC TFs)", "expected in C1 (AC TFs)", "other motifs"),
                      name = NULL, na.translate = FALSE) +
  # symmetric about zero because the quantity is a difference
  scale_x_continuous(limits = function(l) c(-1, 1) * max(abs(l)) * 1.18,
                     expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
  labs(x = "Motif accessibility  (C2 − C1)",
       y = expression(-log[10]~FDR)) +
  theme_cowplot(BASE_PT, font_family = FIG_FAMILY) +
  theme(aspect.ratio = 1,
        axis.text = element_text(size = BASE_PT - 1, colour = "grey20"),
        axis.title.y = element_text(size = PT_SUP),   # carries a subscript
        axis.title.x = element_text(size = BASE_PT),
        axis.line = element_line(linewidth = 0.3, colour = "grey30"),
        axis.ticks = element_line(linewidth = 0.3, colour = "grey30"),
        panel.grid = element_blank(),
        legend.position = "inside", legend.position.inside = c(0.01, 0.99),
        legend.justification = c(0, 1),
        legend.direction = "vertical",
        legend.background = element_blank(), legend.margin = margin(0,0,0,0),
        legend.text = element_text(size = BASE_PT - 1),
        legend.key.size = unit(3.2, "mm"), plot.margin = margin(2, 2, 1, 2))
fwrite(tf, file.path(OUT_S8, sprintf("figS8a_tf_motif_%s.csv", S)))

## ---- Figure 5f: does the accessibility difference reach the RNA? ------------
# Promoter accessibility (TSS +/- 2 kb) against gene expression between the same two clones, in
# the same nuclei. Genes are joined on Ensembl ID and barcodes are shared between the modalities.
ar <- fread(sprintf("%s/atac_rna_clone_%s.tsv", DATA, S))    # reproduce/analysis/multiome/evaluation/atac_to_rna_clone.R
LIN <- data.table(
  gene = c("OLIG1","OLIG2","SOX10","PDGFRA","CSPG4","NKX2-2","SOX8","ASCL1","DLL3","BCAN",
           "GFAP","AQP4","ALDH1L1","S100B","SLC1A2","SLC1A3","NFIA","NFIB","NFIX","SOX9","CLU","ID3"),
  fam  = c(rep("OPC lineage", 10), rep("AC lineage", 12)))
ar <- merge(ar, LIN, by = "gene", all.x = TRUE)
ar[is.na(fam), fam := "other"]
rho_ar  <- cor(ar$d_atac, ar$d_rna, method = "spearman")
conc    <- 100 * ar[abs(d_rna) > 1, mean(sign(d_atac) == sign(d_rna))]
cat(sprintf("panel f: %d genes, Spearman rho = %.3f, %.1f%% concordant among |log2FC RNA| > 1\n",
            nrow(ar), rho_ar, conc))
# Labels: the three strongest genes per lineage plus the two TFs named in the motif panel
# (ASCL1, SOX10), placed by place_labels() (shared/fig_palette.R) against every point in the
# panel, grey included; the corner holding the rho annotation is a keep-out box.
D_LINK <- c("ASCL1", "SOX10")
lab_ar <- ar[fam != "other"][, rk := frank(-abs(d_rna)), by = fam][rk <= 3 | gene %in% D_LINK][order(-abs(d_rna))]
F_XR <- range(ar$d_atac) + c(-1, 1) * 0.10 * diff(range(ar$d_atac))
F_YR <- range(ar$d_rna)  + c(-1, 1) * 0.05 * diff(range(ar$d_rna))
F_MM <- c(47, 47)
LB_ar <- place_labels(ar$d_atac, ar$d_rna, ar$fam != "other", match(lab_ar$gene, ar$gene), lab_ar$gene,
                      F_XR, F_YR, F_MM, pt_mm = fifelse(ar$fam != "other", 1.0, 0.4), text_pt = FS_LGT,
                      keepout = c(0, 13, F_MM[2] - 6, F_MM[2]))
LB_ar <- merge(as.data.table(LB_ar), ar[, .(gene, fam)], by = "gene")

p_rep <- ggplot(ar, aes(d_atac, d_rna)) +
  geom_hline(yintercept = 0, linewidth = 0.22, colour = "grey75") +
  geom_vline(xintercept = 0, linewidth = 0.22, colour = "grey75") +
  geom_point(data = ar[fam == "other"], colour = "grey82", size = 0.35, alpha = 0.45,
             show.legend = FALSE) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
              linewidth = 0.3, linetype = "22", colour = "grey45") +
  geom_segment(data = LB_ar[seg == TRUE], aes(x = ax, y = ay, xend = px, yend = py),
               linewidth = 0.18, colour = "grey60") +
  geom_point(data = ar[fam != "other"], aes(colour = fam), size = 1.15) +
  geom_text(data = LB_ar, aes(lx, ly, label = gene), colour = "grey15", size = FS_LGT/.pt,   # gene names in near-black, as in panel_de.R
            family = FIG_FAMILY, show.legend = FALSE) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.18, vjust = 1.5,
           size = PT_SUP/.pt, fontface = "bold", colour = "grey15", parse = TRUE,
           label = sprintf("rho == %.2f", rho_ar)) +
  # grey here means every other gene in the sample, not "not prespecified" as in the volcano
  scale_colour_manual(values = c(`OPC lineage` = STATE_COLORS[["OPC-like"]],
                                 `AC lineage`  = STATE_COLORS[["AC-like"]],
                                 `other` = "grey55"),
                      breaks = c("OPC lineage", "AC lineage", "other"),
                      labels = c("expected in C2 (OPC lineage)", "expected in C1 (AC lineage)", "all other genes"),
                      name = NULL, na.translate = FALSE) +
  scale_x_continuous(breaks = seq(-6, 6, 2), expand = expansion(mult = 0.10)) +
  scale_y_continuous(expand = expansion(mult = 0.05)) +
  labs(x = expression("Promoter accessibility  " * log[2] * " (C2/C1)"),
       y = expression("Expression  " * log[2] * " (C2/C1)")) +
  theme_cowplot(BASE_PT, font_family = FIG_FAMILY) +
  theme(axis.text = element_text(size = BASE_PT - 1, colour = "grey20"),
        axis.title = element_text(size = PT_SUP),
        axis.line = element_line(linewidth = 0.3, colour = "grey30"),
        axis.ticks = element_line(linewidth = 0.3, colour = "grey30"),
        panel.grid = element_blank(), legend.position = "none",
        aspect.ratio = 1,
        plot.margin = margin(2, 2, 1, 2))
fwrite(ar, file.path(OUT, sprintf("fig5f_atac_vs_rna_%s.csv", S)))

## ---- Supplementary Fig. 8b: what the copy number itself explains ------------
# The subclones are defined by copy number, so the lineage differences above could in principle be
# gene dosage. d_cn is read from the ATAC-only arm so that RNA does not enter both axes (in the
# ATAC + RNA arm the RNA on the y axis also enters the copy-number states on the x axis).
dz <- fread(sprintf("%s/cna_dosage_link_%s_atac.tsv", DATA, S))   # reproduce/analysis/multiome/evaluation/cna_dosage_link.R

# log2((2 + d_cn)/2) is the shift a gene would show if copy number acted at full efficiency with
# no compensation; its maximum over the lineage genes is the dosage ceiling.
dz[, dose_max := abs(log2((2 + d_cn) / 2))]
Lg <- merge(dz, LIN, by = "gene")[!is.na(d_cn)]
CEIL <- max(Lg$dose_max)
cat(sprintf(paste0("Supplementary Fig. 8b: dosage ceiling over %d lineage genes = %.3f log2; ",
                   "largest observed |d_rna| = %.2f (%.0fx); %d of %d copy-neutral\n"),
            nrow(Lg), CEIL, max(abs(Lg$d_rna)), max(abs(Lg$d_rna)) / CEIL,
            Lg[abs(d_cn) <= 0.2, .N], nrow(Lg)))

# Per gene, the grey dot is the shift copy number predicts, log2((2 + dCN)/2), i.e. dosage acting at
# full efficiency; the coloured dot is the observed expression shift; the gap is what dosage does not explain.
Lg[, pred := log2((2 + d_cn) / 2)]                    # signed, unlike dose_max
Lg[, gene := factor(gene, levels = Lg[order(d_rna), gene])]
cat(sprintf("Supplementary Fig. 8b: %d lineage genes; |prediction| <= %.2f, |observed| up to %.2f (%.0fx); %d move against their own copy number\n",
            nrow(Lg), max(abs(Lg$pred)), max(abs(Lg$d_rna)),
            max(abs(Lg$d_rna))/max(abs(Lg$pred)),
            Lg[sign(pred) != sign(d_rna) & abs(pred) > 0.01, .N]))

p_dose <- ggplot(Lg, aes(y = gene)) +
  geom_vline(xintercept = 0, linewidth = 0.22, colour = "grey75") +
  geom_segment(aes(x = pred, xend = d_rna, yend = gene),
               linewidth = 0.3, colour = "grey78") +
  geom_point(aes(x = pred), size = 0.85, colour = "grey50") +
  geom_point(aes(x = d_rna, colour = fam), size = 1.35, show.legend = FALSE) +
  # the two dot types are named directly above the top row instead of in a legend; "observed"
  # is right-aligned because the top row holds the rightmost point of the panel
  annotate("text", x = 0.1, y = nrow(Lg) + 1.6, hjust = 0.5, vjust = 0.5,
           size = PT_SMALL/.pt, colour = "grey40", label = "predicted") +
  geom_text(data = Lg[which.max(d_rna)], aes(x = d_rna, y = nrow(Lg) + 1.6),
            inherit.aes = FALSE, hjust = 1, vjust = 0.5,
            size = PT_SMALL/.pt, colour = "grey30", label = "observed") +
  scale_colour_manual(values = c(`OPC lineage` = STATE_COLORS[["OPC-like"]],
                                 `AC lineage`  = STATE_COLORS[["AC-like"]])) +
  scale_x_continuous(breaks = seq(-2, 4, 2)) +
  scale_y_discrete(expand = expansion(add = c(0.7, 2.6))) +
  labs(x = expression(log[2]*" shift  (C2/C1)"), y = NULL) +
  theme_cowplot(BASE_PT, font_family = FIG_FAMILY) +
  theme(axis.text.y = element_text(size = PT_SMALL, colour = "grey20"),   # 21 gene names in a dense tick row
        axis.text.x = element_text(size = BASE_PT - 1, colour = "grey20"),
        axis.title.x = element_text(size = PT_SUP),
        axis.line.y = element_line(linewidth = 0.3, colour = "grey30"),
        axis.ticks.y = element_blank(),
        axis.line.x = element_line(linewidth = 0.3, colour = "grey30"),
        axis.ticks.x = element_line(linewidth = 0.3, colour = "grey30"),
        panel.grid.major.y = element_line(linewidth = 0.15, colour = "grey93"),
        panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        plot.margin = margin(1, 2, 1, 1))

fwrite(Lg[, .(gene, fam, d_cn, dose_max, d_atac, d_rna)],
       file.path(OUT_S8, sprintf("figS8b_lineage_dosage_%s.csv", S)))

MM <- 1/25.4
if (!exists("FIG5_ASSEMBLING"))
for (x in list(list(p_vol,  "figS8a_motif_volcano.pdf", 70, 62),
               list(p_rep,  "fig5f_atac_rna.pdf", 73, 57),
               list(p_dose, "figS8b_dosage_pergene.pdf", 60, 75))) {
  cairo_pdf(file.path(OUT, x[[2]]), width = x[[3]]*MM, height = x[[4]]*MM, family = FIG_FAMILY)
  print(x[[1]]); dev.off(); cat("wrote", x[[2]], "\n")
}
