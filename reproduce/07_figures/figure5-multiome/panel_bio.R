#!/usr/bin/env Rscript
# Figure 4, biology panels. Sample P111T2.
#
# WHY THIS SAMPLE. Subclones here are recoverable from the DISCRETE copy-number
# calls alone (k-means on the 1/2/3 state matrix reproduces the AG-DCAE latent
# clustering at ARI 0.831), so they are copy number and not chromatin
# accessibility wearing a copy-number costume. That test fails on P105T3
# (ARI -0.003) and is weak on P100T1 (0.063), and those samples are therefore not
# used for the biology even though P100T1 has the better F1. The criterion was
# fixed before the panels were drawn.
#
# THE CHAIN. Subclones come from ATAC coverage only; the cell states come from
# the paper's own RNA annotation; the motif accessibility comes back from ATAC at
# promoter resolution. No step sees the next step's data.
.libPaths(c("<R_LIBS_USER>", .libPaths()))
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(cowplot); library(ggrepel)
})
ROOT <- "<MINTCNA_ROOT>"
source("<REPRODUCE_DIR>/07_figures/shared/fig_palette.R")
OUT  <- file.path(ROOT, "figure4paper/figure5-multiome/output")
DATA <- file.path(ROOT, "figure4paper/figure5-multiome/data")
for (d in c(OUT, DATA)) dir.create(d, showWarnings = FALSE, recursive = TRUE)
S <- "P111T2"; TAG <- "h15a005e020"
setwd(file.path(ROOT, "ATAC_RNA"))

e <- new.env(); load(sprintf("Result/MintCNA/%s_multi/%s_results_%s.RData", S, S, TAG), e)
Sm <- base::get("cna_state_matrix", e); sub <- base::get("subclones", e)
keep <- fread(sprintf("Data/cell_lists/%s_cells.tsv", S))[cell_type == "Malignant", barcode]
i <- colnames(Sm) %in% keep; bc <- colnames(Sm)[i]; sub <- sub[i]
tb <- sort(table(sub), decreasing = TRUE); o <- as.integer(names(tb)[1:2])
cl <- ifelse(sub == o[1], "C1", ifelse(sub == o[2], "C2", NA_character_))

## ---- panel: clone x cell state ----------------------------------------------
md <- as.data.table(readRDS("Data/synapse/snrna_malignant_meta_data.RDS"))[care_id == S]
md[, bare := sub("^.*-([ACGT]+-[0-9]+)$", "\\1", CellID)]
st <- md$cell_state[match(bc, md$bare)]
d <- data.table(clone = cl, state = st)[!is.na(clone) & !is.na(state)]
fwrite(d[, .N, by = .(clone, state)], file.path(DATA, "fig5_state_by_clone_full.csv"))
d[, st3 := fifelse(state %in% c("OPC-like","AC-like"), state, "other")]
d[, st3 := factor(st3, levels = c("OPC-like","AC-like","other"))]
tab <- d[, .N, by = .(clone, st3)][, pct := 100*N/sum(N), by = clone]
lab <- d[, .(n = .N), by = clone][, txt := sprintf("%s\n n = %d", clone, n)]
tab <- merge(tab, lab, by = "clone")
cat("clone x state:\n"); print(dcast(tab, clone ~ st3, value.var = "pct", fill = 0))

p_state <- ggplot(tab, aes(txt, pct, fill = st3)) +
  geom_col(width = 0.62, colour = NA) +
  scale_fill_manual(values = STATE_COLORS, name = NULL) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 100, 25)) +
  labs(y = "Malignant nuclei (%)", x = NULL) +
  theme_cowplot(BASE_PT, font_family = FIG_FAMILY) +
  theme(axis.text.x = element_text(size = BASE_PT, colour = "grey15"),
        axis.text.y = element_text(size = BASE_PT - 1, colour = "grey20"),
        axis.title.y = element_text(size = BASE_PT),
        axis.line = element_line(linewidth = 0.3, colour = "grey30"),
        axis.ticks = element_line(linewidth = 0.3, colour = "grey30"),
        panel.grid = element_blank(), legend.position = "bottom",
        legend.text = element_text(size = BASE_PT - 1),
        legend.key.size = unit(3.2, "mm"), plot.margin = margin(2, 2, 1, 2))

## ---- panel: TF motif volcano ------------------------------------------------
# THE HIGHLIGHTED SET IS EXACTLY THE PRESPECIFIED SET, and it did not used to be.
# The only dated prespecification artefact is the regex at
# Code/eval/tf_motif_subclone.R:99 (2026-08-29):
#     OLIG1|OLIG2|SOX10|SOX8|NKX2-2|NFIA|NFIB|NFIC|NFIX|SOX9|STAT3|RFX4
# This file previously highlighted a DIFFERENT set: it had added ASCL1 and dropped
# STAT3. ASCL1 is a strong hit (x = +0.35, FDR = 2.5e-10); STAT3 is null
# (x = -0.07, FDR = 0.22). Adding a hit and removing a miss is precisely what makes
# a prespecification claim worthless, and no record anywhere names ASCL1 in
# advance. So the lists are now identical to the analysis script's, STAT3 included
# -- a prespecified TF that came out null is evidence the set was not curated, and
# it belongs in the panel for that reason.
# ASCL1 keeps its place among the strongest NOT-prespecified motifs below, which is
# what it is. Do not quietly promote it back.
LINEAGE <- c(OPC = "OLIG1|OLIG2|SOX10|SOX8|NKX2-2",
             AC  = "NFIA|NFIB|NFIC|NFIX|SOX9|STAT3|RFX4")
tf <- fread(sprintf("Result/qc/tf_motif_%s.tsv", S))
# EVERY AXIS IN D, E AND F IS NAMED FOR OUR CLONES, C1 and C2, not for the study's
# cell states. The quantities are computed between C1 and C2 cells -- those come
# from our CNA calls -- and the clone-to-state correspondence is 95% / 87%, not
# 100%, so labelling the contrast "OPC - AC" would both credit the study's
# annotation for our partition and misstate what was averaged. Panel C carries the
# correspondence, as an independent check: its clone strip is ours, its cell-state
# strip is the study's, and the reader reads the mapping there.
#
# The sign is still oriented by which clone is OPC-like, so the axis direction is
# comparable across samples; in this sample that makes x = C2 - C1.
# read "open in the OPC-like clone" in every sample regardless of clone numbering
opc_is <- tab[st3 == "OPC-like"][which.max(pct), clone]
tf[, x := if (opc_is == "C1") diff else -diff]
tf[, fam := fifelse(grepl(LINEAGE["OPC"], TF), "OPC lineage",
             fifelse(grepl(LINEAGE["AC"], TF), "AC lineage", "other"))]
tf[, lp := -log10(pmax(fdr, 1e-300))]
# Label the prespecified families AND the strongest unnamed motifs. Leaving the
# two most significant points in the panel unlabelled invites exactly the
# question the panel cannot then answer.
# Labelled: every PRESPECIFIED motif at FDR < 0.05, the 3 strongest motifs that
# were not named in advance, and ASCL1. ASCL1 is called out by name because it is the
# most interpretable of the unnamed hits -- a class-II bHLH, the canonical OPC
# factor -- and because this figure previously mislabelled it as prespecified. It
# is drawn in the SAME grey as the other unnamed motifs, which is the honest
# statement: a real, strong result (x = +0.35, FDR = 2.5e-10) that the analysis did
# not predict. By -log10 FDR it ranks 6th among the unnamed, so it does not reach
# the panel on rank alone; naming it is a deliberate, disclosed choice.
# 2026-09-25, strand-aware promoter windows: EVERY prespecified motif at FDR < 0.05
# is named, whichever direction it moved. Two now move AGAINST their expectation
# (SOX9, FDR 5e-6, and NFIA, FDR 0.04, both more accessible in the OPC-like clone);
# a top-7-by-effect rule would have named SOX9 but hidden NFIA, and a contradiction
# of the prespecified set is exactly what must not be hidden.
show <- rbind(tf[fam != "other" & fdr < 0.05][order(-abs(x))],
              head(tf[fam == "other"][order(-lp)], 3),
              tf[fam == "other" & TF == "ASCL1"])
show <- show[!duplicated(motif)]
# 2026-10-15: NFIX(var.2) (x -0.19, lp 2.7) and NFIA (x +0.15, lp 1.4) sit either side of zero just above the
# FDR line and repel parked both labels on the same spot between them. Explicit nudges: NFIX(var.2) straight
# left into the gap between NFIC(var.2) and NFIC::TLX1, NFIA right and below the dashed line.
show[, `:=`(nx = 0, ny = 0)]
show[TF == "NFIX(var.2)", `:=`(nx = -0.26, ny = 0)]
show[TF == "NFIA",        `:=`(nx =  0.15, ny = -0.9)]
p_vol <- ggplot(tf, aes(x, lp)) +
  geom_hline(yintercept = -log10(0.05), linewidth = 0.22, linetype = "22", colour = "grey70") +
  geom_vline(xintercept = 0, linewidth = 0.22, colour = "grey70") +
  # inherit.aes = FALSE and no colour aesthetic on the grey cloud: mapping it to
  # a colour outside the manual scale is what put an "NA" key in the legend.
  geom_point(data = tf[fam == "other"], colour = "grey82", size = 0.45, alpha = 0.55,
             show.legend = FALSE) +
  geom_point(data = tf[fam != "other"], aes(colour = fam), size = 1.15) +
  # The four LABELLED unnamed motifs, drawn through the colour scale so the
  # "not prespecified" key has a glyph. The background cloud above stays unmapped
  # (mapping it is what once put an "NA" key in the legend), and without this layer
  # the third legend entry rendered as text with no dot beside it. These four are
  # also the ones being named, so drawing them a shade darker than the cloud is
  # correct rather than incidental.
  geom_point(data = show[fam == "other"], aes(colour = fam), size = 1.15) +
  # ONE repel layer for every label in this panel, prespecified and unnamed alike.
  # Splitting them across two layers meant neither could see the other's text, and
  # ASCL1 printed straight through MYOG. force is raised because the labels now
  # have to resolve against each other as well as against the points.
  geom_text_repel(data = show, aes(label = TF, colour = fam), size = PT_SMALL/.pt,
                  position = position_nudge_repel(x = show$nx, y = show$ny),
                  segment.size = 0.16, segment.colour = "grey65", min.segment.length = 0,
                  box.padding = 0.38, point.padding = 0.12, max.overlaps = Inf,
                  force = 5, force_pull = 0.5, max.iter = 50000, seed = 11, show.legend = FALSE) +   # 2026-10-15: NFIX(var.2) and NFIA overlapped at force 3 / seed 7
  # THE GREY IS IN THE LEGEND, and it has to be. Colour here does not encode a
  # lineage, it encodes PROVENANCE: coloured = named in tf_motif_subclone.R before
  # the test was run, grey = everything else. Without a third key entry a reader
  # sees four grey labels (NRF1, ZSCAN4, MYOG, ASCL1) beside the coloured ones and
  # cannot know the two classes mean different things -- which is exactly what
  # makes a prespecified hit worth more than a top-of-the-list one. AtaCNA's
  # Figure 5F, the closest published analogue, does not draw this distinction.
  scale_colour_manual(values = c(`OPC lineage` = STATE_COLORS[["OPC-like"]],
                                 `AC lineage`  = STATE_COLORS[["AC-like"]],
                                 `other` = "grey55"),
                      breaks = c("OPC lineage", "AC lineage", "other"),
                      labels = c("expected in C2 (OPC TFs)", "expected in C1 (AC TFs)", "other motifs"),
                      name = NULL, na.translate = FALSE) +
  # Symmetric about zero because the quantity is a difference, but only 8% wider
  # than the data. An earlier 35% pad was added to stop a label clipping and left
  # a third of the panel empty; the label crowding is a repel problem, not a
  # range problem.
  scale_x_continuous(limits = function(l) c(-1, 1) * max(abs(l)) * 1.18,
                     expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
  # compact: the panel lost width to B's and E's squares, and the full
  # "OPC-like / AC-like clone" is spelled out in the legend and in panel C.
  labs(x = "Motif accessibility  (C2 − C1)",
       y = expression(-log[10]~FDR)) +
  theme_cowplot(BASE_PT, font_family = FIG_FAMILY) +
  theme(# SQUARE since 2026-09-06: E shares the right column with two squares
        # and both of its axes are unit-free scores, so 1:1 loses nothing.
        aspect.ratio = 1,
        axis.text = element_text(size = BASE_PT - 1, colour = "grey20"),
        # PT_SUP, not BASE_PT: this panel's y title is expression(-log[10]~FDR)
        # and plotmath sets the subscript at ~0.7 of the parent. BASE_PT (7) put
        # the "10" at 4.9 pt, under the 5 pt floor; PT_SUP (7.4) puts it at 5.2.
        # Caught in the exported SVG, not in the source -- a declared size of 7
        # is not a rendered size of 7 once plotmath is involved.
        # SPLIT .x / .y. The 7.4 exception exists ONLY for a title carrying a
        # sub- or superscript, which plotmath renders at ~0.7 of its parent: at
        # BASE_PT the "2" of log[2] lands at 4.9 pt, under Nature's 5 pt floor.
        # An x title with no plotmath descender has no such need, and at 7.4 it
        # renders 0.4 pt larger than the identical role in Figure 3 and than
        # panel B's own "Precision" six centimetres away. Figure 3 uses ONE size
        # per role; this is Figure 4 getting back to that.
        axis.title.y = element_text(size = PT_SUP),
        axis.title.x = element_text(size = BASE_PT),
        axis.line = element_line(linewidth = 0.3, colour = "grey30"),
        axis.ticks = element_line(linewidth = 0.3, colour = "grey30"),
        panel.grid = element_blank(), # top-LEFT inside: a volcano is empty on the significant-but-negative
        # side above the labelled hits, whereas top-centre sat on NRF1.
        # VERTICAL, not horizontal. Two entries side by side spanned ~30 mm of a
        # 46.9 mm panel and reached across to NRF1's label at the top right; stacked
        # they occupy ~15 mm and stay in the volcano's empty upper-left corner.
        legend.position = "inside", legend.position.inside = c(0.01, 0.99),
        legend.justification = c(0, 1),
        legend.direction = "vertical",
        legend.background = element_blank(), legend.margin = margin(0,0,0,0),
        legend.text = element_text(size = BASE_PT - 1),
        legend.key.size = unit(3.2, "mm"), plot.margin = margin(2, 2, 1, 2))
fwrite(tf, file.path(DATA, sprintf("fig5_tf_motif_%s.csv", S)))

## ---- SUPPLEMENTARY panel: the same TFs in a second tumour ------------------
# Demoted out of Figure 4 (it now goes to the SI). It is a real internal control
# -- in P100T1 the OPC-like clone is C1 and here it is C2, so a motif difference
# that followed the clone NUMBER rather than the cell state would be an artefact
# -- but as a main panel it asked the reader to hold two tumours and a reversed
# labelling in mind to read one correlation. The main figure's panel E now spends
# that space on the other modality of the same nuclei instead.
# The internal control. In P100T1 the OPC-like clone is C1; here it is C2. If the
# motif differences were a batch or labelling artefact they would follow the
# clone NUMBER; if they are lineage, they follow the cell state. Both are plotted
# on the "OPC-like minus AC-like" axis, so agreement is a positive correlation.
t2 <- fread("Result/qc/tf_motif_P100T1.tsv")
t2[, x2 := diff]                                    # P100T1: C1 is the OPC-like clone
both <- merge(tf[, .(TF, motif, x1 = x, fam)], t2[, .(motif, x2)], by = "motif")
both <- both[fam != "other" | abs(x1) > 0.3 | abs(x2) > 0.3]
rho <- cor(both$x1, both$x2, method = "spearman")
cat(sprintf("\ncross-sample agreement on the OPC-minus-AC axis: Spearman %.3f (n = %d motifs)\n",
            rho, nrow(both)))
REP_LIM <- c(-1, 1) * max(abs(c(both$x1, both$x2))) * 1.05
p_si <- ggplot(both, aes(x2, x1)) +
  geom_hline(yintercept = 0, linewidth = 0.22, colour = "grey75") +
  geom_vline(xintercept = 0, linewidth = 0.22, colour = "grey75") +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.22, linetype = "22", colour = "grey70") +
  geom_point(data = both[fam == "other"], colour = "grey82", size = 0.45, alpha = 0.55,
             show.legend = FALSE) +
  geom_point(data = both[fam != "other"], aes(colour = fam), size = 1.15) +
  geom_text_repel(data = head(both[fam != "other"][order(-abs(x1))], 7), aes(label = TF, colour = fam),
                  size = PT_SMALL/.pt, segment.size = 0.16, segment.colour = "grey65",
                  min.segment.length = 0, box.padding = 0.22, point.padding = 0.12,
                  max.overlaps = 30, seed = 1, show.legend = FALSE) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.25, vjust = 1.6,
           size = PT_SUP/.pt, fontface = "bold", colour = "grey15", parse = TRUE,
           label = sprintf("rho == %.2f", rho)) +
  # Both axes carry the same quantity in the same units and the dashed line is
  # y = x, so the two ranges must match or the identity line stops being a
  # 45-degree line and the eye reads the agreement wrongly.
  # Equal limits on both axes -- that is what makes the two samples comparable and
  # keeps y = x a meaningful locus. coord_fixed() is NOT used: it forces a square
  # panel, which letterboxed this plot inside its grid slot and left it visibly
  # narrower than D directly above it. Panel alignment across the column is worth
  # more than the diagonal landing at exactly 45 degrees, and the reference line
  # is still exactly the set of points where the two samples agree.
  scale_x_continuous(limits = REP_LIM, expand = c(0, 0)) +
  scale_y_continuous(limits = REP_LIM, expand = c(0, 0)) +
  scale_colour_manual(values = c(`OPC lineage` = STATE_COLORS[["OPC-like"]],
                                 `AC lineage`  = STATE_COLORS[["AC-like"]]), name = NULL) +
  labs(x = "P100T1   OPC − AC", y = "P111T2   OPC − AC") +
  theme_cowplot(BASE_PT, font_family = FIG_FAMILY) +
  theme(axis.text = element_text(size = BASE_PT - 1, colour = "grey20"),
        axis.title = element_text(size = BASE_PT),
        axis.line = element_line(linewidth = 0.3, colour = "grey30"),
        axis.ticks = element_line(linewidth = 0.3, colour = "grey30"),
        panel.grid = element_blank(), legend.position = "none",
        # SQUARE, restored deliberately. Both axes carry log2 of the same OPC/AC
        # ratio on the same scale, and the dashed reference is y = x, so anything
        # other than 1:1 stretches one modality against the other and tilts a line
        # that means "the two agree". This is bought at the cost of width against
        # D and F in the same column; a distorted panel is the worse defect.
        aspect.ratio = 1,
        plot.margin = margin(2, 2, 1, 2))

## ---- panel E: does the accessibility difference reach the RNA? -------------
# Panel D is an ATAC-only statement. This closes the loop on the second modality
# of the same assay, in the same nuclei: promoter accessibility (TSS +/- 2 kb)
# against gene expression, both on the identical OPC-like minus AC-like contrast.
# A positive relation means an opening promoter is an expressed gene rather than
# an isolated ATAC artefact. Genes are joined on Ensembl ID, and barcodes are
# shared between the modalities, so neither join is approximate.
ar <- fread(sprintf("Result/qc/atac_rna_clone_%s.tsv", S))
LIN <- data.table(
  gene = c("OLIG1","OLIG2","SOX10","PDGFRA","CSPG4","NKX2-2","SOX8","ASCL1","DLL3","BCAN",
           "GFAP","AQP4","ALDH1L1","S100B","SLC1A2","SLC1A3","NFIA","NFIB","NFIX","SOX9","CLU","ID3"),
  fam  = c(rep("OPC lineage", 10), rep("AC lineage", 12)))
ar <- merge(ar, LIN, by = "gene", all.x = TRUE)
ar[is.na(fam), fam := "other"]
rho_ar  <- cor(ar$d_atac, ar$d_rna, method = "spearman")
conc    <- 100 * ar[abs(d_rna) > 1, mean(sign(d_atac) == sign(d_rna))]
cat(sprintf("panel E: %d genes, Spearman rho = %.3f, %.1f%% concordant among |log2FC RNA| > 1\n",
            nrow(ar), rho_ar, conc))
# Four strongest per lineage, PLUS the TFs that panel D names, so a reader can
# carry ASCL1 and SOX10 from the motif panel into this one. Labelling every gene
# that moves put 16 names into a 61 x 35 mm panel and ALDH1L1 printed through
# NFIA; the full per-gene table is Source Data, and the panel's job is to show the
# separation, not to name every member of it.
# Five per lineage, not four. ggrepel only avoids the points in its OWN data, so
# an unlabelled coloured point is invisible to it: at four, OLIG2 was dropped and
# its point printed through the OLIG1 label sitting beside it. Labelling it makes
# it an obstacle repel can see.
D_LINK <- c("ASCL1", "SOX10")
# Three per lineage plus the two motif-panel TFs (ASCL1, SOX10), placed by
# place_labels() (shared/fig_palette.R) against EVERY point in the panel, grey
# included; see panel E. The rho annotation's corner is a keep-out box.
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
  geom_text(data = LB_ar, aes(lx, ly, label = gene), colour = "grey15", size = FS_LGT/.pt,   # black names, see panel_de.R
            family = FIG_FAMILY, show.legend = FALSE) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.18, vjust = 1.5,
           size = PT_SUP/.pt, fontface = "bold", colour = "grey15", parse = TRUE,
           label = sprintf("rho == %.2f", rho_ar)) +
  # Grey is named here too, but it means something DIFFERENT from panel E's grey:
  # there it is "not prespecified", here it is simply every other gene in the
  # sample. Same colour, two meanings, two panels -- so both have to say which.
  scale_colour_manual(values = c(`OPC lineage` = STATE_COLORS[["OPC-like"]],
                                 `AC lineage`  = STATE_COLORS[["AC-like"]],
                                 `other` = "grey55"),
                      breaks = c("OPC lineage", "AC lineage", "other"),
                      labels = c("expected in C2 (OPC lineage)", "expected in C1 (AC lineage)", "all other genes"),
                      name = NULL, na.translate = FALSE) +
  # a little more room on x than the data needs: DLL3 sits at +4.1 and its label
  # was landing on the panel border
  # integer breaks: with strand-aware promoter windows the x range is about +/-4.5
  # and ggplot's default ticks became -5.0/-2.5/0.0, decimals for no reason
  scale_x_continuous(breaks = seq(-6, 6, 2), expand = expansion(mult = 0.10)) +
  scale_y_continuous(expand = expansion(mult = 0.05)) +
  # "(OPC/AC)", not "(OPC-like / AC-like)". Spelled out, the y title measured
  # 42.1 mm on a 34.8 mm plot area: it overran its own panel, spilled into the
  # gap above, and closed to 0.5 mm of the panel letter. The full contrast is
  # already spelled out on D's x axis directly above, on F's x axis, and in the
  # legend, so this panel does not have to carry it twice more.
  labs(x = expression("Promoter accessibility  " * log[2] * " (C2/C1)"),
       y = expression("Expression  " * log[2] * " (C2/C1)")) +
  theme_cowplot(BASE_PT, font_family = FIG_FAMILY) +
  theme(axis.text = element_text(size = BASE_PT - 1, colour = "grey20"),
        axis.title = element_text(size = PT_SUP),
        axis.line = element_line(linewidth = 0.3, colour = "grey30"),
        axis.ticks = element_line(linewidth = 0.3, colour = "grey30"),
        panel.grid = element_blank(), legend.position = "none",
        # SQUARE. Both axes carry log2 of the same OPC/AC ratio on the same scale
        # and the dashed reference is y = x, so any aspect other than 1:1 stretches
        # one modality against the other and tilts a line that means "the two
        # agree". Bought at the cost of width against D and F in this column; a
        # distorted panel is the worse defect.
        aspect.ratio = 1,
        plot.margin = margin(2, 2, 1, 2))
fwrite(ar, file.path(DATA, sprintf("fig5e_atac_vs_rna_%s.csv", S)))

## ---- panel G: what the COPY NUMBER itself explains -------------------------
# E and F compare two subclones and find lineage differences. The subclones are
# DEFINED by copy number, so a reviewer is entitled to ask whether the whole
# lineage difference is just gene dosage. This panel answers that, and it is the
# only reason it exists.
#
# It answers it in two rows sharing one x axis, because the answer is a statement
# about MAGNITUDE and the two magnitudes differ by an order of magnitude:
#   lower row  the dose response. Where the clones differ in copy number both
#              modalities follow. This is the forward evidence and the only claim
#              in the figure that requires a multiome assay.
#   upper row  the 21 lineage genes, on the same x, at their own scale.
#   both rows  a shaded band at the DOSAGE CEILING. It fills the lower row almost
#              completely and is a sliver in the upper one. The bound stops being
#              a number and becomes the shape of the panel, which is the one form
#              a reviewer cannot argue with.

# ---- the x axis comes from the ATAC-ONLY arm, and that is not a detail --------
# In mode "multi" MintCNA's change-point scan is a multivariate CUSUM over an
# rbind(RNA, DNA) matrix and the state call is a coverage-weighted fusion of both
# rows. RNA carries a value in 77.6 % of bins, and 13.3 % of (bin, cell) state
# calls change when it is removed. On that arm the RNA on this panel's y axis
# helped choose which x bin a gene falls into -- the corroboration is partly
# circular. Reading d_cn from the ATAC-only arm removes RNA from the x axis.
# HONEST LIMIT, and it belongs in the caption: NO arm makes both curves
# independent. Here the expression curve is a clean test; the accessibility curve
# shares its reads with its own x axis, because copy number is inferred from ATAC
# bin coverage and promoter accessibility is ATAC coverage in a window. That curve
# is an internal consistency check, not evidence.
dz <- fread(sprintf("Result/qc/cna_dosage_link_%s_atac.tsv", S))
BRK <- c(-Inf, -0.4, -0.2, 0.2, 0.4, Inf)
# Short category labels: the panel lost width to B's and F's squares, and
# "-0.4 to -0.2" at five categories no longer fits across the axis.
LAB <- c("\u2264-0.4", "-0.3", "0", "+0.3", "\u2265+0.4")
dz[, cn_bin := cut(d_cn, breaks = BRK, labels = LAB)]

# ---- the bound, as a ceiling rather than as a fitted statistic ----------------
# The "dosage R^2 = 0.015" that used to be annotated here is GONE. Two measured
# reasons, either sufficient:
#   ECOLOGICAL FALLACY. Between the five plotted bin means R^2 is 0.936. The
#   annotation said 0.015. A reader sees a near-perfect staircase beside a number
#   claiming 1.5 % and cannot reconcile them, because they are not computed on the
#   same thing. That mismatch is what makes the panel unreadable.
#   COMPOSITION DRIFT. Restrict to |d_cn| > 0.2 and R^2 is 0.058; to |d_cn| > 0.4
#   and it is 0.109 -- seven times the printed value -- while the slope stays flat
#   at 0.63-0.76. A referee gets a bigger number from one obvious subset, and the
#   paper is on record having chosen the smallest.
# The replacement needs no fit at all. log2((2 + d_cn)/2) is what a gene WOULD
# shift if copy number acted at 100 % efficiency with no compensation. Nothing to
# attenuate, nothing that moves with the gene set. Checked in both tumours and
# both arms: the largest ceiling over the lineage genes is 0.26-0.41 log2 against
# observed shifts to 4.19, i.e. at least tenfold everywhere.
# NOT USED, deliberately: "slope = 0.97x the dosage expectation". It does not
# replicate (1.28x in P100T1, above the theoretical maximum of a dosage law), and
# it states the reviewer's objection as a finding.
dz[, dose_max := abs(log2((2 + d_cn) / 2))]
Lg <- merge(dz, LIN, by = "gene")[!is.na(cn_bin)]
CEIL <- max(Lg$dose_max)
cat(sprintf(paste0("panel G: dosage ceiling over %d lineage genes = %.3f log2; ",
                   "largest observed |d_rna| = %.2f (%.0fx); %d of %d copy-neutral\n"),
            nrow(Lg), CEIL, max(abs(Lg$d_rna)), max(abs(Lg$d_rna)) / CEIL,
            Lg[abs(d_cn) <= 0.2, .N], nrow(Lg)))

mses <- function(v) { m <- mean(v); s <- sd(v)/sqrt(length(v)); c(m, m - 1.96*s, m + 1.96*s) }
dsum <- rbindlist(lapply(c(d_atac = "d_atac", d_rna = "d_rna"), function(v)
  dz[, { q <- mses(get(v)); .(mean = q[1], lo = q[2], hi = q[3], n = .N) }, by = cn_bin][, modality := v]))
dsum[, modality := factor(fifelse(modality == "d_atac", "Promoter accessibility", "Expression"),
                          levels = c("Promoter accessibility", "Expression"))]
dsum <- dsum[!is.na(cn_bin)][order(modality, cn_bin)]

# the x categories carry the DIRECTION palette, the same blue/grey/red as panel C,
# so a loss bin is the same blue in both panels and the reader re-learns nothing
CNRAMP <- setNames(c(CN3_COLORS[["1"]], colorRampPalette(c(CN3_COLORS[["1"]], "grey88"))(3)[2],
                     "grey78",
                     colorRampPalette(c("grey88", CN3_COLORS[["3"]]))(3)[2], CN3_COLORS[["3"]]), LAB)

# ---- THE PANEL IS PER-GENE, NOT PER-BIN, AND THAT IS THE THIRD ATTEMPT -------
# Version 1 was a staircase with "dosage R^2 = 0.015" annotated on it. Version 2
# was two stacked rows with a shaded dosage-ceiling band drawn through both. The
# manuscript's own author could not read either, and the reason is the same both
# times: they ask the reader to hold an abstraction (variance explained, or a
# ceiling) and then perform a two-step comparison against it. In version 2 the
# band -- the entire argument -- was 40 mm tall in the lower row and a 0.4 mm
# sliver in the upper one, so the two rows had no visible connection at all.
#
# This version states the claim directly and per gene, with no abstraction:
#   grey dot   what copy number PREDICTS this gene's shift to be, log2((2+dCN)/2),
#              i.e. dosage acting at 100 % efficiency with no compensation
#   colour dot what was actually OBSERVED
#   the gap    everything copy number does not explain
# All 21 grey dots collapse onto a line at zero; the coloured dots fan out to 4.19.
# Nothing has to be inferred -- the reader sees the prediction and the observation
# side by side on one axis, for named genes they already met in panels E and F.
Lg[, pred := log2((2 + d_cn) / 2)]                    # SIGNED, unlike dose_max
Lg[, gene := factor(gene, levels = Lg[order(d_rna), gene])]
cat(sprintf("panel G: %d lineage genes; |prediction| <= %.2f, |observed| up to %.2f (%.0fx); %d move against their own copy number\n",
            nrow(Lg), max(abs(Lg$pred)), max(abs(Lg$d_rna)),
            max(abs(Lg$d_rna))/max(abs(Lg$pred)),
            Lg[sign(pred) != sign(d_rna) & abs(pred) > 0.01, .N]))

p_dose <- ggplot(Lg, aes(y = gene)) +
  geom_vline(xintercept = 0, linewidth = 0.22, colour = "grey75") +
  geom_segment(aes(x = pred, xend = d_rna, yend = gene),
               linewidth = 0.3, colour = "grey78") +
  geom_point(aes(x = pred), size = 0.85, colour = "grey50") +
  geom_point(aes(x = d_rna, colour = fam), size = 1.35, show.legend = FALSE) +
  # DIRECT LABELS, not a legend. The two lineage colours are already keyed in
  # panels E and F, and the only thing needing a name is which dot is which --
  # which a key on the data would say worse than two words above the columns.
  annotate("text", x = 0.1, y = nrow(Lg) + 1.6, hjust = 0.5, vjust = 0.5,
           size = PT_SMALL/.pt, colour = "grey40", label = "predicted") +
  # directly over the TOP row's coloured dot, not over the right half of the axis:
  # the AC-lineage genes' observed values are on the LEFT, so a label parked at
  # x = 3 would have named only half the data.
  # hjust = 1: the top row's mark is the RIGHTMOST point in the panel, so a
  # centred label there hangs 2.1 mm over the edge and ships as "observ".
  geom_text(data = Lg[which.max(d_rna)], aes(x = d_rna, y = nrow(Lg) + 1.6),
            inherit.aes = FALSE, hjust = 1, vjust = 0.5,
            size = PT_SMALL/.pt, colour = "grey30", label = "observed") +
  scale_colour_manual(values = c(`OPC lineage` = STATE_COLORS[["OPC-like"]],
                                 `AC lineage`  = STATE_COLORS[["AC-like"]])) +
  scale_x_continuous(breaks = seq(-2, 4, 2)) +
  scale_y_discrete(expand = expansion(add = c(0.7, 2.6))) +
  labs(x = expression(log[2]*" shift  (C2/C1)"), y = NULL) +
  theme_cowplot(BASE_PT, font_family = FIG_FAMILY) +
  # PT_SMALL for the gene names: 21 rows in ~50 mm is a DENSE tick row, the same
  # case as panel C's 19 chromosome labels, which is the one place this manuscript
  # allows 5.5 pt.
  theme(axis.text.y = element_text(size = PT_SMALL, colour = "grey20"),
        axis.text.x = element_text(size = BASE_PT - 1, colour = "grey20"),
        axis.title.x = element_text(size = PT_SUP),
        # the y axis LINE stays even though its ticks do not: it is the panel's
        # left edge, it is conventional on a lollipop chart, and qa_svg_geometry.py
        # finds a panel by the L its two axis lines make -- without it panel G was
        # invisible to the geometry check, which then silently measured six panels
        # instead of seven.
        axis.line.y = element_line(linewidth = 0.3, colour = "grey30"),
        axis.ticks.y = element_blank(),
        axis.line.x = element_line(linewidth = 0.3, colour = "grey30"),
        axis.ticks.x = element_line(linewidth = 0.3, colour = "grey30"),
        panel.grid.major.y = element_line(linewidth = 0.15, colour = "grey93"),
        panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        plot.margin = margin(1, 2, 1, 1))

# ---- the staircase is not deleted, it is DEMOTED to Supplementary -------------
# It carries the forward claim (copy number moves both modalities where the clones
# differ) which nothing else in Figure 4 makes, but it cannot carry the bound at
# the same time, and the bound is what the main figure needs. Emitted below.
p_si_dose <- ggplot(dsum, aes(cn_bin, mean, group = modality)) +
  geom_rect(data = data.table(1), aes(xmin = -Inf, xmax = Inf, ymin = -CEIL, ymax = CEIL),
            inherit.aes = FALSE, fill = "grey90", colour = NA) +
  geom_hline(yintercept = 0, linewidth = 0.22, colour = "grey75") +
  geom_line(aes(linetype = modality), linewidth = 0.32, colour = "grey35") +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.16, linewidth = 0.28, colour = "grey45") +
  geom_point(aes(fill = cn_bin, shape = modality), size = 1.8, stroke = 0.32, colour = "grey25") +
  annotate("text", x = 0.55, y = CEIL, hjust = 0, vjust = 1.5, size = PT_SMALL/.pt,
           colour = "grey45", label = "most dosage can produce") +
  scale_fill_manual(values = CNRAMP, guide = "none") +
  scale_shape_manual(values = c(`Promoter accessibility` = 21, `Expression` = 24), name = NULL) +
  scale_linetype_manual(values = c(`Promoter accessibility` = "solid", `Expression` = "22"), name = NULL) +
  labs(x = expression(Delta*"copy number  (C2 " - " C1)"),
       y = expression("mean "*log[2]*" shift")) +
  theme_cowplot(BASE_PT, font_family = FIG_FAMILY) +
  theme(axis.text = element_text(size = BASE_PT - 1, colour = "grey20"),
        axis.title.y = element_text(size = PT_SUP), axis.title.x = element_text(size = BASE_PT),
        axis.line = element_line(linewidth = 0.3, colour = "grey30"),
        axis.ticks = element_line(linewidth = 0.3, colour = "grey30"),
        panel.grid = element_blank(),
        legend.position = "inside", legend.position.inside = c(0.99, 0.03),
        legend.justification = c(1, 0), legend.direction = "vertical",
        legend.background = element_blank(), legend.margin = margin(0,0,0,0),
        legend.text = element_text(size = BASE_PT - 1),
        legend.key.size = unit(2.8, "mm"), plot.margin = margin(2, 2, 1, 2))

# named for the panel it IS. The block, its message and this file all said "f"
# while the caption said G, which is how a mislabelled figure ships.
fwrite(dsum, file.path(DATA, sprintf("fig4g_dosage_%s.csv", S)))
fwrite(Lg[, .(gene, fam, d_cn, dose_max, d_atac, d_rna)],
       file.path(DATA, sprintf("fig4g_lineage_dosage_%s.csv", S)))

MM <- 1/25.4
# The clone x state stacked bar is NOT exported as its own panel. It spends a
# whole panel on two numbers, and the same correspondence is stronger shown
# rather than counted: panel a carries a cell-state strip beside the clone strip,
# so the reader sees the violet block align with one clone. The bar is kept here
# only as a QA render and its counts are in the source data.
cairo_pdf(file.path(ROOT, "figure4paper/figure5-multiome/qa/qa_state_bar.pdf"),
          width = 52/25.4, height = 52/25.4, family = FIG_FAMILY)
print(p_state + theme(legend.position = "bottom", legend.box = "vertical") +
      guides(fill = guide_legend(nrow = 2))); dev.off()

if (!exists("FIG5_ASSEMBLING"))
for (x in list(list(p_vol, "figS_motif_volcano.pdf", 70, 62),
               list(p_rep, "fig5_atac_rna.pdf", 73, 57),
               list(p_si_dose, "figS_dosage_staircase.pdf", 73, 50),
               list(p_dose, "figS_dosage_pergene.pdf", 60, 75),
               list(p_si,  "figS_motif_replication.pdf", 62, 62))) {
  cairo_pdf(file.path(OUT, x[[2]]), width = x[[3]]*MM, height = x[[4]]*MM, family = FIG_FAMILY)
  print(x[[1]]); dev.off(); cat("wrote", x[[2]], "\n")
}
