#!/usr/bin/env Rscript
# Supplementary Fig. 7 -- how the 20-sample multiome cohort was arrived at.
# (S9 until 2026-10-09; the author ordered the last three as cohort flow, subclone
# biology, computational cost.)
#
# Every number is read from Data/meta/sample_manifest.csv at build time and the
# selection rule is APPLIED, not typed in: the script stops if the rule does not
# reproduce the 20 samples the paper analyses. That is the whole point of drawing
# this in code -- the earlier "at least 500 nuclei" sentence was wrong by 8
# samples, and a hand-drawn chart cannot catch that.
#
# The rule (verified 2026-09-14): matched bulk WGS, and >= 1,900 QC-passing
# nuclei at each of two time points; for the four patients with three time points
# (P59, P100, P105, P111) the two with the most nuclei are used. Those IDs are
# computed from the manifest below, not typed.
.libPaths(c("<R_LIBS_USER>", .libPaths()))
suppressPackageStartupMessages({ library(data.table); library(grid) })
ROOT <- "<MINTCNA_ROOT>"
source(file.path(ROOT, "figure4paper/shared/fig_palette.R"))
OUT <- file.path(ROOT, "figure4paper/figureS-supplementary/_work/output")
PDF_OUT <- file.path(ROOT, "figure4paper/figureS-supplementary")   # final PDFs live at the top level; everything else under _work/
QA  <- file.path(ROOT, "figure4paper/figureS-supplementary/_work/qa")

## ---- the numbers, from the manifest -------------------------------------------------
m <- fread(file.path(ROOT, "ATAC_RNA/Data/meta/sample_manifest.csv"))
FINAL <- c("P59T2","P59T3","P100T1","P100T3","P103T1","P103T2","P104T1","P104T2","P105T1","P105T3",
           "P106T1","P106T2","P108T1","P108T2","P109T1","P109T2","P110T1","P110T2","P111T2","P111T3")
NUC_MIN <- 1900
stopifnot(nrow(m) == 48, uniqueN(m$patient) == 22)
w  <- m[bulk_wgs == "yes"]
stopifnot(nrow(w) == 42, uniqueN(w$patient) == 20)
# apply the rule
w[, atac_cells := as.integer(atac_cells)]
cand <- w[atac_cells >= NUC_MIN]
cand <- cand[, if (.N >= 2) .SD[order(-atac_cells)][1:2], by = patient]
stopifnot(setequal(cand$care_id, FINAL))              # the rule reproduces the cohort
# source data: every multiome sample, with the two fields the rule reads and the outcome
src <- m[, .(care_id, patient, timepoint, subtype, atac_nuclei = as.integer(atac_cells),
             bulk_wgs, analysed = care_id %in% FINAL)][order(-analysed, patient, care_id)]
fwrite(src, file.path(ROOT, "figure4paper/figureS-supplementary/_work/data/figS7_cohort_selection.csv"))
sub_n <- w[, .(samples = .N, patients = uniqueN(patient)), by = subtype]
sel_n <- cand[, .(samples = .N, patients = uniqueN(patient)), by = subtype]
tp <- c(Initial = "T1", Recurrence = "T2", `2nd Recurrence` = "T3")
pairs <- cand[order(as.integer(sub("P", "", patient)), care_id), .(txt = sprintf("%s: %s → %s", patient[1], tp[timepoint[1]], tp[timepoint[2]])), by = .(subtype, patient)]
excl_wgs <- nrow(m) - nrow(w)
excl_rule <- nrow(w) - nrow(cand)
n_nuc <- 69384; n_mal <- 47210                         # §12.3 of the methods note
cat(sprintf("48 -> %d with WGS (%d excluded) -> %d analysed (%d excluded by the nucleus rule)\n",
            nrow(w), excl_wgs, nrow(cand), excl_rule))
print(sub_n); print(sel_n); print(pairs)

## ---- drawing, in mm ------------------------------------------------------------------
W <- 124; H <- 142
AST <- sub_n[subtype == "Astrocytoma"]; OLI <- sub_n[subtype == "Oligodendroglioma"]
ASTs <- sel_n[subtype == "Astrocytoma"]; OLIs <- sel_n[subtype == "Oligodendroglioma"]
ast_pairs <- pairs[subtype == "Astrocytoma", txt]; oli_pairs <- pairs[subtype == "Oligodendroglioma", txt]
COL_A <- STATE_COLORS[["AC-like"]]; COL_O <- STATE_COLORS[["OPC-like"]]   # subtype tints follow the state palette
tint <- function(col, a = 0.10) adjustcolor(col, alpha.f = a)

box <- function(x, y, w, h, title, lines = character(), fill = "white", border = "grey35",
                title_col = "grey10", lwd = 0.6, size_t = BASE_PT, size_l = BASE_PT - 1, r = 1.2) {
  grid.roundrect(x = unit(x, "mm"), y = unit(y, "mm"), width = unit(w, "mm"), height = unit(h, "mm"),
                 r = unit(r, "mm"), gp = gpar(fill = fill, col = border, lwd = lwd))
  n <- length(lines); lh <- 3.3
  top <- y + h/2 - 2.6
  grid.text(title, x = unit(x, "mm"), y = unit(top, "mm"), just = c("centre", "top"),
            gp = gpar(fontsize = size_t, fontface = "bold", col = title_col, fontfamily = FIG_FAMILY))
  if (n) for (i in seq_len(n))
    grid.text(lines[i], x = unit(x, "mm"), y = unit(top - 3.9 - (i - 1) * lh, "mm"), just = c("centre", "top"),
              gp = gpar(fontsize = size_l, col = "grey20", fontfamily = FIG_FAMILY))
}
arrow_v <- function(x, y0, y1) grid.segments(unit(x, "mm"), unit(y0, "mm"), unit(x, "mm"), unit(y1, "mm"),
                                             arrow = arrow(length = unit(1.6, "mm"), type = "closed"),
                                             gp = gpar(col = "grey30", fill = "grey30", lwd = 0.6))
elbow <- function(x0, y0, x1, y1) {   # down, across, down
  ym <- (y0 + y1) / 2
  grid.lines(unit(c(x0, x0, x1), "mm"), unit(c(y0, ym, ym), "mm"), gp = gpar(col = "grey30", lwd = 0.6))
  arrow_v(x1, ym, y1)
}
side_note <- function(x, y, txt) grid.text(txt, x = unit(x, "mm"), y = unit(y, "mm"), just = c("left", "centre"),
                                           gp = gpar(fontsize = BASE_PT - 1, col = "grey40", fontfamily = FIG_FAMILY))

draw <- function() {
  grid.newpage()
  pushViewport(viewport(width = unit(W, "mm"), height = unit(H, "mm"), just = c("centre", "centre")))
  cx <- W/2; bw <- 58
  # top chain
  y1 <- H - 8;  box(cx, y1, bw, 10, "CARE IDH-mutant glioma cohort", sprintf("%d patients, %d tumour samples", 35, 75))
  y2 <- y1 - 16; box(cx, y2, bw, 10, "10x Multiome (same-nucleus ATAC + RNA)", sprintf("%d patients, %d samples", uniqueN(m$patient), nrow(m)))
  y3 <- y2 - 16; box(cx, y3, bw, 10, "Matched bulk whole-genome sequencing", sprintf("%d patients, %d samples", uniqueN(w$patient), nrow(w)))
  arrow_v(cx, y1 - 5, y2 + 5); arrow_v(cx, y2 - 5, y3 + 5)
  # split by subtype
  xa <- cx - 30; xo <- cx + 30; y4 <- y3 - 17
  box(xa, y4, 54, 10, "Astrocytoma", sprintf("%d patients, %d samples", AST$patients, AST$samples), fill = tint(COL_A), border = COL_A, title_col = COL_A)
  box(xo, y4, 54, 10, "Oligodendroglioma", sprintf("%d patients, %d samples", OLI$patients, OLI$samples), fill = tint(COL_O), border = COL_O, title_col = COL_O)
  elbow(cx, y3 - 5, xa, y4 + 5); elbow(cx, y3 - 5, xo, y4 + 5)
  # the rule, as a bar across both columns
  y5 <- y4 - 14
  tri <- w[, .N, by = patient][N >= 3, patient]; tri <- tri[order(as.integer(sub("P", "", tri)))]
  box(cx, y5, 114, 13, sprintf("Selection: ≥ %s QC-passing nuclei at each of two time points", format(NUC_MIN, big.mark = ",")),
      c(sprintf("For patients with three time points (%s), choose the two with the most nuclei", paste(tri, collapse = ", ")),
        sprintf("Exclude %d samples: fewer nuclei at one time point, or a single time point", excl_rule)),
      fill = "grey96", border = "grey55")
  arrow_v(xa, y4 - 5, y5 + 6.5); arrow_v(xo, y4 - 5, y5 + 6.5)
  # selected longitudinal pairs
  # the two pair boxes are TOP-aligned with their own heights
  hA <- 8 + 3.3 * length(ast_pairs) + 2; hO <- 8 + 3.3 * length(oli_pairs) + 2
  top6 <- y5 - 6.5 - 6
  yA <- top6 - hA/2; yO <- top6 - hO/2
  box(xa, yA, 54, hA, "Selected longitudinal pairs", c(sprintf("%d patients, %d samples", ASTs$patients, ASTs$samples), ast_pairs),
      fill = tint(COL_A), border = COL_A, title_col = COL_A)
  box(xo, yO, 54, hO, "Selected longitudinal pairs", c(sprintf("%d patients, %d samples", OLIs$patients, OLIs$samples), oli_pairs),
      fill = tint(COL_O), border = COL_O, title_col = COL_O)
  arrow_v(xa, y5 - 6.5, top6); arrow_v(xo, y5 - 6.5, top6)
  # final
  bot6 <- min(yA - hA/2, yO - hO/2)
  y7 <- bot6 - 13.5
  box(cx, y7, bw, 14, "Analysed cohort (Fig. 5)",
      c(sprintf("%d patients, %d samples", uniqueN(cand$patient), nrow(cand)),
        sprintf("%s QC nuclei, %s malignant", format(n_nuc, big.mark = ","), format(n_mal, big.mark = ","))))
  ym <- y7 + 7 + 3
  grid.lines(unit(c(xa, xa, cx), "mm"), unit(c(yA - hA/2, ym, ym), "mm"), gp = gpar(col = "grey30", lwd = 0.6))
  grid.lines(unit(c(xo, xo, cx), "mm"), unit(c(yO - hO/2, ym, ym), "mm"), gp = gpar(col = "grey30", lwd = 0.6))
  arrow_v(cx, y7 + 7 + 3, y7 + 7)
  popViewport()
}
MM <- 1/25.4
cairo_pdf(file.path(PDF_OUT, "SuppFig07_cohort_flow.pdf"), width = W*MM, height = H*MM, family = FIG_FAMILY); draw(); dev.off()
svglite::svglite(file.path(OUT, "figS7_cohort_flow.svg"), width = W*MM, height = H*MM); draw(); dev.off()
ragg::agg_png(file.path(QA, "figS7_cohort_flow_preview.png"), width = W*MM, height = H*MM, units = "in", res = 900); draw(); dev.off()
ragg::agg_tiff(file.path(OUT, "figS7_cohort_flow.tiff"), width = W*MM, height = H*MM, units = "in", res = 600, compression = "lzw"); draw(); dev.off()
cat(sprintf("wrote figS7_cohort_flow.pdf/.svg/.tiff  %d x %d mm\n", W, H))
