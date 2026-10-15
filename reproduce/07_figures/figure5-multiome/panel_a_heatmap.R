#!/usr/bin/env Rscript
# Figure 4a -- single-cell CNA heatmap for one multiome sample, split by subclone.
#
# Grammar is Figure 3a's: direction (blue<->red) is the only blue/red in the
# figure, subclone metadata takes the violet/ochre family, the chromosome ruler
# thins its labels by RENDERED WIDTH rather than by an odd/even rule, and cells
# are ordered by hierarchical clustering inside each subclone block so a shared
# event reads as a solid band instead of scattered streaks.
#
# Two differences from Fig 3a, both forced by the assay:
#   - three states, not five. At 500 kb with ~2 fragments per bin per nucleus the
#     ATAC arm cannot separate one-copy from homozygous loss, so calling five
#     states would be reporting a resolution the data does not have.
#   - no Site strip. Figure 3 had four anatomical sites from one patient; here
#     one sample is one timepoint from one patient, so subclone is the only
#     metadata channel and it is labelled directly on the block.
#
# SAMPLE CHOICE is a stated criterion, not an aesthetic one, and it is two
# criteria in order:
#   (1) the split must not be sequencing depth. Measured as the AUC of nFrags
#       separating the two clones, P103T1 is 0.98 -- depth alone reproduces that
#       split -- against 0.626 for P111T2, 0.621 for P111T3, 0.649 for P100T1.
#   (2) among those, the clones must be LINEAGE-PURE, because this panel is the
#       premise for the TF-motif panels: P111T2 is C1 95.0% AC-like / C2 87.1%
#       OPC-like (chi-square p = 1.3e-185).
# P111T3 has the larger between-clone CNA difference (a 214-bin run against 46)
# and would make the showier heatmap, but its C1 is a MIXTURE -- 44.9% OPC-like,
# 28.0% NPC-like, 15.3% AC-like -- and the prespecified lineage TFs accordingly
# come out weak and sign-inconsistent there (NFIC::TLX1 +0.25 but NFIA -0.11,
# OLIG1/2 n.s.). Showing the sample whose clones are mixed and then testing
# lineage TFs in a different sample would break the figure's own argument, so
# panel a shows P111T2 and the whole lower half of the figure is one tumour.
.libPaths(c("<R_LIBS_USER>", .libPaths()))
suppressPackageStartupMessages({
  library(ComplexHeatmap); library(circlize); library(grid); library(data.table)
})
ROOT <- "<MINTCNA_ROOT>"
source("<REPRODUCE_DIR>/07_figures/shared/fig_palette.R")
OUT <- file.path(ROOT, "figure4paper/figure5-multiome/output")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

S <- Sys.getenv("FIG5_SAMPLE", "P111T2")
TAG <- Sys.getenv("FIG5_MINT_TAG", "h15a005e020")
e <- new.env()
load(sprintf("%s/ATAC_RNA/Result/MintCNA/%s_multi/%s_results_%s.RData", ROOT, S, S, TAG), e)
cna <- get("cna_state_matrix", e); ref <- as.data.frame(get("ref", e))
sub <- get("subclones", e)
stopifnot(all(cna %in% 1:3), ncol(cna) == length(sub))

keep <- fread(sprintf("%s/ATAC_RNA/Data/cell_lists/%s_cells.tsv", ROOT, S))[
          cell_type == "Malignant", barcode]
i <- colnames(cna) %in% keep
cna <- cna[, i, drop = FALSE]; sub <- sub[i]

ord_by_size <- names(sort(table(sub), decreasing = TRUE))
clone <- factor(ifelse(sub == as.integer(ord_by_size[1]), "C1", "C2"), levels = c("C1","C2"))

# The paper's own RNA-derived malignant cell state, matched by barcode. The
# subclones are clustered from the ATAC-trained AG-DCAE latent -- no RNA annotation enters the
# clustering -- so putting the two strips side by side is a test, not a
# tautology: in this sample C2 is 96.4% AC-like against C1's 15.3%
# (chi-square p = 1.5e-189).
md <- as.data.table(readRDS(file.path(ROOT, "ATAC_RNA/Data/synapse/snrna_malignant_meta_data.RDS")))
md <- md[care_id == S]
md[, bare := sub("^.*-([ACGT]+-[0-9]+)$", "\\1", CellID)]
state <- factor(md$cell_state[match(colnames(cna), md$bare)],
                levels = names(CELLSTATE_COLORS))
cat(sprintf("cell state matched for %d / %d nuclei\n", sum(!is.na(state)), length(state)))

ref$chr <- sub("^chr", "", as.character(ref$chr))
kb <- ref$chr %in% as.character(1:22)
ref <- ref[kb, , drop = FALSE]; cna <- cna[kb, , drop = FALSE]
ref$chr <- factor(ref$chr, levels = as.character(1:22))
ob <- order(ref$chr, ref$start); ref <- ref[ob, , drop = FALSE]; cna <- cna[ob, , drop = FALSE]

# Cell order, three nested levels: subclone block -> cell state -> Ward linkage on
# the CNA profile inside each (clone, state) group.
#
# Ward alone inside a clone interleaves the cell states, so the state strip came
# out as speckle and the heatmap body read as noise even where it was not. States
# run in the palette's own order in every clone, so the strip is directly
# comparable between the two blocks, and cells with no RNA annotation go last
# rather than being scattered through the block.
state_ord <- factor(state, levels = names(CELLSTATE_COLORS))
oc <- unlist(lapply(levels(clone), function(k) {
  unlist(lapply(c(levels(state_ord), NA), function(st) {
    idx <- if (is.na(st)) which(clone == k & is.na(state_ord))
           else            which(clone == k & !is.na(state_ord) & state_ord == st)
    if (length(idx) <= 2L) return(idx)
    idx[hclust(dist(t(cna[, idx, drop = FALSE])), method = "ward.D2")$order]
  }), use.names = FALSE)
}), use.names = FALSE)
stopifnot(length(oc) == ncol(cna), !anyDuplicated(oc))
mat <- t(cna[, oc, drop = FALSE]); cl_ord <- clone[oc]; st_ord <- state[oc]

cat(sprintf("%s: %d malignant nuclei, C1 = %d, C2 = %d, %d bins\n",
            S, ncol(cna), sum(clone=="C1"), sum(clone=="C2"), nrow(ref)))
# 1p and 19q are the diagnostic codeletion; report how completely they are called
for (a in list(c("1", 0, 123.4e6, "1p"), c("19", 26.2e6, 58.6e6, "19q"))) {
  j <- which(ref$chr == a[[1]] & ref$start >= as.numeric(a[[2]]) & ref$end <= as.numeric(a[[3]]))
  if (!length(j)) next
  f <- colMeans(cna[j, , drop = FALSE] == 1L)
  cat(sprintf("  %-3s loss in >50%% of its bins: %.1f%% of nuclei (C1 %.1f%%, C2 %.1f%%)\n",
              a[[4]], 100*mean(f > 0.5), 100*mean(f[clone=="C1"] > 0.5), 100*mean(f[clone=="C2"] > 0.5)))
}

## ---- geometry: final size, placed at 100%, never rescaled --------------------
MM <- 1/25.4
W_MM <- as.numeric(Sys.getenv("FIG4A_W", "117"))
H_MM <- as.numeric(Sys.getenv("FIG4A_H", "93"))
FS_CHR <- PT_SMALL; FS_CLONE <- 7; FS_LGT <- 6.5; FS_LGL <- 6   # same constants as Fig 3a; 7 pt = guide maximum (2026-10-09)

ha_left <- rowAnnotation(
  Clone = anno_block(gp = gpar(fill = CLONE_COLORS2[levels(cl_ord)], col = NA),
                     labels = levels(cl_ord),
                     # ink chosen per block: white on the dark one, near-black on
                     # the light one. One fixed "white" put C2 at about 2.7:1.
                     labels_gp = gpar(col = c("white", "grey10")[seq_along(levels(cl_ord))],
                                      fontsize = FS_CLONE, fontface = "bold")),
  State = st_ord,
  col = list(State = CELLSTATE_COLORS),
  na_col = "white",
  annotation_width = unit(c(4, 2.2), "mm"), width = unit(6.2, "mm"),
  show_annotation_name = FALSE, show_legend = FALSE,
  simple_anno_size_adjust = TRUE)

lg_st <- Legend(labels = names(CELLSTATE_COLORS),
                legend_gp = gpar(fill = unname(CELLSTATE_COLORS)),
                title = "Cell state (RNA)", direction = "horizontal", nrow = 1,
                title_position = "leftcenter",
                grid_height = unit(2.4, "mm"), grid_width = unit(2.4, "mm"),
                title_gp = gpar(fontsize = FS_LGT, fontface = "bold"),
                labels_gp = gpar(fontsize = FS_LGL))

chr_lv   <- levels(droplevels(ref$chr))
bins_per <- as.integer(table(droplevels(ref$chr))[chr_lv])
heat_mm  <- W_MM - 6.2 - 8
lbl_mm   <- 2.4
chr_txt  <- chr_lv
chr_txt[bins_per / sum(bins_per) * heat_mm < lbl_mm] <- ""
# Width-based thinning drops 19-22, but chr1 and chr19 carry this figure's claim:
# the 1p/19q codeletion is the diagnostic event the panel exists to show. A label
# the reader needs in order to read the result is not optional, so these two are
# always kept even though chr19 is narrower than its own label. The exception is
# named rather than silently widening the threshold, which would also bring back
# 20-22 and re-create the "19202122" run the rule was written to prevent.
chr_txt[chr_lv %in% c("1","19")] <- chr_lv[chr_lv %in% c("1","19")]
cat("chromosome labels shown:", paste(chr_txt[chr_txt != ""], collapse = " "), "\n")
chr_fill <- rep(c("#C9CDD3", "#5A6069"), length.out = length(chr_lv)); names(chr_fill) <- chr_lv
ha_top <- HeatmapAnnotation(
  chr = anno_block(gp = gpar(fill = unname(chr_fill[chr_lv]), col = NA), labels = chr_txt,
                   labels_gp = gpar(col = rep(c("grey15","white"), length.out = length(chr_lv)),
                                    fontsize = FS_CHR, fontface = "bold")),
  show_annotation_name = FALSE, annotation_height = unit(3.2, "mm"))

lg_cn <- Legend(labels = CN3_LABELS, legend_gp = gpar(fill = unname(CN3_COLORS)),
                title = "Copy number", direction = "horizontal", nrow = 1,
                # Title on the same line as the keys. Stacked title-over-keys for
                # two legends cost ~30% of the panel height; inline titles make
                # each legend one line, and the panel keeps its heatmap.
                title_position = "leftcenter",
                grid_height = unit(2.4, "mm"), grid_width = unit(2.4, "mm"),
                title_gp = gpar(fontsize = FS_LGT, fontface = "bold"),
                labels_gp = gpar(fontsize = FS_LGL))

ht <- Heatmap(mat, name = "CN", col = CN3_COLORS,
              top_annotation = ha_top, left_annotation = ha_left,
              cluster_rows = FALSE, cluster_columns = FALSE,
              show_row_names = FALSE, show_column_names = FALSE,
              row_split = cl_ord, row_title = NULL, row_gap = unit(1, "mm"),
              column_split = droplevels(ref$chr), column_gap = unit(0, "mm"),
              border = FALSE, rect_gp = gpar(col = NA),
              # rasterised: the body is ~1,700 x 5,300 rectangles. 12x device
              # resolution keeps a single 500 kb column from being averaged away.
              # raster_quality is a MULTIPLIER on 72 dpi. At 12 (864 dpi) the widest chromosome
              # slice was 324 px for ~440 bins -- 0.7 px per 500-kb bin, so neighbouring bins
              # were averaged away before the PDF ever saw them. The body is 120 mm for 5,256
              # bins (43.8 bins/mm); two pixels per bin needs 2,225 dpi. 32 x 72 = 2,304 dpi.
              # Rows were already fine at 2.9 px per nucleus. Measured with PyMuPDF on the
              # exported PDF, not assumed.
              use_raster = TRUE, raster_quality = 32, raster_device = "CairoPNG",
              column_title = NULL, show_heatmap_legend = FALSE)

lg_all <- packLegend(lg_cn, lg_st, direction = "vertical", gap = unit(1, "mm"))
# With two stacked heatmaps in panel C the legends are drawn once, under the
# lower one, and the upper one suppresses them. The sample name becomes a title,
# because two heatmaps with no titles is a puzzle.
SHOW_LEG <- Sys.getenv("FIG4C_LEGEND", "TRUE") == "TRUE"
drawit <- function() {
  args <- list(ht, column_title = S,
               column_title_gp = gpar(fontsize = FS_LGT, fontface = "bold"),
               padding = unit(c(1, 1, 1, 1), "mm"))
  if (SHOW_LEG) args <- c(args, list(heatmap_legend_list = list(lg_all),
                                     heatmap_legend_side = "bottom",
                                     legend_grouping = "original"))
  else args <- c(args, list(show_heatmap_legend = FALSE))
  do.call(draw, args)
}

if (!exists("FIG5_ASSEMBLING")) {
  cairo_pdf(file.path(OUT, "fig5a_heatmap.pdf"), width = W_MM*MM, height = H_MM*MM,
            family = FIG_FAMILY); drawit(); dev.off()
  cat(sprintf("wrote fig5a_heatmap.pdf  (%g x %g mm)\n", W_MM, H_MM))
}
