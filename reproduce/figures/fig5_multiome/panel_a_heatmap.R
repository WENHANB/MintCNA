#!/usr/bin/env Rscript
# Figure 5, panel c: per-nucleus copy-number heatmap of one multiome sample (P111T2), malignant nuclei
# split by MintCNA subclone, with the source study's RNA-derived cell state beside each nucleus.
# Input (reproduce/figures/data/fig5_multiome/): P111T2_results_h15a005e020.RData (MintCNA calls and
# subclones; reproduce/analysis/multiome/run_mintcna.sh), P111T2_cells.tsv (cell types;
# reproduce/analysis/multiome/export_cell_lists.R), snrna_malignant_meta_P111T2.tsv (cell states).
# Output: reproduce/figures/output/fig5/fig5a_heatmap.pdf when run on its own; the assembly captures drawit().
# Run from the root of the MintCNA repository: Rscript reproduce/figures/fig5_multiome/panel_a_heatmap.R
# (or through fig5_multiome/assemble_fig5.R).
#
# Three copy-number states rather than five: at 500 kb with about two fragments per bin per nucleus the
# ATAC arm cannot separate a one-copy loss from a homozygous one. P111T2 is shown because its subclone
# split is not explained by sequencing depth and because each subclone is dominated by one cell state
# (AC-like and OPC-like), which is the premise of panels e and f.
suppressPackageStartupMessages({
  library(ComplexHeatmap); library(circlize); library(grid); library(data.table)
})
FIG <- "reproduce/figures"   # run from the root of the MintCNA repository
if (!file.exists(file.path(FIG, "shared/fig_palette.R"))) stop("run this script from the root of the MintCNA repository")
source(file.path(FIG, "shared/fig_palette.R"))
DATA <- file.path(FIG, "data/fig5_multiome")
OUT  <- file.path(FIG, "output/fig5"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

S <- "P111T2"
TAG <- "h15a005e020"
e <- new.env()
load(sprintf("%s/%s_results_%s.RData", DATA, S, TAG), e)
cna <- get("cna_state_matrix", e); ref <- as.data.frame(get("ref", e))
sub <- get("subclones", e)
stopifnot(all(cna %in% 1:3), ncol(cna) == length(sub))

keep <- fread(sprintf("%s/%s_cells.tsv", DATA, S))[
          cell_type == "Malignant", barcode]
i <- colnames(cna) %in% keep
cna <- cna[, i, drop = FALSE]; sub <- sub[i]

ord_by_size <- names(sort(table(sub), decreasing = TRUE))
clone <- factor(ifelse(sub == as.integer(ord_by_size[1]), "C1", "C2"), levels = c("C1","C2"))

# The source study's RNA-derived cell state, matched by barcode. The subclones are clustered from the
# ATAC-trained AG-DCAE latent and no RNA annotation enters the clustering, so the agreement between the
# clone blocks and the state strip is a check rather than a tautology.
md <- fread(sprintf("%s/snrna_malignant_meta_%s.tsv", DATA, S))
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

# Cell order, three nested levels: subclone block, then cell state in the palette's order, then Ward
# linkage on the CNA profile inside each (clone, state) group. Ward alone inside a clone interleaves
# the states and turns the state strip into speckle; nuclei with no RNA annotation go last.
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

## ---- geometry: final size, drawn at 100% and never rescaled ------------------
# When sourced by the assembly the panel takes its slot size from grid_fig5.R; on its own it is
# written at 117 x 93 mm. W_MM also sets the heatmap body width from which the chromosome labels are thinned.
MM <- 1/25.4
if (exists("FIG5_ASSEMBLING")) { W_MM <- LEFT_W; H_MM <- C_HM_H } else { W_MM <- 117; H_MM <- 93 }
FS_CHR <- PT_SMALL; FS_CLONE <- 7; FS_LGT <- 6.5; FS_LGL <- 6   # text sizes, pt

ha_left <- rowAnnotation(
  Clone = anno_block(gp = gpar(fill = CLONE_COLORS2[levels(cl_ord)], col = NA),
                     labels = levels(cl_ord),
                     # label ink per block: white on the dark block, near-black on the light one
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

# Chromosome labels are dropped where the rendered slice is narrower than the label (2.4 mm), except
# chr1 and chr19, which are always labelled because the 1p/19q codeletion is the event the panel shows.
chr_lv   <- levels(droplevels(ref$chr))
bins_per <- as.integer(table(droplevels(ref$chr))[chr_lv])
heat_mm  <- W_MM - 6.2 - 8
lbl_mm   <- 2.4
chr_txt  <- chr_lv
chr_txt[bins_per / sum(bins_per) * heat_mm < lbl_mm] <- ""
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
                title_position = "leftcenter",   # inline title keeps each legend to one line
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
              # The body is about 1,700 x 5,300 cells, so it is rasterised. raster_quality multiplies
              # 72 dpi; 32 (2,304 dpi) gives about two pixels per 500 kb bin at the drawn size, so a
              # single-bin column is not averaged away.
              use_raster = TRUE, raster_quality = 32, raster_device = "CairoPNG",
              column_title = NULL, show_heatmap_legend = FALSE)

# Both legends are drawn under the heatmap; the sample name is the column title.
lg_all <- packLegend(lg_cn, lg_st, direction = "vertical", gap = unit(1, "mm"))
drawit <- function() {
  args <- list(ht, column_title = S,
               column_title_gp = gpar(fontsize = FS_LGT, fontface = "bold"),
               padding = unit(c(1, 1, 1, 1), "mm"))
  args <- c(args, list(heatmap_legend_list = list(lg_all),
                       heatmap_legend_side = "bottom",
                       legend_grouping = "original"))
  do.call(draw, args)
}

if (!exists("FIG5_ASSEMBLING")) {
  cairo_pdf(file.path(OUT, "fig5a_heatmap.pdf"), width = W_MM*MM, height = H_MM*MM,
            family = FIG_FAMILY); drawit(); dev.off()
  cat(sprintf("wrote fig5a_heatmap.pdf  (%g x %g mm)\n", W_MM, H_MM))
}
