library(data.table)
cell_names <- colnames(cna_multicna_ae)

label <- data.frame(
  SRR = cell_names,
  Site = sub("^CRC\\d+_([A-Za-z]+)[0-9]*_\\d+$", "\\1", cell_names),
  stringsAsFactors = FALSE
)

ord <- order(label$Site)

load("<MINTCNA_ROOT>/Result/CRC/Benchmark/MultiCNA_AE_197_d250.RData")

ref_multicna_ae$chr <- sub("^chr", "", ref_multicna_ae$chr)

## Plot
library(ComplexHeatmap)
library(circlize)
library(grid)

# ---------------------------
# match labels
# ---------------------------
label <- label[match(colnames(cna_multicna_ae), label$SRR), ]
label$Site <- factor(label$Site, levels = c("LN", "PT", "ML", "MP"))
ord_cell <- order(label$Site)

# ---------------------------
# clean chromosome variable
# ---------------------------
ref_multicna_ae$chr <- sub("^chr", "", as.character(ref_multicna_ae$chr))
ref_multicna_ae$chr[!ref_multicna_ae$chr %in% as.character(1:22)] <- NA

chr_levels <- as.character(1:22)
keep_bin <- !is.na(ref_multicna_ae$chr)

ref_multicna_ae <- ref_multicna_ae[keep_bin, , drop = FALSE]
cna_multicna_ae <- cna_multicna_ae[keep_bin, , drop = FALSE]

ref_multicna_ae$chr <- factor(ref_multicna_ae$chr, levels = chr_levels)

# ---------------------------
# order bins
# ---------------------------
ord_bin <- order(ref_multicna_ae$chr, ref_multicna_ae$start, ref_multicna_ae$end)

mat_plot <- t(cna_multicna_ae[, ord_cell, drop = FALSE])
mat_plot <- mat_plot[, ord_bin, drop = FALSE]
ref_ord <- ref_multicna_ae[ord_bin, , drop = FALSE]
label_ord <- label[ord_cell, , drop = FALSE]

# ---------------------------
# split by chromosome
# ---------------------------
col_split <- droplevels(ref_ord$chr)
slice_levels <- levels(col_split)

# ---------------------------
# left site blocks
# ---------------------------
site_col <- c(
  LN = "#4DAF4A",
  PT = "#377EB8",
  ML = "#E41A1C",
  MP = "#984EA3"
)

row_split <- label_ord$Site

ha_left <- rowAnnotation(
  Site = anno_block(
    gp = gpar(fill = site_col[levels(row_split)], col = NA),
    labels = levels(row_split),
    labels_gp = gpar(col = "white", fontsize = 10, fontface = "bold")
  ),
  width = unit(5, "mm"),
  show_annotation_name = FALSE
)

# ---------------------------
# top chromosome blocks
# ---------------------------
chr_fill <- rep(c("#BDBDBD", "#000000"), 11)
names(chr_fill) <- chr_levels

ha_top <- HeatmapAnnotation(
  chr = anno_block(
    gp = gpar(fill = unname(chr_fill[slice_levels]), col = NA),
    labels = slice_levels,
    labels_gp = gpar(
      col = ifelse(chr_fill[slice_levels] == "#000000", "white", "black"),
      fontsize = 9,
      fontface = "bold"
    )
  ),
  show_annotation_name = FALSE,
  annotation_height = unit(4, "mm")
)

# ---------------------------
# 5-state colors, matched to your example
# 1 = Del.d, 2 = Del.s, 3 = Norm, 4 = Dup.s, 5 = Dup.d
# ---------------------------
cn_cols <- c(
  "1" = "#2C7BB6",  # dark blue
  "2" = "#43B5E8",  # light blue
  "3" = "#F3ECDD",  # norm
  "4" = "#F4A582",  # light salmon
  "5" = "#B2182B"   # dark red
)

## discrete color mapping (not continuous)
cn_col_discrete <- structure(
  c("#2C7BB6", "#43B5E8", "#F3ECDD", "#F4A582", "#B2182B"),
  names = as.character(1:5)
)

## custom discrete legend
cn_legend <- Legend(
  labels = c("Del", "Loss", "Nor", "Gain", "Amp"),
  legend_gp = gpar(fill = cn_col_discrete),
  title = "CN",
  direction = "horizontal",
  nrow = 1
)

# ---------------------------
# heatmap
# ---------------------------
ht <- Heatmap(
  mat_plot,
  name = "CN",
  col = cn_col_discrete,
  top_annotation = ha_top,
  left_annotation = ha_left,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  show_column_dend = FALSE,
  show_row_dend = FALSE,
  row_split = row_split,
  row_title = NULL,
  row_gap = unit(0, "mm"),
  column_split = col_split,
  column_gap = unit(0, "mm"),
  border = FALSE,
  rect_gp = gpar(col = NA),
  use_raster = TRUE,
  raster_quality = 2,
  column_title = NULL,
  show_heatmap_legend = FALSE
)

pdf("<MINTCNA_ROOT>/Result/CRC/Segmentation/Evaluation/MintCNA.pdf",
    width = 13, height = 10)
draw(ht,
     heatmap_legend_list = list(cn_legend),
     heatmap_legend_side = "bottom")

## draw thin vertical lines at chromosome boundaries (no horizontal lines)
n_row_slices <- length(levels(row_split))
n_col_slices <- length(slice_levels)
for (ri in seq_len(n_row_slices)) {
  for (ci in seq_len(n_col_slices)) {
    decorate_heatmap_body("CN", slice = ri, column_slice = ci, {
      ## left border
      grid.lines(x = unit(c(0, 0), "npc"), y = unit(c(0, 1), "npc"),
                 gp = gpar(col = "grey40", lwd = 0.3))
      ## right border
      grid.lines(x = unit(c(1, 1), "npc"), y = unit(c(0, 1), "npc"),
                 gp = gpar(col = "grey40", lwd = 0.3))
    })
  }
}
dev.off()
