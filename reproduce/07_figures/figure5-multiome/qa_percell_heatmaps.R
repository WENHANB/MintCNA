#!/usr/bin/env Rscript
# QA, not a deliverable: what the consensus step actually sees, per method.
#
# Each panel is the SAME matrix the F1 is computed from -- every method's
# per-cell output discretised by AtaCNV's published rule (copy ratio > 1.25 gain,
# < 0.75 loss) -- so a reader can check by eye why inferCNV and epiAneufinder
# reach no modal call: not because the threshold is wrong for them, but because
# their calls never line up across cells.
.libPaths(c("<R_LIBS_USER>", .libPaths()))
suppressPackageStartupMessages({
  library(ComplexHeatmap); library(grid); library(data.table); library(circlize)
})
ROOT <- "<MINTCNA_ROOT>"
source(file.path(ROOT, "figure4paper/shared/fig_palette.R"))
setwd(file.path(ROOT, "ATAC_RNA")); source("Code/eval/pearson_helpers.R")
OUT <- file.path(ROOT, "figure4paper/figure5-multiome/qa")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

S <- "P111T3"; CHRS <- paste0("chr", 1:22); LO <- log2(0.75); HI <- log2(1.25)
truth <- fread(sprintf("Result/truth/%s_truth_100kb.tsv.gz", S))[chr %in% CHRS]
keep  <- fread(sprintf("Data/cell_lists/%s_cells.tsv", S))[cell_type == "Malignant", barcode]
M <- percell_matrices(S, truth, keep)

# display grid: 500 kb, so 28,760 columns become ~5,300 and single-bin noise is
# not lost to rasterisation
g   <- paste(truth$chr, truth$start %/% 5e5)
grp <- split(seq_along(g), factor(g, levels = unique(g)))
key <- do.call(rbind, strsplit(names(grp), " "))
chrf <- factor(sub("^chr", "", key[, 1]), levels = as.character(1:22))
ord  <- order(chrf, as.numeric(key[, 2]))

# Every method, ordered by pseudobulk Pearson so the panel reads as a descent.
# The point of showing all of them together is that the consensus step can only
# take what is vertically coherent: a method whose calls do not line up across
# cells has nothing for the mode to find, however much it calls per cell.
SHOW <- c("MintCNA", "MintCNA(ATAC)", "AtaCNV", "Copy-scAT",
          "RIDDLER", "epiAneufinder", "inferCNV", "CONGAS+")
SHOW <- SHOW[SHOW %in% names(M)]
pl <- list()
# bulk WGS first, as the reference every panel below is scored against. One row,
# repeated so it renders at a visible height rather than as a hairline.
{
  tvv <- truth$cnv_call
  tb <- vapply(grp, function(j) {
    z <- tvv[j]; z <- z[z %in% c("-", "0", "+")]
    if (!length(z)) NA_integer_ else c("-" = 1L, "0" = 2L, "+" = 3L)[names(which.max(table(z)))]
  }, integer(1))
  B <- matrix(rep(tb[ord], each = 40), nrow = 40)
  pl[["bulk WGS"]] <- list(B = B, name = "bulk WGS", stat = "reference")
  cat("bulk WGS reference row\n")
}
for (m in SHOW) {
  X <- M[[m]]; if (is.null(X)) next
  D <- ifelse(X < LO, 1L, ifelse(X > HI, 3L, 2L)); D[!is.finite(X)] <- NA_integer_
  # majority state per 500 kb block, per cell
  B <- vapply(grp, function(j) {
    z <- D[, j, drop = FALSE]
    n1 <- rowSums(z == 1L, na.rm = TRUE); n3 <- rowSums(z == 3L, na.rm = TRUE)
    n2 <- rowSums(z == 2L, na.rm = TRUE)
    as.integer(ifelse(n1 > n2 & n1 >= n3, 1L, ifelse(n3 > n2 & n3 > n1, 3L, 2L)))
  }, integer(nrow(D)))
  B <- B[, ord, drop = FALSE]
  # order cells by how much they call, so a sparse caller reads as sparse rather
  # than as an arbitrary shuffle
  B <- B[order(rowMeans(B != 2L, na.rm = TRUE), decreasing = TRUE), , drop = FALSE]
  frac <- mean(B != 2L, na.rm = TRUE)
  agree <- apply(B, 2, function(z) max(mean(z == 1L, na.rm=TRUE), mean(z == 3L, na.rm=TRUE)))
  # Short row title, horizontal. The statistics go in a light annotation ABOVE
  # each block: as a rotated row_title they were nine long strings stacked in a
  # 20 mm gutter and overlapped into an unreadable band.
  pl[[m]] <- list(B = B, name = m,
                  stat = sprintf("%d cells   %.1f%% of cell-bins called   max agreement %.2f",
                                 nrow(B), 100*frac, max(agree, na.rm = TRUE)))
  cat(sprintf("%-14s %s\n", m, pl[[m]]$stat))
}

chr_lv <- levels(droplevels(chrf[ord]))
fill   <- rep(c("#C9CDD3", "#5A6069"), length.out = length(chr_lv)); names(fill) <- chr_lv
MM <- 1/25.4
cairo_pdf(file.path(OUT, sprintf("qa_percell_%s.pdf", S)),
          width = 180*MM, height = 26*MM*length(pl), family = FIG_FAMILY)
hl <- NULL
for (m in names(pl)) {
  ha <- HeatmapAnnotation(chr = anno_block(
          gp = gpar(fill = unname(fill[chr_lv]), col = NA), labels = chr_lv,
          labels_gp = gpar(col = rep(c("grey15","white"), length.out = length(chr_lv)),
                           fontsize = 5, fontface = "bold")),
        show_annotation_name = FALSE, annotation_height = unit(2.6, "mm"))
  ann_txt <- HeatmapAnnotation(
    which = "column", show_annotation_name = FALSE,
    stat = anno_text(c(pl[[m]]$stat, rep("", ncol(pl[[m]]$B) - 1)),
                     gp = gpar(fontsize = 5.2, col = "grey35"),
                     just = "left", location = unit(0, "npc"),
                     height = unit(2.6, "mm")))
  h <- Heatmap(pl[[m]]$B, col = CN3_COLORS, na_col = "white",
               top_annotation = if (m == names(pl)[1]) ha else ann_txt,
               cluster_rows = FALSE, cluster_columns = FALSE,
               show_row_names = FALSE, show_column_names = FALSE,
               column_split = droplevels(chrf[ord]), column_gap = unit(0, "mm"),
               row_title = pl[[m]]$name, row_title_gp = gpar(fontsize = 7, fontface = "bold"),
               row_title_rot = 0, row_title_side = "left",
               border = FALSE, rect_gp = gpar(col = NA),
               use_raster = TRUE, raster_quality = 8, raster_device = "CairoPNG",
               column_title = NULL, show_heatmap_legend = (m == names(pl)[1]),
               heatmap_legend_param = list(title = "State", at = 1:3,
                 labels = CN3_LABELS, direction = "horizontal",
                 title_gp = gpar(fontsize = 6, fontface = "bold"),
                 labels_gp = gpar(fontsize = 6)))
  hl <- if (is.null(hl)) h else hl %v% h
}
draw(hl, heatmap_legend_side = "bottom")
dev.off()
cat("wrote", file.path(OUT, sprintf("qa_percell_%s.pdf", S)), "\n")
