#!/usr/bin/env Rscript
# Figure 4a: single-cell copy-number heatmap of the CRC01 scTrio-seq2 cells, split by MintCNA subclone.
# Input (reproduce/figures/data/fig4_crc/): MultiCNA_AE_197_st5_d250.RData, the MintCNA 5-state calls
#   with their bin coordinates, and latent_z_dna.csv, the AG-DCAE latent embedding the clones are
#   clustered on.
# Output: reproduce/figures/output/fig4/fig4a_heatmap.pdf (117 x 93 mm).
# Run from the root of the MintCNA repository:
#   Rscript reproduce/figures/fig4_crc/panel_a_heatmap.R (or through fig4_crc/assemble_fig4.R)

suppressPackageStartupMessages({
  library(ComplexHeatmap); library(circlize); library(grid); library(data.table)
  library(mclust); library(cluster)
})
FIG <- "reproduce/figures"   # run from the root of the MintCNA repository
if (!file.exists(file.path(FIG, "shared/fig_palette.R"))) stop("run this script from the root of the MintCNA repository")
source(file.path(FIG, "fig4_crc/utils.R"))
source("MintCNA/call.R")                      # cluster_latent_kmeans()
source(file.path(FIG, "shared/fig_palette.R"))
OUT <- file.path(FIG, "output/fig4"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# The calls are the AG-DCAE latent_dim 250 run (197 cells), with CRC01_ML3_085 removed by the
# per-cell QC rule described in utils.R.
e <- new.env(); load(file.path(BENCH_DIR, "MultiCNA_AE_197_st5_d250.RData"), envir = e)
cna <- e$cna_multicna_ae; ref <- e$ref_multicna_ae
EXCL <- "CRC01_ML3_085"   # per-cell QC exclusion, see utils.R
if (length(EXCL)) cna <- cna[, !colnames(cna) %in% EXCL, drop = FALSE]
NEU <- 3L; stopifnot(all(cna %in% 1:5))

# Clone labels are k-means with K = 2 on the unscaled 200-dimensional AG-DCAE latent in
# latent_z_dna.csv, restricted to the cells present in the call matrix. scale_Z = FALSE because the
# latent dimensions come out of one encoder and already share a scale; standardising them inflates
# the near-dead dimensions to unit variance and the PCA then spends components on them.
LATENT_FILE <- "latent_z_dna.csv"
latent <- as.matrix(read.csv(file.path(BENCH_DIR, LATENT_FILE),
                             row.names = 1, check.names = FALSE))
latent <- latent[intersect(rownames(latent), colnames(cna)), , drop = FALSE]
r1 <- cluster_latent_kmeans(latent, seed = 1, K_min = 2, K_max = 2, scale_Z = FALSE)
cat(sprintf("published clustering: %d cells of the %d-dim latent, K = 2 fixed, unscaled; sizes %s\n", nrow(latent), ncol(latent), paste(table(r1$clusters), collapse = "/")))
cl <- r1$clusters
ord_by_size <- names(sort(table(cl), decreasing = TRUE))
clone <- factor(ifelse(cl == as.integer(ord_by_size[1]), "C1", "C2"), levels = c("C1", "C2"))
names(clone) <- names(cl)

cells <- intersect(colnames(cna), names(clone))
clone <- droplevels(clone[cells])
site  <- factor(sub("^CRC\\d+_([A-Za-z]+)[0-9]*_\\d+$", "\\1", cells), levels = SITE_ORDER)

ref$chr <- sub("^chr", "", as.character(ref$chr))
kb <- ref$chr %in% as.character(1:22)
ref <- ref[kb, , drop = FALSE]; cna <- cna[kb, cells, drop = FALSE]
ref$chr <- factor(ref$chr, levels = as.character(1:22))
ob <- order(ref$chr, ref$start); ref <- ref[ob, , drop = FALSE]; cna <- cna[ob, , drop = FALSE]
# Cell order: clone block, then site block, then Ward clustering of the 5-state profiles inside each
# clone x site block, so that a shared event appears as a contiguous band rather than as scattered
# streaks. row_split is on clone only, so the heatmap has two slices.
oc <- unlist(lapply(levels(clone), function(k)
  unlist(lapply(levels(site), function(s) {
    idx <- which(clone == k & site == s)
    if (length(idx) <= 2L) return(idx)
    idx[hclust(dist(t(cna[, idx, drop = FALSE])), method = "ward.D2")$order]
  }))), use.names = FALSE)
stopifnot(length(oc) == ncol(cna), !anyDuplicated(oc))
mat <- t(cna[, oc, drop = FALSE]); cl_ord <- clone[oc]; st_ord <- site[oc]
cat("cells:", ncol(cna), "  C1 =", sum(clone == "C1"), " C2 =", sum(clone == "C2"), "\n")
# Agreement between the clone labels and the chr5 loss they are expected to capture.
{
  i5 <- which(ref$chr == "5")
  c5 <- colMeans(cna[i5, , drop = FALSE] < NEU)
  cat(sprintf("chr5 lost (>50%% of the chromosome) in %d of %d cells;  C2 %d/%d, C1 %d/%d\n",
              sum(c5 > 0.5), length(c5),
              sum(c5[clone == "C2"] > 0.5), sum(clone == "C2"),
              sum(c5[clone == "C1"] > 0.5), sum(clone == "C1")))
}

## ---- geometry: drawn at final size, never rescaled ---------------------------
MM <- 1/25.4
W_MM <- 117; H_MM <- 93          # left column of the top row of Fig. 4 (117 + 66 = 183 mm)
FS_CHR <- PT_SMALL; FS_CLONE <- 7; FS_LGT <- 6.5; FS_LGL <- 6; FS_ANN <- 6

# Copy-number direction keeps blue/red to itself; the Clone and Site strips use a separate
# neutral/violet family, so no categorical colour repeats a copy-number colour.
ha_left <- rowAnnotation(
  Clone = anno_block(gp = gpar(fill = CLONE_COLORS2[levels(cl_ord)], col = NA),
                     labels = levels(cl_ord),
                     # label ink per block: white on the dark clone colour, near-black on the light one
                     labels_gp = gpar(col = c("white", "grey10")[seq_along(levels(cl_ord))],
                                      fontsize = FS_CLONE, fontface = "bold")),
  Site  = anno_simple(as.character(st_ord), col = SITE_COLORS2, width = unit(3, "mm")),
  # Annotation names are not drawn; the legends below already name the strips.
  width = unit(9, "mm"), show_annotation_name = FALSE, gap = unit(0.8, "mm"))

chr_lv <- levels(droplevels(ref$chr))
# Chromosome labels are thinned by rendered width. Each block is as wide as that chromosome's share
# of the bins; a block narrower than its own two-digit label gets a blank label and is identified by
# the alternating block shade instead.
bins_per <- as.integer(table(droplevels(ref$chr))[chr_lv])
heat_mm  <- W_MM - 9 - 8                       # panel width minus row annotation and margins
lbl_mm   <- 2.4                                # two digits at FS_CHR plus breathing room
chr_txt  <- chr_lv
chr_txt[bins_per / sum(bins_per) * heat_mm < lbl_mm] <- ""
cat("chromosome labels shown:", paste(chr_txt[chr_txt != ""], collapse = " "), "\n")
chr_fill <- rep(c("#C9CDD3", "#5A6069"), length.out = length(chr_lv)); names(chr_fill) <- chr_lv
ha_top <- HeatmapAnnotation(
  chr = anno_block(gp = gpar(fill = unname(chr_fill[chr_lv]), col = NA), labels = chr_txt,
                   labels_gp = gpar(col = rep(c("grey15", "white"), length.out = length(chr_lv)),
                                    fontsize = FS_CHR, fontface = "bold")),
  show_annotation_name = FALSE, annotation_height = unit(3.2, "mm"))

LGT <- gpar(fontsize = FS_LGT, fontface = "bold"); LGL <- gpar(fontsize = FS_LGL)
# Both legends run horizontally with the title inline with the keys; the tight item gaps keep the
# legend row inside the 117 mm panel. There is no subclone legend: C1 and C2 are written directly
# on their blocks.
lg_cn <- Legend(labels = CN5_LABELS, legend_gp = gpar(fill = unname(CN5_COLORS)),
                title = "Copy number", direction = "horizontal", nrow = 1,
                title_position = "leftcenter",
                grid_height = unit(2.4, "mm"), grid_width = unit(2.4, "mm"),
                column_gap = unit(1.4, "mm"),
                title_gp = LGT, labels_gp = LGL)
lg_st <- Legend(labels = names(SITE_COLORS2), legend_gp = gpar(fill = unname(SITE_COLORS2)),
                title = "Site", direction = "horizontal", nrow = 1,
                title_position = "leftcenter",
                grid_height = unit(2.4, "mm"), grid_width = unit(2.4, "mm"),
                column_gap = unit(1.4, "mm"),
                title_gp = LGT, labels_gp = LGL)
ht_opt$legend_gap <- unit(2.5, "mm")   # gap between the Copy number and Site legends (default 4)

ht <- Heatmap(mat, name = "CN", col = CN5_COLORS,
              top_annotation = ha_top, left_annotation = ha_left,
              cluster_rows = FALSE, cluster_columns = FALSE,
              show_row_names = FALSE, show_column_names = FALSE,
              row_split = cl_ord, row_title = NULL, row_gap = unit(1, "mm"),
              column_split = droplevels(ref$chr), column_gap = unit(0, "mm"),
              border = FALSE, rect_gp = gpar(col = NA),
              # The body (196 cells x 25,262 bins) is rasterised; raster_quality 12 keeps the raster
              # well above 600 dpi at this panel width, so single-bin events are not averaged away.
              use_raster = TRUE, raster_quality = 12, raster_device = "CairoPNG",
              column_title = NULL, show_heatmap_legend = FALSE)

drawit <- function() draw(ht, heatmap_legend_list = list(lg_cn, lg_st),
                          heatmap_legend_side = "bottom", legend_grouping = "original")

if (!exists("FIG4_ASSEMBLING")) {
# The font family is set on the device: ComplexHeatmap's gpar() calls carry no fontfamily and
# inherit the device default.
cairo_pdf(file.path(OUT, "fig4a_heatmap.pdf"), width = W_MM*MM, height = H_MM*MM,
          family = FIG_FAMILY); drawit(); dev.off()
cat(sprintf("wrote fig4a_heatmap.pdf  (%d x %d mm)\n", W_MM, H_MM))
}
