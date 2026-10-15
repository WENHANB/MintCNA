#!/usr/bin/env Rscript
# Figure 3a -- single-cell CNA heatmap, cells split by MintCNA subclone.
# REVISION of 10_subclone_heatmap.R.  Only the colour assignment changes; the
# clustering, ordering, CN encoding and geometry are the original's.
#
# What was wrong: three colour channels were competing inside one 180 mm panel.
#   heatmap body   blue -> red     Del .. Amp          (direction)
#   Site strip     blue/green/orange/red               (categorical metadata)
#   Clone blocks   green #1B9E77 / orange #D95F02      (categorical metadata)
# The Site strip sits 2 mm from the CN ramp and reuses BOTH of its poles, so a
# blue Site cell and a blue deletion are the same blue with different meanings.
# Worse, the two Clone hues are the two method hues from panel b: #1B9E77 is
# MintCNA there, and #D95F02 is one JND from rcCAE's #D55E00. Inside one figure
# the same green meant "our method" and "subclone C1".
#
# Fix: direction keeps blue/red exclusively. All categorical metadata moves to a
# violet/ochre family that appears nowhere else in the figure. Site becomes one
# ordered violet ramp, darkest at PT, so the strip reads as a progression from
# the primary rather than as four unrelated flags.

.libPaths(c("<R_LIBS_USER>", .libPaths()))
suppressPackageStartupMessages({
  library(ComplexHeatmap); library(circlize); library(grid); library(data.table)
  library(mclust); library(cluster)
})
ROOT <- "<MINTCNA_ROOT>"
source(file.path(ROOT, "Result/CRC/figures/scripts/utils.R"))
source(file.path(ROOT, "MintCNA_git/MintCNA/call.R"))   # cluster_latent_kmeans moved here in the repo refactor
source(file.path(ROOT, "figure4paper/shared/fig_palette.R"))
OUT <- file.path(ROOT, "figure4paper/figure4-CRC/output")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

SFX <- Sys.getenv("CRC_SUFFIX", "_d250")          # default "_d250" = 2026-10-14 re-run (AG-DCAE latent_dim 250, 197 cells); "" = the pre-2026-10-14 files
e <- new.env(); load(file.path(BENCH_DIR, paste0("MultiCNA_AE_197_st5", SFX, ".RData")), envir = e)
cna <- e$cna_multicna_ae; ref <- e$ref_multicna_ae
EXCL <- strsplit(Sys.getenv("CRC_EXCLUDE", "CRC01_ML3_085"), ",")[[1]]; EXCL <- EXCL[nzchar(EXCL)]   # per-cell QC exclusion, see utils.R
if (length(EXCL)) cna <- cna[, !colnames(cna) %in% EXCL, drop = FALSE]
NEU <- 3L; stopifnot(all(cna %in% 1:5))

LATENT_FILE <- if (SFX == "_d250") "latent_z_dna_d250_197.csv" else "latent_z_dna.csv"
latent <- as.matrix(read.csv(file.path(ROOT, "Data/CRC/DNAmeth/processed_l2r", LATENT_FILE),
                             row.names = 1, check.names = FALSE))
# scale_Z = FALSE. The packaged default z-scores every latent dimension before
# PCA, which is wrong for an autoencoder embedding: all 200 dimensions come out of
# the same encoder and already share a scale, so their variance spread (17x between
# the largest and smallest here, with the 50 quietest dimensions holding only 12%
# of total variance) IS the signal-to-noise ordering. Standardising inflates the
# near-dead dimensions to unit variance and PCA then spends components on them.
# Measured against the chr5-loss subclone defined from the CNA data itself:
#   scale_Z = TRUE   ARI 0.460, precision 0.55, small-cluster silhouette -0.201
#   scale_Z = FALSE  ARI 0.626, precision 0.76, small-cluster silhouette -0.104
# Strictly better on every axis, including K = 2 stability (0.781 vs 0.755).
if (Sys.getenv("CRC_CLUSTER", "documented") == "documented") {
  # Supplementary Methods procedure, no outlier removal: standardised latent, PCA >= 95 %
  # variance capped at 20 PCs, K by subsampling stability (200 x 80 %, mean pairwise ARI,
  # one-SE rule), final k-means with 200 starts -- the package defaults.
  latent <- latent[intersect(rownames(latent), colnames(cna)), , drop = FALSE]
  # K by the stability rule. With CRC01_ML3_085 removed by the per-cell QC rule (CRC_EXCLUDE,
  # see utils.R) it selects K = 2 on the 196 cells; the guard below makes any other outcome fail
  # loudly instead of silently drawing a different number of clones.
  r1 <- cluster_latent_kmeans(latent, seed = 1)
  cat(sprintf("documented procedure: K_hat = %d, sizes %s\n", r1$K_hat, paste(table(r1$clusters), collapse = "/")))
  if (r1$K_hat != 2) stop("documented clustering selected K = ", r1$K_hat, ", not 2 -- report before changing the figure")
  cl <- r1$clusters
} else {
r0 <- cluster_latent_kmeans(latent, seed = 1, K_min = 2, K_max = 2, scale_Z = FALSE)
sz <- table(r0$clusters)
outlier <- names(r0$clusters)[r0$clusters == as.integer(names(sz)[which.min(sz)])]
keep <- setdiff(rownames(latent), outlier)
r1 <- cluster_latent_kmeans(latent[keep, , drop = FALSE], seed = 1, K_min = 2, K_max = 2,
                            scale_Z = FALSE)
cl <- r1$clusters
}
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
# Cell order: clone block -> site block -> hierarchical clustering INSIDE each
# clone x site block. The site grouping the reader relies on is untouched; what
# changes is the order within a block, which was previously the arbitrary original
# column order and made shared events read as scattered streaks instead of solid
# bands. Ward linkage on the 5-state profile puts similar cells side by side, so a
# real event shows up as a contiguous block and a one-cell artefact stays a single
# line. row_split stays on clone only, so there are still just two slices.
oc <- unlist(lapply(levels(clone), function(k)
  unlist(lapply(levels(site), function(s) {
    idx <- which(clone == k & site == s)
    if (length(idx) <= 2L) return(idx)
    idx[hclust(dist(t(cna[, idx, drop = FALSE])), method = "ward.D2")$order]
  }))), use.names = FALSE)
stopifnot(length(oc) == ncol(cna), !anyDuplicated(oc))
mat <- t(cna[, oc, drop = FALSE]); cl_ord <- clone[oc]; st_ord <- site[oc]
cat("cells:", ncol(cna), "  C1 =", sum(clone == "C1"), " C2 =", sum(clone == "C2"), "\n")
# report how well the clusters line up with the chr5 event they are meant to capture
{
  i5 <- which(ref$chr == "5")
  c5 <- colMeans(cna[i5, , drop = FALSE] < NEU)
  cat(sprintf("chr5 lost (>50%% of the chromosome) in %d of %d cells;  C2 %d/%d, C1 %d/%d\n",
              sum(c5 > 0.5), length(c5),
              sum(c5[clone == "C2"] > 0.5), sum(clone == "C2"),
              sum(c5[clone == "C1"] > 0.5), sum(clone == "C1")))
}

## ---- geometry: written at FINAL size, placed at 100%, never rescaled ---------
MM <- 1/25.4
W_MM <- 117; H_MM <- 93          # left column of Fig 3 row 1 (117 + 66 = 183 mm)
FS_CHR <- PT_SMALL; FS_CLONE <- 7; FS_LGT <- 6.5; FS_LGL <- 6; FS_ANN <- 6

ha_left <- rowAnnotation(
  Clone = anno_block(gp = gpar(fill = CLONE_COLORS2[levels(cl_ord)], col = NA),
                     labels = levels(cl_ord),
                     # ink per block: white on the dark one, near-black on the
                     # light one. Matches Fig 4a. A single "white" put the second
                     # clone label at about 2.7:1 once CLONE_COLORS2 became two
                     # neutrals, which it had to become because the old violet C1
                     # collided with this panel's violet Site ramp.
                     labels_gp = gpar(col = c("white", "grey10")[seq_along(levels(cl_ord))],
                                      fontsize = FS_CLONE, fontface = "bold")),
  Site  = anno_simple(as.character(st_ord), col = SITE_COLORS2, width = unit(3, "mm")),
  # No annotation names: they render at the bottom of the strip and collide with
  # the legend row directly beneath, and the legends already say Subclone / Site.
  width = unit(9, "mm"), show_annotation_name = FALSE, gap = unit(0.8, "mm"))

chr_lv <- levels(droplevels(ref$chr))
# Chromosome labels are thinned by RENDERED WIDTH, not by an odd/even rule. Each
# block is as wide as that chromosome's share of the 25,262 bins; at 117 mm the
# short arms 19-22 are narrower than their own two-digit label, which is what
# makes them run together as "19202122". Anything too narrow gets a blank label
# and is identified by its alternating block instead.
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
lg_cn <- Legend(labels = CN5_LABELS, legend_gp = gpar(fill = unname(CN5_COLORS)),
                title = "Copy number", direction = "horizontal", nrow = 1,
                # title inline with the keys, matching Fig 4a. Stacked titles cost
                # a whole text line per legend, which Fig 4a could not afford, and
                # two figures should not spell a legend two ways.
                title_position = "leftcenter",
                grid_height = unit(2.4, "mm"), grid_width = unit(2.4, "mm"),
                # column_gap 2 -> 1.4 mm (2026-10-09): "Nor" became "Neutral" and the
                # row overflowed the 117 mm panel, clipping "MP"; tighter item gaps
                # and a 2.5 mm gap between the two legends bring it back inside.
                column_gap = unit(1.4, "mm"),
                title_gp = LGT, labels_gp = LGL)
# No Subclone legend: C1 and C2 are written directly on their blocks, and the
# skill's Nature sample notes that legends are dropped where direct labelling works.
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
              # raster_quality 12: the body is 196 cells x 25,262 bins = 5M rectangles, far too
              # many to emit as vector paths, so it is rasterised. 12x the device resolution
              # puts the raster at ~7,200 dpi-equivalent for this panel width, which is well
              # past the 600 dpi combination-artwork requirement and keeps single-bin events
              # (one 100 kb column) from being averaged away.
              use_raster = TRUE, raster_quality = 12, raster_device = "CairoPNG",
              column_title = NULL, show_heatmap_legend = FALSE)

drawit <- function() draw(ht, heatmap_legend_list = list(lg_cn, lg_st),
                          heatmap_legend_side = "bottom", legend_grouping = "original")

if (!exists("FIG4_ASSEMBLING")) {
# family = on the device: ComplexHeatmap's gpar() calls carry no fontfamily and
# inherit the device default. Without this, panel a rendered in DejaVu Sans while
# b and c rendered in Nimbus Sans, and both got embedded in the same PDF.
cairo_pdf(file.path(OUT, "fig4a_heatmap.pdf"), width = W_MM*MM, height = H_MM*MM,
          family = FIG_FAMILY); drawit(); dev.off()
cat(sprintf("wrote fig4a_heatmap.pdf  (%d x %d mm)\n", W_MM, H_MM))
}
