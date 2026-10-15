#!/usr/bin/env Rscript
# QA: does tuning rescue epiAneufinder or inferCNV?
#
# Both produce almost nothing at their published defaults -- recall 0.002 and
# 0.012 against 49,000 truly altered bins -- and the question is whether that is
# a parameter choice or the method on this data. The tuned runs already exist on
# disk (epiAneufinder w1mb_mf5000, inferCNV cutoff 0.1 + denoise), so this draws
# both configurations of each, same sample, same discretisation, same scale.
.libPaths(c("<R_LIBS_USER>", .libPaths()))
suppressPackageStartupMessages({
  library(ComplexHeatmap); library(grid); library(data.table)
})
ROOT <- "<MINTCNA_ROOT>"
source("<REPRODUCE_DIR>/07_figures/shared/fig_palette.R")
setwd(file.path(ROOT, "ATAC_RNA")); source("<REPRODUCE_DIR>/lib/eval/pearson_helpers.R")
OUT <- file.path(ROOT, "figure4paper/figure5-multiome/qa")
S <- "P111T3"; CHRS <- paste0("chr", 1:22); LO <- log2(0.75); HI <- log2(1.25)
truth <- fread(sprintf("Result/truth/%s_truth_100kb.tsv.gz", S))[chr %in% CHRS]
keep  <- fread(sprintf("Data/cell_lists/%s_cells.tsv", S))[cell_type == "Malignant", barcode]
g   <- paste(truth$chr, truth$start %/% 5e5)
grp <- split(seq_along(g), factor(g, levels = unique(g)))
key <- do.call(rbind, strsplit(names(grp), " "))
chrf <- factor(sub("^chr", "", key[, 1]), levels = as.character(1:22)); ord <- order(chrf, as.numeric(key[,2]))

CFG <- list(c("epiAneufinder","EPI_TAG","w100kb_mf20000","as published"),
            c("epiAneufinder","EPI_TAG","w1mb_mf5000","tuned  (1 Mb window, minFrags 5000)"),
            c("inferCNV","INFERCNV_TAG","default","as published"),
            c("inferCNV","INFERCNV_TAG","cut01_denoiseT","tuned  (cutoff 0.1, denoise)"))
pl <- list()
tvv <- truth$cnv_call
tb <- vapply(grp, function(j){ z <- tvv[j]; z <- z[z %in% c("-","0","+")]
        if(!length(z)) NA_integer_ else c("-"=1L,"0"=2L,"+"=3L)[names(which.max(table(z)))] }, integer(1))
pl[["bulk WGS"]] <- list(B = matrix(rep(tb[ord], each = 40), nrow = 40),
                         name = "bulk WGS", stat = "reference")
for (cf in CFG) {
  do.call(Sys.setenv, setNames(list(cf[3]), cf[2]))
  M <- percell_matrices(S, truth, keep); X <- M[[cf[1]]]
  if (is.null(X)) { cat(cf[1], cf[3], ": absent\n"); next }
  D <- ifelse(X < LO, 1L, ifelse(X > HI, 3L, 2L)); D[!is.finite(X)] <- NA_integer_
  B <- vapply(grp, function(j){ z <- D[, j, drop=FALSE]
    n1 <- rowSums(z==1L, na.rm=TRUE); n3 <- rowSums(z==3L, na.rm=TRUE); n2 <- rowSums(z==2L, na.rm=TRUE)
    as.integer(ifelse(n1 > n2 & n1 >= n3, 1L, ifelse(n3 > n2 & n3 > n1, 3L, 2L))) }, integer(nrow(D)))
  B <- B[, ord, drop = FALSE]
  B <- B[order(rowMeans(B != 2L, na.rm=TRUE), decreasing = TRUE), , drop = FALSE]
  agree <- apply(B, 2, function(z) max(mean(z==1L,na.rm=TRUE), mean(z==3L,na.rm=TRUE)))
  k <- sprintf("%s %s", cf[1], cf[4])
  pl[[k]] <- list(B = B, name = sprintf("%s\n%s", cf[1], cf[4]),
    stat = sprintf("%d cells   %.1f%% of cell-bins called   max agreement %.2f",
                   nrow(B), 100*mean(B != 2L, na.rm=TRUE), max(agree, na.rm=TRUE)))
  cat(sprintf("%-14s %-38s %s\n", cf[1], cf[4], pl[[k]]$stat))
}
Sys.setenv(EPI_TAG = "w100kb_mf20000", INFERCNV_TAG = "default")

chr_lv <- levels(droplevels(chrf[ord]))
fill <- rep(c("#C9CDD3","#5A6069"), length.out = length(chr_lv)); names(fill) <- chr_lv
MM <- 1/25.4
cairo_pdf(file.path(OUT, sprintf("qa_tuned_%s.pdf", S)),
          width = 180*MM, height = 30*MM*length(pl), family = FIG_FAMILY)
hl <- NULL
for (k in names(pl)) {
  ha <- HeatmapAnnotation(chr = anno_block(gp = gpar(fill = unname(fill[chr_lv]), col = NA),
         labels = chr_lv, labels_gp = gpar(col = rep(c("grey15","white"), length.out=length(chr_lv)),
         fontsize = 5, fontface = "bold")), show_annotation_name = FALSE,
         annotation_height = unit(2.6, "mm"))
  at <- HeatmapAnnotation(which = "column", show_annotation_name = FALSE,
         stat = anno_text(c(pl[[k]]$stat, rep("", ncol(pl[[k]]$B) - 1)),
                          gp = gpar(fontsize = 5.2, col = "grey35"), just = "left",
                          location = unit(0, "npc"), height = unit(2.6, "mm")))
  h <- Heatmap(pl[[k]]$B, col = CN3_COLORS, na_col = "white",
               top_annotation = if (k == names(pl)[1]) ha else at,
               cluster_rows = FALSE, cluster_columns = FALSE,
               show_row_names = FALSE, show_column_names = FALSE,
               column_split = droplevels(chrf[ord]), column_gap = unit(0, "mm"),
               row_title = pl[[k]]$name, row_title_rot = 0, row_title_side = "left",
               row_title_gp = gpar(fontsize = 6.5), border = FALSE, rect_gp = gpar(col = NA),
               use_raster = TRUE, raster_quality = 8, raster_device = "CairoPNG",
               column_title = NULL, show_heatmap_legend = (k == names(pl)[1]),
               heatmap_legend_param = list(title = "State", at = 1:3, labels = CN3_LABELS,
                 direction = "horizontal", title_gp = gpar(fontsize = 6, fontface = "bold"),
                 labels_gp = gpar(fontsize = 6)))
  hl <- if (is.null(hl)) h else hl %v% h
}
draw(hl, heatmap_legend_side = "bottom")
dev.off()
cat("wrote", file.path(OUT, sprintf("qa_tuned_%s.pdf", S)), "\n")
