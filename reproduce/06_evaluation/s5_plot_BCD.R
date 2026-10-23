#!/usr/bin/env Rscript
# S5 panels B/C/D: precision-recall by CN gain state (CN3/CN4/CN5) for MintCNA (DNA+RNA),
# MintCNA(DNA) and MintCNA(RNA), at coverage 0.05X (B), 0.1X (C), 0.2X (D).
# Same design as Simulation/code/Evaulation/status_F1_1.R (iso-F1 contours, triangle points,
# facet label carries dF1 = DNA+RNA - DNA) but with the third mode added and the axis window
# widened to [0.1,1] so the RNA points (min precision 0.173) are not clipped.
# Data: eval/s5_cn_pr_3cov.csv (verified; 0.1X reproduces the published numbers exactly).
.libPaths(c("<R_LIBS_USER>", .libPaths()))
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(ggh4x); library(patchwork) })
RB  <- "<MINTCNA_ROOT>/Simulation/cnasim/c3_N100_100kb_l1e4/RNA_benchmark"
FIG <- file.path(RB,"figures")

d <- fread(file.path(RB,"eval/s5_cn_pr_3cov.csv"))
d[, method_lab := fifelse(method=="DNA+RNA","MintCNA",
                   fifelse(method=="DNA-only","MintCNA(DNA)","MintCNA(RNA)"))]
d[, coverage := factor(coverage, levels=c("0.05X","0.1X","0.2X"))]

## facet labels: dF1 = MintCNA - MintCNA(DNA), as in status_F1_1.R.
## The F1 column in the csv is rounded to 3dp; subtracting two ROUNDED F1s drifts the label by up
## to 0.001 from the published panels. status_F1_1.R differences the UNROUNDED F1s, so recompute
## F1 at full precision from the raw TP/FP/FN. This reproduces the published labels exactly for
## 0.05X and 0.2X (7 of 9 facets); the two that still differ are 0.1X CN4/CN5 -- see NOTE below.
d[, f1_exact := { p <- TP/(TP+FP); r <- TP/(TP+FN); 2*p*r/(p+r) }]
lab <- dcast(d, coverage + cn ~ method, value.var="f1_exact")
lab[, facet_lab := sprintf("%s (△F1 = %+0.03f)", cn, `DNA+RNA` - `DNA-only`)]
## NOTE (0.1X only): published CN4 +0.024 / CN5 +0.009 vs recomputed +0.034 / +0.013. The baseline
## MintCNA_h5{,_single}.RData were re-run 2026-03-19 12:43-12:54, a day AFTER panel C was plotted
## 03-18 16:14, so panel C's source data was overwritten and its published numbers are not
## reproducible. The plotted POINTS differ too, not just the label -- so the label here is kept
## consistent with the points actually shown rather than hard-coded to the paper's value.
d <- merge(d, lab[, .(coverage, cn, facet_lab)], by=c("coverage","cn"))
setorder(d, coverage, cn)

## iso-F1 contour field + labels — placed exactly as in status_F1_1.R (anchor at max precision,
## hjust=1/vjust=0), to match the published S5 panels.
grid <- CJ(precision=seq(0,1,length.out=400), recall=seq(0,1,length.out=400))
grid[, f1 := fifelse(precision+recall==0, 0, 2*precision*recall/(precision+recall))]
f1_levels <- c(0.2,0.4,0.6,0.8)
f1_text <- rbindlist(lapply(f1_levels, function(v){ tmp <- grid[abs(f1-v) < 0.005]
  if(!nrow(tmp)) return(NULL); tmp[which.max(precision), .(precision, recall, label=paste0("F1=",v))] }))
pal <- c("MintCNA"="#1B9E77", "MintCNA(DNA)"="#D95F02", "MintCNA(RNA)"="#7570B3")
# legend display labels: keep the (DNA)/(RNA) keys internally but render DNA/RNA as superscripts.
method_labels <- parse(text = c("MintCNA", "MintCNA^DNA", "MintCNA^RNA"))
LIM <- c(0.1, 1)

mk_panel <- function(cov, title_txt, show_x=TRUE, show_legend=FALSE){
  dd <- d[coverage==cov]
  ord <- unique(dd[order(cn)]$facet_lab); dd[, facet_lab := factor(facet_lab, levels=ord)]
  p <- ggplot(dd, aes(precision, recall)) +
    geom_contour(data=grid, aes(x=precision, y=recall, z=f1), breaks=f1_levels,
                 color="grey85", linewidth=0.3, inherit.aes=FALSE) +
    geom_text(data=f1_text, aes(x=precision, y=recall, label=label),
              color="grey60", size=3, vjust=0, hjust=1, inherit.aes=FALSE) +
    geom_point(aes(color=method_lab), shape=17, size=5) +
    scale_color_manual(values=pal, breaks=names(pal), labels=method_labels) +
    coord_fixed(xlim=LIM, ylim=LIM) +
    ggh4x::facet_wrap2(~facet_lab, nrow=1, axes="all") +
    labs(x=if(show_x) "Precision" else NULL, y="Recall", color=NULL, title=title_txt) +
    theme_classic(base_size=11) +
    theme(legend.position=if(show_legend) "bottom" else "none",
          plot.title=element_text(face="bold", size=12))
  if(!show_x) p <- p + theme(axis.title.x=element_blank())
  p
}

## ---- individual panels: NO per-panel legend, NO hard-coded panel letter ----
## (the legend ships once as S5_legend.pdf, mirroring Result/legend_method_F1.pdf; add your own
##  (B)/(C)/(D) letters at assembly so they don't duplicate)
specs <- list(c("0.05X","S5_B_cov0.05.pdf"), c("0.1X","S5_C_cov0.1.pdf"), c("0.2X","S5_D_cov0.2.pdf"))
for(s in specs){
  p <- mk_panel(s[1], sprintf("coverage %s", s[1]), show_x=TRUE, show_legend=FALSE)
  cairo_pdf(file.path(FIG, s[2]), width=9, height=5.0); print(p); dev.off()
  cat("wrote", s[2], "\n")
}

## ---- standalone shared legend (place once) ----
gl <- ggplotGrob(mk_panel("0.1X", "", show_legend=TRUE))
gi <- grep("guide-box", vapply(gl$grobs, function(x) x$name, character(1)))
gi <- gi[vapply(gi, function(i) !inherits(gl$grobs[[i]], "zeroGrob"), logical(1))]
stopifnot(length(gi) >= 1)
cairo_pdf(file.path(FIG,"S5_legend.pdf"), width=5.2, height=0.45)
grid::grid.newpage(); grid::grid.draw(gl$grobs[[gi[1]]]); dev.off()
cat("wrote S5_legend.pdf\n")

## ---- combined B/C/D with ONE shared legend ----
pB <- mk_panel("0.05X","coverage 0.05X", show_x=FALSE, show_legend=TRUE)
pC <- mk_panel("0.1X", "coverage 0.1X",  show_x=FALSE, show_legend=TRUE)
pD <- mk_panel("0.2X", "coverage 0.2X",  show_x=TRUE,  show_legend=TRUE)
comb <- (pB / pC / pD) + plot_layout(guides="collect") & theme(legend.position="bottom")
cairo_pdf(file.path(FIG,"S5_BCD_combined.pdf"), width=9, height=13.0); print(comb); dev.off()
cat("wrote S5_BCD_combined.pdf\n")

cat("\npoint ranges (sanity): precision", sprintf("%.3f-%.3f", min(d$precision), max(d$precision)),
    " recall", sprintf("%.3f-%.3f", min(d$recall), max(d$recall)), " | axis window", LIM[1], "-", LIM[2], "\n")
cat("all points inside window:", all(d$precision>=LIM[1] & d$precision<=LIM[2] & d$recall>=LIM[1] & d$recall<=LIM[2]), "\n")
