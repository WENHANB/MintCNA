library(rtracklayer)
library(data.table)
library(dplyr)
setwd("<MINTCNA_ROOT>/")

## region of interest
hg19_gene_info <- readGFF("Data/CRC/hg19/gencode.v19.annotation.gtf.gz")
hg19_gene_info <- hg19_gene_info[hg19_gene_info$type == "gene",]
cytoband <- fread("Data/CRC/hg19/cytoBand.txt", header = FALSE, stringsAsFactors = FALSE)
colnames(cytoband) <- c("chrom", "chromStart", "chromEnd", "name", "gieStain")
region_cnv <- fread("Data/CRC/region_cnv.csv",header = T)
region_cnv <- merge(region_cnv, cytoband[,1:4], by = c("chrom","name"))
region_cnv <- merge(region_cnv, hg19_gene_info[,c("start","end","gene_name")], by = c("gene_name"))
colnames(region_cnv)[colnames(region_cnv) == "chromStart"] <- "cytoband_start";colnames(region_cnv)[colnames(region_cnv) == "chromEnd"] <- "cytoband_end"
colnames(region_cnv)[colnames(region_cnv) == "start"] <- "gene_start";colnames(region_cnv)[colnames(region_cnv) == "end"] <- "gene_end"

# get_region_bins <- function(chr, s, e, ref) {
#   which(ref$chr == chr & ref$start < e & ref$end > s)
# }
# 
# region_cnv[, n_bins_cytoband := sapply(seq_len(nrow(region_cnv)), function(i) {
#   chr <- region_cnv$chrom[i]
#   s   <- region_cnv$cytoband_start[i]
#   e   <- region_cnv$cytoband_end[i]
#   length(get_region_bins(chr, s, e, ref))
# })]

## Comparing
get_region_bins <- function(chr, s, e, ref) {
  which(ref$chr == chr & ref$start < e & ref$end > s)
}

# Works for 1–5 coding or raw integer CN
to_signs <- function(mat){
  z <- matrix(0L, nrow(mat), ncol(mat))
  if (max(mat, na.rm=TRUE) == 5) { z[mat <=2] <- -1L; z[mat >=4] <-  1L }
  else                           { z[mat < 2] <- -1L; z[mat > 2] <-  1L }
  z
}

region_percent <- function(mat, ref, regions,
                           rule = c("fraction","any"),
                           alpha = 1,
                           method_name = NULL) {
  rule  <- match.arg(rule)
  signs <- to_signs(mat)                 # [bins x cells]
  n_cells_total <- ncol(mat)
  
  out <- data.table::rbindlist(lapply(seq_len(nrow(regions)), function(i) {
    b <- get_region_bins(regions$chrom[i], regions$gene_start[i], regions$gene_end[i], ref)
    #b <- get_region_bins(regions$chrom[i], regions$cytoband_start[i], regions$cytoband_end[i], ref)
    if (!length(b)) return(NULL)
    sub <- signs[b,, drop = FALSE]
    
    # per-cell summaries
    frac_gain <- colMeans(sub ==  1)
    frac_loss <- colMeans(sub == -1)
    any_gain  <- colSums(sub ==  1) > 0
    any_loss  <- colSums(sub == -1) > 0
    any_alt   <- colSums(sub !=  0) > 0
    
    # cell-level positive according to rule + panel
    status <- tolower(regions$cn_status[i])
    pos <- switch(status,
                  "gain" = if (rule=="fraction") frac_gain >= alpha else any_gain,
                  "dup"  = if (rule=="fraction") frac_gain >= alpha else any_gain,
                  "amp"  = if (rule=="fraction") frac_gain >= alpha else any_gain,
                  "loss" = if (rule=="fraction") frac_loss >= alpha else any_loss,
                  "del"  = if (rule=="fraction") frac_loss >= alpha else any_loss,
                  "nor"  = if (rule=="fraction") (pmax(frac_gain, frac_loss) >= alpha) else any_alt,
                  "normal"=if (rule=="fraction") (pmax(frac_gain, frac_loss) >= alpha) else any_alt,
                  { warning("Unknown cn_status: ", status); rep(FALSE, n_cells_total) }
    )
    
    cna_cells   <- sum(pos, na.rm = TRUE)
    no_cna_cells<- n_cells_total - cna_cells
    
    data.table::data.table(
      gene_name  = regions$gene_name[i],
      chrom      = regions$chrom[i],
      cn_status  = regions$cn_status[i],
      n_bins     = length(b),
      n_cells    = n_cells_total,
      CNA_cell   = cna_cells,
      No_CNA_cell= no_cna_cells,
      pct_cells  = 100 * cna_cells / n_cells_total
    )
  }), use.names = TRUE)
  
  if (!is.null(method_name) && nrow(out)) out[, method := method_name]
  out
}

# Run for both methods
load("<MINTCNA_ROOT>/Result/CRC/Benchmark/SCOPE_197.RData")
load("<MINTCNA_ROOT>/Result/CRC/Benchmark/MultiCNA_197.RData")
load("<MINTCNA_ROOT>/Result/CRC/Benchmark/MultiCNA_AE_197_d250.RData")
load("<MINTCNA_ROOT>/Result/CRC/Benchmark/HMMcopy_197.RData")
load("<MINTCNA_ROOT>/Result/CRC/Benchmark/DNAcopy.RData")

colnames(ref_scope)[1] <- "chr";
cna_scope <- cna_scope[,colnames(cna_scope) %in% colnames(cna_multicna)]
cna_scope <- cna_scope + 1L; cna_scope[cna_scope >= 5] <- 5L 
cna_hmmcopy <- cna_hmmcopy[,colnames(cna_hmmcopy) %in% colnames(cna_multicna)]

regions     <- as.data.table(region_cnv)

## Result and Merge together
res_S  <- region_percent(cna_multicna,       ref_multicna,       regions, rule="any")[,  method := "MultiCNA"]
res_AE  <- region_percent(cna_multicna_ae,       ref_multicna_ae,       regions, rule="any")[,  method := "MultiCNA-AE"]
res_HC <- region_percent(cna_hmmcopy, ref_hmmcopy,       regions, rule="any")[,  method := "HMMcopy"]
res_SC <- region_percent(cna_scope,     ref_scope, regions, rule="any")[,  method := "SCOPE"]
# res_S  <- region_percent(cna_multicna,       ref_multicna,       regions, rule="fraction")[,  method := "MultiCNA"]
# res_HC <- region_percent(cna_hmmcopy, ref_hmmcopy,       regions, rule="fraction")[,  method := "HMMcopy"]
# res_SC <- region_percent(cna_scope,     ref_scope, regions, rule="fraction")[,  method := "SCOPE"]
# res_FL <- region_percent(cna_flcna[,colnames(cna_multicna)],
#                          ref_flcna, regions, rule="fraction")[,  method := "FLCNA"]
res_all <- data.table::rbindlist(list(res_S, res_AE,res_HC, res_SC), use.names=TRUE, fill=TRUE)


library(ggplot2)
library(viridisLite)

## --- ensure MultiCNA is always on top in legend and bars ---
methods <- c("MultiCNA", setdiff(unique(res_all$method), "MultiCNA"))
res_all$method <- factor(res_all$method,
                         levels = c("MultiCNA", "MultiCNA-AE","HMMcopy", "SCOPE"))

## --- build viridis palette matching the order ---
pal <- setNames(viridis(length(methods), option = "D", begin = 0, end = 1), methods)

## --- plotting ---
res_all$cn_status <- factor(res_all$cn_status, levels = c("gain","loss","nor"),
                            labels = c("DUP","DEL","Nor"))

ggplot(res_all, aes(x = gene_name, y = pct_cells, fill = method)) +
  geom_col(position = position_dodge(width = 0.8, reverse = FALSE), width = 0.75) +
  facet_wrap(~ cn_status, ncol = 1, scales = "free_x") +
  scale_fill_manual(values = pal, guide = guide_legend(reverse = FALSE)) +
  xlab("Canonical region") + ylab("% of cells with CNA")  +
  theme_bw(base_size = 13) +
  theme(
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(face = "bold"),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 11),
    axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("<MINTCNA_ROOT>/Result/CRC/Segmentation/Evaluation/pseudo/CNA_region_panel.png",
       width=7, height=10, dpi=320)


########################## Recall,Specificity ##########################
## ===================== Add-on: Recall & Specificity =====================
library(data.table)

## --- Function to compute 2×2 confusion metrics per method ---
get_perf_table <- function(res_all, method_name) {
  ra <- copy(res_all)[method == method_name]
  ra[, status := tolower(trimws(cn_status))]
  ra[, truth := fifelse(status %in% c("gain","dup","amp","loss","del"), "cna",
                        fifelse(status %in% c("nor","normal"), "nor", NA_character_))]
  ra <- ra[!is.na(truth)]
  
  per_region <- ra[, .(
    TP = ifelse(truth=="cna", CNA_cell,       0L),
    FN = ifelse(truth=="cna", No_CNA_cell,    0L),
    FP = ifelse(truth=="nor", CNA_cell,       0L),
    TN = ifelse(truth=="nor", No_CNA_cell,    0L)
  ), by = .(gene_name, truth)]
  
  cm <- per_region[, .(
    TP = sum(TP), FN = sum(FN), FP = sum(FP), TN = sum(TN)
  )]
  
  cm[, `:=`(
    method           = method_name,
    Total            = TP + FN + FP + TN,
    Sensitivity      = TP / (TP + FN),
    Specificity      = TN / (TN + FP),
    Precision  = TP / (TP + FP),
    Accuracy         = (TP + TN) / (TP + TN + FP + FN),
    Balanced_Accuracy = (TP / (TP + FN) + TN / (TN + FP)) / 2
  )]
  cm[, `:=`(
    F1_score = 2 * (Precision * Sensitivity) / (Precision + Sensitivity)
  )]
  cm[, .(method, TP, FN, FP, TN, Total, Sensitivity, Specificity, Precision, F1_score,
         Accuracy, Balanced_Accuracy)]
}


## --- Run for all methods and combine ---
perf_plot <- rbindlist(list(
  get_perf_table(res_all, "MultiCNA"),
  get_perf_table(res_all, "MultiCNA-AE"),
  get_perf_table(res_all, "HMMcopy"),
  get_perf_table(res_all, "SCOPE")
), use.names = TRUE, fill = TRUE)


# Visualization: method-level Recall vs Specificity
library(ggplot2)
library(ggrepel)
perf_plot <- within(perf_plot, {
  Sensitivity <- 100 * Sensitivity
  Specificity <- 100 * Specificity
  Accuracy <-  100 * Accuracy
})

## ensure MultiCNA first and palette matches legend
# methods <- c("MultiCNA", setdiff(perf_plot$method, "MultiCNA"))
# perf_plot$method <- factor(perf_plot$method, levels = methods)
# pal <- setNames(viridis(length(methods), option = "D", begin = 0, end = 1), methods)
# 
# ## plot
# ggplot(perf_plot, aes(x = Specificity, y = Sensitivity, color = method)) +
#   geom_point(size = 4) +
#   geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.4, color = "grey50") +
#   scale_color_manual(values = pal, name = "Method") +
#   scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
#   scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
#   xlab("Specificity on Normal Regions (%)") +
#   ylab("Sensitivity on DUP/DEL Regions (%)") +
#   ggtitle("Comparison of Sensitivity and Specificity Across CNA Calling Methods") +
#   theme_bw(base_size = 14) +
#   theme(
#     plot.title = element_text(face = "bold", hjust = 0.5),
#     legend.position = "right",
#     legend.title = element_text(face = "bold"),
#     axis.title = element_text(face = "bold")
#   )

# ggsave("<MINTCNA_ROOT>/Result/CRC/Segmentation/Evaluation/pseudo/CNA_method_sens_spec.png",
#        width=10, height=6, dpi=320)

## plot -- sensitivity, specificity, and accuracy
# perf_plot has: method, TP, FN, FP, TN, Total, Sensitivity, Specificity, Accuracy, Balanced_Accuracy

# 1) order + colors (MultiCNA first)
methods <- c("MultiCNA","MultiCNA-AE","HMMcopy","SCOPE")
perf_plot[, method := factor(as.character(method), levels = methods)]
pal <- setNames(viridis(length(methods), option = "D"), methods)

# 2) long format for 3 metrics
perf_long <- melt(
  perf_plot[, .(method, Sensitivity, Specificity, Accuracy)],
  id.vars = "method",
  variable.name = "metric",
  value.name = "value"
)

# 3) plot: three panels (bars)
ggplot(perf_long, aes(x = method, y = value, fill = method)) +
  geom_col(width = 0.7) +
  facet_wrap(~ metric, nrow = 1) +
  scale_fill_manual(values = pal, name = "Method") +
  scale_y_continuous(breaks = seq(70, 100, 5)) +   # just set breaks
  coord_cartesian(ylim = c(70, 100)) +             # zoom, no data dropped
  xlab(NULL) + ylab("Percentage (%)") +
  ggtitle("CNA Caller Performance by Metric") +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    strip.background = element_rect(fill = "grey92"),
    strip.text = element_text(face = "bold"),
    panel.grid.major.y = element_line(color = "grey85", linewidth = 0.3),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  )
ggsave("<MINTCNA_ROOT>/Result/CRC/Segmentation/Evaluation/pseudo/CNA_Metrics.png",
       width=10, height=6, dpi=320)

# F1-score vs recall vs precision
# perf_plot has columns: method, TP, FN, FP, TN, Total, Sensitivity, Specificity, Precision, F1_score
perf_plot <- within(perf_plot, {
  Sensitivity <- Sensitivity/100
  Specificity <- Specificity/100
  Accuracy <-  Accuracy/100
})
perf_plot$Recall <- perf_plot$Sensitivity

# methods order and palette (MultiCNA first)
methods <- c("MultiCNA", setdiff(perf_plot$method, "MultiCNA"))
perf_plot$method <- factor(perf_plot$method, levels = methods)
pal <- setNames(viridis(length(methods), option = "D"), methods)

# --- iso-F1 background (dashed contours) ---
f_levels <- c(0.9, 0.8, 0.7)
grid <- expand.grid(
  Precision = seq(0.5, 1, by = 0.002),
  Recall    = seq(0.5, 1, by = 0.002)
)
grid$F1 <- with(grid, 2 * Precision * Recall / (Precision + Recall))

# label positions for the iso-F1 curves (near the right edge)
P_lab <- 0.98
lab_df <- data.frame(
  F1 = f_levels,
  Precision = P_lab,
  Recall = (f_levels * P_lab) / (2 * P_lab - f_levels),
  label = paste0("f = ", f_levels)
)

ggplot(perf_plot, aes(x = Precision, y = Recall, color = method)) +
  stat_contour(data = grid, aes(z = F1), breaks = f_levels,
               color = "grey50", linetype = "dashed", linewidth = 0.5) +
  geom_text(data = lab_df, aes(x = Precision, y = Recall,label = label), inherit.aes = FALSE,
            hjust = -0.1, vjust = 0.5, size = 3.8, color = "grey35") +
  geom_point(size = 5) +
  scale_color_manual(values = pal, name = "Method") +
  coord_cartesian(xlim = c(0.6, 1.0), ylim = c(0.6, 1.0)) +
  xlab("Precision = TP / (TP + FP)") +
  ylab("Recall = TP / (TP + FN)") +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    plot.margin = margin(10, 40, 10, 10) # room for f-labels on right
  )

ggsave("<MINTCNA_ROOT>/Result/CRC/Segmentation/Evaluation/pseudo/F1_score.png",
       width=8, height=6, dpi=600)

