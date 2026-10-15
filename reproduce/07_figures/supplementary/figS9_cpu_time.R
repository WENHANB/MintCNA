#!/usr/bin/env Rscript
# Supplementary Fig. 9 (S8 until 2026-10-09) -- computational cost of every method in the paper.
#
# WHY TWO PANELS
# The manuscript benchmarks MintCNA twice: against seven scDNA callers on
# simulations and the colorectal scTrio-seq2 dataset (Figs 2-4), and against six
# scATAC callers on the 20-tumour glioma multiome cohort (Fig 5). One runtime
# figure has to cover both, and the two arms cannot share an axis of "per
# dataset" cost because the multiome arm was run per tumour. So:
#   a  scDNA arm    one number per method per dataset, from the timing sheet
#   b  multiome arm one number per method per tumour (20 points), median marked
#
# WHAT IS PLOTTED IN b, AND WHAT IS NOT
# CPU time (TotalCPU from SLURM accounting, core-minutes), not wall time. The
# methods were sized from their own measured peaks and ran on 2-32 cores, so wall
# time would reward the method given the most cores; core-minutes is the quantity
# that does not. Only the REPORTED configuration of each method is drawn -- the
# parameter sweeps that also sit in the accounting log are excluded by job name
# and array task. MintCNA's number is the whole chain from fragments to calls
# (binning, preprocessing, AG-DCAE, alignment, scan, call) as one job per tumour;
# 17 tumours ran that way, the three pilot tumours ran in stages and are left out.
# The ATAC-only ablation shares the denoising stage with the full model and is not
# a separate bar; its own scan + call cost is printed for the legend.
#
# THE TWO CAVEATS THAT MUST REACH THE LEGEND
#   - RIDDLER's window matrix and Copy-scAT's binning were re-implemented in
#     data.table (verified identical) because the published input steps could not
#     read a 10x Multiome fragments file; their times are therefore FASTER than
#     the published code would give.
#   - core-minutes still favour a method that parallelises badly, because idle
#     cores cost nothing in TotalCPU.
#
# STYLE
# The author's original S7: dodged bars on a log axis, one colour per dataset,
# methods on x with rotated labels, legend below. Panel B keeps exactly that
# grammar with a third dataset colour; its bars are the per-tumour median and the
# whisker is the interquartile range over the 20 tumours.

.libPaths(c("<R_LIBS_USER>", .libPaths()))
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(cowplot)
})
ROOT <- "<MINTCNA_ROOT>"
source(file.path(ROOT, "figure4paper/shared/fig_palette.R"))
DATA <- file.path(ROOT, "figure4paper/figureS-supplementary/_work/data")
OUT  <- file.path(ROOT, "figure4paper/figureS-supplementary/_work/output")
PDF_OUT <- file.path(ROOT, "figure4paper/figureS-supplementary")   # final PDFs live at the top level; everything else under _work/
QA   <- file.path(ROOT, "figure4paper/figureS-supplementary/_work/qa")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
dir.create(QA,  showWarnings = FALSE, recursive = TRUE)

## ---- data ----------------------------------------------------------------------
dna <- fread(file.path(DATA, "figS9_dna_cpu_time.csv"))
ORD_A <- c("MintCNA", "MintCNA(DNA)", "SCOPE", "HMMcopy", "SCONCE", "SCCNV", "rcCAE", "SCICoNE")
DS_A  <- c("Simulation (N = 100 cells)", "Colorectal cancer (scTrio-seq2)")
stopifnot(all(dna$method %in% ORD_A), all(dna$dataset %in% DS_A))
PM_A <- c("MintCNA" = "MintCNA", "MintCNA(DNA)" = "MintCNA^DNA", "SCOPE" = "SCOPE", "HMMcopy" = "HMMcopy",
          "SCONCE" = "SCONCE", "SCCNV" = "SCCNV", "rcCAE" = "rcCAE", "SCICoNE" = "SCICoNE")
lab_by_name <- function(pm) function(x) parse(text = pm[as.character(x)])

sa <- fread(file.path(DATA, "figS9_atac_runtime_sacct.tsv"))
sa <- sa[state == "COMPLETED" & !is.na(cpu_min)]
sa[, array_id := sub("_.*", "", jobid)]
# the REPORTED run of each method, by job name and, where an array mixed
# configurations, by task range. Sample = array task (all arrays share one
# sample order); a sample re-run keeps its latest completed task.
pick <- rbindlist(list(
  sa[jobname == "mint_all",                            .(method = "MintCNA",       sample = paste0(array_id, "_", array_task), cpu_min, req_cpus, peak_rss_gb)],
  sa[jobname == "atacnv",                              .(method = "AtaCNV",        sample = array_task, cpu_min, req_cpus, peak_rss_gb, array_id)][order(-as.numeric(array_id))][!duplicated(sample)][, array_id := NULL],
  sa[jobname == "copyscat",                            .(method = "Copy-scAT",     sample = array_task, cpu_min, req_cpus, peak_rss_gb)],
  sa[jobname == "riddler" & array_id == "39637608",    .(method = "RIDDLER",       sample = array_task, cpu_min, req_cpus, peak_rss_gb)],
  # AS PUBLISHED (100-kb window, minFrags 20000): array 39335320 tasks 1-20. That is the
  # configuration Fig 5A/B scores (pearson_vs_truth.R, EPI_TAG default). The adjusted
  # 1-Mb run (39388639 tasks 61-80, ~43 min) is the sensitivity analysis, not the headline.
  sa[jobname == "epiAneu" & array_id == "39335320" & as.integer(array_task) <= 20,
                                                       .(method = "epiAneufinder", sample = array_task, cpu_min, req_cpus, peak_rss_gb)],
  sa[jobname == "icnv_default",                        .(method = "inferCNV",      sample = array_task, cpu_min, req_cpus, peak_rss_gb)],
  sa[jobname == "congas" & array_id == "39388738",     .(method = "CONGAS+",       sample = array_task, cpu_min, req_cpus, peak_rss_gb)]))
cnt <- pick[, .(n = .N, median_min = round(median(cpu_min), 1), q1 = quantile(cpu_min, .25), q3 = quantile(cpu_min, .75),
                min = round(min(cpu_min), 1), max = round(max(cpu_min), 1), cores = as.integer(median(req_cpus)),
                peak_gb = round(max(peak_rss_gb, na.rm = TRUE), 1)), by = method]
cat("panel b, per-tumour CPU minutes:\n"); print(cnt)
stopifnot(cnt[method != "MintCNA", all(n == 20)], cnt[method == "MintCNA", n == 17])
fwrite(cnt[, .(method, n, median_min, min, max, cores, peak_gb)], file.path(DATA, "figS9_atac_runtime_summary.csv"))
fwrite(pick, file.path(DATA, "figS9_atac_runtime_per_tumour.csv"))
abl <- sa[jobname == "atac_only", .(sample = array_task, scan = cpu_min)]
abl <- merge(abl, sa[jobname == "atac_call" & array_id == "39938586", .(sample = array_task, call = cpu_min)], by = "sample")
cat(sprintf("MintCNA(ATAC) scan+call stages: median %.1f min over %d tumours (scan %.1f, call %.1f)\n",
            median(abl$scan + abl$call), nrow(abl), median(abl$scan), median(abl$call)))

ORD_B <- c("MintCNA", "AtaCNV", "Copy-scAT", "RIDDLER", "epiAneufinder", "inferCNV", "CONGAS+")
PM_B  <- c("MintCNA" = "MintCNA", "AtaCNV" = "'AtaCNA'", "Copy-scAT" = "'Copy-scAT'", "RIDDLER" = "RIDDLER",
           "epiAneufinder" = "epiAneufinder", "inferCNV" = "inferCNV", "CONGAS+" = "'CONGAS+'")
DS_B <- "Glioma multiome (20 tumours, per tumour)"

## ---- one grammar for both panels ---------------------------------------------------
DS_COL <- c("#3C5488", "#E64B35", "#00A087")   # the author's navy / red, plus one teal
names(DS_COL) <- c(DS_A, DS_B)
DS_LAB <- c("Simulation (N = 100 cells)", "Colorectal cancer (CRC01)", DS_B)   # cell counts differ by method (197/198/205), stated in the legend (2026-10-14)
Y_BREAKS <- c(1, 2, 5, 10, 20, 50, 100, 200, 500, 1000)
bar_theme <- theme_bw(BASE_PT, base_family = FIG_FAMILY) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = PT_SUP, colour = "black"),
        axis.text.y = element_text(size = BASE_PT - 1, colour = "black"),
        axis.title.y = element_text(size = BASE_PT),
        panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(linewidth = 0.25, colour = "grey88"),
        panel.border = element_rect(linewidth = 0.4, colour = "grey30"),
        legend.position = "bottom", legend.title = element_text(size = BASE_PT),
        legend.text = element_text(size = BASE_PT - 1), legend.key.size = unit(3.2, "mm"),
        # room above the panel for the letter; without it "A" sat on the 1000 tick
        plot.margin = margin(8, 3, 1, 2))

dna[, method := factor(method, levels = ORD_A)]
dna[, dataset := factor(dataset, levels = names(DS_COL))]
p_a <- ggplot(dna, aes(method, cpu_min, fill = dataset)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_y_log10(breaks = Y_BREAKS, labels = Y_BREAKS, limits = c(1, 1000), expand = expansion(mult = c(0, 0.02))) +
  scale_x_discrete(labels = lab_by_name(PM_A)) +
  scale_fill_manual(values = DS_COL, breaks = names(DS_COL), labels = DS_LAB, name = "Dataset", drop = FALSE) +
  labs(x = NULL, y = "CPU time (minutes, log scale)") +
  bar_theme

cnt[, method := factor(method, levels = ORD_B)]
cnt[, dataset := factor(DS_B, levels = names(DS_COL))]
p_b <- ggplot(cnt, aes(method, median_min, fill = dataset)) +
  geom_col(width = 0.7) +
  # interquartile range over the 20 tumours
  geom_errorbar(aes(ymin = q1, ymax = q3), width = 0.25, linewidth = 0.35, colour = "grey25") +
  scale_y_log10(breaks = Y_BREAKS, labels = Y_BREAKS, limits = c(1, 1000), expand = expansion(mult = c(0, 0.02))) +
  scale_x_discrete(labels = lab_by_name(PM_B)) +
  scale_fill_manual(values = DS_COL, breaks = names(DS_COL), labels = DS_LAB, name = "Dataset", drop = FALSE) +
  labs(x = NULL, y = NULL) +
  bar_theme

## ---- assemble: one shared legend below both panels --------------------------------
# the legend is drawn from a dummy that holds all three datasets: taken from p_a
# alone the teal key came out empty
dummy <- ggplot(data.table(x = names(DS_COL), y = 1, dataset = factor(names(DS_COL), levels = names(DS_COL))),
                aes(x, y, fill = dataset)) + geom_col() +
  scale_fill_manual(values = DS_COL, breaks = names(DS_COL), labels = DS_LAB, name = "Dataset") +
  bar_theme + guides(fill = guide_legend(nrow = 1))
leg <- get_plot_component(dummy, "guide-box-bottom", return_all = TRUE)
if (is.list(leg) && !inherits(leg, "gtable")) leg <- leg[[1]]
p_a <- p_a + theme(legend.position = "none"); p_b <- p_b + theme(legend.position = "none")
W_MM <- 165; H_MM <- 76;  mm <- 1/25.4     # 165 = SI page text width (was 180)
A_W <- 0.55; LEG_H <- 0.10
body <- plot_grid(p_a, p_b, ncol = 2, rel_widths = c(A_W, 1 - A_W), align = "h", axis = "tb",
                  labels = c("a", "b"), label_size = 8, label_fontfamily = FIG_FAMILY, label_fontface = "bold",
                  label_x = c(0.005, 0.005), label_y = 1.0, vjust = 1.1)
fig <- plot_grid(body, leg, ncol = 1, rel_heights = c(1 - LEG_H, LEG_H))

cairo_pdf(file.path(PDF_OUT, "SuppFig09_cpu_time.pdf"), width = W_MM*mm, height = H_MM*mm, family = FIG_FAMILY)
print(fig); dev.off()
svglite::svglite(file.path(OUT, "figS9_cpu_time.svg"), width = W_MM*mm, height = H_MM*mm)
print(fig); dev.off()
ragg::agg_png(file.path(QA, "figS9_cpu_time_preview.png"), width = W_MM*mm, height = H_MM*mm, units = "in", res = 900)
print(fig); dev.off()
cat(sprintf("\nwrote figS9_cpu_time.pdf/.svg  %d x %d mm\n", W_MM, H_MM))
