load("<MINTCNA_ROOT>/Result/CRC/Benchmark/SCOPE.RData")
#load("<MINTCNA_ROOT>/Result/CRC/Benchmark/MultiCNA.RData")
load("<MINTCNA_ROOT>/Result/CRC/Benchmark/MultiCNA_AE_197_d250.RData")
load("<MINTCNA_ROOT>/Result/CRC/Benchmark/HMMcopy.RData")
load("<MINTCNA_ROOT>/Result/CRC/Benchmark/DNAcopy.RData")
load("<MINTCNA_ROOT>/Result/CRC/Benchmark/SCCNV.RData")
load("<MINTCNA_ROOT>/Result/CRC/Benchmark/SCONCE.RData")
load("<MINTCNA_ROOT>/Result/CRC/Benchmark/rcCAE.RData")
load("<MINTCNA_ROOT>/Result/CRC/Benchmark/SCICoNE.RData")

colnames(ref_scope)[1] <- "chr";
cna_scope <- cna_scope[,colnames(cna_scope) %in% colnames(cna_multicna_ae)]
cna_scope <- cna_scope + 1L; cna_scope[cna_scope >= 5] <- 5L 
cna_hmmcopy <- cna_hmmcopy[,colnames(cna_hmmcopy) %in% colnames(cna_multicna_ae)]
cna_sconce <- cna_sconce[,colnames(cna_sconce) %in% colnames(cna_multicna_ae)]
cna_sccnv <- cna_sccnv[,colnames(cna_sccnv) %in% colnames(cna_multicna_ae)]
cna_rccae <- cna_rccae[,colnames(cna_rccae) %in% colnames(cna_multicna_ae)]
cna_scicone <- cna_scicone[,colnames(cna_scicone) %in% colnames(cna_multicna_ae)]

library(data.table)
library(viridis)

## --- helpers from before (keep as-is) ---
coerce_123 <- function(mat) {
  v <- as.integer(mat)
  # If already 1/2/3 only, keep as-is
  if (all(v %in% c(1L,2L,3L, NA_integer_))) {
    return(matrix(v, nrow = nrow(mat), ncol = ncol(mat),
                  dimnames = dimnames(mat)))
  }
  # Otherwise map 1..5 -> 1/2/3
  out <- rep(NA_integer_, length(v))
  out[v %in% c(1L,2L)] <- 1L  # deletions
  out[v %in% 3L]        <- 2L  # neutral
  out[v %in% c(4L,5L)] <- 3L  # duplications
  matrix(out, nrow = nrow(mat), ncol = ncol(mat),
         dimnames = dimnames(mat))
}

cna_hmmcopy  = coerce_123(cna_hmmcopy)
#cna_multicna = coerce_123(cna_multicna)
cna_multicna_ae = coerce_123(cna_multicna_ae)
cna_scope    = coerce_123(cna_scope)
cna_sconce    = coerce_123(cna_sconce)
cna_sccnv    = coerce_123(cna_sccnv)
cna_rccae    = coerce_123(cna_rccae)
cna_scicone    = coerce_123(cna_scicone)

library(data.table)
library(ggplot2)

plot_scna_freq <- function(cna_mat, ref_df, method_name = "Method") {
  # cna_mat : bins × cells, values {1 = del, 2 = normal, 3 = amp}
  # ref_df  : data.frame with column 'chr' (one row per bin)
  
  dt_ref <- as.data.table(ref_df)
  dt_ref[, bin_id := .I]
  
  states  <- as.matrix(cna_mat)
  n_cells <- ncol(states)
  
  # per-bin frequencies
  freq_dt <- data.table(
    bin_id   = dt_ref$bin_id,
    chr_raw  = as.character(dt_ref$chr),
    freq_amp = rowSums(states == 3, na.rm = TRUE) / n_cells,
    freq_del = rowSums(states == 1, na.rm = TRUE) / n_cells
  )
  
  # keep chr1–22 and order them
  freq_dt[, chr := sub("^chr", "", chr_raw)]
  freq_dt <- freq_dt[chr %in% as.character(1:22)]
  freq_dt[, chr := factor(chr, levels = as.character(1:22))]
  setorder(freq_dt, chr, bin_id)
  
  # simple genomic index
  freq_dt[, gpos := .I]
  chr_info <- freq_dt[, .(start = min(gpos), end = max(gpos)), by = chr]
  chr_info[, xmid := (start + end) / 2]
  
  # long + signed for plotting
  freq_long <- melt(
    freq_dt,
    id.vars      = c("bin_id", "chr", "gpos"),
    measure.vars = c("freq_amp", "freq_del"),
    variable.name = "type",
    value.name   = "freq"
  )
  
  freq_long[, CNA_type :=
              ifelse(type == "freq_amp", "Amplification", "Deletion")]
  freq_long[, freq_signed :=
              ifelse(CNA_type == "Deletion", -freq, freq)]
  freq_long[, CNA_type := factor(CNA_type,
                                 levels = c("Amplification","Deletion"))]
  
  # plot
  p <- ggplot(freq_long, aes(x = gpos, y = freq_signed, fill = CNA_type)) +
    geom_col(width = 1) +
    geom_hline(yintercept = 0, colour = "black", size = 0.3) +
    geom_vline(xintercept = chr_info$start, colour = "grey80", size = 0.2) +
    scale_x_continuous(
      breaks = chr_info$xmid,
      labels = chr_info$chr,
      expand = c(0, 0)
    ) +
    scale_fill_manual(
      name   = "CNA type",
      values = c("Amplification" = "#E41A1C",
                 "Deletion"      = "#377EB8")
    ) +
    coord_cartesian(ylim = c(-1, 1)) +
    labs(
      title = method_name,
      x     = "Chromosome",
      y     = "SCNA frequency"
    ) +
    theme_bw(base_size = 10) +
    theme(
      panel.grid   = element_blank(),
      axis.text.x  = element_text(angle = 0),
      legend.position = "none"
    )
  
  return(p)
}
p_wgs    <- plot_scna_freq(cna_dnacopy,    ref_dnacopy,    "WGS(Ground Truth)")
#p_multicna    <- plot_scna_freq(cna_multicna,    ref_multicna,    "MultiCNA")
p_multicna_ae <- plot_scna_freq(cna_multicna_ae, ref_multicna_ae, "MintCNA")
p_rccae <- plot_scna_freq(cna_rccae, ref_multicna_ae, "rcCAE")
p_hmmcopy     <- plot_scna_freq(cna_hmmcopy,     ref_hmmcopy,     "HMMcopy")
p_scope       <- plot_scna_freq(cna_scope,       ref_scope,       "SCOPE")
p_sconce       <- plot_scna_freq(cna_sconce,       ref_sconce,       "SCONCE")
p_sccnv       <- plot_scna_freq(cna_sccnv,       ref_sccnv,       "SCCNV")
p_scicone       <- plot_scna_freq(cna_scicone,       ref_multicna_ae,       "SCICoNE")

library(ggplot2)

outdir <- "<MINTCNA_ROOT>/Result/CRC/Segmentation/Evaluation/frequency"

save_freq_plot <- function(p, fname, h = 2.2) {
  ggsave(
    filename = file.path(outdir, fname),
    plot     = p,
    width    = 8.5,
    height   = h,
    dpi = 600
  )
}

save_freq_plot(p_wgs,        "CRC01_freq_WGS.pdf")
#save_freq_plot(p_multicna,   "CRC01_freq_MultiCNA.png")
save_freq_plot(p_multicna_ae,"CRC01_freq_MintCNA.pdf")
save_freq_plot(p_hmmcopy,    "CRC01_freq_HMMcopy.pdf")
save_freq_plot(p_scope,      "CRC01_freq_SCOPE.pdf")
save_freq_plot(p_sconce,      "CRC01_freq_SCONCE.pdf")
save_freq_plot(p_sccnv,      "CRC01_freq_SCCNV.pdf")
save_freq_plot(p_rccae,      "CRC01_freq_rcCAE.pdf")
save_freq_plot(p_scicone,      "CRC01_freq_SCICoNE.pdf")


# library(patchwork)
# 
# # combine: 5 rows × 1 column, shared legend
# p_all <- (
#   p_wgs /
#     p_multicna_ae /
#     p_hmmcopy /
#     p_scope /
#     p_sconce /
#     p_sccnv
# ) + plot_layout(guides = "collect")
# 
# # draw in RStudio
# p_all
# 
# # save combined figure
# ggsave("<MINTCNA_ROOT>/Result/CRC/Segmentation/Evaluation/frequency/CRC01_SCNA_freq.pdf",
#        plot   = p_all,
#        width  = 8.5,   # adjust as you like
#        height = 10,
#        dpi = 600)  # taller for 5 panels

