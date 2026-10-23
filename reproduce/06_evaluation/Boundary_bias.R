library(data.table)
library(IRanges)
library(stringr)

# pred events for one method matrix pred_mat
pred_seg_one <- function(pred_mat) {
  summarize_cna_all_cells(
    CNmat_sub = pred_mat,
    baseline_cn = 2L,
    keep_neutral = FALSE,
    smooth = TRUE
  )
}

match_tp_and_bb_one_cell_state <- function(T, P, overlap_frac = 0.20) {
  # T, P: data.table for one cell + one event type (DEL or DUP), may have multiple chrs
  if (nrow(T) == 0L || nrow(P) == 0L) return(data.table())
  
  T <- copy(T); P <- copy(P)
  T[, tid := .I]; P[, pid := .I]
  
  cand <- rbindlist(lapply(intersect(unique(T$chr), unique(P$chr)), function(cc) {
    Tcc <- T[chr == cc]
    Pcc <- P[chr == cc]
    if (nrow(Tcc) == 0L || nrow(Pcc) == 0L) return(NULL)
    
    hits <- findOverlaps(IRanges(Tcc$start, Tcc$end),
                         IRanges(Pcc$start, Pcc$end))
    if (length(hits) == 0L) return(NULL)
    
    ti <- queryHits(hits)
    pj <- subjectHits(hits)
    
    ov_start <- pmax(Tcc$start[ti], Pcc$start[pj])
    ov_end   <- pmin(Tcc$end[ti],   Pcc$end[pj])
    
    # same convention as your F1 code (no +1)
    ov_len <- pmax(0, ov_end - ov_start)
    
    truth_len <- pmax(1, Tcc$len_bp[ti])  # avoid 0
    prop <- ov_len / truth_len
    
    data.table(
      tid = Tcc$tid[ti],
      pid = Pcc$pid[pj],
      chr = cc,
      t_start = Tcc$start[ti], t_end = Tcc$end[ti],
      p_start = Pcc$start[pj], p_end = Pcc$end[pj],
      ov_len = ov_len,
      prop = prop
    )[prop >= overlap_frac]
  }), use.names = TRUE, fill = TRUE)
  
  if (is.null(cand) || nrow(cand) == 0L) return(data.table())
  
  cand[, BB := abs(t_start - p_start) + abs(t_end - p_end)]
  
  # greedy one-to-one: maximize prop, then overlap, then smaller BB
  setorder(cand, -prop, -ov_len, BB)
  
  used_t <- rep(FALSE, nrow(T))
  used_p <- rep(FALSE, nrow(P))
  
  out <- vector("list", nrow(cand))
  k <- 0L
  for (i in seq_len(nrow(cand))) {
    tt <- cand$tid[i]
    pp <- cand$pid[i]
    if (!used_t[tt] && !used_p[pp]) {
      used_t[tt] <- TRUE
      used_p[pp] <- TRUE
      k <- k + 1L
      out[[k]] <- cand[i]
    }
  }
  if (k == 0L) return(data.table())
  rbindlist(out[seq_len(k)])
}

bb_per_method_from_segments <- function(truth_seg, pred_seg, overlap_frac = 0.20,
                                        summary = c("mean","median")) {
  summary <- match.arg(summary)
  
  cells <- intersect(unique(truth_seg$cell), unique(pred_seg$cell))
  if (length(cells) == 0L) stop("No overlapping cells between truth and pred segments.")
  
  res <- rbindlist(lapply(cells, function(cell_id) {
    T0 <- truth_seg[cell == cell_id]
    P0 <- pred_seg [cell == cell_id]
    
    # split by event type (DEL/DUP) to match like your F1
    out_pairs <- rbindlist(lapply(intersect(unique(T0$event), unique(P0$event)), function(ev) {
      match_tp_and_bb_one_cell_state(
        T = T0[event == ev],
        P = P0[event == ev],
        overlap_frac = overlap_frac
      )[, event := ev]
    }), use.names = TRUE, fill = TRUE)
    
    if (is.null(out_pairs) || nrow(out_pairs) == 0L) {
      return(data.table(cell = cell_id, BB_cell = NA_real_, n_TP = 0L))
    }
    
    bb_cell <- if (summary == "mean") mean(out_pairs$BB) else median(out_pairs$BB)
    data.table(cell = cell_id, BB_cell = bb_cell, n_TP = nrow(out_pairs))
  }))
  
  res[, n_cells_used := length(cells)]
  res[]
}

# build truth segments once
truth_seg <- summarize_cna_all_cells(truth_mat, baseline_cn = 2L, keep_neutral = FALSE, smooth = TRUE)

bb_tab <- rbindlist(lapply(names(methods), function(m) {
  print(m)
  pred_seg <- summarize_cna_all_cells(methods[[m]], baseline_cn = 2L, keep_neutral = FALSE, smooth = TRUE)
  
  dt <- bb_per_method_from_segments(truth_seg, pred_seg,
                                    overlap_frac = 0.20,
                                    summary = "mean")
  dt[, method := m]
  dt
}))

bb_tab[, method_disp := name_map[method]]
bb_tab[, method_disp := factor(method_disp, levels = method_order)]
bb_tab[, BB_cell_mb := BB_cell / 1e6]

bb_dt <- as.data.table(bb_tab)
bb_summary <- bb_dt[, .(
  n_cells_used = unique(n_cells_used)[1],
  TP_total     = sum(n_TP, na.rm = TRUE),
  BB_mean_bp   = mean(BB_cell, na.rm = TRUE),
  BB_mean_mb   = mean(BB_cell_mb, na.rm = TRUE)
), by = .(method = method_disp)][order(method)]


# ---- 6) boxplot (BB per cell) ----
p_bb <- ggplot(bb_tab[!is.na(bb_tab$BB_cell_mb), ],
               aes(x = method_disp, y = BB_cell_mb, fill = method_disp)) +
  geom_boxplot(width = 0.6, outlier.size = 0.35, linewidth = 0.3) +
  scale_fill_manual(values = color_map, guide = "none") +  
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    
    # x-axis
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, size = 11),
    
    # y-axis: make font larger
    axis.title.y = element_text(size = 14),
    axis.text.y  = element_text(size = 12),
    
    # if you also feel y ticks are thin/small
    axis.ticks.y = element_line(linewidth = 0.4)
  ) + theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank())
p_bb
ggsave(file.path(CNA_DIR,"Result/Boundary Bias.pdf"), plot = p_bb, width = 6.5, height = 6, units = "in", dpi = 600)

library(ggpubr)

leg_df <- data.frame(method_disp = factor(method_order, levels = method_order),
                     y = 1)

p_leg <- ggplot(leg_df, aes(x = method_disp, y = y, fill = method_disp)) +
  geom_col() +
  scale_fill_manual(name = "Method:", values = color_map) +
  theme_void(base_size = 11) +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 12),
    legend.key.size = unit(0.45, "cm")
  ) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

leg <- get_legend(p_leg)

ggsave(file.path(CNA_DIR,"Result/legend_method.pdf"), as_ggplot(leg), width = 10, height = 0.8)
