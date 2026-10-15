library(data.table)
library(IRanges)
library(stringr)
library(ggplot2)

# -----------------------------
# 0) methods + site-specific truth (same mapping as RMSE)
# -----------------------------
methods <- list(
  "MintCNA"      = cna_multicna_ae,
  "SCOPE"        = cna_scope,
  "HMMcopy"      = cna_hmmcopy,
  "SCONCE"       = cna_sconce,
  "SCCNV"        = cna_sccnv,
  "rcCAE"        = cna_rccae
  # "SCICoNE"    = cna_scicone
)

methods <- lapply(methods, as.matrix)
truth_base <- as.matrix(cna_dnacopy)   # cols: PT, LN1, ML1, MP1

cell_to_truthcol <- function(cell_names) {
  site <- sub("^.*_(PT|LN|ML|MP).*?$", "\\1", cell_names)
  site[!site %in% c("PT","LN","ML","MP")] <- NA_character_
  
  truth_col <- fifelse(site == "PT", "PT",
                       fifelse(site == "LN", "LN1",
                               fifelse(site == "ML", "ML1",
                                       fifelse(site == "MP", "MP1", NA_character_))))
  data.table(cell = cell_names, site = site, truth_col = truth_col)
}

summarize_cna_all_cells <- function(CNmat_sub, baseline_cn = 2L, keep_neutral = FALSE,
                                    smooth = TRUE,
                                    max_short_bins = 1L, min_neighbor_bins = 5L, max_step = 1L) {
  
  bins_dt <- .parse_bins(rownames(CNmat_sub))  # must contain idx, chr, start, end
  
  res <- rbindlist(lapply(seq_len(ncol(CNmat_sub)), function(j) {
    cell <- colnames(CNmat_sub)[j]
    seg  <- summarize_cna_one_cell(CNmat_sub[, j], bins_dt,
                                   baseline_cn = baseline_cn,
                                   keep_neutral = keep_neutral,
                                   smooth = smooth,
                                   max_short_bins = max_short_bins,
                                   min_neighbor_bins = min_neighbor_bins,
                                   max_step = max_step)
    if (nrow(seg) == 0L) return(NULL)
    seg[, cell := cell]
    setcolorder(seg, c("cell","chr","start","end","len_bp","n_bins","state","event"))
    seg
  }), use.names = TRUE, fill = TRUE)
  
  if (is.null(res)) return(data.table())
  res[]
}

smooth_cn_by_runs <- function(dt, baseline_cn = 2L,
                              max_short_bins = 1L,
                              min_neighbor_bins = 5L) {
  setorder(dt, chr, idx)
  dt[, run_id := cumsum(is.na(data.table::shift(chr)) |
                          chr != data.table::shift(chr) |
                          cn  != data.table::shift(cn))]
  
  runs <- dt[, .(chr=dplyr::first(chr), n=.N, cn=dplyr::first(cn)), by=run_id]
  
  runs[, ev := sign(cn - baseline_cn)]
  runs[, `:=`(
    prev_cn = data.table::shift(cn), next_cn = data.table::shift(cn, type="lead"),
    prev_n  = data.table::shift(n),  next_n  = data.table::shift(n, type="lead"),
    prev_ev = data.table::shift(ev), next_ev = data.table::shift(ev, type="lead"),
    prev_chr = data.table::shift(chr), next_chr = data.table::shift(chr, type="lead")
  )]
  
  runs[, new_cn := cn]
  
  runs[n <= max_short_bins &
         !is.na(prev_ev) & !is.na(next_ev) &
         chr == prev_chr & chr == next_chr &
         prev_ev != 0L & prev_ev == next_ev &
         pmax(prev_n,0L) >= min_neighbor_bins & pmax(next_n,0L) >= min_neighbor_bins,
       new_cn := fifelse(prev_n >= next_n, prev_cn, next_cn)]
  
  runs[n <= max_short_bins &
         !is.na(next_ev) & chr == next_chr &
         next_ev != 0L & ev == next_ev &
         next_n >= min_neighbor_bins,
       new_cn := next_cn]
  
  runs[n <= max_short_bins &
         !is.na(prev_ev) & chr == prev_chr &
         prev_ev != 0L & ev == prev_ev &
         prev_n >= min_neighbor_bins,
       new_cn := prev_cn]
  
  dt[runs[, .(run_id, new_cn)], on="run_id", cn := i.new_cn]
  dt[, run_id := NULL]
  dt[]
}

summarize_cna_one_cell <- function(cn_vec, bins_dt, baseline_cn = 2L, keep_neutral = FALSE,
                                   smooth = TRUE,
                                   max_short_bins = 1L, min_neighbor_bins = 5L, max_step = 1L) {
  
  dt <- cbind(bins_dt, cn = as.integer(round(cn_vec)))
  setorder(dt, idx)  # preserve CNmat_sub adjacency
  
  if (isTRUE(smooth)) {
    dt <- smooth_cn_by_runs(dt, baseline_cn = baseline_cn,
                            max_short_bins = max_short_bins,
                            min_neighbor_bins = min_neighbor_bins)
  }
  
  prev_chr <- data.table::shift(dt$chr)
  prev_cn  <- data.table::shift(dt$cn)
  new_seg  <- is.na(prev_chr) | dt$chr != prev_chr | dt$cn != prev_cn
  dt[, seg_id := cumsum(new_seg)]
  
  seg <- dt[, .(
    chr    = dplyr::first(chr),
    start  = dplyr::first(start),
    end    = data.table::last(end),
    n_bins = .N,
    state  = dplyr::first(cn),
    len_bp = data.table::last(end) - dplyr::first(start)
  ), by = seg_id][, seg_id := NULL]
  
  seg[, event := fifelse(state < baseline_cn, "DEL",
                         fifelse(state > baseline_cn, "DUP", "NEU"))]
  if (!keep_neutral) seg <- seg[event != "NEU"]
  seg[]
}

# Build a truth matrix whose columns are the predicted cells,
# but each column pulls from its site truth column (PT/LN1/ML1/MP1).
build_truth_for_pred <- function(truth_base, pred_mat) {
  pred_mat <- as.matrix(pred_mat)
  
  # align rows (bins)
  if (!is.null(rownames(truth_base)) && !is.null(rownames(pred_mat))) {
    common_bins <- intersect(rownames(truth_base), rownames(pred_mat))
    if (length(common_bins) == 0L) stop("No overlapping rownames (bins) between truth and prediction.")
    Tbin <- truth_base[common_bins, , drop = FALSE]
    Pbin <- pred_mat[common_bins, , drop = FALSE]
  } else {
    stopifnot(nrow(truth_base) == nrow(pred_mat))
    Tbin <- truth_base
    Pbin <- pred_mat
  }
  
  map <- cell_to_truthcol(colnames(Pbin))
  map <- map[!is.na(truth_col) & truth_col %in% colnames(Tbin)]
  if (nrow(map) == 0L) stop("No cells in pred_mat could be mapped to PT/LN/ML/MP truth columns.")
  
  Tcell <- sapply(map$truth_col, function(tc) Tbin[, tc], simplify = "matrix")
  colnames(Tcell) <- map$cell
  rownames(Tcell) <- rownames(Pbin)
  
  list(truth_cell = Tcell, pred_aligned = Pbin[, map$cell, drop = FALSE], map = map)
}

# -----------------------------
# 1) your existing matching functions (unchanged)
# -----------------------------
match_tp_and_bb_one_cell_state <- function(T, P, overlap_frac = 0.20) {
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
    ov_len   <- pmax(0, ov_end - ov_start)
    
    truth_len <- pmax(1, Tcc$len_bp[ti])
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

.parse_bins <- function(bin_ids) {
  m <- str_match(bin_ids, "^(chr[^:]+):(\\d+)-(\\d+)$")
  if (any(is.na(m[,1]))) stop("Bad bin_id: ", bin_ids[which(is.na(m[,1]))[1]])
  data.table(
    idx    = seq_along(bin_ids),   # <-- original order in CNmat_sub
    bin_id = bin_ids,
    chr    = m[,2],
    start  = as.integer(m[,3]),
    end    = as.integer(m[,4])
  )
}

# -----------------------------
# 2) REAL-DATA: compute BB using site-specific truth per cell
# -----------------------------
bb_tab <- rbindlist(lapply(names(methods), function(m) {
  message("BB: ", m)
  
  A <- build_truth_for_pred(truth_base, methods[[m]])
  truth_mat_cell <- A$truth_cell
  pred_mat_cell  <- A$pred_aligned
  
  truth_seg <- summarize_cna_all_cells(truth_mat_cell, baseline_cn = 2L,
                                       keep_neutral = FALSE, smooth = TRUE)
  pred_seg  <- summarize_cna_all_cells(pred_mat_cell, baseline_cn = 2L,
                                       keep_neutral = FALSE, smooth = TRUE)
  
  dt <- bb_per_method_from_segments(truth_seg, pred_seg,
                                    overlap_frac = 0.20,
                                    summary = "mean")
  dt[, method := m]
  dt
}))

bb_tab[, method := factor(method, levels = names(color_map))]
bb_tab[, BB_cell_mb := BB_cell / 1e6]

# summary table (optional)
bb_summary <- bb_tab[, .(
  n_cells_used = unique(n_cells_used)[1],
  TP_total     = sum(n_TP, na.rm = TRUE),
  BB_mean_bp   = mean(BB_cell, na.rm = TRUE),
  BB_mean_mb   = mean(BB_cell_mb, na.rm = TRUE)
), by = .(method)][order(method)]

print(bb_summary)

# plot
p_bb <- ggplot(bb_tab[!is.na(BB_cell_mb)],
               aes(x = method, y = BB_cell_mb, fill = method)) +
  geom_boxplot(width = 0.6, outlier.size = 0.35, linewidth = 0.3) +
  scale_fill_manual(values = color_map, guide = "none") +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 12),
    axis.ticks.y = element_line(linewidth = 0.4),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank()
  ) + ylim(0,35)
p_bb
ggsave("<MINTCNA_ROOT>/Result/CRC/Segmentation/Evaluation/wgs/Boundary_Bias.pdf", 
       plot = p_bb, width = 6.5, height = 6, units = "in", dpi = 600)
