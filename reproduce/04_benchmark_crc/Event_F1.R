load("<MINTCNA_ROOT>/Result/CRC/Benchmark/SCOPE.RData")
load("<MINTCNA_ROOT>/Result/CRC/Benchmark/MultiCNA_AE_197_d250.RData")
load("<MINTCNA_ROOT>/Result/CRC/Benchmark/HMMcopy.RData")
load("<MINTCNA_ROOT>/Result/CRC/Benchmark/DNAcopy.RData")
load("<MINTCNA_ROOT>/Result/CRC/Benchmark/SCCNV.RData")
load("<MINTCNA_ROOT>/Result/CRC/Benchmark/SCONCE.RData")
load("<MINTCNA_ROOT>/Result/CRC/Benchmark/rcCAE.RData")
load("<MINTCNA_ROOT>/Result/CRC/Benchmark/MultiCNA_AE_197_single_d250.RData")

colnames(ref_scope)[1] <- "chr";
cna_scope <- cna_scope[,colnames(cna_scope) %in% colnames(cna_multicna_ae)]
cna_scope <- cna_scope + 1L; cna_scope[cna_scope >= 5] <- 5L 
cna_hmmcopy <- cna_hmmcopy[,colnames(cna_hmmcopy) %in% colnames(cna_multicna_ae)]
cna_sconce <- cna_sconce[,colnames(cna_sconce) %in% colnames(cna_multicna_ae)]
cna_sccnv <- cna_sccnv[,colnames(cna_sccnv) %in% colnames(cna_multicna_ae)]
cna_rccae <- cna_rccae[,colnames(cna_rccae) %in% colnames(cna_multicna_ae)]

# build bin IDs from the reference you want to use (dnacopy here)
bin_id <- with(ref_dnacopy, paste0(chr, ":", start, "-", end))

# enforce identical rownames across all matrices
rownames(cna_dnacopy)      <- bin_id
rownames(cna_multicna_ae)  <- bin_id
rownames(cna_multicna_ae_single)  <- bin_id
rownames(cna_scope)        <- bin_id
rownames(cna_hmmcopy)      <- bin_id
rownames(cna_sconce)       <- bin_id
rownames(cna_sccnv)        <- bin_id
rownames(cna_rccae)        <- bin_id
# rownames(cna_scicone)    <- bin_id   # if applicable

library(data.table)
library(viridis)
library(ggplot2)

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
cna_multicna_ae = coerce_123(cna_multicna_ae)
cna_multicna_ae_single = coerce_123(cna_multicna_ae_single)
cna_scope    = coerce_123(cna_scope)
cna_sconce    = coerce_123(cna_sconce)
cna_sccnv    = coerce_123(cna_sccnv)
cna_rccae    = coerce_123(cna_rccae)

######## Event F1 ############
library(data.table)
library(IRanges)
library(stringr)

# -----------------------------
# A) map cell -> site truth column (same as RMSE/BB)
# -----------------------------
cell_to_truthcol <- function(cell_names) {
  site <- sub("^.*_(PT|LN|ML|MP).*?$", "\\1", cell_names)
  site[!site %in% c("PT","LN","ML","MP")] <- NA_character_
  
  truth_col <- fifelse(site == "PT", "PT",
                       fifelse(site == "LN", "LN1",
                               fifelse(site == "ML", "ML1",
                                       fifelse(site == "MP", "MP1", NA_character_))))
  data.table(cell = cell_names, site = site, truth_col = truth_col)
}

build_truth_for_pred <- function(truth_base, pred_mat) {
  pred_mat <- as.matrix(pred_mat)
  
  # align bins (rows)
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
  
  # truth per cell = pull the appropriate truth column
  Tcell <- sapply(map$truth_col, function(tc) Tbin[, tc], simplify = "matrix")
  colnames(Tcell) <- map$cell
  rownames(Tcell) <- rownames(Pbin)
  
  list(
    truth_cell = Tcell,
    pred_aligned = Pbin[, map$cell, drop = FALSE],
    map = map
  )
}

# -----------------------------
# B) event-level TP/FP/FN per cell+state (your function)
# -----------------------------
cnvwise_counts_nostrata <- function(truth_seg, pred_seg, overlap_frac = 0.20) {
  truth_seg <- copy(truth_seg)
  pred_seg  <- copy(pred_seg)
  
  # label states like your old code (DEL/DUP only)
  truth_seg[, state_lab := fifelse(state == 1L, "del", "dup")]
  pred_seg[,  state_lab := fifelse(state == 1L, "del", "dup")]
  
  keys <- unique(rbind(
    truth_seg[, .(cell, state_lab)],
    pred_seg[,  .(cell, state_lab)]
  ), use.names = TRUE, fill = TRUE)
  
  res_list <- vector("list", nrow(keys))
  
  for (k in seq_len(nrow(keys))) {
    cell_k <- keys$cell[k]
    st_k   <- keys$state_lab[k]
    
    T <- truth_seg[cell == cell_k & state_lab == st_k]
    P <- pred_seg [cell == cell_k & state_lab == st_k]
    
    TP <- 0L
    FN <- nrow(T)
    FP <- nrow(P)
    
    if (nrow(T) > 0L && nrow(P) > 0L) {
      T[, tid := .I]
      P[, pid := .I]
      
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
        ov_len   <- pmax(0, ov_end - ov_start)  # same convention as your code
        
        truth_len <- pmax(1, Tcc$len_bp[ti])
        prop <- ov_len / truth_len
        
        data.table(tid = Tcc$tid[ti], pid = Pcc$pid[pj], prop = prop, ov_len = ov_len)[prop > 0]
      }), use.names = TRUE, fill = TRUE)
      
      if (!is.null(cand) && nrow(cand) > 0L) {
        cand <- cand[prop >= overlap_frac]
        if (nrow(cand) > 0L) {
          setorder(cand, -prop, -ov_len)
          used_t <- rep(FALSE, nrow(T))
          used_p <- rep(FALSE, nrow(P))
          m <- 0L
          for (ii in seq_len(nrow(cand))) {
            tt <- cand$tid[ii]
            pp <- cand$pid[ii]
            if (!used_t[tt] && !used_p[pp]) {
              used_t[tt] <- TRUE
              used_p[pp] <- TRUE
              m <- m + 1L
            }
          }
          TP <- m
          FN <- nrow(T) - TP
          FP <- nrow(P) - TP
        }
      }
    }
    
    res_list[[k]] <- data.table(cell = cell_k, state = st_k, TP = TP, FP = FP, FN = FN)
  }
  
  res <- rbindlist(res_list, use.names = TRUE, fill = TRUE)
  res[, precision := fifelse((TP + FP) == 0, NA_real_, TP/(TP + FP))]
  res[, recall    := fifelse((TP + FN) == 0, NA_real_, TP/(TP + FN))]
  res[, F1        := fifelse(is.na(precision) | is.na(recall) | (precision + recall) == 0,
                             NA_real_, 2*precision*recall/(precision + recall))]
  res[]
}

# -----------------------------
# C) run per method on real data
# -----------------------------
truth_base <- as.matrix(cna_dnacopy)   # cols: PT, LN1, ML1, MP1

all_counts <- rbindlist(lapply(names(methods), function(m) {
  message("EventF1: ", m)
  
  A <- build_truth_for_pred(truth_base, methods[[m]])
  truth_mat_cell <- A$truth_cell
  pred_mat_cell  <- A$pred_aligned
  
  # segments
  truth_seg <- summarize_cna_all_cells(truth_mat_cell, baseline_cn = 2L, keep_neutral = FALSE, smooth = TRUE)
  pred_seg  <- summarize_cna_all_cells(pred_mat_cell,  baseline_cn = 2L, keep_neutral = FALSE, smooth = TRUE)
  
  # counts + label
  cnt <- cnvwise_counts_nostrata(truth_seg, pred_seg, overlap_frac = 0.20)
  cnt[, method := m]
  cnt
}), use.names = TRUE, fill = TRUE)

# -----------------------------
# D) summarize to method x state table (same as simulation)
# -----------------------------
all_counts[, method := factor(method, levels = names(color_map))]
all_counts[, state  := factor(state, levels = c("del","dup"))]

summary_all <- all_counts[, .(
  TP = sum(TP, na.rm = TRUE),
  FP = sum(FP, na.rm = TRUE),
  FN = sum(FN, na.rm = TRUE),
  n_cells = uniqueN(cell)
), by = .(method, state)][order(method, state)]

summary_all[, precision := fifelse((TP + FP) == 0, NA_real_, TP/(TP + FP))]
summary_all[, recall    := fifelse((TP + FN) == 0, NA_real_, TP/(TP + FN))]
summary_all[, F1        := fifelse(is.na(precision) | is.na(recall) | (precision + recall) == 0,
                                   NA_real_, 2*precision*recall/(precision + recall))]

# ordering (reuse your method_order if defined)
summary_all[, method := factor(method, levels = method_order)]
summary_all[, state  := factor(state, levels = c("del","dup"))]

summary_one <- summary_all[, .(
  TP = sum(TP, na.rm=TRUE),
  FP = sum(FP, na.rm=TRUE),
  FN = sum(FN, na.rm=TRUE),
  n_cells = max(n_cells, na.rm=TRUE)  # just to keep a value around
), by = .(method)]

summary_one[, precision := fifelse((TP + FP) == 0, NA_real_, TP/(TP + FP))]
summary_one[, recall    := fifelse((TP + FN) == 0, NA_real_, TP/(TP + FN))]
summary_one[, F1        := fifelse(is.na(precision) | is.na(recall) | (precision + recall) == 0,
                                   NA_real_, 2*precision*recall/(precision + recall))]
summary_one[, method := factor(method, levels = method_order)]

## Plot
grid_data <- CJ(precision = seq(0, 1, length.out = 300),
                recall    = seq(0, 1, length.out = 300))
grid_data[, f1 := fifelse(precision + recall == 0, 0,
                          2 * precision * recall / (precision + recall))]

f1_levels <- c(0.2, 0.4, 0.6, 0.8)

label_positions <- rbindlist(lapply(f1_levels, function(f1v) {
  tmp <- grid_data[abs(f1 - f1v) < 0.005]
  if (nrow(tmp) == 0) return(NULL)
  # put label near upper-right portion of each contour
  tmp[which.max(precision)][, .(precision, recall, label = paste0("F1=", f1v))]
}), use.names = TRUE, fill = TRUE)

shape_map <- setNames(rep(21, length(method_order)), method_order)  # filled circles
shape_map[1:2] <- 24

p_pr_one <- ggplot(summary_one, aes(x = precision, y = recall)) +
  geom_contour(data = grid_dup,
               aes(x = precision, y = recall, z = f1),
               breaks = f1_levels, color = "grey80", linewidth = 0.3,
               inherit.aes = FALSE) +
  geom_text(data = label_dup,
            aes(x = precision, y = recall, label = label),
            color = "grey60", size = 3, vjust = 0, hjust = 1,
            inherit.aes = FALSE) +
  geom_point(aes(color = method, fill = method, shape = method),
             size = 6.5, stroke = 0.4) +
  scale_shape_manual(values = shape_map, guide = "none") +
  scale_color_manual(values = color_map, guide = "none") +
  scale_fill_manual(values = color_map, guide = "none") +
  coord_fixed(ratio = 1, xlim = c(0,1), ylim = c(0,1)) +
  labs(x = "Precision", y = "Recall") +
  theme_classic(base_size = 11) +
  theme(legend.position = "none")

print(p_pr_one)
ggsave("<MINTCNA_ROOT>/Result/CRC/Segmentation/Evaluation/wgs/PR_overall_onepoint.pdf",
       plot = p_pr_one, width = 6.5, height = 5.2, units = "in", dpi = 600)

## =========================================================
## Plot B: split by DEL vs DUP (still no length stratification)
##         -> two panels, each method one point per panel
## =========================================================
# --------------------------
# DUP
# --------------------------
summary_dup <- subset(summary_all, state == "dup")

grid_dup  <- if ("state" %in% names(grid_data)) subset(grid_data, state == "dup") else grid_data
label_dup <- if ("state" %in% names(label_positions)) subset(label_positions, state == "dup") else label_positions

p_pr_dup <- ggplot(summary_dup, aes(x = precision, y = recall)) +
  geom_contour(data = grid_dup,
               aes(x = precision, y = recall, z = f1),
               breaks = f1_levels, color = "grey80", linewidth = 0.3,
               inherit.aes = FALSE) +
  geom_text(data = label_dup,
            aes(x = precision, y = recall, label = label),
            color = "grey60", size = 3, vjust = 0, hjust = 1,
            inherit.aes = FALSE) +
  geom_point(aes(color = method, fill = method, shape = method),
             size = 6.5, stroke = 0.4) +
  scale_shape_manual(values = shape_map, guide = "none") +
  scale_color_manual(values = color_map, guide = "none") +
  scale_fill_manual(values = color_map, guide = "none") +
  coord_fixed(ratio = 1, xlim = c(0,1), ylim = c(0,1)) +
  labs(x = "Precision", y = "Recall") +
  theme_classic(base_size = 11) +
  theme(legend.position = "none")

print(p_pr_dup)

ggsave("<MINTCNA_ROOT>/Result/CRC/Segmentation/Evaluation/wgs/PR_overall_dup.pdf",
       plot = p_pr_dup, width = 6.5, height = 5.2, units = "in", dpi = 600)

# --------------------------
# DEL
# --------------------------
summary_del <- subset(summary_all, state == "del")

grid_del  <- if ("state" %in% names(grid_data)) subset(grid_data, state == "del") else grid_data
label_del <- if ("state" %in% names(label_positions)) subset(label_positions, state == "del") else label_positions

p_pr_del <- ggplot(summary_del, aes(x = precision, y = recall)) +
  geom_contour(data = grid_del,
               aes(x = precision, y = recall, z = f1),
               breaks = f1_levels, color = "grey80", linewidth = 0.3,
               inherit.aes = FALSE) +
  geom_text(data = label_del,
            aes(x = precision, y = recall, label = label),
            color = "grey60", size = 3, vjust = 0, hjust = 1,
            inherit.aes = FALSE) +
  geom_point(aes(color = method, fill = method, shape = method),
             size = 6.5, stroke = 0.4) +
  scale_shape_manual(values = shape_map, guide = "none") +
  scale_color_manual(values = color_map, guide = "none") +
  scale_fill_manual(values = color_map, guide = "none") +
  coord_fixed(ratio = 1, xlim = c(0,1), ylim = c(0,1)) +
  labs(x = "Precision", y = "Recall") +
  theme_classic(base_size = 11) +
  theme(legend.position = "none")

print(p_pr_del)

ggsave("<MINTCNA_ROOT>/Result/CRC/Segmentation/Evaluation/wgs/PR_overall_del.pdf",
       plot = p_pr_dup, width = 6.5, height = 5.2, units = "in", dpi = 600)
