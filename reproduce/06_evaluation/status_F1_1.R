truth_mat <- CNmat_sub
truth_mat[truth_mat >= 5] <- 5L 

#truth_seg <- summarize_cna_all_cells(truth_mat, baseline_cn = 2, keep_neutral = FALSE)

library(data.table)

# ---------------------------
# Parse bin IDs -> chr/start/end
# ---------------------------
parse_bins <- function(bin_ids) {
  b <- gsub("^chr", "", bin_ids, ignore.case = TRUE)
  
  if (any(grepl(":", b))) {
    chr   <- sub(":.*$", "", b)
    pos   <- sub("^.*?:", "", b)
    start <- as.integer(sub("-.*$", "", pos))
    end   <- as.integer(sub("^.*?-", "", pos))
  } else {
    sp <- strsplit(b, "_", fixed = TRUE)
    chr   <- vapply(sp, `[[`, "", 1)
    start <- as.integer(vapply(sp, `[[`, "", 2))
    end   <- as.integer(vapply(sp, `[[`, "", 3))
  }
  data.table(chr = chr, start = start, end = end)
}

# ---------------------------
# Build segments by RLE within chr (Strategy 1)
# ---------------------------
make_segments <- function(states, bins_dt) {
  stopifnot(length(states) == nrow(bins_dt))
  dt <- copy(bins_dt)
  dt[, state := states]
  
  segs <- dt[, {
    r <- rle(state)
    ends   <- cumsum(r$lengths)
    starts <- ends - r$lengths + 1L
    data.table(
      state = r$values,
      start = start[starts],
      end   = end[ends]
    )
  }, by = chr]
  
  segs
}

# truth CN -> direction
truth_dir <- function(cn) {
  fifelse(cn %in% c(0L, 1L), "DEL",
          fifelse(cn == 2L, "NEU", "DUP"))  # 3/4 -> DUP
}

# pred 1/2/3 -> direction
pred_dir <- function(x) {
  fifelse(x == 1L, "DEL",
          fifelse(x == 2L, "NEU", "DUP"))
}

# ---------------------------
# CNk directional CNV-wise counts for ONE cell
# ---------------------------
cnk_counts_one_cell <- function(truth_vec, pred_vec, bins_dt, cnk) {
  truth_vec <- as.integer(truth_vec)   # keep 0..5
  pred_vec  <- as.integer(pred_vec)
  
  truth_segs_cn <- make_segments(truth_vec, bins_dt)            # state: 0..5
  pred_segs_dir <- make_segments(pred_dir(pred_vec), bins_dt)   # state: DEL/NEU/DUP
  
  dir_k <- truth_dir(cnk)
  
  # predicted segments in correct direction (same for CN3/4/5)
  pred_k <- pred_segs_dir[state == dir_k]
  pred_k[, pid := .I]
  
  # truth segments in the SAME DIRECTION SET (same for CN3/4/5)
  truth_dir_pos <- truth_segs_cn[
    (dir_k == "DUP" & state %in% c(3L,4L,5L)) |
      (dir_k == "DEL" & state %in% c(0L,1L))
  ]
  
  # ---- FP: computed regardless of whether CNk exists in this cell
  FP <- 0L
  if (nrow(pred_k) > 0L) {
    if (nrow(truth_dir_pos) == 0L) {
      FP <- nrow(pred_k)
    } else {
      setkey(pred_k,        chr, start, end)
      setkey(truth_dir_pos, chr, start, end)
      ov2 <- foverlaps(pred_k, truth_dir_pos, type="any", nomatch=0L)
      covered <- unique(ov2$pid)
      FP <- nrow(pred_k) - length(covered)
    }
  }
  
  # ---- TP/FN: only depends on truth CNk segments
  truth_pos <- truth_segs_cn[state == cnk]
  if (nrow(truth_pos) == 0L) {
    return(data.table(cn=paste0("CN", cnk),
                      TP=0L, FP=as.integer(FP), FN=0L))
  }
  truth_pos[, id := .I]
  
  if (nrow(pred_k) == 0L) {
    TP <- 0L
    FN <- nrow(truth_pos)
  } else {
    setkey(truth_pos, chr, start, end)
    setkey(pred_k,   chr, start, end)
    ov <- foverlaps(truth_pos, pred_k, type="any", nomatch=0L)
    hit_ids <- unique(ov$id)
    TP <- length(hit_ids)
    FN <- nrow(truth_pos) - TP
  }
  
  data.table(cn=paste0("CN", cnk),
             TP=as.integer(TP), FP=as.integer(FP), FN=as.integer(FN))
}

# ---------------------------
# Evaluate one method matrix (CNV-wise), merged over length
# ---------------------------
cnk_eval <- function(truth_mat, pred_mat, cnks = c(0L,1L,3L,4L,5L)) {
  stopifnot(all(dim(truth_mat) == dim(pred_mat)))
  stopifnot(identical(rownames(truth_mat), rownames(pred_mat)))
  stopifnot(identical(colnames(truth_mat), colnames(pred_mat)))
  
  bins_dt <- parse_bins(rownames(truth_mat))
  cells <- colnames(truth_mat)
  
  per_cell <- rbindlist(lapply(cells, function(cell) {
    tv <- truth_mat[, cell]
    pv <- pred_mat[, cell]
    rbindlist(lapply(cnks, function(k) cnk_counts_one_cell(tv, pv, bins_dt, k)),
              use.names=TRUE, fill=TRUE)[, cell := cell]
  }), use.names=TRUE, fill=TRUE)
  
  # aggregate across cells to get one number per (cn)
  out <- per_cell[, .(
    TP = sum(TP, na.rm=TRUE),
    FP = sum(FP, na.rm=TRUE),
    FN = sum(FN, na.rm=TRUE),
    n_cells = uniqueN(cell)
  ), by = .(cn)]
  
  out[, precision := fifelse((TP+FP)==0, NA_real_, TP/(TP+FP))]
  out[, recall    := fifelse((TP+FN)==0, NA_real_, TP/(TP+FN))]
  out[, F1        := fifelse(is.na(precision) | is.na(recall) | (precision+recall)==0,
                             NA_real_, 2*precision*recall/(precision+recall))]
  out[]
}

pred_list <- list(
  MintCNA = cna_multicna,
  MintCNA_single = cna_multicna_single
)

summary_cnk <- rbindlist(lapply(names(pred_list), function(m) {
  tmp <- cnk_eval(truth_mat, pred_list[[m]], cnks=c(0L,1L,3L,4L,5L))
  tmp[, method := m]
  tmp
}), use.names=TRUE, fill=TRUE)

library(data.table); library(ggplot2); library(ggh4x)

df <- as.data.table(summary_cnk)
df[method == "MintCNA_single", method := "MintCNA(DNA)"]
df <- df[cn %in% c("CN3", "CN4", "CN5")]
df <- df[method %in% c("MintCNA", "MintCNA(DNA)")]

# facet labels with delta F1
lab_df <- dcast(df, cn ~ method, value.var = "F1")
lab_df[, facet_lab := sprintf("%s (\u25B3F1 = %+0.03f)", 
                              cn, MintCNA - `MintCNA(DNA)`)]

df <- merge(df, lab_df[, .(cn, facet_lab)], by = "cn")

grid <- CJ(
  precision = seq(0, 1, length.out = 300),
  recall    = seq(0, 1, length.out = 300)
)
grid[, f1 := fifelse(precision + recall == 0, 0,
                     2 * precision * recall / (precision + recall))]
f1_levels <- c(0.2, 0.4, 0.6, 0.8)

f1_text <- rbindlist(lapply(f1_levels, function(v) {
  tmp <- grid[abs(f1 - v) < 0.005]
  if (!nrow(tmp)) return(NULL)
  tmp[which.max(precision), .(precision, recall, label = paste0("F1=", v))]
}))

p <- ggplot(df, aes(precision, recall)) +
  geom_contour(
    data = grid,
    aes(x = precision, y = recall, z = f1),
    breaks = f1_levels, color = "grey85", linewidth = 0.3,
    inherit.aes = FALSE
  ) +
  geom_text(
    data = f1_text,
    aes(x = precision, y = recall, label = label),
    color = "grey60", size = 3, vjust = 0, hjust = 1,
    inherit.aes = FALSE
  ) +
  geom_point(aes(color = method), shape = 17, size = 6) +
  scale_color_manual(
    values = c(
      "MintCNA" = "#1B9E77",
      "MintCNA(DNA)" = "#D95F02"
    ),
    breaks = c("MintCNA", "MintCNA(DNA)")
  ) +
  coord_fixed(xlim = c(0.4, 1), ylim = c(0.4, 1)) +
  ggh4x::facet_wrap2(~facet_lab, nrow = 1, axes = "all") +
  labs(x = "Precision", y = "Recall", color = NULL) +
  theme_classic(base_size = 11) +
  theme(legend.position = "bottom")

# print(p)
# ggsave(file.path(CNA_DIR, "Result/Status_F1.pdf"),
#        plot = p, width = 9, height = 5.2, units = "in", dpi = 600)

cairo_pdf(file.path(CNA_DIR, "Result/Status_F1.pdf"), width = 9, height = 5.2)
print(p)
dev.off()
