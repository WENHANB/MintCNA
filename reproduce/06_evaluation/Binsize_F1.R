BASE   <- "c3_N100_100kb_l1e4"
CNA_DIR <- file.path("<MINTCNA_ROOT>/Simulation/cnasim", BASE)  

load(file.path(CNA_DIR,"Calling/MintCNA_h5.RData"))
load(file.path(CNA_DIR,"Calling/MintCNA_h5_single.RData"))

library(data.table)
library(viridis)

coerce_123 <- function(mat) {
  v <- as.integer(mat)
  
  # If already 1/2/3 only, keep as-is
  if (all(v %in% c(1L, 2L, 3L, NA_integer_))) {
    return(matrix(v, nrow=nrow(mat), ncol=ncol(mat), dimnames=dimnames(mat)))
  }
  
  out <- rep(NA_integer_, length(v))
  out[v %in% c(1L, 2L)] <- 1L   # deletions
  out[v == 3L]          <- 2L   # neutral
  out[v >= 4L]          <- 3L   # duplications: 4,5,... -> 3
  
  matrix(out, nrow=nrow(mat), ncol=ncol(mat), dimnames=dimnames(mat))
}

#---------------------------
# Helper: length stratum from bp length
#---------------------------
len_stratum <- function(len_bp) {
  fifelse(len_bp < 5e6, "<5Mb",
          fifelse(len_bp <= 10e6, "5-10Mb", ">10Mb"))
}

#---------------------------
# Convert bin-level states to CNV segments (events) per cell
# Returns: data.table(cell, chr, start, end, len_bp, state)
# Only keeps state %in% keep_states (e.g., c(1,3) for del/dup)
#---------------------------
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

summarize_cna_all_cells <- function(CNmat_sub, baseline_cn = 2L, keep_neutral = FALSE,
                                    smooth = TRUE,
                                    max_short_bins = 1L, min_neighbor_bins = 5L, max_step = 1L) {
  
  bins_dt <- .parse_bins(rownames(CNmat_sub))
  
  res <- rbindlist(lapply(seq_len(ncol(CNmat_sub)), function(j) {
    cell <- colnames(CNmat_sub)[j]
    seg  <- summarize_cna_one_cell(CNmat_sub[, j], bins_dt,
                                   baseline_cn = baseline_cn,
                                   keep_neutral = keep_neutral,
                                   smooth = smooth,
                                   max_short_bins = max_short_bins,
                                   min_neighbor_bins = min_neighbor_bins,
                                   max_step = max_step)
    if (nrow(seg) == 0) return(NULL)
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
  # dt must have: idx, chr, cn (already in CNmat_sub row order)
  setorder(dt, chr, idx)
  dt[, run_id := cumsum(is.na(data.table::shift(chr)) | chr != data.table::shift(chr) | cn != data.table::shift(cn))]
  
  runs <- dt[, .(chr=dplyr::first(chr), n=.N, cn=dplyr::first(cn)), by=run_id]
  
  # event direction: -1 DEL, 0 NEU, +1 DUP
  runs[, ev := sign(cn - baseline_cn)]
  runs[, `:=`(
    prev_cn = data.table::shift(cn), next_cn = data.table::shift(cn, type="lead"),
    prev_n  = data.table::shift(n),  next_n  = data.table::shift(n, type="lead"),
    prev_ev = data.table::shift(ev), next_ev = data.table::shift(ev, type="lead"),
    prev_chr = data.table::shift(chr), next_chr = data.table::shift(chr, type="lead")
  )]
  
  runs[, new_cn := cn]
  
  # --- sandwiched short run: merge if both sides share same event direction (DUP/DUP or DEL/DEL)
  # choose the *longer* neighbor's CN as the target (handles 4 vs 6 vs 7 automatically)
  runs[n <= max_short_bins &
         !is.na(prev_ev) & !is.na(next_ev) &
         chr == prev_chr & chr == next_chr &
         prev_ev != 0L & prev_ev == next_ev &
         pmax(prev_n,0L) >= min_neighbor_bins & pmax(next_n,0L) >= min_neighbor_bins,
       new_cn := fifelse(prev_n >= next_n, prev_cn, next_cn)]
  
  # --- edge short run: absorb into adjacent run if same event direction and neighbor long enough
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
  
  # write back to bins
  dt[runs[, .(run_id, new_cn)], on="run_id", cn := i.new_cn]
  dt[, run_id := NULL]
  dt[]
}


summarize_cna_one_cell <- function(cn_vec, bins_dt, baseline_cn = 2L, keep_neutral = FALSE,
                                   smooth = TRUE,
                                   max_short_bins = 1L, min_neighbor_bins = 5L, max_step = 1L) {
  dt <- cbind(bins_dt, cn = as.integer(round(cn_vec)))
  setorder(dt, idx)  # keep CNmat_sub adjacency
  dt <- smooth_cn_by_runs(dt, baseline_cn = baseline_cn,
                          max_short_bins = 1L, min_neighbor_bins = 5L)
  
  
  prev_chr <- data.table::shift(dt$chr)
  prev_cn  <- data.table::shift(dt$cn)
  
  new_seg <- is.na(prev_chr) | dt$chr != prev_chr | dt$cn != prev_cn
  dt[, seg_id := cumsum(new_seg)]
  
  seg <- dt[, .(
    chr       = dplyr::first(chr),
    start     = dplyr::first(start),
    end       = last(end),
    n_bins    = .N,
    state     = dplyr::first(cn),
    len_bp = last(end) - dplyr::first(start)
  ), by = seg_id][, seg_id := NULL]
  
  seg[, event := fifelse(state < baseline_cn, "DEL",
                         fifelse(state > baseline_cn, "DUP", "NEU"))]
  if (!keep_neutral) seg <- seg[event != "NEU"]
  seg[]
}

#---------------------------
# CNV-wise one-to-one matching within (cell, state)
# TP rule: overlap / truth_len >= overlap_frac (default 0.2)
# Stratification:
#   - truth_stratum: based on truth segment length
#   - pred_stratum : based on predicted segment length (for FP counting)
#
# Returns: TP/FP/FN per cell x state x stratum
#---------------------------
cnvwise_counts <- function(truth_seg, pred_seg, overlap_frac = 0.20) {
  # add strata
  truth_seg <- copy(truth_seg)
  pred_seg  <- copy(pred_seg)
  truth_seg[, stratum := len_stratum(len_bp)]
  pred_seg[,  stratum := len_stratum(len_bp)]
  
  # states: 1=del, 3=dup
  truth_seg[, state_lab := fifelse(state == 1L, "del", "dup")]
  pred_seg[,  state_lab := fifelse(state == 1L, "del", "dup")]
  
  # unique keys to iterate
  keys <- unique(rbind(
    truth_seg[, .(cell, state_lab)],
    pred_seg[,  .(cell, state_lab)]
  ), use.names = TRUE, fill = TRUE)
  
  strata_levels <- c("<5Mb", "5-10Mb", ">10Mb")
  
  res_list <- vector("list", nrow(keys))
  for (k in seq_len(nrow(keys))) {
    cell_k <- keys$cell[k]
    st_k   <- keys$state_lab[k]
    
    T0 <- truth_seg[cell == cell_k & state_lab == st_k]
    P0 <- pred_seg [cell == cell_k & state_lab == st_k]
    
    # If nothing at all, still return 0 rows across strata
    if (nrow(T0) == 0L) T0 <- T0[0]
    if (nrow(P0) == 0L) P0 <- P0[0]
    
    # Prepare per-stratum counts
    cell_state_rows <- vector("list", length(strata_levels))
    
    for (s in seq_along(strata_levels)) {
      str_s <- strata_levels[s]
      T <- T0[stratum == str_s]
      P <- P0[stratum == str_s]  # FP counted within same stratum by predicted length
      
      # default counts
      TP <- 0L
      FN <- nrow(T)
      FP <- nrow(P)
      
      if (nrow(T) > 0L && nrow(P) > 0L) {
        # assign ids
        T[, tid := .I]
        P[, pid := .I]
        
        # candidate pairs on same chr with positive overlap
        # build overlaps by chr to keep it efficient
        cand <- rbindlist(lapply(intersect(unique(T$chr), unique(P$chr)), function(cc) {
          Tcc <- T[chr == cc]
          Pcc <- P[chr == cc]
          if (nrow(Tcc) == 0L || nrow(Pcc) == 0L) return(NULL)
          
          hits <- findOverlaps(
            IRanges(Tcc$start, Tcc$end),
            IRanges(Pcc$start, Pcc$end)
          )
          if (length(hits) == 0L) return(NULL)
          
          ti <- queryHits(hits)
          pj <- subjectHits(hits)
          
          ov_start <- pmax(Tcc$start[ti], Pcc$start[pj])
          ov_end   <- pmin(Tcc$end[ti],   Pcc$end[pj])
          ov_len   <- pmax(0, ov_end - ov_start)
          
          data.table(
            tid = Tcc$tid[ti],
            pid = Pcc$pid[pj],
            ov_len = ov_len,
            truth_len = Tcc$len_bp[ti]
          )[, prop := ov_len / truth_len][prop > 0]
        }), use.names = TRUE, fill = TRUE)
        
        if (!is.null(cand) && nrow(cand) > 0L) {
          # keep those meeting the TP criterion
          cand <- cand[prop >= overlap_frac]
          if (nrow(cand) > 0L) {
            # greedy maximum-prop one-to-one matching
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
      
      cell_state_rows[[s]] <- data.table(
        cell = cell_k,
        state = st_k,
        stratum = str_s,
        TP = TP, FP = FP, FN = FN
      )
    }
    
    res_list[[k]] <- rbindlist(cell_state_rows)
  }
  
  res <- rbindlist(res_list, use.names = TRUE, fill = TRUE)
  
  # add F1 (event-wise) if you want it directly
  res[, precision := fifelse((TP + FP) == 0, NA_real_, TP / (TP + FP))]
  res[, recall    := fifelse((TP + FN) == 0, NA_real_, TP / (TP + FN))]
  res[, F1 := fifelse(is.na(precision) | is.na(recall) | (precision + recall) == 0,
                      NA_real_, 2 * precision * recall / (precision + recall))]
  
  res[]
}

summary_all_bin <- summary_mss_bin <- summary_one_bin <- c()

for (bin in c(50,100,200,500)) {
  BASE   <- paste0("c3_N100_",bin,"kb_l1e4")
  CNA_DIR <- file.path("<MINTCNA_ROOT>/Simulation/cnasim", BASE) 
  print(CNA_DIR)
  
  load(file.path(CNA_DIR,"Calling/MintCNA_h5.RData"))
  load(file.path(CNA_DIR,"Calling/MintCNA_h5_single.RData"))

##################### Loading Ground Truth #####################

# ---------------------------
# 1) Read CN bin truth (profiles.tsv)
# ---------------------------
profiles <- fread(file.path(CNA_DIR, "profiles.tsv"))
colnames(profiles)[5] <- "cn_states"
setnames(profiles, tolower(names(profiles)))

# cn_total from "1.1" (two haplotypes split by ".")
profiles[, cn_total := {
  x <- trimws(as.character(cn_states))
  ifelse(grepl("\\.", x),
         rowSums(matrix(as.numeric(unlist(tstrsplit(x, "\\.", perl=TRUE))), ncol=2, byrow=FALSE)),
         round(as.numeric(x))) |> as.integer()
}]

# standardize chrom and make ids
profiles[!grepl("^chr", chrom), chrom := paste0("chr", chrom)]
profiles[, bin_id := paste0(chrom, ":", start, "-", end)]

# numeric keys for proper ordering (cell1,cell2,... and chr1,chr2,...)
profiles[, cell_num := as.integer(sub("^cell", "", cell))]
profiles[, chr_num  := as.integer(sub("^chr",  "", chrom))]

# keep only real cells if ancestors are present (optional; comment out if not needed)
profiles <- profiles[grepl("^cell\\d+$", cell)]

# order rows (clean)
setorder(profiles, chr_num, start, end, cell_num)

# correct cell order
cells <- profiles[, unique(.SD[order(cell_num)]$cell), .SDcols = c("cell","cell_num")]

# bin x cell CN matrix with correct column order
CNdt  <- dcast(profiles, bin_id ~ cell, value.var = "cn_total")
bin_ids <- CNdt$bin_id
CNmat <- as.matrix(CNdt[, ..cells])
rownames(CNmat) <- bin_ids

# enforce row order by chr/start
ord_bin <- order(
  as.integer(sub("^chr", "", sub(":.*$", "", rownames(CNmat)))),
  as.integer(sub("^.*:(\\d+)-.*$", "\\1", rownames(CNmat)))
)
CNmat <- CNmat[ord_bin, , drop=FALSE]

cn_bins <- rownames(CNmat)

chr  <- sub("^chr", "", sub(":.*$", "", cn_bins))          # "1"
pos  <- sub("^.*:", "", cn_bins)                           # "0-100000"
st   <- as.integer(sub("-.*$", "", pos))                   # 0
en   <- as.integer(sub("^.*-", "", pos))                   # 100000

key_CN <- paste(chr, st, en, sep = "_")                    # "1_0_100000"
names(key_CN) <- cn_bins                                   # map key -> original rowname

# 2) Keys needed from cna_multicna (already in that style)
key_need <- rownames(cna_multicna)

# 3) Match and subset CNmat to those rows (same order as cna_multicna)
idx <- match(key_need, key_CN)

# Sanity check for missing bins
if (anyNA(idx)) {
  missing <- key_need[is.na(idx)]
  stop("Some bins in cna_multicna not found in CNmat. Example: ", missing[1])
}

label <- fread(file.path(CNA_DIR,"cell_types.tsv"))
nc_samples <- label[which(label$normal == "normal"), root]
nor_idx <- which(colnames(CNmat) %in% nc_samples)

CNmat_sub <- CNmat[idx, , drop = FALSE]
CNmat_sub <- CNmat_sub[,-nor_idx]

# set all method rownames identical to truth_mat rownames
set_truth_rownames <- function(mat, truth_mat) {
  stopifnot(nrow(mat) == nrow(truth_mat))
  rownames(mat) <- rownames(truth_mat)
  mat
}
###############################################################
truth_mat <- coerce_123(CNmat_sub+1)

cna_multicna        <- set_truth_rownames(coerce_123(cna_multicna),        truth_mat)
cna_multicna_single <- set_truth_rownames(coerce_123(cna_multicna_single), truth_mat)

cna_multicna <- coerce_123(cna_multicna)
cna_multicna_single <- coerce_123(cna_multicna_single)

methods <- list(mintcna = cna_multicna,mintcna_single = cna_multicna_single)

name_map     <- c(mintcna="MintCNA", mintcna_single="MintCNA(DNA)")

method_order <- c("MintCNA", "MintCNA(DNA)")
truth_seg <- summarize_cna_all_cells(truth_mat, baseline_cn = 2, keep_neutral = FALSE)

# run each method
align_to_truth_by_colname_fillNA <- function(mat, truth_mat) {
  stopifnot(nrow(mat) == nrow(truth_mat))
  
  rownames(mat) <- rownames(truth_mat)
  
  out <- matrix(NA_integer_, nrow = nrow(truth_mat), ncol = ncol(truth_mat),
                dimnames = dimnames(truth_mat))
  
  common <- intersect(colnames(truth_mat), colnames(mat))
  if (length(common) == 0L) stop("No overlapping colnames between mat and truth_mat.")
  
  out[, common] <- mat[, common, drop = FALSE]
  out
}

all_counts <- rbindlist(lapply(names(methods), function(m) {
  print(m)
  pred_mat <- coerce_123(methods[[m]])
  
  # choose ONE:
  #pred_mat <- align_to_truth_by_colname(pred_mat, truth_mat)
  pred_mat <- align_to_truth_by_colname_fillNA(pred_mat, truth_mat)
  
  pred_seg <- summarize_cna_all_cells(pred_mat, baseline_cn = 2, keep_neutral = FALSE)
  
  cnt <- cnvwise_counts(truth_seg, pred_seg, overlap_frac = 0.20)
  cnt[, method := name_map[m]]
  cnt
}), use.names = TRUE, fill = TRUE)


# ordering for plots/tables
all_counts[, method := factor(method, levels = method_order)]
all_counts[, state  := factor(state, levels = c("del","dup"))]
all_counts[, stratum := factor(stratum, levels = c("<5Mb","5-10Mb",">10Mb"))]

## Across cells
# all_counts has: method, cell, state, stratum, TP, FP, FN
summary_mss <- all_counts[, .(
  TP = sum(TP, na.rm = TRUE),
  FP = sum(FP, na.rm = TRUE),
  FN = sum(FN, na.rm = TRUE),
  n_cells = uniqueN(cell)
), by = .(method, state, stratum)]

summary_mss[, precision := fifelse((TP + FP) == 0, NA_real_, TP/(TP + FP))]
summary_mss[, recall    := fifelse((TP + FN) == 0, NA_real_, TP/(TP + FN))]
summary_mss[, F1        := fifelse(is.na(precision) | is.na(recall) | (precision + recall) == 0,
                                   NA_real_, 2*precision*recall/(precision + recall))]

# optional: enforce plotting order
summary_mss[, method  := factor(method, levels = method_order)]
summary_mss[, state   := factor(state, levels = c("del","dup"))]
summary_mss[, stratum := factor(stratum, levels = c("<5Mb","5-10Mb",">10Mb"))]
summary_mss$Bin <- bin

summary_all <- all_counts[, .(
  TP = sum(TP, na.rm = TRUE),
  FP = sum(FP, na.rm = TRUE),
  FN = sum(FN, na.rm = TRUE),
  n_cells = uniqueN(cell)
), by = .(method, state)]

summary_all[, precision := fifelse((TP + FP) == 0, NA_real_, TP/(TP + FP))]
summary_all[, recall    := fifelse((TP + FN) == 0, NA_real_, TP/(TP + FN))]
summary_all[, F1        := fifelse(is.na(precision) | is.na(recall) | (precision + recall) == 0,
                                   NA_real_, 2*precision*recall/(precision + recall))]

# ordering (reuse your method_order if defined)
summary_all[, method := factor(method, levels = method_order)]
summary_all[, state  := factor(state, levels = c("del","dup"))]
summary_all$Bin <- bin

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
summary_one$Bin <- bin

summary_all_bin <- rbind(summary_all_bin,summary_all)
summary_mss_bin <- rbind(summary_mss_bin,summary_mss)
summary_one_bin <- rbind(summary_one_bin,summary_one)
}

################################ Plot ################################
keep_bins <- c(50, 100, 200, 500) 

df_plot <- summary_mss_bin %>%
  filter(stratum %in% c("5-10Mb", ">10Mb"),
         method %in% c("MintCNA", "MintCNA(DNA)"),
         state  %in% c("del", "dup"),
         Bin %in% keep_bins) %>%
  mutate(
    Bin = factor(Bin, levels = keep_bins),  # 离散轴 => 等距
    method = factor(method, levels = c("MintCNA", "MintCNA(DNA)")),
    state  = factor(state,  levels = c("del", "dup")),
    stratum = factor(stratum, levels = c("5-10Mb", ">10Mb"))
  )

p <- ggplot(df_plot,
            aes(x = Bin, y = F1,
                group = interaction(method, state),
            color = method,
            shape = state)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2.6) +
    facet_wrap(~ stratum, nrow = 1) +
    scale_shape_manual(values = c(del = 17, dup = 15)) +  # del=triangle, dup=square
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
        legend.position = "bottom") +
    labs(x = "Bin (kb)", y = "F1 score", color = "Method", shape = "State")

print(p)
ggsave("<MINTCNA_ROOT>/Simulation/results/Bin_F1.png", plot = p, width = 8, height = 6, units = "in", dpi = 600)



