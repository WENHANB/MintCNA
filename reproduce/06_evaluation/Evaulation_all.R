args <- commandArgs(trailingOnly = TRUE)
BASE <- if (length(args) >= 1) args[1] else "c3_N100_100kb_l1e4"
print(BASE)
CNA_DIR <- file.path("<MINTCNA_ROOT>/Simulation/cnasim", BASE)

RES_DIR <- file.path(CNA_DIR, "Result")
dir.create(RES_DIR, recursive = TRUE, showWarnings = FALSE)

load("<MINTCNA_ROOT>/Simulation/results/cnasim/c3_N100_100kb_l1e4/Calling/MintCNA_h5.RData")
load(file.path(CNA_DIR,"Calling/MintCNA_h5.RData"))
#load(file.path(CNA_DIR,"Calling/MintCNA_h5_single.RData"))
load(file.path(CNA_DIR,"Calling/MintCNA_h5_single.RData"))
load(file.path(CNA_DIR,"Calling/HMMcopy.RData"))
load(file.path(CNA_DIR,"Calling/SCCNV.RData"))
load(file.path(CNA_DIR,"Calling/SCOPE.RData"))
load(file.path(CNA_DIR,"Calling/SCONCE.RData"))
load(file.path(CNA_DIR,"Calling/rcCAE.RData"))
load(file.path(CNA_DIR,"Calling/SCICoNE.RData"))

cna_hmmcopy <- cna_hmmcopy[,colnames(cna_hmmcopy) %in% colnames(cna_multicna)]
cna_hmmcopy[cna_hmmcopy >= 5] <- 5L 
cna_scope <- cna_scope[,colnames(cna_scope) %in% colnames(cna_multicna)]
cna_scope <- cna_scope + 1L; cna_scope[cna_scope >= 5] <- 5L 

cna_sconce <- cna_sconce[, colnames(cna_sconce) %in% colnames(cna_multicna), drop = FALSE]
cell_num <- suppressWarnings(as.integer(gsub(".*?(\\d+).*", "\\1", colnames(cna_sconce))))
ord <- order(cell_num, colnames(cna_sconce), na.last = TRUE)  # tie-break by name; NAs go last
cna_sconce <- cna_sconce[, ord, drop = FALSE]

cna_sccnv <- cna_sccnv[,colnames(cna_sccnv) %in% colnames(cna_multicna)]
cna_rccae <- cna_rccae[,colnames(cna_rccae) %in% colnames(cna_multicna)]
cna_scicone <- cna_scicone[,colnames(cna_scicone) %in% colnames(cna_multicna)]

library(data.table)
library(viridis)
library(ggpubr)

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
#cna_multicna_ae        <- set_truth_rownames(coerce_123(cna_multicna_ae),        truth_mat)
cna_multicna_single <- set_truth_rownames(coerce_123(cna_multicna_single), truth_mat)
cna_hmmcopy         <- set_truth_rownames(coerce_123(cna_hmmcopy),         truth_mat)
cna_scope           <- set_truth_rownames(coerce_123(cna_scope),           truth_mat)
cna_sconce          <- set_truth_rownames(coerce_123(cna_sconce),          truth_mat)
cna_sccnv           <- set_truth_rownames(coerce_123(cna_sccnv),           truth_mat)
cna_rccae           <- set_truth_rownames(coerce_123(cna_rccae),           truth_mat)
cna_scicone           <- set_truth_rownames(coerce_123(cna_scicone),           truth_mat)

cna_multicna <- coerce_123(cna_multicna)
cna_multicna_single <- coerce_123(cna_multicna_single)
cna_hmmcopy <- coerce_123(cna_hmmcopy)
cna_scope    = coerce_123(cna_scope)
cna_sconce    = coerce_123(cna_sconce)
cna_sccnv    = coerce_123(cna_sccnv)
cna_rccae    = coerce_123(cna_rccae)
cna_scicone    = coerce_123(cna_scicone)

################################################### F1 by length#################################################################

library(data.table);library(IRanges);library(ggplot2); library(stringr); library(dplyr);
library(tidyverse);
library(ggh4x)
#---------------------------
# Helper: length stratum from bp length
#---------------------------
len_stratum <- function(len_bp) {
  fifelse(len_bp < 3e6, "<3Mb",
          fifelse(len_bp <= 10e6, "3-10Mb", ">10Mb"))
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
  
  strata_levels <- c("<3Mb", "3-10Mb", ">10Mb")
  
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

#===========================
# F1-score
#===========================
#BASE   <- "c2_N50_100kb_l1e4"
#CNA_DIR <- file.path("<MINTCNA_ROOT>/Simulation/cnasim", BASE)  
# load(file.path(CNA_DIR,"Calling/MintCNA_h2.RData"))
# load(file.path(CNA_DIR,"Calling/MintCNA_h2_single.RData"))
#load(file.path(CNA_DIR,"Calling/MintCNA_h10_AE.RData"))
#cna_multicna        <- set_truth_rownames(coerce_123(cna_multicna),        truth_mat)
#cna_multicna_ae        <- set_truth_rownames(coerce_123(cna_multicna_ae),        truth_mat)
#cna_multicna_single <- set_truth_rownames(coerce_123(cna_multicna_single), truth_mat)

# cell_counts is what you asked for:
# columns: cell, state (del/dup), stratum (<3Mb / 3-10Mb / >10Mb), TP, FP, FN, precision, recall, F1
# truth segments once
methods <- list(
  hmmcopy  = cna_hmmcopy,
  mintcna = cna_multicna,
  mintcna_single = cna_multicna_single,
  # mintcna_ae = cna_multicna_ae,
  scope    = cna_scope,
  sconce    = cna_sconce,
  sccnv    = cna_sccnv,
  rccae = cna_rccae
  ,scicone = cna_scicone
)

name_map     <- c(hmmcopy="HMMcopy", mintcna="MintCNA",scope="SCOPE", 
                  sconce ="SCONCE",sccnv="SCCNV", rccae="rcCAE",
                  mintcna_single="MintCNA(DNA)", scicone="SCICoNE"
                  #mintcna_ae = "MintCNA(AE)"
)

method_order <- c("MintCNA",
                  "MintCNA(DNA)",
                  #"MintCNA(AE)",
                  "SCOPE","HMMcopy","SCONCE","SCCNV","rcCAE","SCICoNE")
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
  
  cnt <- cnvwise_counts(truth_seg, pred_seg, overlap_frac = 0.2)
  cnt[, method := name_map[m]]
  cnt
}), use.names = TRUE, fill = TRUE)


# ordering for plots/tables
all_counts[, method := factor(method, levels = method_order)]
all_counts[, state  := factor(state, levels = c("del","dup"))]
all_counts[, stratum := factor(stratum, levels = c("<3Mb","3-10Mb",">10Mb"))]

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
summary_mss[, stratum := factor(stratum, levels = c("<3Mb","3-10Mb",">10Mb"))]

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

# -----------------------
# aesthetics (simple + robust)
# -----------------------
# If you already have color_map/shape_map, reuse them.
# Otherwise, keep shapes uniform and let color distinguish methods.
shape_map <- setNames(rep(21, length(method_order)), method_order)  # filled circles
shape_map[1:2] <- 24

color_map <- c(
  "MintCNA"      = "#1B9E77",
  "MintCNA(DNA)" = "#D95F02",
  "SCOPE"        = "#7570B3",
  "HMMcopy"      = "#E7298A",
  "SCONCE"       = "#4DBBD5",  # changed (was #66A61E)
  "SCCNV"        = "#E6AB02",
  "rcCAE"        = "#A6761D",
  "SCICoNE"      = "#7A7A7A"
)

p_pr <- ggplot(summary_mss, aes(x = precision, y = recall)) +
  geom_contour(data = grid_data,
               aes(x = precision, y = recall, z = f1),
               breaks = f1_levels, color = "grey80", linewidth = 0.3,
               inherit.aes = FALSE) +
  geom_text(data = label_positions,
            aes(x = precision, y = recall, label = label),
            color = "grey60", size = 3, vjust = 0, hjust = 1,
            inherit.aes = FALSE) +
  geom_point(aes(color = method, fill = method, shape = method), size = 4) +
  scale_shape_manual(values = shape_map) +
  scale_color_manual(values = color_map) +
  scale_fill_manual(values = color_map) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) +
  ggh4x::facet_grid2(state ~ stratum, axes = "all") +   # <-- key line
  labs(x = "Precision", y = "Recall") +
  theme_classic(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "white", color = NA),
    legend.position = "none"
  ) +
  guides(color = guide_legend(nrow = 2),
         fill  = guide_legend(nrow = 2),
         shape = guide_legend(nrow = 2))

print(p_pr)
ggsave(file.path(CNA_DIR,"Result/F1_CNV.pdf"), plot = p_pr, width = 12, height = 6, units = "in", dpi = 600)

################################################### F1 totall (no length)################################################################

cnvwise_counts_nostrata <- function(truth_seg, pred_seg, overlap_frac = 0.20) {
  truth_seg <- copy(truth_seg)
  pred_seg  <- copy(pred_seg)

  # label states like your old code
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
        ov_len   <- pmax(0, ov_end - ov_start)  # keep same convention as your code

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

truth_seg <- summarize_cna_all_cells(truth_mat, baseline_cn = 2, keep_neutral = FALSE)

all_counts <- rbindlist(lapply(names(methods), function(m) {
  message(m)
  pred_mat <- coerce_123(methods[[m]])
  pred_mat <- align_to_truth_by_colname_fillNA(pred_mat, truth_mat)

  pred_seg <- summarize_cna_all_cells(pred_mat, baseline_cn = 2, keep_neutral = FALSE)

  cnt <- cnvwise_counts_nostrata(truth_seg, pred_seg, overlap_frac = 0.2)
  cnt[, method := name_map[m]]
  cnt
}), use.names = TRUE, fill = TRUE)

all_counts[, method := factor(method, levels = method_order)]
all_counts[, state  := factor(state, levels = c("del","dup"))]

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

p_pr_one <- ggplot(summary_one, aes(x = precision, y = recall)) +
  geom_contour(data = grid_data,
               aes(x = precision, y = recall, z = f1),
               breaks = f1_levels, color = "grey80", linewidth = 0.3,
               inherit.aes = FALSE) +
  geom_text(data = label_positions,
            aes(x = precision, y = recall, label = label),
            color = "grey60", size = 3, vjust = 0, hjust = 1,
            inherit.aes = FALSE) +
  geom_point(aes(color = method, fill = method, shape = method), size = 2.8) +
  scale_shape_manual(values = shape_map) +
  scale_color_manual(values = color_map) +
  scale_fill_manual(values = color_map) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) +
  labs(x = "Precision", y = "Recall") +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal"
  ) +
  guides(color = guide_legend(nrow = 2),
         fill  = guide_legend(nrow = 2),
         shape = guide_legend(nrow = 2))

print(p_pr_one)
ggsave(file.path(CNA_DIR, "Result/PR_overall_onepoint1.pdf"),
       plot = p_pr_one, width = 6.5, height = 5.2, units = "in", dpi = 600)

## =========================================================
## Plot B: split by DEL vs DUP (still no length stratification)
##         -> two panels, each method one point per panel
## =========================================================
p_pr_state <- ggplot(summary_all, aes(x = precision, y = recall)) +
  geom_contour(data = grid_data,
               aes(x = precision, y = recall, z = f1),
               breaks = f1_levels, color = "grey80", linewidth = 0.3,
               inherit.aes = FALSE) +
  geom_text(data = label_positions,
            aes(x = precision, y = recall, label = label),
            color = "grey60", size = 3, vjust = 0, hjust = 1,
            inherit.aes = FALSE) +
  geom_point(aes(color = method, fill = method, shape = method), size = 2.6) +
  scale_shape_manual(values = shape_map) +
  scale_color_manual(values = color_map) +
  scale_fill_manual(values = color_map) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) +
  ggh4x::facet_wrap2(~ state, nrow = 1, axes = "all") +   # <-- key line
  labs(x = "Precision", y = "Recall") +
  theme_classic(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "white", color = NA),
    legend.position = "bottom",
    legend.direction = "horizontal"
  ) +
  guides(color = guide_legend(nrow = 2),
         fill  = guide_legend(nrow = 2),
         shape = guide_legend(nrow = 2))

print(p_pr_state)
ggsave(file.path(CNA_DIR, "Result/PR_overall_by_state.pdf"),
       plot = p_pr_state, width = 10.5, height = 5.2, units = "in", dpi = 600)

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

ggsave(file.path(CNA_DIR, "Result/PR_overall_dup.pdf"),
       plot = p_pr_dup, width = 5.2, height = 5.2, units = "in", dpi = 600)

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

ggsave(file.path(CNA_DIR, "Result/PR_overall_del.pdf"),
       plot = p_pr_del, width = 5.2, height = 5.2, units = "in", dpi = 600)

# desired order
leg_df <- data.frame(
  method_disp = factor(method_order, levels = method_order),
  y = 1
)

# solid shapes (no outline)
shape_map <- setNames(rep(19, length(method_order)), method_order)  # solid circle
shape_map[method_order[1:2]] <- 17                                  # solid triangle

p_leg <- ggplot(leg_df, aes(x = method_disp, y = y,
                            colour = method_disp, shape = method_disp)) +
  geom_point(size = 5) +
  scale_colour_manual(name = "Method:", values = color_map, breaks = method_order) +
  scale_shape_manual(values = shape_map, breaks = method_order) +
  theme_void(base_size = 11) +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 12),
    legend.key.size = unit(0.45, "cm")
  ) +
  guides(
    colour = guide_legend(
      nrow = 1, byrow = TRUE,
      override.aes = list(
        shape = unname(shape_map[method_order]),
        size  = 5
      )
    ),
    shape = "none"
  )

leg <- cowplot::get_legend(p_leg)
leg_plot <- ggpubr::as_ggplot(leg)
ggsave(file.path(CNA_DIR,"Result/legend_method_F1.pdf"), leg_plot, width = 10, height = 0.8)

########################### RMSE #######################################

library(data.table)
library(ggplot2)

# per-cell RMSE across bins
rmse_per_cell <- function(truth_mat, pred_mat) {
  stopifnot(nrow(truth_mat) == nrow(pred_mat))
  stopifnot(!is.null(colnames(truth_mat)), !is.null(colnames(pred_mat)))

  keep_cells <- intersect(colnames(truth_mat), colnames(pred_mat))
  if (length(keep_cells) == 0L) stop("No overlapping cell names.")

  Tsub <- truth_mat[, keep_cells, drop = FALSE]
  Psub <- pred_mat[,  keep_cells, drop = FALSE]

  rmse_i <- sqrt(colMeans((Tsub - Psub)^2, na.rm = TRUE))
  data.table(cell = keep_cells, rmse = as.numeric(rmse_i), n_cells_used = length(keep_cells))
}

# compute for all methods
rmse_tab <- rbindlist(lapply(names(methods), function(m) {
  dt <- rmse_per_cell(truth_mat, methods[[m]])
  dt[, method := m]
  dt
}))

# map names/order
rmse_tab[, method_disp := name_map[method]]
rmse_tab[, method_disp := factor(method_disp, levels = method_order)]

# boxplot
p_rmse <- ggplot(rmse_tab, aes(x = method_disp, y = rmse, fill = method_disp)) +
  geom_boxplot(width = 0.6, outlier.size = 0.35, linewidth = 0.3) +
  scale_fill_manual(values = color_map, guide = "none") +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, size = 11),
    axis.text.y = element_text(size = 12),
    axis.ticks.y = element_line(linewidth = 0.4)
  ) +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  coord_cartesian(ylim = c(0, 1))   # better than ylim() (doesn't drop points)

print(p_rmse)
ggsave(file.path(CNA_DIR,"Result/RMSE.pdf"), plot = p_rmse, width = 6.5, height = 6, units = "in", dpi = 600)

########################### Boundary Bias #######################################

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
# p_bb <- ggplot(bb_tab[!is.na(BB_cell_mb)],
#                aes(x = method_disp, y = BB_cell_mb, fill = method_disp)) +
#   geom_boxplot(width = 0.6, outlier.size = 0.4) +
#   scale_fill_manual(values = color_map) +
#   labs(x = NULL, y = "Mean Boundary Bias(Mb)") +
#   theme_bw(base_size = 11) +
#   theme(legend.position = "none",
#         axis.text.x = element_text(angle = 35, hjust = 1))+
#   ylim(0,1.5)
# print(p_bb)
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
    axis.ticks.x = element_blank())+
  coord_cartesian(ylim = c(0, 20))
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
