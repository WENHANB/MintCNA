args <- commandArgs(trailingOnly = TRUE)
BASE <- if (length(args) >= 1) args[1] else "c3_N100_100kb_l1e4"
print(BASE)
CNA_DIR <- file.path("<MINTCNA_ROOT>/Simulation/cnasim", BASE)

OUT_DIR <- file.path("<MINTCNA_ROOT>/Simulation/results/tables", BASE)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

###############################################################
# Load method results
###############################################################
load(file.path(CNA_DIR,"Calling/MintCNA_h5.RData"))
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
ord <- order(cell_num, colnames(cna_sconce), na.last = TRUE)
cna_sconce <- cna_sconce[, ord, drop = FALSE]

cna_sccnv <- cna_sccnv[,colnames(cna_sccnv) %in% colnames(cna_multicna)]
cna_rccae <- cna_rccae[,colnames(cna_rccae) %in% colnames(cna_multicna)]
cna_scicone <- cna_scicone[,colnames(cna_scicone) %in% colnames(cna_multicna)]

library(data.table)
library(IRanges)
library(stringr)
library(dplyr)

###############################################################
# Helper functions
###############################################################
coerce_123 <- function(mat) {
  v <- as.integer(mat)
  if (all(v %in% c(1L, 2L, 3L, NA_integer_))) {
    return(matrix(v, nrow=nrow(mat), ncol=ncol(mat), dimnames=dimnames(mat)))
  }
  out <- rep(NA_integer_, length(v))
  out[v %in% c(1L, 2L)] <- 1L
  out[v == 3L]          <- 2L
  out[v >= 4L]          <- 3L
  matrix(out, nrow=nrow(mat), ncol=ncol(mat), dimnames=dimnames(mat))
}

set_truth_rownames <- function(mat, truth_mat) {
  stopifnot(nrow(mat) == nrow(truth_mat))
  rownames(mat) <- rownames(truth_mat)
  mat
}

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

len_stratum <- function(len_bp) {
  fifelse(len_bp < 3e6, "<3Mb",
          fifelse(len_bp <= 10e6, "3-10Mb", ">10Mb"))
}

.parse_bins <- function(bin_ids) {
  m <- str_match(bin_ids, "^(chr[^:]+):(\\d+)-(\\d+)$")
  if (any(is.na(m[,1]))) stop("Bad bin_id: ", bin_ids[which(is.na(m[,1]))[1]])
  data.table(
    idx    = seq_along(bin_ids),
    bin_id = bin_ids,
    chr    = m[,2],
    start  = as.integer(m[,3]),
    end    = as.integer(m[,4])
  )
}

smooth_cn_by_runs <- function(dt, baseline_cn = 2L,
                              max_short_bins = 1L,
                              min_neighbor_bins = 5L) {
  setorder(dt, chr, idx)
  dt[, run_id := cumsum(is.na(data.table::shift(chr)) | chr != data.table::shift(chr) | cn != data.table::shift(cn))]
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
  setorder(dt, idx)
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

cnvwise_counts <- function(truth_seg, pred_seg, overlap_frac = 0.20) {
  truth_seg <- copy(truth_seg)
  pred_seg  <- copy(pred_seg)
  truth_seg[, stratum := len_stratum(len_bp)]
  pred_seg[,  stratum := len_stratum(len_bp)]
  truth_seg[, state_lab := fifelse(state == 1L, "del", "dup")]
  pred_seg[,  state_lab := fifelse(state == 1L, "del", "dup")]
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
    if (nrow(T0) == 0L) T0 <- T0[0]
    if (nrow(P0) == 0L) P0 <- P0[0]
    cell_state_rows <- vector("list", length(strata_levels))
    for (s in seq_along(strata_levels)) {
      str_s <- strata_levels[s]
      T <- T0[stratum == str_s]
      P <- P0[stratum == str_s]
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
          hits <- findOverlaps(IRanges(Tcc$start, Tcc$end), IRanges(Pcc$start, Pcc$end))
          if (length(hits) == 0L) return(NULL)
          ti <- queryHits(hits)
          pj <- subjectHits(hits)
          ov_start <- pmax(Tcc$start[ti], Pcc$start[pj])
          ov_end   <- pmin(Tcc$end[ti],   Pcc$end[pj])
          ov_len   <- pmax(0, ov_end - ov_start)
          data.table(tid = Tcc$tid[ti], pid = Pcc$pid[pj], ov_len = ov_len,
                     truth_len = Tcc$len_bp[ti])[, prop := ov_len / truth_len][prop > 0]
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
      cell_state_rows[[s]] <- data.table(cell = cell_k, state = st_k, stratum = str_s,
                                         TP = TP, FP = FP, FN = FN)
    }
    res_list[[k]] <- rbindlist(cell_state_rows)
  }
  res <- rbindlist(res_list, use.names = TRUE, fill = TRUE)
  res[, precision := fifelse((TP + FP) == 0, NA_real_, TP / (TP + FP))]
  res[, recall    := fifelse((TP + FN) == 0, NA_real_, TP / (TP + FN))]
  res[, F1 := fifelse(is.na(precision) | is.na(recall) | (precision + recall) == 0,
                      NA_real_, 2 * precision * recall / (precision + recall))]
  res[]
}

cnvwise_counts_nostrata <- function(truth_seg, pred_seg, overlap_frac = 0.20) {
  truth_seg <- copy(truth_seg)
  pred_seg  <- copy(pred_seg)
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
        hits <- findOverlaps(IRanges(Tcc$start, Tcc$end), IRanges(Pcc$start, Pcc$end))
        if (length(hits) == 0L) return(NULL)
        ti <- queryHits(hits)
        pj <- subjectHits(hits)
        ov_start <- pmax(Tcc$start[ti], Pcc$start[pj])
        ov_end   <- pmin(Tcc$end[ti],   Pcc$end[pj])
        ov_len   <- pmax(0, ov_end - ov_start)
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

match_tp_and_bb_one_cell_state <- function(T, P, overlap_frac = 0.20) {
  if (nrow(T) == 0L || nrow(P) == 0L) return(data.table())
  T <- copy(T); P <- copy(P)
  T[, tid := .I]; P[, pid := .I]
  cand <- rbindlist(lapply(intersect(unique(T$chr), unique(P$chr)), function(cc) {
    Tcc <- T[chr == cc]
    Pcc <- P[chr == cc]
    if (nrow(Tcc) == 0L || nrow(Pcc) == 0L) return(NULL)
    hits <- findOverlaps(IRanges(Tcc$start, Tcc$end), IRanges(Pcc$start, Pcc$end))
    if (length(hits) == 0L) return(NULL)
    ti <- queryHits(hits)
    pj <- subjectHits(hits)
    ov_start <- pmax(Tcc$start[ti], Pcc$start[pj])
    ov_end   <- pmin(Tcc$end[ti],   Pcc$end[pj])
    ov_len <- pmax(0, ov_end - ov_start)
    truth_len <- pmax(1, Tcc$len_bp[ti])
    prop <- ov_len / truth_len
    data.table(tid = Tcc$tid[ti], pid = Pcc$pid[pj], chr = cc,
               t_start = Tcc$start[ti], t_end = Tcc$end[ti],
               p_start = Pcc$start[pj], p_end = Pcc$end[pj],
               ov_len = ov_len, prop = prop)[prop >= overlap_frac]
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

###############################################################
# Load ground truth
###############################################################
profiles <- fread(file.path(CNA_DIR, "profiles.tsv"))
colnames(profiles)[5] <- "cn_states"
setnames(profiles, tolower(names(profiles)))

profiles[, cn_total := {
  x <- trimws(as.character(cn_states))
  ifelse(grepl("\\.", x),
         rowSums(matrix(as.numeric(unlist(tstrsplit(x, "\\.", perl=TRUE))), ncol=2, byrow=FALSE)),
         round(as.numeric(x))) |> as.integer()
}]

profiles[!grepl("^chr", chrom), chrom := paste0("chr", chrom)]
profiles[, bin_id := paste0(chrom, ":", start, "-", end)]
profiles[, cell_num := as.integer(sub("^cell", "", cell))]
profiles[, chr_num  := as.integer(sub("^chr",  "", chrom))]
profiles <- profiles[grepl("^cell\\d+$", cell)]
setorder(profiles, chr_num, start, end, cell_num)

cells <- profiles[, unique(.SD[order(cell_num)]$cell), .SDcols = c("cell","cell_num")]
CNdt  <- dcast(profiles, bin_id ~ cell, value.var = "cn_total")
bin_ids <- CNdt$bin_id
CNmat <- as.matrix(CNdt[, ..cells])
rownames(CNmat) <- bin_ids

ord_bin <- order(
  as.integer(sub("^chr", "", sub(":.*$", "", rownames(CNmat)))),
  as.integer(sub("^.*:(\\d+)-.*$", "\\1", rownames(CNmat)))
)
CNmat <- CNmat[ord_bin, , drop=FALSE]

cn_bins <- rownames(CNmat)
chr  <- sub("^chr", "", sub(":.*$", "", cn_bins))
pos  <- sub("^.*:", "", cn_bins)
st   <- as.integer(sub("-.*$", "", pos))
en   <- as.integer(sub("^.*-", "", pos))
key_CN <- paste(chr, st, en, sep = "_")
names(key_CN) <- cn_bins

key_need <- rownames(cna_multicna)
idx <- match(key_need, key_CN)
if (anyNA(idx)) {
  missing <- key_need[is.na(idx)]
  stop("Some bins in cna_multicna not found in CNmat. Example: ", missing[1])
}

label <- fread(file.path(CNA_DIR,"cell_types.tsv"))
nc_samples <- label[which(label$normal == "normal"), root]
nor_idx <- which(colnames(CNmat) %in% nc_samples)

CNmat_sub <- CNmat[idx, , drop = FALSE]
CNmat_sub <- CNmat_sub[,-nor_idx]

###############################################################
# Coerce all matrices to 1/2/3 states
###############################################################
truth_mat <- coerce_123(CNmat_sub+1)

cna_multicna        <- set_truth_rownames(coerce_123(cna_multicna),        truth_mat)
cna_multicna_single <- set_truth_rownames(coerce_123(cna_multicna_single), truth_mat)
cna_hmmcopy         <- set_truth_rownames(coerce_123(cna_hmmcopy),         truth_mat)
cna_scope           <- set_truth_rownames(coerce_123(cna_scope),           truth_mat)
cna_sconce          <- set_truth_rownames(coerce_123(cna_sconce),          truth_mat)
cna_sccnv           <- set_truth_rownames(coerce_123(cna_sccnv),          truth_mat)
cna_rccae           <- set_truth_rownames(coerce_123(cna_rccae),          truth_mat)
cna_scicone         <- set_truth_rownames(coerce_123(cna_scicone),        truth_mat)

cna_multicna        <- coerce_123(cna_multicna)
cna_multicna_single <- coerce_123(cna_multicna_single)
cna_hmmcopy         <- coerce_123(cna_hmmcopy)
cna_scope           <- coerce_123(cna_scope)
cna_sconce          <- coerce_123(cna_sconce)
cna_sccnv           <- coerce_123(cna_sccnv)
cna_rccae           <- coerce_123(cna_rccae)
cna_scicone         <- coerce_123(cna_scicone)

###############################################################
# Method list and name mapping
###############################################################
methods <- list(
  hmmcopy  = cna_hmmcopy,
  mintcna = cna_multicna,
  mintcna_single = cna_multicna_single,
  scope    = cna_scope,
  sconce   = cna_sconce,
  sccnv    = cna_sccnv,
  rccae    = cna_rccae,
  scicone  = cna_scicone
)

name_map <- c(hmmcopy="HMMcopy", mintcna="MintCNA", scope="SCOPE",
              sconce="SCONCE", sccnv="SCCNV", rccae="rcCAE",
              mintcna_single="MintCNA(DNA)", scicone="SCICoNE")

method_order <- c("MintCNA", "MintCNA(DNA)", "SCOPE", "HMMcopy",
                  "SCONCE", "SCCNV", "rcCAE", "SCICoNE")

###############################################################
# Table 1: summary_mss (F1 by method x state x length stratum)
###############################################################
message("Computing summary_mss ...")
truth_seg <- summarize_cna_all_cells(truth_mat, baseline_cn = 2, keep_neutral = FALSE)

all_counts_strata <- rbindlist(lapply(names(methods), function(m) {
  message(m)
  pred_mat <- coerce_123(methods[[m]])
  pred_mat <- align_to_truth_by_colname_fillNA(pred_mat, truth_mat)
  pred_seg <- summarize_cna_all_cells(pred_mat, baseline_cn = 2, keep_neutral = FALSE)
  cnt <- cnvwise_counts(truth_seg, pred_seg, overlap_frac = 0.2)
  cnt[, method := name_map[m]]
  cnt
}), use.names = TRUE, fill = TRUE)

all_counts_strata[, method  := factor(method, levels = method_order)]
all_counts_strata[, state   := factor(state, levels = c("del","dup"))]
all_counts_strata[, stratum := factor(stratum, levels = c("<3Mb","3-10Mb",">10Mb"))]

summary_mss <- all_counts_strata[, .(
  TP = sum(TP, na.rm = TRUE),
  FP = sum(FP, na.rm = TRUE),
  FN = sum(FN, na.rm = TRUE),
  n_cells = uniqueN(cell)
), by = .(method, state, stratum)]

summary_mss[, precision := fifelse((TP + FP) == 0, NA_real_, TP/(TP + FP))]
summary_mss[, recall    := fifelse((TP + FN) == 0, NA_real_, TP/(TP + FN))]
summary_mss[, F1        := fifelse(is.na(precision) | is.na(recall) | (precision + recall) == 0,
                                   NA_real_, 2*precision*recall/(precision + recall))]
summary_mss[, method  := factor(method, levels = method_order)]
summary_mss[, state   := factor(state, levels = c("del","dup"))]
summary_mss[, stratum := factor(stratum, levels = c("<3Mb","3-10Mb",">10Mb"))]

###############################################################
# Tables 2-4: summary_del, summary_dup, summary_one (no strata)
###############################################################
message("Computing summary_del, summary_dup, summary_one ...")
truth_seg <- summarize_cna_all_cells(truth_mat, baseline_cn = 2, keep_neutral = FALSE)

all_counts_nostrata <- rbindlist(lapply(names(methods), function(m) {
  message(m)
  pred_mat <- coerce_123(methods[[m]])
  pred_mat <- align_to_truth_by_colname_fillNA(pred_mat, truth_mat)
  pred_seg <- summarize_cna_all_cells(pred_mat, baseline_cn = 2, keep_neutral = FALSE)
  cnt <- cnvwise_counts_nostrata(truth_seg, pred_seg, overlap_frac = 0.2)
  cnt[, method := name_map[m]]
  cnt
}), use.names = TRUE, fill = TRUE)

all_counts_nostrata[, method := factor(method, levels = method_order)]
all_counts_nostrata[, state  := factor(state, levels = c("del","dup"))]

summary_all <- all_counts_nostrata[, .(
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

summary_del <- subset(summary_all, state == "del")
summary_dup <- subset(summary_all, state == "dup")

summary_one <- summary_all[, .(
  TP = sum(TP, na.rm=TRUE),
  FP = sum(FP, na.rm=TRUE),
  FN = sum(FN, na.rm=TRUE),
  n_cells = max(n_cells, na.rm=TRUE)
), by = .(method)]

summary_one[, precision := fifelse((TP + FP) == 0, NA_real_, TP/(TP + FP))]
summary_one[, recall    := fifelse((TP + FN) == 0, NA_real_, TP/(TP + FN))]
summary_one[, F1        := fifelse(is.na(precision) | is.na(recall) | (precision + recall) == 0,
                                   NA_real_, 2*precision*recall/(precision + recall))]
summary_one[, method := factor(method, levels = method_order)]

###############################################################
# Table 5: bb_summary (boundary bias)
###############################################################
message("Computing bb_summary ...")
truth_seg <- summarize_cna_all_cells(truth_mat, baseline_cn = 2L, keep_neutral = FALSE, smooth = TRUE)

bb_tab <- rbindlist(lapply(names(methods), function(m) {
  message(m)
  pred_seg <- summarize_cna_all_cells(methods[[m]], baseline_cn = 2L, keep_neutral = FALSE, smooth = TRUE)
  dt <- bb_per_method_from_segments(truth_seg, pred_seg, overlap_frac = 0.20, summary = "mean")
  dt[, method := m]
  dt
}))

bb_tab[, method_disp := name_map[method]]
bb_tab[, method_disp := factor(method_disp, levels = method_order)]
bb_tab[, BB_cell_mb := BB_cell / 1e6]

bb_summary <- as.data.table(bb_tab)[, .(
  n_cells_used = unique(n_cells_used)[1],
  TP_total     = sum(n_TP, na.rm = TRUE),
  BB_mean_bp   = mean(BB_cell, na.rm = TRUE),
  BB_mean_mb   = mean(BB_cell_mb, na.rm = TRUE)
), by = .(method = method_disp)][order(method)]

###############################################################
# Table 7: RMSE per cell
###############################################################
message("Computing RMSE per cell ...")
rmse_tab <- rbindlist(lapply(names(methods), function(m) {
  message(m)
  dt <- rmse_per_cell(truth_mat, methods[[m]])
  dt[, method := m]
  dt
}))
rmse_tab[, method_disp := name_map[method]]
rmse_tab[, method_disp := factor(method_disp, levels = method_order)]

###############################################################
# Save all tables
###############################################################
message("Saving tables to ", OUT_DIR)
fwrite(summary_del, file.path(OUT_DIR, "summary_del.csv"))
fwrite(summary_dup, file.path(OUT_DIR, "summary_dup.csv"))
fwrite(summary_mss, file.path(OUT_DIR, "summary_mss.csv"))
fwrite(summary_one, file.path(OUT_DIR, "summary_one.csv"))
fwrite(bb_summary,  file.path(OUT_DIR, "bb_summary.csv"))
fwrite(bb_tab[, .(cell, BB_cell, BB_cell_mb, n_TP, n_cells_used, method = method_disp)],
       file.path(OUT_DIR, "bb_percell.csv"))
fwrite(rmse_tab[, .(cell, rmse, n_cells_used, method = method_disp)],
       file.path(OUT_DIR, "RMSE_percell.csv"))

message("Done: ", BASE)
