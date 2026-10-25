#!/usr/bin/env Rscript
# Accuracy of every method on one simulation setting, against the CNAsim truth:
# event-level precision / recall / F1 (an event is recovered when a called segment of the same
# direction overlaps >= 20 % of it, greedy one-to-one matching, pooled over cells), by direction
# and by event length; boundary bias of the recovered events; per-cell RMSE of the state matrix.
#   Rscript reproduce/analysis/simulation/evaluate.R <setting>
# Input : <setting>/profiles.tsv (truth), cell_types.tsv, and the call sets in <setting>/Calling/
#         (MintCNA_h5.RData, MintCNA_h5_single.RData, HMMcopy, SCCNV, SCOPE, SCONCE, rcCAE, SCICoNE).
#         To score a run of the released package instead of the deposited MintCNA calls, set
#         MINTCNA_RESULTS and MINTCNA_DNA_RESULTS to its mintcna_results.RData files (run_mintcna.sh).
# Output: Simulation/results/tables/<setting>/{summary_one,summary_del,summary_dup,summary_mss,
#         bb_summary,bb_percell,RMSE_percell}.csv  (merged over settings by merge_tables.R)
# Run from the root of the MintCNA repository with MINTCNA_DATA set (see reproduce/README.md).
suppressPackageStartupMessages({ library(data.table); library(IRanges); library(stringr); library(dplyr) })
ROOT <- Sys.getenv("MINTCNA_DATA"); if (!nzchar(ROOT)) stop("set MINTCNA_DATA (see reproduce/README.md)")
args <- commandArgs(trailingOnly = TRUE)
BASE <- if (length(args) >= 1) args[1] else "c3_N100_100kb_l1e4"
CNA_DIR <- file.path(ROOT, "Simulation/cnasim", BASE)
OUT_DIR <- file.path(ROOT, "Simulation/results/tables", BASE)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- call sets -------------------------------------------------------------------------------
# MintCNA: the deposited calls (development version, Calling/MintCNA_h5*.RData, objects cna_multicna /
# cna_multicna_single) or a run of the released package (cna_state_matrix, same 1/2/3 coding).
read_mintcna <- function(env_var, deposited, object) {
  f <- Sys.getenv(env_var)
  if (nzchar(f)) { e <- new.env(); load(f, e); return(get("cna_state_matrix", e)) }
  e <- new.env(); load(file.path(CNA_DIR, "Calling", deposited), e); get(object, e)
}
cna_multicna        <- read_mintcna("MINTCNA_RESULTS",     "MintCNA_h5.RData",        "cna_multicna")
cna_multicna_single <- read_mintcna("MINTCNA_DNA_RESULTS", "MintCNA_h5_single.RData", "cna_multicna_single")
for (f in c("HMMcopy", "SCCNV", "SCOPE", "SCONCE", "rcCAE", "SCICoNE")) load(file.path(CNA_DIR, "Calling", paste0(f, ".RData")))

tumour <- colnames(cna_multicna)
cna_hmmcopy <- cna_hmmcopy[, colnames(cna_hmmcopy) %in% tumour]; cna_hmmcopy[cna_hmmcopy >= 5] <- 5L
cna_scope   <- cna_scope[, colnames(cna_scope) %in% tumour] + 1L; cna_scope[cna_scope >= 5] <- 5L
cna_sconce  <- cna_sconce[, colnames(cna_sconce) %in% tumour, drop = FALSE]
cna_sconce  <- cna_sconce[, order(suppressWarnings(as.integer(gsub(".*?(\\d+).*", "\\1", colnames(cna_sconce)))),
                                  colnames(cna_sconce), na.last = TRUE), drop = FALSE]
cna_sccnv   <- cna_sccnv[, colnames(cna_sccnv) %in% tumour]
cna_rccae   <- cna_rccae[, colnames(cna_rccae) %in% tumour]
cna_scicone <- cna_scicone[, colnames(cna_scicone) %in% tumour]

# ---- helpers ---------------------------------------------------------------------------------
# collapse any state coding to 1 (loss) / 2 (neutral) / 3 (gain)
coerce_123 <- function(mat) {
  v <- as.integer(mat)
  if (all(v %in% c(1L, 2L, 3L, NA_integer_))) return(matrix(v, nrow(mat), ncol(mat), dimnames = dimnames(mat)))
  out <- rep(NA_integer_, length(v))
  out[v %in% c(1L, 2L)] <- 1L; out[v == 3L] <- 2L; out[v >= 4L] <- 3L
  matrix(out, nrow(mat), ncol(mat), dimnames = dimnames(mat))
}
set_truth_rownames <- function(mat, truth_mat) { stopifnot(nrow(mat) == nrow(truth_mat)); rownames(mat) <- rownames(truth_mat); mat }
align_to_truth_by_colname_fillNA <- function(mat, truth_mat) {
  stopifnot(nrow(mat) == nrow(truth_mat)); rownames(mat) <- rownames(truth_mat)
  out <- matrix(NA_integer_, nrow(truth_mat), ncol(truth_mat), dimnames = dimnames(truth_mat))
  common <- intersect(colnames(truth_mat), colnames(mat))
  if (length(common) == 0L) stop("no overlapping cells between a call set and the truth")
  out[, common] <- mat[, common, drop = FALSE]
  out
}
len_stratum <- function(len_bp) fifelse(len_bp < 3e6, "<3Mb", fifelse(len_bp <= 10e6, "3-10Mb", ">10Mb"))
.parse_bins <- function(bin_ids) {
  m <- str_match(bin_ids, "^(chr[^:]+):(\\d+)-(\\d+)$")
  if (any(is.na(m[, 1]))) stop("bad bin id: ", bin_ids[which(is.na(m[, 1]))[1]])
  data.table(idx = seq_along(bin_ids), bin_id = bin_ids, chr = m[, 2], start = as.integer(m[, 3]), end = as.integer(m[, 4]))
}

# single-bin runs sandwiched between long runs of one altered state are absorbed into their neighbours
smooth_cn_by_runs <- function(dt, baseline_cn = 2L, max_short_bins = 1L, min_neighbor_bins = 5L) {
  setorder(dt, chr, idx)
  dt[, run_id := cumsum(is.na(data.table::shift(chr)) | chr != data.table::shift(chr) | cn != data.table::shift(cn))]
  runs <- dt[, .(chr = dplyr::first(chr), n = .N, cn = dplyr::first(cn)), by = run_id]
  runs[, ev := sign(cn - baseline_cn)]
  runs[, `:=`(prev_cn = data.table::shift(cn), next_cn = data.table::shift(cn, type = "lead"),
              prev_n = data.table::shift(n), next_n = data.table::shift(n, type = "lead"),
              prev_ev = data.table::shift(ev), next_ev = data.table::shift(ev, type = "lead"),
              prev_chr = data.table::shift(chr), next_chr = data.table::shift(chr, type = "lead"))]
  runs[, new_cn := cn]
  runs[n <= max_short_bins & !is.na(prev_ev) & !is.na(next_ev) & chr == prev_chr & chr == next_chr &
         prev_ev != 0L & prev_ev == next_ev & pmax(prev_n, 0L) >= min_neighbor_bins & pmax(next_n, 0L) >= min_neighbor_bins,
       new_cn := fifelse(prev_n >= next_n, prev_cn, next_cn)]
  runs[n <= max_short_bins & !is.na(next_ev) & chr == next_chr & next_ev != 0L & ev == next_ev & next_n >= min_neighbor_bins,
       new_cn := next_cn]
  runs[n <= max_short_bins & !is.na(prev_ev) & chr == prev_chr & prev_ev != 0L & ev == prev_ev & prev_n >= min_neighbor_bins,
       new_cn := prev_cn]
  dt[runs[, .(run_id, new_cn)], on = "run_id", cn := i.new_cn]
  dt[, run_id := NULL]
  dt[]
}

# state vector of one cell -> table of altered segments (events)
summarize_cna_one_cell <- function(cn_vec, bins_dt, baseline_cn = 2L, keep_neutral = FALSE) {
  dt <- cbind(bins_dt, cn = as.integer(round(cn_vec)))
  setorder(dt, idx)
  dt <- smooth_cn_by_runs(dt, baseline_cn = baseline_cn)
  new_seg <- is.na(data.table::shift(dt$chr)) | dt$chr != data.table::shift(dt$chr) | dt$cn != data.table::shift(dt$cn)
  dt[, seg_id := cumsum(new_seg)]
  seg <- dt[, .(chr = dplyr::first(chr), start = dplyr::first(start), end = dplyr::last(end), n_bins = .N, state = dplyr::first(cn),
                len_bp = dplyr::last(end) - dplyr::first(start)), by = seg_id][, seg_id := NULL]
  seg[, event := fifelse(state < baseline_cn, "DEL", fifelse(state > baseline_cn, "DUP", "NEU"))]
  if (!keep_neutral) seg <- seg[event != "NEU"]
  seg[]
}
summarize_cna_all_cells <- function(CNmat_sub, baseline_cn = 2L, keep_neutral = FALSE) {
  bins_dt <- .parse_bins(rownames(CNmat_sub))
  res <- rbindlist(lapply(seq_len(ncol(CNmat_sub)), function(j) {
    seg <- summarize_cna_one_cell(CNmat_sub[, j], bins_dt, baseline_cn = baseline_cn, keep_neutral = keep_neutral)
    if (nrow(seg) == 0) return(NULL)
    seg[, cell := colnames(CNmat_sub)[j]]
    setcolorder(seg, c("cell", "chr", "start", "end", "len_bp", "n_bins", "state", "event"))
    seg
  }), use.names = TRUE, fill = TRUE)
  if (is.null(res)) return(data.table())
  res[]
}

# greedy one-to-one matching of truth and called events of one cell and direction:
# candidate pairs overlap >= overlap_frac of the truth event; best (overlap fraction, overlap length) first,
# ties broken by boundary bias when the matches are used to measure it
match_events <- function(T, P, overlap_frac, by_bb = FALSE) {
  T[, tid := .I]; P[, pid := .I]
  cand <- rbindlist(lapply(intersect(unique(T$chr), unique(P$chr)), function(cc) {
    Tcc <- T[chr == cc]; Pcc <- P[chr == cc]
    if (nrow(Tcc) == 0L || nrow(Pcc) == 0L) return(NULL)
    hits <- findOverlaps(IRanges(Tcc$start, Tcc$end), IRanges(Pcc$start, Pcc$end))
    if (length(hits) == 0L) return(NULL)
    ti <- queryHits(hits); pj <- subjectHits(hits)
    ov_len <- pmax(0, pmin(Tcc$end[ti], Pcc$end[pj]) - pmax(Tcc$start[ti], Pcc$start[pj]))
    data.table(tid = Tcc$tid[ti], pid = Pcc$pid[pj], ov_len = ov_len, chr = cc,
               t_start = Tcc$start[ti], t_end = Tcc$end[ti], p_start = Pcc$start[pj], p_end = Pcc$end[pj],
               prop = ov_len / pmax(1, Tcc$len_bp[ti]))[prop > 0]
  }), use.names = TRUE, fill = TRUE)
  if (is.null(cand) || nrow(cand) == 0L) return(data.table())
  cand <- cand[prop >= overlap_frac]
  if (nrow(cand) == 0L) return(data.table())
  cand[, BB := abs(t_start - p_start) + abs(t_end - p_end)]
  if (by_bb) setorder(cand, -prop, -ov_len, BB) else setorder(cand, -prop, -ov_len)
  used_t <- rep(FALSE, nrow(T)); used_p <- rep(FALSE, nrow(P)); keep <- logical(nrow(cand))
  for (i in seq_len(nrow(cand))) {
    tt <- cand$tid[i]; pp <- cand$pid[i]
    if (!used_t[tt] && !used_p[pp]) { used_t[tt] <- TRUE; used_p[pp] <- TRUE; keep[i] <- TRUE }
  }
  cand[keep]
}

# TP / FP / FN per cell and direction (optionally by truth-event length stratum)
cnvwise_counts <- function(truth_seg, pred_seg, overlap_frac = 0.20, strata = TRUE) {
  truth_seg <- copy(truth_seg); pred_seg <- copy(pred_seg)
  truth_seg[, state_lab := fifelse(state == 1L, "del", "dup")]
  pred_seg[,  state_lab := fifelse(state == 1L, "del", "dup")]
  strata_levels <- if (strata) c("<3Mb", "3-10Mb", ">10Mb") else "all"
  if (strata) { truth_seg[, stratum := len_stratum(len_bp)]; pred_seg[, stratum := len_stratum(len_bp)] } else {
    truth_seg[, stratum := "all"]; pred_seg[, stratum := "all"] }
  keys <- unique(rbind(truth_seg[, .(cell, state_lab)], pred_seg[, .(cell, state_lab)]))
  res <- rbindlist(lapply(seq_len(nrow(keys)), function(k) {
    rbindlist(lapply(strata_levels, function(str_s) {
      T <- truth_seg[cell == keys$cell[k] & state_lab == keys$state_lab[k] & stratum == str_s]
      P <- pred_seg [cell == keys$cell[k] & state_lab == keys$state_lab[k] & stratum == str_s]
      TP <- if (nrow(T) && nrow(P)) nrow(match_events(T, P, overlap_frac)) else 0L
      data.table(cell = keys$cell[k], state = keys$state_lab[k], stratum = str_s, TP = TP, FP = nrow(P) - TP, FN = nrow(T) - TP)
    }))
  }))
  if (!strata) res[, stratum := NULL]
  res[]
}
add_prf <- function(d) {
  d[, precision := fifelse((TP + FP) == 0, NA_real_, TP / (TP + FP))]
  d[, recall    := fifelse((TP + FN) == 0, NA_real_, TP / (TP + FN))]
  d[, F1 := fifelse(is.na(precision) | is.na(recall) | (precision + recall) == 0, NA_real_, 2 * precision * recall / (precision + recall))]
  d[]
}

# boundary bias: |start difference| + |end difference| of each recovered event, averaged per cell
bb_per_method_from_segments <- function(truth_seg, pred_seg, overlap_frac = 0.20) {
  cells <- intersect(unique(truth_seg$cell), unique(pred_seg$cell))
  if (length(cells) == 0L) stop("no overlapping cells between truth and called segments")
  res <- rbindlist(lapply(cells, function(cell_id) {
    T0 <- truth_seg[cell == cell_id]; P0 <- pred_seg[cell == cell_id]
    pairs <- rbindlist(lapply(intersect(unique(T0$event), unique(P0$event)), function(ev)
      match_events(copy(T0[event == ev]), copy(P0[event == ev]), overlap_frac, by_bb = TRUE)), use.names = TRUE, fill = TRUE)
    if (is.null(pairs) || nrow(pairs) == 0L) return(data.table(cell = cell_id, BB_cell = NA_real_, n_TP = 0L))
    data.table(cell = cell_id, BB_cell = mean(pairs$BB), n_TP = nrow(pairs))
  }))
  res[, n_cells_used := length(cells)]
  res[]
}
rmse_per_cell <- function(truth_mat, pred_mat) {
  keep_cells <- intersect(colnames(truth_mat), colnames(pred_mat))
  if (length(keep_cells) == 0L) stop("no overlapping cell names")
  rmse <- sqrt(colMeans((truth_mat[, keep_cells, drop = FALSE] - pred_mat[, keep_cells, drop = FALSE])^2, na.rm = TRUE))
  data.table(cell = keep_cells, rmse = as.numeric(rmse), n_cells_used = length(keep_cells))
}

# ---- truth: total copy number per bin from profiles.tsv, tumour cells, MintCNA's bins -----------
profiles <- fread(file.path(CNA_DIR, "profiles.tsv"))
colnames(profiles)[5] <- "cn_states"
setnames(profiles, tolower(names(profiles)))
profiles[, cn_total := { x <- trimws(as.character(cn_states))
  ifelse(grepl("\\.", x), rowSums(matrix(as.numeric(unlist(tstrsplit(x, "\\.", perl = TRUE))), ncol = 2, byrow = FALSE)),
         round(as.numeric(x))) |> as.integer() }]
profiles[!grepl("^chr", chrom), chrom := paste0("chr", chrom)]
profiles[, bin_id := paste0(chrom, ":", start, "-", end)]
profiles[, cell_num := as.integer(sub("^cell", "", cell))]
profiles[, chr_num := as.integer(sub("^chr", "", chrom))]
profiles <- profiles[grepl("^cell\\d+$", cell)]
setorder(profiles, chr_num, start, end, cell_num)
cells <- profiles[, unique(.SD[order(cell_num)]$cell), .SDcols = c("cell", "cell_num")]
CNdt <- dcast(profiles, bin_id ~ cell, value.var = "cn_total")
CNmat <- as.matrix(CNdt[, ..cells]); rownames(CNmat) <- CNdt$bin_id
CNmat <- CNmat[order(as.integer(sub("^chr", "", sub(":.*$", "", rownames(CNmat)))),
                     as.integer(sub("^.*:(\\d+)-.*$", "\\1", rownames(CNmat)))), , drop = FALSE]

cn_bins <- rownames(CNmat)
key_CN <- paste(sub("^chr", "", sub(":.*$", "", cn_bins)), as.integer(sub("-.*$", "", sub("^.*:", "", cn_bins))),
                as.integer(sub("^.*-", "", sub("^.*:", "", cn_bins))), sep = "_")
idx <- match(rownames(cna_multicna), key_CN)
if (anyNA(idx)) stop("bins of the MintCNA calls missing from the truth, e.g. ", rownames(cna_multicna)[is.na(idx)][1])
label <- fread(file.path(CNA_DIR, "cell_types.tsv"))
nc_samples <- label[normal == "normal", root]
CNmat_sub <- CNmat[idx, !colnames(CNmat) %in% nc_samples, drop = FALSE]
truth_mat <- coerce_123(CNmat_sub + 1)

methods <- list(hmmcopy = cna_hmmcopy, mintcna = cna_multicna, mintcna_single = cna_multicna_single,
                scope = cna_scope, sconce = cna_sconce, sccnv = cna_sccnv, rccae = cna_rccae, scicone = cna_scicone)
methods <- lapply(methods, function(m) set_truth_rownames(coerce_123(m), truth_mat))
name_map <- c(hmmcopy = "HMMcopy", mintcna = "MintCNA", mintcna_single = "MintCNA(DNA)", scope = "SCOPE",
              sconce = "SCONCE", sccnv = "SCCNV", rccae = "rcCAE", scicone = "SCICoNE")
method_order <- c("MintCNA", "MintCNA(DNA)", "SCOPE", "HMMcopy", "SCONCE", "SCCNV", "rcCAE", "SCICoNE")

# ---- event-level precision / recall / F1 ------------------------------------------------------
truth_seg <- summarize_cna_all_cells(truth_mat)
pred_segs <- lapply(methods, function(m) summarize_cna_all_cells(align_to_truth_by_colname_fillNA(m, truth_mat)))

counts <- function(strata) rbindlist(lapply(names(methods), function(m) {
  message(m); cnvwise_counts(truth_seg, pred_segs[[m]], overlap_frac = 0.2, strata = strata)[, method := name_map[m]]
}))
message("by length stratum ...")
all_strata <- counts(TRUE)
summary_mss <- add_prf(all_strata[, .(TP = sum(TP), FP = sum(FP), FN = sum(FN), n_cells = uniqueN(cell)), by = .(method, state, stratum)])
summary_mss[, `:=`(method = factor(method, levels = method_order), state = factor(state, levels = c("del", "dup")),
                   stratum = factor(stratum, levels = c("<3Mb", "3-10Mb", ">10Mb")))]
message("all lengths ...")
all_nostrata <- counts(FALSE)
summary_all <- add_prf(all_nostrata[, .(TP = sum(TP), FP = sum(FP), FN = sum(FN), n_cells = uniqueN(cell)), by = .(method, state)])
summary_all[, `:=`(method = factor(method, levels = method_order), state = factor(state, levels = c("del", "dup")))]
summary_del <- summary_all[state == "del"]
summary_dup <- summary_all[state == "dup"]
summary_one <- add_prf(summary_all[, .(TP = sum(TP), FP = sum(FP), FN = sum(FN), n_cells = max(n_cells)), by = method])

# ---- boundary bias and RMSE -----------------------------------------------------------------
message("boundary bias ...")
bb_tab <- rbindlist(lapply(names(methods), function(m)
  bb_per_method_from_segments(truth_seg, summarize_cna_all_cells(methods[[m]]))[, method := m]))
bb_tab[, method_disp := factor(name_map[method], levels = method_order)]
bb_tab[, BB_cell_mb := BB_cell / 1e6]
bb_summary <- bb_tab[, .(n_cells_used = unique(n_cells_used)[1], TP_total = sum(n_TP, na.rm = TRUE),
                         BB_mean_bp = mean(BB_cell, na.rm = TRUE), BB_mean_mb = mean(BB_cell_mb, na.rm = TRUE)),
                     by = .(method = method_disp)][order(method)]
message("RMSE ...")
rmse_tab <- rbindlist(lapply(names(methods), function(m) rmse_per_cell(truth_mat, methods[[m]])[, method := m]))
rmse_tab[, method_disp := factor(name_map[method], levels = method_order)]

# ---- write ------------------------------------------------------------------------------------
fwrite(summary_del, file.path(OUT_DIR, "summary_del.csv"))
fwrite(summary_dup, file.path(OUT_DIR, "summary_dup.csv"))
fwrite(summary_mss, file.path(OUT_DIR, "summary_mss.csv"))
fwrite(summary_one, file.path(OUT_DIR, "summary_one.csv"))
fwrite(bb_summary,  file.path(OUT_DIR, "bb_summary.csv"))
fwrite(bb_tab[, .(cell, BB_cell, BB_cell_mb, n_TP, n_cells_used, method = method_disp)], file.path(OUT_DIR, "bb_percell.csv"))
fwrite(rmse_tab[, .(cell, rmse, n_cells_used, method = method_disp)], file.path(OUT_DIR, "RMSE_percell.csv"))
message("wrote ", OUT_DIR)
