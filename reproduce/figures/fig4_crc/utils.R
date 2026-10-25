# Shared helpers for Figure 4 (CRC01 scTrio-seq2 benchmark): load the call set of every method and
# the bulk WGS truth, and define the method, site, state and chromosome constants of the panels.
# Input : reproduce/figures/data/fig4_crc/*.RData, one file per method (MintCNA calls in
#         MultiCNA_AE_197_d250.RData and MultiCNA_AE_197_single_d250.RData) plus DNAcopy_qc.RData,
#         the post-QC bulk WGS truth (reproduce/analysis/crc01/wgs_truth_postqc.R)
# Sourced by fig4_crc/panel_a_heatmap.R, panel_b_accuracy.R and panel_c_tracks.R, which run
# through assemble_fig4.R from the root of the MintCNA repository.

suppressPackageStartupMessages({
  library(data.table)
})

BENCH_DIR <- "reproduce/figures/data/fig4_crc"

METHOD_ORDER <- c("MintCNA", "MintCNA(DNA)", "SCOPE", "HMMcopy",
                  "SCONCE", "SCCNV", "rcCAE", "SCICoNE")

# Method colours: two greens for the MintCNA pair, six non-green hues for the competitors.
METHOD_COLORS <- c(
  "MintCNA"      = "#1B9E77",
  "MintCNA(DNA)" = "#046B44",
  "SCOPE"        = "#0072B2",
  "HMMcopy"      = "#E69F00",
  "SCONCE"       = "#984EA3",
  "SCCNV"        = "#56B4E9",
  "rcCAE"        = "#D55E00",
  "SCICoNE"      = "#CC79A7"
)

# Display labels keyed by method name, not by position, so a panel ordered by performance cannot
# mislabel a box; the data key "MintCNA(DNA)" is rendered as MintCNA with a superscript DNA.
# shared/fig_palette.R, sourced after this file by the panel scripts, defines the same names for
# all figures.
METHOD_PLOTMATH <- c(
  "MintCNA"      = "MintCNA",
  "MintCNA(DNA)" = "MintCNA^DNA",
  "SCOPE"        = "SCOPE",
  "HMMcopy"      = "HMMcopy",
  "SCONCE"       = "SCONCE",
  "SCCNV"        = "SCCNV",
  "rcCAE"        = "rcCAE",
  "SCICoNE"      = "SCICoNE"
)
method_labels <- function(ms) {
  stopifnot(all(ms %in% names(METHOD_PLOTMATH)))
  parse(text = METHOD_PLOTMATH[ms])
}
METHOD_LABELS <- method_labels(METHOD_ORDER)

# Display order for the benchmark panels. The two MintCNA modes lead, in that order, as the
# methods under study (a layout choice, not a result); the competitors follow ranked by the
# median of the metric the panel shows, with `lower_is_better` flipping the direction for error
# metrics. The one vector returned serves as factor levels, scale limits and summary-table row
# order, and two panels sharing a legend must be given the same order.
rank_methods <- function(dt, value, lower_is_better = FALSE,
                         pin = c("MintCNA", "MintCNA(DNA)")) {
  s  <- dt[, .(m = median(get(value), na.rm = TRUE)), by = method]
  s  <- s[order(if (lower_is_better) m else -m)]
  ms <- as.character(s$method)
  c(intersect(pin, ms), setdiff(ms, pin))
}

SITE_ORDER  <- c("PT", "LN", "ML", "MP")
SITE_COLORS <- c(PT = "#4C72B0", LN = "#55A868", ML = "#DD8452", MP = "#C44E52")

CNA_5_COLORS <- c(
  "1" = "#2C7BB6", "2" = "#43B5E8", "3" = "#F3ECDD",
  "4" = "#F4A582", "5" = "#B2182B"
)
CNA_5_LABELS <- c("Del", "Loss", "Nor", "Gain", "Amp")

CHR_ORDER <- paste0("chr", 1:22)
CHR17_CENTROMERE_MB <- 24

CRC_DRIVERS <- data.frame(
  symbol = c("APC", "KRAS", "TP53",  "NF1",  "SMAD4", "DCC",  "MYC",  "ERBB2", "BRAF", "PIK3CA"),
  chr    = c("chr5","chr12","chr17", "chr17","chr18", "chr18","chr8", "chr17", "chr7", "chr3"),
  mb     = c(112.1, 25.4,   7.67,    29.42,  48.6,    50.7,   128.7,  39.72,   140.8,  178.9),
  role   = c("TSG", "Onc",  "TSG",   "TSG",  "TSG",   "TSG",  "Onc",  "Onc",   "Onc",  "Onc")
)

# Collapse a five-state matrix (1 Del, 2 Loss, 3 Nor, 4 Gain, 5 Amp) to three states (1 loss,
# 2 neutral, 3 gain); a matrix already coded 1-3 is returned unchanged.
coerce_123 <- function(mat) {
  v <- as.integer(mat)
  if (all(v %in% c(1L, 2L, 3L, NA_integer_))) {
    return(matrix(v, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat)))
  }
  out <- rep(NA_integer_, length(v))
  out[v %in% c(1L, 2L)] <- 1L
  out[v %in% 3L]        <- 2L
  out[v %in% c(4L, 5L)] <- 3L
  matrix(out, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
}

# truth = "qc" (default): the post-QC bulk WGS truth of reproduce/analysis/crc01/wgs_truth_postqc.R,
# in which CBS segments shorter than 6 Mb are cleared to neutral; "raw" is the unfiltered DNAcopy
# output. sfx selects the MintCNA call set: "_d250" is the AG-DCAE run with latent_dim 250 (197 cells).
load_all_methods_raw <- function(truth = c("qc", "raw"), sfx = "_d250") {
  truth <- match.arg(truth)
  e <- new.env()
  load(file.path(BENCH_DIR, paste0("MultiCNA_AE_197", sfx, ".RData")),        envir = e)
  load(file.path(BENCH_DIR, paste0("MultiCNA_AE_197_single", sfx, ".RData")), envir = e)
  load(file.path(BENCH_DIR, "SCOPE.RData"),                  envir = e)
  load(file.path(BENCH_DIR, "HMMcopy.RData"),                envir = e)
  if (truth == "qc" && file.exists(file.path(BENCH_DIR, "DNAcopy_qc.RData"))) {
    load(file.path(BENCH_DIR, "DNAcopy_qc.RData"), envir = e)
  } else {
    if (truth == "qc") warning("DNAcopy_qc.RData missing -- falling back to the raw truth")
    load(file.path(BENCH_DIR, "DNAcopy.RData"), envir = e)
    truth <- "raw"
  }
  message("bulk WGS truth: ", truth,
          if (truth == "qc") sprintf(" (%d of %d altered bins kept)",
                                     sum(e$cna_dnacopy != 2L), sum(e$cna_dnacopy_raw != 2L)) else "")
  load(file.path(BENCH_DIR, "SCONCE.RData"),                 envir = e)
  load(file.path(BENCH_DIR, "SCCNV.RData"),                  envir = e)
  load(file.path(BENCH_DIR, "rcCAE.RData"),                  envir = e)
  load(file.path(BENCH_DIR, "SCICoNE.RData"),                envir = e)

  common_cells <- colnames(e$cna_multicna_ae)
  sub <- function(m) m[, colnames(m) %in% common_cells, drop = FALSE]

  # SCOPE reports integer copy number: shift it onto the five-state scale (3 = neutral), capped at 5
  e$cna_scope <- e$cna_scope + 1L
  e$cna_scope[e$cna_scope >= 5] <- 5L

  raw <- list(
    "MintCNA"      = e$cna_multicna_ae,
    "MintCNA(DNA)" = e$cna_multicna_ae_single,
    "SCOPE"        = sub(e$cna_scope),
    "HMMcopy"      = sub(e$cna_hmmcopy),
    "SCONCE"       = sub(e$cna_sconce),
    "SCCNV"        = sub(e$cna_sccnv),
    "rcCAE"        = sub(e$cna_rccae),
    "SCICoNE"      = sub(e$cna_scicone)
  )
  list(
    methods   = lapply(raw, function(m) coerce_123(as.matrix(m))),
    truth_mat = as.matrix(e$cna_dnacopy),
    truth_src = truth,
    ref       = e$ref_multicna_ae
  )
}

# The MintCNA calls alone, bins ordered by chromosome and start, with the site of each cell parsed
# from its name.
load_mintcna <- function() {
  e <- new.env()
  load(file.path(BENCH_DIR, "MultiCNA_AE_197_d250.RData"), envir = e)
  mat <- e$cna_multicna_ae
  ref <- e$ref_multicna_ae
  ref$chr <- factor(ref$chr, levels = CHR_ORDER)
  ord <- order(ref$chr, ref$start)
  mat <- mat[ord, , drop = FALSE]
  ref <- ref[ord, , drop = FALSE]
  rownames(ref) <- rownames(mat)

  site <- sub("^CRC01_([A-Z]+)[0-9]*_.*", "\\1", colnames(mat))
  site <- factor(site, levels = SITE_ORDER)

  list(mat = mat, ref = ref, site = site)
}

# The bulk WGS truth has one column per site: PT, LN1, ML1 and MP1.
cell_to_truthcol <- function(cell_names) {
  site <- sub("^.*_(PT|LN|ML|MP).*?$", "\\1", cell_names)
  site[!site %in% SITE_ORDER] <- NA_character_
  truth_col <- fifelse(site == "PT", "PT",
                       fifelse(site == "LN", "LN1",
                               fifelse(site == "ML", "ML1",
                                       fifelse(site == "MP", "MP1", NA_character_))))
  data.table(cell = cell_names, site = site, truth_col = truth_col)
}

chr_bin_midpoints <- function(ref) {
  r <- as.data.frame(ref)
  r$bin_idx <- seq_len(nrow(r))
  do.call(rbind, lapply(split(r, r$chr), function(g) {
    if (nrow(g) == 0) return(NULL)
    data.frame(chr = unique(g$chr),
               mid = mean(g$bin_idx),
               start_idx = min(g$bin_idx),
               end_idx = max(g$bin_idx))
  }))
}

classify_arm17 <- function(ref) {
  mid_mb <- (ref$start + ref$end) / 2e6
  ifelse(ref$chr != "chr17", NA_character_,
         ifelse(mid_mb <  CHR17_CENTROMERE_MB, "17p", "17q"))
}

cell_arm17_loss_frac <- function(mat, ref, arm_label) {
  arm <- classify_arm17(ref)
  idx <- which(!is.na(arm) & arm == arm_label)
  colMeans(mat[idx, , drop = FALSE] == 1)
}

sig_label <- function(p) {
  ifelse(is.na(p), "ns",
         ifelse(p < 1e-4, "****",
                ifelse(p < 1e-3, "***",
                       ifelse(p < 1e-2, "**",
                              ifelse(p < 5e-2, "*", "ns")))))
}

# Per-cell QC exclusion applied to every caller at load time. CRC01_ML3_085 is the one cell whose
# smoothed log2-ratio SD exceeds 3x the cohort median (2.54 against a median of 0.51; whole
# chromosomes 4, 11 and 19 sit near log2 -6).
load_all_methods <- function(truth = c("qc", "raw"), sfx = "_d250",
                             exclude = "CRC01_ML3_085") {
  d <- load_all_methods_raw(truth, sfx)
  exclude <- exclude[nzchar(exclude)]
  if (length(exclude)) {
    d$methods <- lapply(d$methods, function(m) { k <- colnames(m) %in% exclude; if (any(k)) m[, !k, drop = FALSE] else m })
    message("dropped ", paste(exclude, collapse = ","), " from every caller")
  }
  d
}
