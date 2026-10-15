suppressPackageStartupMessages({
  library(data.table)
})

BENCH_DIR <- "<MINTCNA_ROOT>/Result/CRC/Benchmark"
ROOT_DIR  <- "<MINTCNA_ROOT>/Result/CRC/figures"
PLOTS_DIR <- file.path(ROOT_DIR, "plots")
TBLS_DIR  <- file.path(ROOT_DIR, "tables")

METHOD_ORDER <- c("MintCNA", "MintCNA(DNA)", "SCOPE", "HMMcopy",
                  "SCONCE", "SCCNV", "rcCAE", "SCICoNE")

# Two colour families (validated colorblind-safe, matching main-text Fig 2/3):
# ours = two greens (MintCNA hero + MintCNA(DNA) dark-green sibling); competitors = six non-green hues.
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

# Display labels: keep "MintCNA(DNA)" as the internal data key but render it as
# MintCNA with a superscript DNA (plotmath via parse()). Keyed BY METHOD, not by
# position, so a figure ordered by performance cannot silently mislabel a box --
# which a positional vector aligned to METHOD_ORDER would do the moment the order
# stops being METHOD_ORDER.
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

# Display ranking for the benchmark panels. The two MintCNA modes lead, in that
# order, as the methods under study; the competitors after them are ranked
# strictly by the metric that panel actually shows, so the axis order is never a
# claim the data does not support. `lower_is_better` flips the direction for
# error metrics. Note the leading pair is a layout choice, not a result: MintCNA
# and rcCAE are not separable on per-cell Pearson (paired Wilcoxon p = 0.24).
# Returns one character vector to use as factor levels, scale limits and
# summary-table row order at once, which is what keeps the three from drifting
# apart. A two-panel figure sharing one legend must pass a SINGLE order for both
# panels (see 02_03_rmse_pearson.R) or a method lands at two different x.
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

# truth = "qc"  -> post-QC bulk WGS (Code/CRC/Evaluation/WGS_postQC.R): CBS segments
#                  shorter than 6 Mb are cleared to neutral. This is the default.
#        = "raw" -> the unfiltered DNAcopy output, for before/after comparisons.
# CRC_SUFFIX (env) selects an alternative MintCNA call set, e.g. "_d250" for the 2026-10-14
# re-run (AG-DCAE latent_dim 250, 197 cells); "" (default) is the published set.
load_all_methods_raw <- function(truth = c("qc", "raw"), sfx = Sys.getenv("CRC_SUFFIX", "_d250")) {
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

load_mintcna <- function() {
  e <- new.env()
  load(file.path(BENCH_DIR, paste0("MultiCNA_AE_197", Sys.getenv("CRC_SUFFIX", "_d250"), ".RData")), envir = e)
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

# CRC_EXCLUDE (env, comma-separated cell IDs): per-cell QC exclusion applied to every caller at
# load time. 2026-10-14: "CRC01_ML3_085" -- its smoothed log2-ratio SD is 2.54 (cohort median 0.51,
# 5x; whole chromosomes 4, 11 and 19 sit at log2 ~ -6), the only cell above 3x the median.
load_all_methods <- function(truth = c("qc", "raw"), sfx = Sys.getenv("CRC_SUFFIX", "_d250"),
                             exclude = strsplit(Sys.getenv("CRC_EXCLUDE", "CRC01_ML3_085"), ",")[[1]]) {
  d <- load_all_methods_raw(truth, sfx)
  exclude <- exclude[nzchar(exclude)]
  if (length(exclude)) {
    d$methods <- lapply(d$methods, function(m) { k <- colnames(m) %in% exclude; if (any(k)) m[, !k, drop = FALSE] else m })
    message("CRC_EXCLUDE: dropped ", paste(exclude, collapse = ","), " from every caller")
  }
  d
}
