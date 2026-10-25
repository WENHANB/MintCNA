# MintCNA - Preprocessing: GC/mappability, DNA / RNA / ATAC log2-ratios, alignment, QN, genome smoothing
# Source-able R functions; sourced by scripts/run_mintcna.R.


# ======================= GC and mappability =======================

# GC content and mappability attachment
#
# The user supplies read counts and bin locations only. MintCNA ships a 5 kb
# gc/mappability track per genome build and aggregates it onto the user's bins.

#' Load a bundled 5 kb gc/mappability reference track.
#'
#' @param genome Genome build. "hg19" and "hg38" are both bundled.
#' @param ref_dir Directory holding "<genome>_gcmapp_5kb.csv.gz".
#' @return Data frame with columns chr, start, end, gc, mapp (+ n_valid weight).
.load_gcmapp_track <- function(genome, ref_dir) {
  f <- file.path(ref_dir, paste0(genome, "_gcmapp_5kb.csv.gz"))
  if (!file.exists(f)) {
    stop("Bundled gc/mappability track not found: ", f,
         "\n  Supported builds: hg19, hg38.",
         "\n  Either set input$genome to a supported build, supply gc/mapp",
         "\n  columns in the bin annotation, or add a track to MintCNA/ref/.")
  }
  read.csv(f, stringsAsFactors = FALSE)
}

#' Aggregate a fine track onto arbitrary bins (length-weighted overlap).
#'
#' GC is weighted by the number of valid (non-N) bases in each overlapping tile
#' (near-exact); mappability is weighted by overlap length (standard). Bins with
#' no overlapping track tile get NA.
#'
#' @param bins Data frame with columns chr, start, end (0/1-based both fine at
#'   this resolution).
#' @param track Fine track from .load_gcmapp_track (regular grid of width `tile`).
#' @param tile Track bin width in bp (default 5000).
#' @return Data frame with columns gc, mapp, frac_valid (one row per input bin).
#'   frac_valid is the fraction of the bin's length covered by track tiles that
#'   carry a GC value; a bin on a centromere or assembly gap has GC/mappability
#'   defined on only a sliver of itself (used by the ATAC arm's valid_min filter).
aggregate_gcmapp <- function(bins, track, tile = 5000L) {
  bins  <- as.data.frame(bins, stringsAsFactors = FALSE)
  track <- as.data.frame(track, stringsAsFactors = FALSE)
  has_w <- "n_valid" %in% names(track)

  # Per-chromosome lookups keyed by tile index (= start %/% tile).
  rows_by_chr <- split(seq_len(nrow(track)), track$chr)
  lut <- lapply(rows_by_chr, function(ii) {
    ti <- track$start[ii] %/% tile
    w  <- if (has_w) ifelse(is.na(track$n_valid[ii]), 0, track$n_valid[ii]) else rep(tile, length(ii))
    list(gc = setNames(track$gc[ii],   ti),
         mp = setNames(track$mapp[ii], ti),
         w  = setNames(w,              ti))
  })

  gc_out <- mp_out <- fv_out <- rep(NA_real_, nrow(bins))
  for (i in seq_len(nrow(bins))) {
    L <- lut[[ as.character(bins$chr[i]) ]]
    if (is.null(L)) next
    s <- bins$start[i]; e <- bins$end[i]
    if (is.na(s) || is.na(e) || e <= s) next

    ti  <- seq.int(s %/% tile, (e - 1L) %/% tile)
    key <- as.character(ti)
    ts  <- ti * tile; te <- ts + tile
    ov  <- pmin(te, e) - pmax(ts, s)          # overlap length per tile

    g <- L$gc[key]; m <- L$mp[key]; wv <- L$w[key]
    wv[is.na(wv)] <- 0
    fv_out[i] <- sum(ov[!is.na(g)]) / sum(ov)   # fraction of the bin with a GC value

    wgc <- ov * wv / tile                     # GC weight: valid bases in overlap
    okg <- !is.na(g) & wgc > 0
    if (any(okg)) gc_out[i] <- sum(g[okg] * wgc[okg]) / sum(wgc[okg])

    okm <- !is.na(m) & ov > 0                 # mappability weight: overlap length
    if (any(okm)) mp_out[i] <- sum(m[okm] * ov[okm]) / sum(ov[okm])
  }
  data.frame(gc = gc_out, mapp = mp_out, frac_valid = fv_out)
}

#' Parse bin locations from count-matrix row names.
#'
#' Accepts IDs like "chr1:900000-1000000" or "chr1_900000_1000000" so a user can
#' supply only a count matrix (bin locations encoded in its row names) with no
#' separate annotation file.
#'
#' @param ids Character vector of bin identifiers.
#' @return Data frame with columns chr, start, end.
parse_bins_from_ids <- function(ids) {
  m  <- regmatches(ids, regexec("^(chr[0-9XYM]+)[:_]([0-9]+)[-_]([0-9]+)$", ids))
  ok <- lengths(m) == 4L
  if (!all(ok)) {
    stop(sum(!ok), " count row name(s) are not parseable bin locations; ",
         "provide input$bin_annotation with chr,start,end columns. ",
         "First unparseable: '", ids[which(!ok)[1]], "'")
  }
  data.frame(chr   = vapply(m, `[`, "", 2L),
             start = as.numeric(vapply(m, `[`, "", 3L)),
             end   = as.numeric(vapply(m, `[`, "", 4L)),
             stringsAsFactors = FALSE)
}

#' Attach gc/mapp columns to user bins from the bundled genome track.
#'
#' @param bins Data frame with chr, start, end (extra columns preserved).
#' @param genome Genome build (default "hg19").
#' @param ref_dir Directory of bundled tracks. Defaults to option
#'   "mintcna.ref_dir" or "MintCNA/ref".
#' @param tile Track resolution in bp (default 5000).
#' @return `bins` with gc, mapp and frac_valid columns appended.
attach_gc_mapp <- function(bins, genome = "hg19", ref_dir = NULL, tile = 5000L) {
  if (is.null(ref_dir)) ref_dir <- getOption("mintcna.ref_dir", "MintCNA/ref")
  track <- .load_gcmapp_track(genome, ref_dir)
  gm <- aggregate_gcmapp(bins, track, tile = tile)
  out <- as.data.frame(bins, stringsAsFactors = FALSE)
  out$gc <- gm$gc; out$mapp <- gm$mapp; out$frac_valid <- gm$frac_valid
  n_na <- sum(is.na(out$gc) | is.na(out$mapp))
  if (n_na > 0L) {
    message(sprintf("[attach_gc_mapp] %d/%d bins had no track overlap (gc/mapp set to NA)",
                    n_na, nrow(out)))
  }
  out
}

# ======================= DNA preprocessing =======================

# MintCNA - DNA Preprocessing
# L2R computation against normal controls + smoothing

#' Compute log2 ratios for DNA read counts against normal controls
#'
#' Takes GC/MAP-normalized and library-size-normalized read counts,
#' computes log2 ratios against the median of normal control cells,
#' applies median shift normalization.
#'
#' @param rc_norm Normalized read count matrix (bins x cells).
#'   Output from GC.MAP.normalization() + libsize_normalize().
#' @param ref Data frame with bin annotation (chr, start, end, gc, mapp).
#' @param normal_cells Character vector of cell IDs that are normal controls.
#' @param eps Small constant to avoid log(0). Default: 0.01.
#' @return List with:
#'   $l2r: bins x cells log2-ratio matrix (tumor cells only);
#'   $ref: bin reference (reordered by genome position).
preprocess_dna <- function(rc_norm, ref, normal_cells,
                           eps = 0.01) {
  rc_norm <- as.matrix(rc_norm)

  # --- L2R against normal control median ---
  message("[MintCNA DNA] Computing log2 ratios against normal controls...")
  nc_in_data <- intersect(normal_cells, colnames(rc_norm))
  if (length(nc_in_data) == 0) stop("No normal control cells found in the data.")
  nc_median <- apply(rc_norm[, nc_in_data, drop = FALSE], 1, median)

  tumor_cells <- setdiff(colnames(rc_norm), normal_cells)
  l2r <- apply(rc_norm[, tumor_cells, drop = FALSE], 2, function(col) {
    log2((col + eps) / (nc_median + eps))
  })
  rownames(l2r) <- rownames(rc_norm)
  colnames(l2r) <- tumor_cells

  # --- Median shift normalization ---
  message("[MintCNA DNA] Median shift normalization...")
  l2r <- apply(l2r, 2, function(col) col - median(col, na.rm = TRUE))

  # --- Order by genome position ---
  ord <- order(as.numeric(sub("chr", "", ref$chr)), ref$start)
  l2r <- l2r[ord, , drop = FALSE]
  ref <- ref[ord, ]

  message(sprintf("[MintCNA DNA] Done. Output: %d bins x %d cells.", nrow(l2r), ncol(l2r)))

  list(l2r = l2r, ref = ref)
}

#' Full DNA preprocessing pipeline
#'
#' Runs the complete DNA pipeline: QC filtering, GC/MAP normalization,
#' library size normalization, L2R computation, and smoothing.
#'
#' @param Y Raw read count matrix (bins x cells).
#' @param ref Data frame with columns: chr, start, end, gc, mapp.
#' @param normal_cells Character vector of normal control cell IDs.
#' @param gc_min Minimum GC for bin filtering. Default: 0.2.
#' @param gc_max Maximum GC for bin filtering. Default: 0.8.
#' @param mapp_min Minimum mappability for bin filtering. Default: 0.9.
#' @return List with:
#'   $l2r: bins x cells log2-ratio matrix (tumor cells only);
#'   $ref: filtered and ordered bin reference;
#'   $rc_norm: normalized read counts (before L2R).
run_dna_preprocessing <- function(Y, ref, normal_cells,
                                  gc_min = 0.2, gc_max = 0.8, mapp_min = 0.9) {
  # Step 1: QC filter bins
  filt <- filter_bins(ref, Y, gc_min = gc_min, gc_max = gc_max, mapp_min = mapp_min)
  Y_filt   <- filt$Y
  ref_filt <- filt$ref

  # Step 2: GC + mappability normalization
  message("[MintCNA DNA] GC/mappability normalization...")
  norm_result <- GC.MAP.normalization(Y_filt, ref_filt)

  # Step 3: Library size normalization
  message("[MintCNA DNA] Library size normalization...")
  lib_result <- libsize_normalize(norm_result$RC_norm, method = "cpm")

  # Step 4: L2R
  dna_result <- preprocess_dna(lib_result$X_norm, ref_filt, normal_cells)

  dna_result$rc_norm <- lib_result$X_norm
  dna_result
}

# ======================= Normalization =======================

# MintCNA - Normalization Functions
# GC-content and mappability bias correction for read count data

#' GC-content bias correction via median scaling
#'
#' @param RC Numeric vector of read counts for one cell.
#' @param GCContent Numeric vector of GC content per bin (same length as RC).
#' @param step Step size for GC binning.
#' @return List with $RCNorm (corrected counts), $Median, $StepGC.
CorrectGC <- function(RC, GCContent, step) {
  stepseq <- seq(min(GCContent), max(GCContent), by = step)
  MasterMedian <- median(RC, na.rm = TRUE)
  MedianGC <- rep(0, length(stepseq) - 1)
  RCNormMedian <- RC
  for (i in 1:(length(stepseq) - 1)) {
    if (i == 1) {
      ind <- which(GCContent >= stepseq[i] & GCContent <= stepseq[i + 1])
    } else {
      ind <- which(GCContent > stepseq[i] & GCContent <= stepseq[i + 1])
    }
    if (length(ind) > 0) {
      m <- median(RC[ind], na.rm = TRUE)
      if (m > 0) {
        MedianGC[i] <- m
        RCNormMedian[ind] <- RC[ind] * MasterMedian / m
      }
    }
  }
  list(Median = MedianGC, StepGC = stepseq[1:(length(stepseq) - 1)], RCNorm = RCNormMedian)
}

#' Mappability bias correction via median scaling
#'
#' @param RC Numeric vector of read counts for one cell.
#' @param MAPContent Numeric vector of mappability per bin.
#' @param step Step size for mappability binning.
#' @return List with $RCNorm (corrected counts), $Median, $StepMAP.
CorrectMAP <- function(RC, MAPContent, step) {
  stepseq <- seq(min(MAPContent), max(MAPContent), by = step)
  MasterMedian <- median(RC, na.rm = TRUE)
  MedianMAP <- rep(0, length(stepseq) - 1)
  RCNormMedian <- RC
  for (i in 1:(length(stepseq) - 1)) {
    if (i == 1) {
      ind <- which(MAPContent >= stepseq[i] & MAPContent <= stepseq[i + 1])
    } else {
      ind <- which(MAPContent > stepseq[i] & MAPContent <= stepseq[i + 1])
    }
    if (length(ind) > 0) {
      m <- median(RC[ind], na.rm = TRUE)
      if (m > 0) {
        MedianMAP[i] <- m
        RCNormMedian[ind] <- RC[ind] * MasterMedian / m
      }
    }
  }
  list(Median = MedianMAP, StepMAP = stepseq[1:(length(stepseq) - 1)], RCNorm = RCNormMedian)
}

#' Combined GC and mappability normalization
#'
#' Applies mappability correction then GC correction per cell,
#' then computes log2 ratios against per-cell median.
#'
#' @param Y Raw read count matrix (bins x cells).
#' @param ref_qc Data frame with columns: start, end, gc, mapp.
#' @return List with $RC_norm (normalized counts) and $log2Rdata (log2 ratios).
GC.MAP.normalization <- function(Y, ref_qc) {
  gc  <- ref_qc$gc
  map <- ref_qc$mapp

  RC_norm <- matrix(NA, nrow = nrow(Y), ncol = ncol(Y),
                    dimnames = list(rownames(Y), colnames(Y)))
  for (sub in 1:ncol(Y)) {
    RCMAPNorm <- CorrectMAP(Y[, sub], MAPContent = map, step = 0.01)$RCNorm
    RCGCNorm  <- CorrectGC(RCMAPNorm, GCContent = gc, step = 0.05)$RCNorm
    RC_norm[, sub] <- RCGCNorm
  }

  lrr <- matrix(NA, nrow = nrow(RC_norm), ncol = ncol(RC_norm),
                dimnames = list(rownames(RC_norm), colnames(RC_norm)))
  for (i in 1:ncol(RC_norm)) {
    med <- median(RC_norm[, i])
    lrr[, i] <- log2((RC_norm[, i] + 0.01) / med)
  }

  list(RC_norm = RC_norm, log2Rdata = lrr)
}

#' Library size normalization
#'
#' @param X Raw count matrix (features x cells).
#' @param method "cpm" (counts per million) or "prop" (proportions).
#' @param log1p If TRUE, apply log1p transformation after normalization.
#' @param scale Scaling factor (default 1e6 for CPM).
#' @param eps Floor for library size to avoid division by zero.
#' @return List with $X_norm and $libsize.
libsize_normalize <- function(X, method = c("cpm", "prop"),
                              log1p = FALSE, scale = 1e6, eps = 1e-8) {
  method <- match.arg(method)
  X <- as.matrix(X)
  lib <- colSums(X, na.rm = TRUE)
  lib_safe <- pmax(lib, eps)
  if (method == "prop") scale <- 1
  Xn <- sweep(X, 2, lib_safe, "/") * scale
  if (log1p) Xn <- log1p(Xn)
  list(X_norm = Xn, libsize = lib)
}

# ======================= QC filtering =======================

# MintCNA - QC Filtering
# Bin-level quality control based on GC content and mappability

#' Filter bins by GC content and mappability
#'
#' Removes bins with extreme GC content or low mappability.
#'
#' @param ref Data frame with bin annotation (must have columns: gc, mapp).
#' @param Y Optional count matrix (bins x cells) to subset in parallel.
#' @param gc_min Minimum GC content (default 0.2).
#' @param gc_max Maximum GC content (default 0.8).
#' @param mapp_min Minimum mappability (default 0.9).
#' @return List with $ref (filtered ref), $Y (filtered matrix, if provided),
#'   $keep_idx (logical vector of kept bins), $n_removed (count removed).
filter_bins <- function(ref, Y = NULL, gc_min = 0.2, gc_max = 0.8, mapp_min = 0.9) {
  keep <- (ref$gc >= gc_min) & (ref$gc <= gc_max) & (ref$mapp >= mapp_min)
  keep[is.na(keep)] <- FALSE

  n_removed <- sum(!keep)
  message(sprintf("[MintCNA] QC filter: keeping %d / %d bins (removed %d)",
                  sum(keep), length(keep), n_removed))

  out <- list(ref = ref[keep, ], keep_idx = keep, n_removed = n_removed)
  if (!is.null(Y)) out$Y <- Y[keep, , drop = FALSE]
  out
}

# ======================= RNA preprocessing =======================

# MintCNA - RNA Preprocessing (CopyKAT-style)
# CPM normalization, filtering, log-Freeman-Tukey, DLM smoothing, L2R

library(edgeR)
library(dlm)

#' Preprocess scRNA-seq counts for CNA inference
#'
#' Applies the CopyKAT-style pipeline: CPM normalization, cell/gene QC filtering,
#' log-Freeman-Tukey transformation, DLM smoothing, normal control subtraction,
#' and median shift normalization.
#'
#' @param counts Gene x cell raw count matrix (integer counts).
#' @param gene_ref Data frame with gene annotation. Must have column 'gene_id'
#'   matching rownames of counts, plus chr, start, end.
#' @param normal_cells Character vector of cell IDs that are normal controls.
#'   Used as reference for L2R computation.
#' @param min_genes_per_cell Minimum genes detected per cell. Default: 200.
#' @param min_detection_rate Minimum fraction of cells a gene must be detected in. Default: 0.05.
#' @param dlm_dV Observation variance for DLM smoothing. Default: 0.16.
#' @param dlm_dW System variance for DLM smoothing. Default: 0.001.
#' @param mc_cores Number of cores for parallel DLM smoothing. Default: 1.
#' @return List with:
#'   $l2r: gene x cell log2-ratio matrix (tumor cells only, after all processing);
#'   $ref: filtered gene reference (aligned to rows of l2r);
#'   $filtered_cells: cell IDs retained after QC;
#'   $filtered_genes: gene IDs retained after QC.
preprocess_rna <- function(counts, gene_ref, normal_cells,
                           min_genes_per_cell = 200,
                           min_detection_rate = 0.05,
                           dlm_dV = 0.16, dlm_dW = 0.001,
                           mc_cores = 1) {
  counts <- as.matrix(counts)

  # --- Step 1: CPM normalization ---
  message("[MintCNA RNA] Step 1: CPM normalization...")
  dge <- DGEList(counts = counts)
  cpm_mat <- cpm(dge)

  # --- Step 2: Cell QC filter (min genes per cell) ---
  message("[MintCNA RNA] Step 2: Cell QC filtering...")
  genes_per_cell <- colSums(cpm_mat > 0)
  keep_cells <- genes_per_cell >= min_genes_per_cell
  if (sum(keep_cells) == 0) stop("No cells have more than ", min_genes_per_cell, " genes.")
  if (sum(!keep_cells) > 0) {
    message(sprintf("  Filtered out %d cells with < %d genes; %d cells remain.",
                    sum(!keep_cells), min_genes_per_cell, sum(keep_cells)))
  }
  cpm_mat <- cpm_mat[, keep_cells, drop = FALSE]

  # --- Step 3: Gene QC filter (detection rate) ---
  message("[MintCNA RNA] Step 3: Gene QC filtering...")
  detection_rate <- rowSums(cpm_mat > 0) / ncol(cpm_mat)
  keep_genes <- detection_rate > min_detection_rate
  if (sum(keep_genes) == 0) stop("No genes pass detection rate filter.")
  message(sprintf("  %d genes preserved after %.0f%% detection rate filter.",
                  sum(keep_genes), min_detection_rate * 100))
  cpm_mat <- cpm_mat[keep_genes, , drop = FALSE]
  ref_filtered <- gene_ref[keep_genes, ]

  if (nrow(cpm_mat) < 7000) {
    warning("Low data quality: fewer than 7000 genes after filtering.")
  }

  # --- Step 4: Log-Freeman-Tukey transformation ---
  message("[MintCNA RNA] Step 4: Log-Freeman-Tukey transformation...")
  norm_mat <- log(sqrt(cpm_mat) + sqrt(cpm_mat + 1))

  # --- Step 5: DLM smoothing ---
  message("[MintCNA RNA] Step 5: DLM smoothing...")
  dlm_smooth_col <- function(c) {
    model <- dlmModPoly(order = 1, dV = dlm_dV, dW = dlm_dW)
    x <- dlmSmooth(norm_mat[, c], model)$s
    x <- x[2:length(x)]
    x - mean(x)
  }

  smoothed_list <- parallel::mclapply(1:ncol(norm_mat), dlm_smooth_col, mc.cores = mc_cores)
  norm_mat_smooth <- matrix(unlist(smoothed_list), ncol = ncol(norm_mat), byrow = FALSE)
  colnames(norm_mat_smooth) <- colnames(norm_mat)
  rownames(norm_mat_smooth) <- rownames(cpm_mat)

  # --- Step 6: Normal control subtraction ---
  message("[MintCNA RNA] Step 6: Normal control subtraction...")
  nc_in_data <- intersect(normal_cells, colnames(norm_mat_smooth))
  if (length(nc_in_data) == 0) stop("No normal control cells found in the data.")
  nc_median <- apply(norm_mat_smooth[, nc_in_data, drop = FALSE], 1, median)
  norm_mat_smooth <- norm_mat_smooth - nc_median

  # Remove normal cells from output
  tumor_cells <- setdiff(colnames(norm_mat_smooth), normal_cells)
  l2r <- norm_mat_smooth[, tumor_cells, drop = FALSE]

  # --- Step 7: Median shift normalization ---
  message("[MintCNA RNA] Step 7: Median shift normalization...")
  l2r <- apply(l2r, 2, function(col) col - median(col, na.rm = TRUE))

  message(sprintf("[MintCNA RNA] Done. Output: %d genes x %d cells.", nrow(l2r), ncol(l2r)))

  list(l2r = l2r, ref = ref_filtered,
       filtered_cells = colnames(cpm_mat),
       filtered_genes = rownames(cpm_mat))
}

# ======================= ATAC preprocessing =======================

# MintCNA - scATAC Preprocessing
# Turns single-cell ATAC fragment counts into the same L2R matrix the DNA arm
# emits, so everything downstream runs unchanged. Because ATAC coverage also
# reflects chromatin accessibility and is sparse per bin, this arm (i) filters
# bins at the resolution they were built at, plus the ENCODE blacklist, (ii)
# pools the raw counts of the reference nuclei instead of taking a per-cell
# median, and (iii) skips per-cell GC/mappability correction, since the shared
# bias cancels in the ratio to the reference. The ratio is formed in count space
# with a Poisson pseudocount.

`%||%` <- function(x, y) if (is.null(x) || (length(x) == 1 && is.na(x))) y else x

#' Log2 ratios for scATAC bin counts against a pooled reference profile
#'
#' @param Y Raw fragment count matrix (bins x cells), one count per fragment.
#' @param ref Bin annotation; needs chr, start, end and (post-filter) nothing else.
#' @param normal_cells Character vector of reference (diploid) cell IDs.
#' @param pseudocount Poisson pseudocount added to observed and expected counts.
#' @param keep_normal If TRUE the reference nuclei stay in the output as a
#'   built-in negative control; the DNA arm drops them, and so does the default.
#' @return List with $l2r (bins x cells), $ref, $p_ref (reference profile),
#'   $libsize, and $diagnostics.
preprocess_atac <- function(Y, ref, normal_cells, pseudocount = 1,
                            keep_normal = FALSE) {
  Y <- as.matrix(Y)
  nc <- intersect(normal_cells, colnames(Y))
  if (!length(nc)) stop("No reference cells found in the count matrix.")

  # --- pooled reference profile ------------------------------------------
  # Sum of raw counts: this weights each nucleus by its depth.
  ref_tot <- rowSums(Y[, nc, drop = FALSE])
  p_ref   <- ref_tot / sum(ref_tot)
  message(sprintf("[MintCNA ATAC] Reference: %d nuclei, %s fragments, median %.0f per bin.",
                  length(nc), format(sum(ref_tot), big.mark = ","), median(ref_tot)))

  out_cells <- if (keep_normal) colnames(Y) else setdiff(colnames(Y), normal_cells)
  X   <- Y[, out_cells, drop = FALSE]
  lib <- colSums(X)

  # --- expected counts, then the ratio in count space ---------------------
  E <- outer(p_ref, lib)
  l2r <- log2((X + pseudocount) / (E + pseudocount))
  dimnames(l2r) <- dimnames(X)

  # --- median shift, as in the DNA arm ------------------------------------
  l2r <- sweep(l2r, 2, apply(l2r, 2, median, na.rm = TRUE))

  ord <- order(as.numeric(sub("chr", "", ref$chr)), ref$start)
  lam <- median(lib) * p_ref
  diagnostics <- list(
    n_reference   = length(nc),
    lambda_median = median(lam[lam > 0]),
    noise_log2    = median(1.443 / sqrt(lam[lam > 0])),
    frac_zero     = mean(X == 0))
  message(sprintf("[MintCNA ATAC] lambda = %.2f in the median bin -> per-bin log2 noise %.2f; %.0f%% of entries are zero.",
                  diagnostics$lambda_median, diagnostics$noise_log2, 100 * diagnostics$frac_zero))

  list(l2r = l2r[ord, , drop = FALSE], ref = ref[ord, ],
       p_ref = p_ref[ord], libsize = lib, diagnostics = diagnostics)
}

#' Bin-level QC for scATAC, at the resolution the bins were actually built at
#'
#' @param ref Bin annotation with gc, mapp and -- if available -- frac_valid and
#'   frac_blacklist.
#' @param Y Optional count matrix (bins x cells) to subset alongside.
#' @param mapp_min Minimum mappability. Default 0.5.
#' @param blacklist_max Maximum fraction of the bin inside the ENCODE blacklist.
#' @param valid_min Minimum fraction of the bin carrying a GC/mappability value.
#' @param max_zero_frac Drop bins that are zero in more than this fraction of
#'   cells (epiAneufinder's default sparse-bin filter). Set NULL to skip.
filter_bins_atac <- function(ref, Y = NULL, gc_min = 0.2, gc_max = 0.8,
                             mapp_min = 0.5, blacklist_max = 0.9, valid_min = 0.5,
                             max_zero_frac = 0.85) {
  keep <- (ref$gc >= gc_min) & (ref$gc <= gc_max) & (ref$mapp >= mapp_min)
  if (!is.null(ref$frac_blacklist)) keep <- keep & (ref$frac_blacklist <= blacklist_max)
  if (!is.null(ref$frac_valid))     keep <- keep & (ref$frac_valid     >= valid_min)
  keep[is.na(keep)] <- FALSE
  n_annot <- sum(!keep)

  if (!is.null(Y) && !is.null(max_zero_frac)) {
    zf <- rowMeans(as.matrix(Y) == 0)
    keep <- keep & (zf <= max_zero_frac)
  }

  message(sprintf("[MintCNA ATAC] Bin QC: keeping %d / %d bins (%d on annotation, %d on sparsity).",
                  sum(keep), length(keep), n_annot, sum(!keep) - n_annot))
  out <- list(ref = ref[keep, ], keep_idx = keep, n_removed = sum(!keep))
  if (!is.null(Y)) out$Y <- Y[keep, , drop = FALSE]
  out
}

#' Full scATAC preprocessing pipeline
#'
#' @inheritParams preprocess_atac
#' @inheritParams filter_bins_atac
#' @return List with $l2r, $ref, $p_ref, $libsize, $diagnostics.
run_atac_preprocessing <- function(Y, ref, normal_cells,
                                   gc_min = 0.2, gc_max = 0.8, mapp_min = 0.5,
                                   blacklist_max = 0.9, valid_min = 0.5,
                                   max_zero_frac = 0.85,
                                   pseudocount = 1, keep_normal = FALSE) {
  f <- filter_bins_atac(ref, Y, gc_min = gc_min, gc_max = gc_max, mapp_min = mapp_min,
                        blacklist_max = blacklist_max, valid_min = valid_min,
                        max_zero_frac = max_zero_frac)
  preprocess_atac(f$Y, f$ref, normal_cells, pseudocount = pseudocount,
                  keep_normal = keep_normal)
}

# ======================= Multi-omic alignment =======================

# MintCNA - Multi-omic alignment
# Puts the RNA arm onto the DNA arm's bin grid so both modalities share one
# coordinate axis, which is what the multivariate CUSUM scan requires.

#' Align RNA gene-level L2R onto the DNA bin grid
#'
#' Each gene is assigned to every bin it overlaps, gene values are collapsed to
#' the per-bin MEDIAN within a cell, and the result is placed on the full DNA bin
#' axis. Bins containing no gene stay NA -- the CUSUM (`cusum.vector`) and the
#' fusion weights both count only non-NA entries, so an RNA-empty bin correctly
#' contributes nothing rather than a zero.
#'
#' @param dna_l2r Bins x cells DNA (or ATAC) log2-ratio matrix.
#' @param dna_ref Bin annotation with chr, start, end.
#' @param rna_l2r Genes x cells RNA log2-ratio matrix.
#' @param rna_ref Gene annotation with chr, start, end.
#' @return List with $combined_list (per cell, a 2 x n_bins matrix with rows
#'   RNA and DNA), $ref (the ordered bin annotation) and $shared_samples.
align_omics <- function(dna_l2r, dna_ref, rna_l2r, rna_ref) {
  stopifnot(requireNamespace("data.table", quietly = TRUE))
  dna_ref$bin_id <- paste(dna_ref$chr, dna_ref$start, dna_ref$end, sep = "_")

  rna_gr <- GenomicRanges::GRanges(as.character(rna_ref$chr),
                                   IRanges::IRanges(rna_ref$start, rna_ref$end))
  dna_gr <- GenomicRanges::GRanges(as.character(dna_ref$chr),
                                   IRanges::IRanges(dna_ref$start, dna_ref$end))
  GenomeInfoDb::seqlevelsStyle(dna_gr) <- GenomeInfoDb::seqlevelsStyle(rna_gr)

  shared_samples <- intersect(colnames(rna_l2r), colnames(dna_l2r))
  if (length(shared_samples) == 0) stop("No shared cell IDs between DNA and RNA matrices.")
  message(sprintf("[MintCNA] Found %d shared cells between DNA and RNA.", length(shared_samples)))

  # gene -> bin, keeping every overlap (a gene spanning a boundary hits two bins)
  hits <- GenomicRanges::findOverlaps(rna_gr, dna_gr)
  g2b  <- data.table::data.table(gene_idx = S4Vectors::queryHits(hits),
                                 bin_idx  = S4Vectors::subjectHits(hits))
  message(sprintf("[MintCNA] %d of %d genes map to a bin (%d gene-bin pairs).",
                  data.table::uniqueN(g2b$gene_idx), nrow(rna_l2r), nrow(g2b)))

  # Per-bin median over the genes in that bin, for all cells at once.
  R <- as.matrix(rna_l2r[, shared_samples, drop = FALSE])
  rna_binned <- matrix(NA_real_, nrow = nrow(dna_ref), ncol = length(shared_samples),
                       dimnames = list(dna_ref$bin_id, shared_samples))
  genes_in_bin <- split(g2b$gene_idx, g2b$bin_idx)
  for (b in names(genes_in_bin)) {
    gi <- genes_in_bin[[b]]
    rna_binned[as.integer(b), ] <- if (length(gi) == 1L) R[gi, ] else
      matrixStats::colMedians(R[gi, , drop = FALSE], na.rm = TRUE)
  }
  rna_binned[is.nan(rna_binned)] <- NA_real_      # a bin whose genes are all NA

  # Order both arms by genome position
  ord <- order(as.numeric(sub("chr", "", dna_ref$chr)), dna_ref$start)
  dna_l2r_ordered <- dna_l2r[ord, shared_samples, drop = FALSE]
  dna_ref_ordered <- dna_ref[ord, ]
  rna_binned      <- rna_binned[ord, , drop = FALSE]
  rownames(dna_l2r_ordered) <- dna_ref_ordered$bin_id

  combined_list <- lapply(seq_along(shared_samples), function(j)
    rbind(RNA = rna_binned[, j], DNA = dna_l2r_ordered[, j]))
  names(combined_list) <- shared_samples

  ref_out <- data.frame(chr = dna_ref_ordered$chr, start = dna_ref_ordered$start,
                        end = dna_ref_ordered$end, stringsAsFactors = FALSE)
  list(combined_list = combined_list, ref = ref_out, shared_samples = shared_samples)
}

# ======================= Quantile-normal transformation =======================

# MintCNA - Quantile-Normal Transformation
# Transforms each modality row to N(0,1) via rank-based quantile mapping

#' Quantile-normal transformation per modality
#'
#' Maps each row of a matrix to N(mean_target, sd_target) via rank-based
#' quantile normalization. Optionally applies per-chromosome.
#'
#' @param X Numeric matrix (modalities x bins). Rows = omics, cols = bins.
#' @param chr_ranges Optional list of integer vectors (column indices per chromosome).
#'   If NULL, transformation is applied globally per row.
#' @param mean_target Target mean (scalar or vector of length nrow(X)).
#' @param sd_target Target SD (scalar or vector of length nrow(X)).
#' @return Transformed matrix of the same dimensions.
qn_modal <- function(X, chr_ranges = NULL, mean_target = 0, sd_target = 1) {
  X <- as.matrix(X)
  out <- X
  R <- nrow(X)
  if (length(sd_target) == 1) sd_target <- rep(sd_target, R)
  if (length(mean_target) == 1) mean_target <- rep(mean_target, R)

  if (is.null(chr_ranges)) {
    for (r in 1:R) {
      ok <- is.finite(X[r, ])
      m  <- sum(ok)
      if (m) {
        rk <- rank(X[r, ok], ties.method = "average")
        u  <- (rk - 0.5) / m
        z  <- qnorm(u)
        out[r, ok] <- mean_target[r] + sd_target[r] * z
      }
    }
  } else {
    for (r in 1:R) {
      for (rg in chr_ranges) {
        ok <- is.finite(X[r, rg])
        m  <- sum(ok)
        if (m) {
          rk <- rank(X[r, rg][ok], ties.method = "average")
          u  <- (rk - 0.5) / m
          z  <- qnorm(u)
          out[r, rg][ok] <- mean_target[r] + sd_target[r] * z
        }
      }
    }
  }
  out
}

# ======================= Genome-direction smoothing =======================

# MintCNA - Genome-direction smoothing
#
# AG-DCAE shares information across cells but does not know that adjacent bins
# are adjacent. A local-level state space model (random walk + noise) recovers
# the piecewise-constant structure along the genome, one cell at a time, within
# a chromosome; the Kalman smoother's posterior mean is the denoised track. Only
# the ratio dV/dW matters, so `ratio` is the single parameter. Smoothing also
# lowers the amplitude at real breakpoints, which is why it is exposed rather
# than fixed. It is meant for coverage assays that are sparse per bin (ATAC); the
# RNA arm is already smooth after per-bin aggregation. Off unless a configuration
# turns it on.

#' Local-level Kalman smoothing along the genome, one track
#'
#' @param y Numeric vector, one cell's signal in genomic order.
#' @param idx List of integer vectors, one per chromosome: which positions of `y`
#'   belong to it. Smoothing never crosses a chromosome boundary.
#' @param mod A dlm model object from `dlm::dlmModPoly`.
#' @return A vector of the same length; NA positions stay NA.
smooth_track <- function(y, idx, mod) {
  for (i in idx) {
    z  <- y[i]
    ok <- is.finite(z)
    # A chromosome with fewer than two observations has nothing to smooth, and
    # dlmSmooth() errors on an empty series rather than returning it unchanged.
    if (sum(ok) < 2L) next
    # dlmSmooth returns n+1 states, the t=0 prior first; drop it.
    y[i[ok]] <- dlm::dlmSmooth(z[ok], mod)$s[-1]
  }
  y
}

#' Smooth selected modalities of a MintCNA cell list
#'
#' @param cell_list Named list, one entry per cell, each a (modalities x bins)
#'   matrix whose rownames are the modality names.
#' @param chr Chromosome of each bin, length = ncol of a cell's matrix.
#' @param ratios Named numeric vector, modality -> dV/dW. Modalities absent from
#'   it, or given 0, are left untouched. Larger is smoother.
#' @param mc_cores Cells are independent; split them over workers.
#' @return The cell list, same shape.
smooth_cell_list <- function(cell_list, chr, ratios, mc_cores = 1) {
  ratios <- ratios[is.finite(ratios) & ratios > 0]
  if (!length(ratios)) return(cell_list)
  if (!requireNamespace("dlm", quietly = TRUE))
    stop("smooth_cell_list() needs the 'dlm' package")

  present <- rownames(cell_list[[1]])
  unknown <- setdiff(names(ratios), present)
  if (length(unknown))
    stop("smooth_genome names a modality this run does not have: ",
         paste(unknown, collapse = ", "), " (have: ", paste(present, collapse = ", "), ")")
  stopifnot(length(chr) == ncol(cell_list[[1]]))

  idx  <- split(seq_along(chr), factor(chr, levels = unique(chr)))
  mods <- lapply(ratios, function(r) dlm::dlmModPoly(order = 1, dV = 1, dW = 1 / r))

  one <- function(m) {
    for (nm in names(ratios)) m[nm, ] <- smooth_track(m[nm, ], idx, mods[[nm]])
    m
  }
  if (mc_cores > 1L) {
    out <- parallel::mclapply(cell_list, one, mc.cores = mc_cores)
    bad <- vapply(out, function(x) is.null(x) || inherits(x, "try-error"), logical(1))
    if (any(bad)) stop("smooth_cell_list(): ", sum(bad), " mclapply workers failed")
    names(out) <- names(cell_list)
    out
  } else lapply(cell_list, one)
}
