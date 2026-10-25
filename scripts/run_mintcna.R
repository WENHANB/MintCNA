#!/usr/bin/env Rscript
# MintCNA - R Pipeline Entry Point
#
# Runs the R portion of the MintCNA pipeline. Stages can be run individually
# or as a single end-to-end call via --stage all_post_denoise.
#
# Usage:
#   Rscript run_mintcna.R --config <yaml> --stage <preprocess|align|segment|all>
#
# Stages:
#   preprocess        — DNA (and RNA if multi) preprocessing → l2r CSV
#   scan              — QN + empirical null + CUSUM scan → scan_h<H>.rds  [slow]
#   call              — threshold a cached scan + TGMM → results          [fast]
#   align             — Combine denoised DNA + processed RNA into multi-omic list
#   segment           — Segmentation + CNA calling → output matrix
#   all_post_denoise  — Run align + segment (after Python denoising is done)
#   all               — Run preprocess + align + segment (does NOT run Python)

suppressPackageStartupMessages({
  library(yaml)
  library(optparse)
})

# -------------------- Argument parsing --------------------
option_list <- list(
  make_option(c("-c", "--config"),  type = "character", default = NULL,
              help = "Path to YAML config file"),
  make_option(c("-s", "--stage"),   type = "character", default = "all",
              help = "Stage: preprocess | align | segment | all_post_denoise | all"),
  make_option(c("-r", "--repo"),    type = "character", default = NULL,
              help = "Path to MintCNA_git repo root (auto-detected if NULL)"),
  # Overrides for the calling stage, so a threshold sweep is one shell command
  # off a cached scan rather than a config edit and a full re-run.
  make_option("--alpha", type = "double",    default = NA, help = "Override segmentation.alpha"),
  make_option("--st",    type = "integer",   default = NA, help = "Override calling.n_states"),
  make_option("--eps",   type = "double",    default = NA, help = "Neutral band, both sides (st=3)"),
  make_option("--eps0",  type = "double",    default = NA, help = "Neutral band, LOSS side"),
  make_option("--eps1",  type = "double",    default = NA, help = "Neutral band, GAIN side"),
  make_option("--h",     type = "integer",   default = NA, help = "Override segmentation.h (scan half-window, in bins)"),
  make_option("--smooth",     type = "double", default = NA,
              help = "Genome-direction smoothing strength (dV/dW) for the DNA/ATAC arm. 0 disables."),
  make_option("--smooth-rna", type = "double", default = NA, dest = "smooth_rna",
              help = "Same, for the RNA arm. Measured unnecessary on multiome (see MintCNA/preprocess.R)."),
  make_option("--tag",   type = "character", default = "", help = "Suffix for the call output file")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$config)) stop("--config <yaml> is required.")

# Auto-detect repo root from this script's location
if (is.null(opt$repo)) {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", args, value = TRUE)
  if (length(file_arg)) {
    script_path <- normalizePath(sub("--file=", "", file_arg))
    opt$repo <- dirname(dirname(script_path))
  } else {
    opt$repo <- normalizePath("..")
  }
}
message("[MintCNA] Repo root: ", opt$repo)

# -------------------- Source all R modules --------------------
r_files <- list.files(file.path(opt$repo, "MintCNA"),
                      pattern = "\\.R$", recursive = FALSE, full.names = TRUE)
# Source utils.R first so other files can use it
utils_first <- function(f) c(grep("utils\\.R$", f, value = TRUE),
                              grep("utils\\.R$", f, value = TRUE, invert = TRUE))
r_files <- utils_first(r_files)
invisible(lapply(r_files, source))

# -------------------- Load config --------------------
cfg <- yaml::read_yaml(opt$config)
out_dir <- cfg$output$dir
prefix  <- if (is.null(cfg$output$prefix) || !nzchar(cfg$output$prefix)) "" else paste0(cfg$output$prefix, "_")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Auto-detect mode if not set
if (is.null(cfg$mode) || is.na(cfg$mode)) {
  cfg$mode <- if (!is.null(cfg$input$rna_counts) && nzchar(cfg$input$rna_counts)) "multi" else "dna"
}
message("[MintCNA] Mode: ", cfg$mode)

# Helper paths (intermediate files)
path_l2r_dna  <- file.path(out_dir, paste0(prefix, "l2r_dna.csv"))
path_l2r_rna  <- file.path(out_dir, paste0(prefix, "l2r_rna.csv"))
path_ref_dna  <- file.path(out_dir, paste0(prefix, "ref_dna.csv"))
path_ref_rna  <- file.path(out_dir, paste0(prefix, "ref_rna.csv"))
path_denoised <- file.path(out_dir, paste0(prefix, "ae_dna.csv"))
path_combined <- file.path(out_dir, paste0(prefix, "combined.RData"))
path_results  <- file.path(out_dir, paste0(prefix, "mintcna_results.RData"))

# Two checkpoints: one after the scan (independent of alpha/st/bounds) and one
# after calling. The scan checkpoint is keyed by h, so h is resolved in one place.
resolve_h <- function() if (!is.na(opt$h)) opt$h else cfg$segmentation$h
# Smoothing strength PER MODALITY. The config may give a scalar (applies to the
# DNA/ATAC arm, the only one that needs it) or a mapping, e.g.
#   denoise:
#     smooth_genome: {DNA: 100, RNA: 0}
resolve_smooth <- function() {
  cf <- cfg$denoise$smooth_genome
  r <- if (is.null(cf)) c(DNA = 0) else if (is.list(cf)) unlist(cf) else c(DNA = as.numeric(cf))
  if (!is.na(opt$smooth))     r["DNA"] <- opt$smooth
  if (!is.na(opt$smooth_rna)) r["RNA"] <- opt$smooth_rna
  r <- r[is.finite(r) & r > 0]
  r
}
# Smoothing is applied inside the scan stage, so its strength is part of the checkpoint key.
smooth_key <- function(r) if (!length(r)) "" else
  paste0("_sm", paste(sprintf("%s%g", names(r), r), collapse = ""))
scan_file <- function(h, sm = resolve_smooth())
  file.path(out_dir, sprintf("%sscan_h%s%s.rds", prefix, h, smooth_key(sm)))
call_file <- function(tag) file.path(out_dir, sprintf("%scall%s.rds", prefix,
                                     if (nzchar(tag)) paste0("_", tag) else ""))

# -------------------- Shared pieces --------------------
# Shared by scan/call and the one-shot segment stage.

# Quantile-normal transform, per modality. sd_target is a cohort median, so this
# is joint across cells and has to happen before the scan.
qn_all <- function(combined_list) {
  if (cfg$mode == "multi") {
    sd_targets <- sapply(seq_len(nrow(combined_list[[1]])), function(r) {
      median(sapply(combined_list, function(m) sd(m[r, ], na.rm = TRUE)))
    })
    lapply(combined_list, function(m) qn_modal(m, sd_target = sd_targets))
  } else {
    sd_target <- median(sapply(combined_list, function(m) sd(m[1, ], na.rm = TRUE)))
    lapply(combined_list, function(m) qn_modal(m, sd_target = sd_target))
  }
}

# Segment-level CN, the bins x cells CN matrix, subclone clustering, saving and
# `seg` is the output of call_cells() / run_segmentation().
assemble_results <- function(seg, combined_list = NULL, ref = NULL) {
  if (is.null(ref)) ref <- seg$ref
  if (is.null(combined_list)) {
    sf <- scan_file(resolve_h())
    combined_list <- readRDS(sf)$cell_list
  }
  cna_state_matrix <- seg$cna_matrix
  cna_df_all       <- seg$cna_df_all
  converged_cells  <- seg$converged_cells
  res_null         <- seg$res_null

  # CN = round(2 * 2^mean_L2R) on the DNA row, capped at [0, 10].
  message("[MintCNA] Computing CN values per segment...")
  cna_df_all$cn <- NA_integer_
  for (i in seq_len(nrow(cna_df_all))) {
    cell <- cna_df_all$cell[i]
    if (!cell %in% names(combined_list)) next
    dna_row  <- combined_list[[cell]]["DNA", , drop = TRUE]
    m <- mean(dna_row[cna_df_all$start_idx[i]:cna_df_all$end_idx[i]], na.rm = TRUE)
    cna_df_all$cn[i] <- if (is.finite(m)) max(0L, min(10L, as.integer(round(2 * 2^m)))) else NA_integer_
  }

  cna_cn_matrix <- matrix(2L, nrow = nrow(ref), ncol = length(converged_cells),
                          dimnames = list(rownames(cna_state_matrix), converged_cells))
  for (i in seq_len(nrow(cna_df_all))) {
    cell <- cna_df_all$cell[i]
    if (!cell %in% colnames(cna_cn_matrix)) next
    if (is.na(cna_df_all$cn[i])) next
    cna_cn_matrix[cna_df_all$start_idx[i]:cna_df_all$end_idx[i], cell] <- cna_df_all$cn[i]
  }

  subclones <- cluster_info <- NULL
  if (isTRUE(cfg$clustering$enabled)) {
    latent_path <- file.path(out_dir, paste0(prefix, "latent_z.csv"))
    if (file.exists(latent_path)) {
      message("[MintCNA] Running subclone clustering on latent embeddings...")
      Z <- as.matrix(read.csv(latent_path, row.names = 1, check.names = FALSE))
      Z_conv <- Z[rownames(Z) %in% converged_cells, , drop = FALSE]
      if (nrow(Z_conv) >= 4) {
        cluster_info <- tryCatch(cluster_latent_kmeans(Z_conv),
                                 error = function(e) { message("Clustering failed: ", conditionMessage(e)); NULL })
        if (!is.null(cluster_info)) subclones <- cluster_info$clusters
      } else message("[MintCNA] Skipping clustering: too few converged cells.")
    } else message("[MintCNA] Skipping clustering: ", latent_path, " not found.")
  }

  save(cna_df_all, cna_state_matrix, cna_cn_matrix, subclones, cluster_info,
       converged_cells, res_null, ref, file = path_results)
  message("[MintCNA] Results saved: ", path_results)
  message("    cna_df_all       — ", nrow(cna_df_all), " CNA segments")
  message("    cna_state_matrix — ", nrow(cna_state_matrix), " x ", ncol(cna_state_matrix))
  message("    cna_cn_matrix    — ", nrow(cna_cn_matrix),    " x ", ncol(cna_cn_matrix))

  invisible(list(cna_df_all = cna_df_all, cna_state_matrix = cna_state_matrix,
                 cna_cn_matrix = cna_cn_matrix, subclones = subclones))
}

# -------------------- Stage: preprocess --------------------
stage_preprocess <- function() {
  message("\n[MintCNA] === Stage: preprocess ===")

  # Read normal cell labels
  cell_types <- read.table(cfg$input$cell_types, header = TRUE, sep = "\t",
                           stringsAsFactors = FALSE)
  normal_cells <- cell_types[cell_types$normal == "normal", 1]
  message("[MintCNA] ", length(normal_cells), " normal control cells.")

  # --- DNA preprocessing ---
  dna_data_path <- cfg$input$dna_counts
  if (grepl("\\.RData$|\\.Rdata$", dna_data_path)) {
    e <- new.env()
    load(dna_data_path, envir = e)
    if ("Y" %in% ls(e) && "gr" %in% ls(e)) {
      Y <- e$Y
      gr <- e$gr
      ref_dna <- data.frame(
        chr = paste0("chr", as.character(GenomicRanges::seqnames(gr))),
        start = GenomicRanges::start(gr),
        end = GenomicRanges::end(gr),
        gc = gr$gc,
        mapp = gr$mapp,
        stringsAsFactors = FALSE
      )
      ref_dna$chr <- sub("^chrchr", "chr", ref_dna$chr)
    } else {
      stop("Expected 'Y' and 'gr' objects in DNA RData file.")
    }
  } else {
    Y <- as.matrix(read.csv(dna_data_path, row.names = 1, check.names = FALSE))

    # Bin locations: from the annotation CSV if given, else parsed from the
    # count-matrix row names (e.g. "chr1:900000-1000000").
    if (!is.null(cfg$input$bin_annotation)) {
      ref_dna <- read.csv(cfg$input$bin_annotation, stringsAsFactors = FALSE)
    } else {
      ref_dna <- parse_bins_from_ids(rownames(Y))
    }

    # GC/mappability are genome properties, not sample data: MintCNA supplies
    # them from the bundled track if the user did not include them.
    if (!all(c("gc", "mapp") %in% names(ref_dna))) {
      message("[MintCNA] Attaching GC/mappability from bundled ",
              cfg$input$genome, " track (user supplied counts + bin locations only)...")
      ref_dna <- attach_gc_mapp(ref_dna, genome = cfg$input$genome,
                                ref_dir = file.path(opt$repo, "MintCNA", "ref"))
    }
  }

  # The DNA arm accepts scATAC bin counts in place of read counts; the ATAC
  # normalisation differs (see MintCNA/preprocess.R).
  modality <- if (is.null(cfg$preprocessing$modality)) "dna" else cfg$preprocessing$modality
  message("[MintCNA] DNA-arm modality: ", modality)

  dna_result <- if (modality == "atac") {
    run_atac_preprocessing(
      Y = Y, ref = ref_dna, normal_cells = normal_cells,
      gc_min = cfg$preprocessing$gc_filter_min,
      gc_max = cfg$preprocessing$gc_filter_max,
      mapp_min = cfg$preprocessing$mappability_min,
      blacklist_max = cfg$preprocessing$blacklist_max %||% 0.9,
      valid_min     = cfg$preprocessing$valid_min     %||% 0.5,
      max_zero_frac = cfg$preprocessing$max_zero_frac %||% 0.85,
      pseudocount   = cfg$preprocessing$pseudocount   %||% 1,
      keep_normal   = isTRUE(cfg$preprocessing$keep_normal))
  } else {
    run_dna_preprocessing(
      Y = Y, ref = ref_dna, normal_cells = normal_cells,
      gc_min = cfg$preprocessing$gc_filter_min,
      gc_max = cfg$preprocessing$gc_filter_max,
      mapp_min = cfg$preprocessing$mappability_min)
  }

  write.csv(dna_result$l2r, path_l2r_dna)
  write.csv(dna_result$ref, path_ref_dna, row.names = FALSE)
  message("[MintCNA] Wrote ", path_l2r_dna)

  # --- RNA preprocessing (multi-omics only) ---
  if (cfg$mode == "multi") {
    rna_counts <- as.matrix(read.table(cfg$input$rna_counts, header = TRUE,
                                       sep = "\t", row.names = 1, check.names = FALSE))
    gene_ref <- read.table(cfg$input$gene_annotation, header = TRUE,
                           sep = "\t", stringsAsFactors = FALSE)
    rownames(gene_ref) <- gene_ref$gene_id
    gene_ref <- gene_ref[rownames(rna_counts), ]

    rna_result <- preprocess_rna(
      counts = rna_counts, gene_ref = gene_ref, normal_cells = normal_cells
    )

    write.csv(rna_result$l2r, path_l2r_rna)
    write.csv(rna_result$ref, path_ref_rna, row.names = FALSE)
    message("[MintCNA] Wrote ", path_l2r_rna)
  }

  message("[MintCNA] preprocess stage complete.")
}

# -------------------- Stage: align --------------------
stage_align <- function() {
  message("\n[MintCNA] === Stage: align ===")

  # Load denoised DNA (post-AG-DCAE) — fall back to processed L2R if not present
  if (file.exists(path_denoised)) {
    dna_l2r <- as.matrix(read.csv(path_denoised, row.names = 1, check.names = FALSE))
    message("[MintCNA] Using denoised DNA: ", path_denoised)
  } else {
    dna_l2r <- as.matrix(read.csv(path_l2r_dna, row.names = 1, check.names = FALSE))
    message("[MintCNA] WARNING: denoised DNA not found, using ", path_l2r_dna)
  }
  dna_ref <- read.csv(path_ref_dna, stringsAsFactors = FALSE)

  if (cfg$mode == "dna") {
    # DNA-only: build per-cell vector list directly
    combined_list <- lapply(seq_len(ncol(dna_l2r)), function(j) {
      matrix(dna_l2r[, j], nrow = 1, dimnames = list("DNA", rownames(dna_l2r)))
    })
    names(combined_list) <- colnames(dna_l2r)
    ref <- dna_ref[, c("chr", "start", "end")]
    save(combined_list, ref, file = path_combined)
    message("[MintCNA] DNA-only combined data saved: ", path_combined)
    return(invisible())
  }

  # Multi-omics: align RNA to DNA bins
  rna_l2r <- as.matrix(read.csv(path_l2r_rna, row.names = 1, check.names = FALSE))
  rna_ref <- read.csv(path_ref_rna, stringsAsFactors = FALSE)
  if (!"gene_id" %in% names(rna_ref)) rna_ref$gene_id <- rownames(rna_l2r)

  align_result <- align_omics(dna_l2r = dna_l2r, dna_ref = dna_ref,
                              rna_l2r = rna_l2r, rna_ref = rna_ref)
  combined_list <- align_result$combined_list
  ref <- align_result$ref
  save(combined_list, ref, file = path_combined)
  message("[MintCNA] Multi-omics combined data saved: ", path_combined)
}

# -------------------- Stage: segment --------------------
stage_segment <- function() {
  message("\n[MintCNA] === Stage: segment ===")
  load(path_combined)                       # combined_list, ref
  combined_list <- qn_all(combined_list)

  seg <- run_segmentation(
    cell_list = combined_list, ref = ref, mode = cfg$mode,
    h        = cfg$segmentation$h,
    penalty  = cfg$segmentation$penalty,
    alpha    = cfg$segmentation$alpha,
    simT     = cfg$segmentation$simT %||% 1,
    seed     = cfg$segmentation$seed %||% 1,
    L        = cfg$fusion$em_iterations,
    st       = cfg$calling$n_states,
    scan_file = scan_file(resolve_h())
  )
  assemble_results(seg, combined_list, ref)
  invisible(NULL)
}

# -------------------- Stage: scan (checkpoint A) --------------------
# Everything up to and including the CUSUM scan. Expensive, and independent of
# alpha, st and the truncation bounds -- so it is done once per h and cached.
stage_scan <- function() {
  message("\n[MintCNA] === Stage: scan ===")
  load(path_combined)                       # combined_list, ref
  # Genome-direction smoothing, after AG-DCAE and before the quantile transform
  # (qn_all rescales each bin to a target SD, which smoothing would otherwise change).
  sm <- resolve_smooth()
  if (length(sm)) {
    message(sprintf("[MintCNA] genome smoothing: %s",
                    paste(sprintf("%s dV/dW=%g", names(sm), sm), collapse = ", ")))
    mc <- if (is.null(cfg$segmentation$mc_cores)) 1 else cfg$segmentation$mc_cores
    combined_list <- smooth_cell_list(combined_list, chr = ref$chr,
                                      ratios = sm, mc_cores = mc)
  }
  combined_list <- qn_all(combined_list)
  h <- resolve_h()
  scan_cells(cell_list = combined_list, ref = ref, mode = cfg$mode,
             h = h, penalty = cfg$segmentation$penalty,
             simT = cfg$segmentation$simT %||% 1,
             seed = cfg$segmentation$seed %||% 1,
             mc_cores = if (is.null(cfg$segmentation$mc_cores)) 1 else cfg$segmentation$mc_cores,
             file = scan_file(h))
  message("[MintCNA] scan complete. Now sweep thresholds off it, e.g.")
  message("    --stage call --alpha 0.05 --st 3 --eps 0.30 --tag a05")
  invisible(NULL)
}

# -------------------- Stage: call (checkpoint B) --------------------
# Threshold the cached scan and run TGMM. Seconds, so sweep freely:
#   --stage call --alpha 0.05 --st 3 --eps 0.30 --tag a05
stage_call <- function() {
  message("\n[MintCNA] === Stage: call ===")
  h  <- resolve_h()
  sf <- scan_file(h)
  if (!file.exists(sf)) stop("No scan checkpoint at ", sf, " — run --stage scan first.")

  alpha <- if (!is.na(opt$alpha)) opt$alpha else cfg$segmentation$alpha
  st    <- if (!is.na(opt$st))    opt$st    else cfg$calling$n_states
  # The neutral band may be asymmetric: --eps0 / --eps1 set the loss and gain side separately.
  eps   <- if (!is.na(opt$eps))  opt$eps  else (cfg$calling$eps %||% 0.30)
  eps0  <- if (!is.na(opt$eps0)) opt$eps0 else eps
  eps1  <- if (!is.na(opt$eps1)) opt$eps1 else eps

  sc  <- readRDS(sf)
  message(sprintf("[MintCNA] call: neutral band (-%.3f, %.3f]", eps0, eps1))
  res <- call_cells(sc, alpha = alpha, st = st, L = cfg$fusion$em_iterations,
                    eps0 = eps0, eps1 = eps1,
                    mc_cores = cfg$segmentation$mc_cores %||% 1,
                    file = call_file(opt$tag))
  assemble_results(res, sc$cell_list, sc$ref)
  invisible(NULL)
}

# -------------------- Dispatch --------------------
switch(opt$stage,
       preprocess        = stage_preprocess(),
       align             = stage_align(),
       scan              = stage_scan(),
       call              = stage_call(),
       segment           = stage_segment(),
       all_post_denoise  = { stage_align(); stage_segment() },
       all               = { stage_preprocess(); stage_align(); stage_segment() },
       stop("Unknown stage: ", opt$stage))

message("\n[MintCNA] Done.")
