#!/usr/bin/env Rscript
# Pseudobulk Pearson against the bulk WGS truth, for every method.
#   Usage: Rscript pearson_vs_truth.R [care_id ...]      (default: all 20)
#
# ---------------------------------------------------------------------------
# HOW THE FIELD MAKES INCOMPARABLE OUTPUTS COMPARABLE
#
# The six callers emit six different things: discrete states, linear copy
# ratios, arm-level pseudo copy numbers, multipliers of a per-cell baseline,
# clone-level integers, and residual expression. epiAneufinder's own benchmark
# solved this with two moves, and both are load-bearing:
#
#   1. AVERAGE OVER CELLS FIRST. Mapping its three states to 1/2/3 and taking
#      the mean across cells turns a discrete per-cell caller into a continuous
#      pseudobulk vector -- no invented amplitude required, the continuity comes
#      from the cell population. epiAneufinder reports Pearson on exactly this.
#   2. USE PEARSON, NOT RMSE. Pearson is invariant to affine rescaling, so
#      whatever arbitrary amplitude a state map implies cancels. epiAneufinder's
#      code computes MSE for every method pair and discards it before plotting;
#      AtaCNV states the reason outright, lifting epiAneufinder into its
#      continuous track with a declared map (gain 1.5 / neutral 1 / loss 0.5)
#      but refusing the same trick for squared error, "because it only classifies
#      genomic regions ... but does not give exact copy numbers".
#
# So: average in each method's NATIVE units, convert once to log2 ratio using
# that method's own neutral value, and report Pearson. RMSE is reported
# separately and only for methods that emit a real copy ratio.
#
# Cells: malignant only. Bulk WGS of a tumour specimen contains normal cells
# too, but Pearson is scale-invariant, so the dilution changes the amplitude of
# the truth vector and not its shape -- the comparison stays valid, and using
# tumour cells is what every precedent does.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({library(data.table)})
ROOT <- "<MINTCNA_ROOT>/ATAC_RNA"
setwd(ROOT)
CHRS <- paste0("chr", 1:22)

# EVALUATION GRID = 500 kb, MintCNA's own bin width.
#
# The truth is stored at 100 kb, but scoring there is misleading: the callers'
# native widths are 100 kb (epiAneufinder), 1 Mb (AtaCNV), 1 Mb / chromosome arm
# (Copy-scAT) and whole arms (CONGAS+), so a 100 kb grid replicates one AtaCNV
# value across 10 bins and one Copy-scAT arm across several hundred. That inflates
# n without adding information and lets the correlation be driven by how many bins
# a coarse call happens to span.
#
# Measured across three samples, moving 100 kb -> 500 kb leaves the RANKING
# unchanged, lifts every method by ~0.02-0.03 (bin noise averaged away), and
# shrinks the own-subset vs common-subset discrepancy by ~35%. 500 kb also means
# MintCNA enters this table with no projection at all.
GRID <- 5e5

# Which configuration of each method is scored. Every one of these is the
# author-default / as-published setting; the alternatives exist on disk and are
# reported separately, never in the headline table.
EPI_TAG      <- Sys.getenv("EPI_TAG",      "w100kb_mf20000")   # as published
RIDDLER_DIR  <- Sys.getenv("RIDDLER_DIR",  "1MB")              # frag_min = 10000
INFERCNV_TAG <- Sys.getenv("INFERCNV_TAG", "default")          # cutoff 1, denoise FALSE
# The configuration the paper reports. h = 15 beat h = 10 and h = 5 on all three
# samples with matched WGS (mean Pearson 0.814 / 0.805 / 0.740, no reversal), with
# alpha = 0.10 and eps = 0.20. Chosen ON those three samples, so it is a choice
# made on the evaluation set -- fix it here and report the other 17 without
# re-tuning, or the number is optimistic.
CHOSEN_MINT  <- "h15a005e020"
MINT_TAG     <- Sys.getenv("MINT_TAG",     CHOSEN_MINT)

# Which Malignant nuclei every method is scored on.
#   own    : each method on the Malignant nuclei IT retained. What the field does,
#            but the retained fractions differ a lot -- epiAneufinder at its own
#            default keeps 51% of them, and those are the DEEPEST 51% (median
#            32,569 fragments vs 20,580 for the full Malignant set), so the
#            methods are not being graded on the same cells.
#   common : the intersection across all methods. Fewer cells, identical exam.
CELLSET <- Sys.getenv("CELLSET", "own")
stopifnot(CELLSET %in% c("own", "common"))

args <- commandArgs(trailingOnly = TRUE)
SAMPLES <- if (length(args)) args else
  c("P59T2","P59T3","P104T1","P104T2","P106T1","P106T2","P100T1","P100T3",
    "P103T1","P103T2","P105T1","P105T3","P108T1","P108T2","P109T1","P109T2",
    "P110T1","P110T2","P111T2","P111T3")

# Project a method's (chr, start, end, value) profile onto the truth's 100 kb
# grid. Methods are coarser than the truth, so each of their bins covers many
# truth bins; assign by midpoint containment.
project <- function(truth, chr, from, to, val) {
  ok <- is.finite(val) & chr %in% CHRS
  if (!any(ok)) return(rep(NA_real_, nrow(truth)))
  # `from` is the roll key, so after the join it holds the QUERY position and the
  # method bin's own start is gone. Keep a copy under another name first.
  d <- data.table(chr = as.character(chr[ok]), from = as.numeric(from[ok]),
                  bin_from = as.numeric(from[ok]), bin_to = as.numeric(to[ok]),
                  val = as.numeric(val[ok]))
  setkey(d, chr, from)
  q <- data.table(chr = as.character(truth$chr),
                  from = (as.numeric(truth$start) + as.numeric(truth$end)) / 2,
                  i = seq_len(nrow(truth)))
  setkey(q, chr, from)
  # Rolling join: for each truth midpoint take the last method bin starting at or
  # before it, then require the midpoint to actually fall inside that bin -- a
  # plain roll would otherwise carry the last bin of a chromosome across its
  # centromere gap and past its end.
  j <- d[q, roll = TRUE]
  j[is.na(bin_to) | from < bin_from | from > bin_to, val := NA_real_]
  out <- rep(NA_real_, nrow(truth))
  out[j$i] <- j$val
  out
}

# --- per-method pseudobulk in log2 ratio -----------------------------------
get_epi <- function(s, truth, keep) {
  # AS PUBLISHED: windowSize 1e5 (README/vignette) and minFrags 20000 (package
  # default). The w1mb_mf5000 directory is a tuned configuration found by a
  # parameter search and is NOT what this table reports -- benchmarking a method
  # on settings we chose for it is not benchmarking the method.
  f <- file.path("Result/benchmark/epiAneufinder", EPI_TAG, s,
                 "epiAneufinder_results/results_table.tsv")
  if (!file.exists(f)) return(NULL)
  ef <- fread(f); cc <- grep("^cell-", names(ef), value = TRUE)
  bc <- sub("^cell-", "", cc); use <- cc[bc %in% keep]
  if (!length(use)) return(NULL)
  # states 0/1/2 -> copy number 1/2/3, epiAneufinder's own map; mean over cells
  m <- rowMeans(as.matrix(ef[, ..use]) + 1)
  project(truth, ef$seq, ef$start, ef$end, log2(m / 2))
}
get_atacnv <- function(s, truth, keep) {
  f <- file.path("Result/benchmark/AtaCNV", s, "atacnv_result.rds")
  if (!file.exists(f)) return(NULL)
  A <- readRDS(f)$cnv$copy_ratio
  i <- rownames(A) %in% keep; if (!any(i)) return(NULL)
  m  <- colMeans(A[i, , drop = FALSE])              # neutral = 1
  cn <- do.call(rbind, strsplit(colnames(A), "_"))
  project(truth, cn[, 1], as.numeric(cn[, 2]), as.numeric(cn[, 3]), log2(m))
}
get_copyscat <- function(s, truth, keep, arms) {
  f <- file.path("Result/benchmark/CopyscAT", s, paste0(s, "clean_cnv_cnv_scores.csv"))
  if (!file.exists(f)) return(NULL)
  cs <- read.csv(f, row.names = 1, check.names = FALSE)
  i <- rownames(cs) %in% keep; if (!any(i)) return(NULL)
  # arms it did not call were tested and found neutral -> 2, not missing
  tested <- as.character(readRDS(file.path("Result/benchmark/CopyscAT", s,
                                           "copyscat_result.rds"))$collapse[[1]])
  m <- setNames(rep(2, length(tested)), tested)
  cm <- colMeans(cs[i, , drop = FALSE])
  m[names(cm)] <- cm                                # neutral = 2
  a <- arms[label %in% names(m)]
  project(truth, a$chr, a$from, a$to, log2(m[a$label] / 2))
}
get_riddler <- function(s, truth, keep) {
  f <- file.path("Result/benchmark/RIDDLER", RIDDLER_DIR, s, "riddler_result.rds")
  if (!file.exists(f)) return(NULL)
  rd <- readRDS(f)$riddler
  # riddler_fit() reads its matrix with read.table(header = TRUE), which turns
  # barcodes into syntactic names: "AAACAGCC-1" comes back as "AAACAGCC.1".
  # Matching against the whitelist without undoing that silently selects no cells.
  cn_bc <- sub("\\.(\\d+)$", "-\\1", colnames(rd$CNV))
  i <- cn_bc %in% keep; if (!any(i)) return(NULL)
  m  <- rowMeans(rd$CNV[, i, drop = FALSE])         # neutral = 1, a multiplier
  w  <- rd$window; mi <- match(rownames(rd$CNV), as.character(w[[4]]))
  project(truth, as.character(w[[1]])[mi], w[[2]][mi], w[[3]][mi], log2(m))
}
get_congas <- function(s, truth, keep) {
  f <- file.path("Result/benchmark/CONGASplus", s, "congas_result.rds")
  if (!file.exists(f)) return(NULL)
  bf <- readRDS(f)$fit$best_fit
  cna <- as.data.table(if (!is.null(bf$CNA_real)) bf$CNA_real else bf$CNA)
  ass <- as.data.table(bf$cluster_assignments)
  ass[, cell := sub("-(RNA|ATAC)$", "", cell)]
  ass <- unique(ass[, .(cell, cluster)])[cell %in% keep]
  if (!nrow(ass)) return(NULL)
  w <- ass[, .(n = .N), by = cluster]
  p <- merge(cna, w, by = "cluster")[, .(cn = sum(value * n) / sum(n)), by = segment_id]
  p[, c("c", "f", "t") := tstrsplit(segment_id, ":", fixed = TRUE)]
  project(truth, p$c, as.numeric(p$f), as.numeric(p$t), log2(p$cn / 2))   # neutral = 2
}
# The ATAC-only arm of the ablation, entering the table as its own method so it
# can be read against the six benchmarks rather than only against its own
# multi-omic counterpart. eps is 0.25 there, not 0.20: the fused segment mean is
# wider without RNA averaged into it, and 0.25 is the value at which the two arms
# call the same share of the genome (25.5% against 24.5%). See
# Notes/parameter_sweeps.md 3c.
MINT_ATAC_TAG <- Sys.getenv("MINT_ATAC_TAG", "h15a005e025")
get_mintcna_atac <- function(s, truth, keep)
  .mint_track(sprintf("Result/MintCNA/%s_atac/%s_results_%s.RData", s, s, MINT_ATAC_TAG),
              truth, keep)

get_mintcna <- function(s, truth, keep)
  .mint_track(sprintf("Result/MintCNA/%s_multi/%s_results_%s.RData", s, s, MINT_TAG),
              truth, keep)

.mint_track <- function(f, truth, keep) {
  if (!file.exists(f)) return(NULL)
  e <- new.env(); load(f, e)
  Sm <- get("cna_state_matrix", e); rf <- get("ref", e)
  i <- colnames(Sm) %in% keep; if (!any(i)) return(NULL)
  # discrete states -> continuous by averaging the codes, then ONE log2 through
  # the method's own neutral code. Same convention as epiAneufinder above.
  neutral <- if (grepl("st5", basename(f))) 3 else 2
  m <- rowMeans(Sm[, i, drop = FALSE])
  project(truth, rf$chr, rf$start, rf$end, log2(m / neutral))
}
get_infercnv <- function(s, truth, keep) {
  f <- file.path("Result/benchmark/inferCNV", INFERCNV_TAG, s, "infercnv_export.rds")
  if (!file.exists(f)) return(NULL)
  iv <- readRDS(f)
  i <- colnames(iv$expr) %in% keep; if (!any(i)) return(NULL)
  m  <- rowMeans(iv$expr[, i, drop = FALSE])        # neutral = 1
  go <- iv$gene_order
  project(truth, as.character(go$chr), go$start, go$stop, log2(m))
}

# Which barcodes each method actually reports. Only read when CELLSET = "common".
# Kept next to the getters above so the two cannot describe different cell sets.
cells_of <- list(
  epiAneufinder = function(s) {
    f <- file.path("Result/benchmark/epiAneufinder", EPI_TAG, s,
                   "epiAneufinder_results/results_table.tsv")
    if (!file.exists(f)) return(NULL)
    sub("^cell-", "", grep("^cell-", names(fread(f, nrows = 0)), value = TRUE))
  },
  AtaCNV = function(s) {
    f <- file.path("Result/benchmark/AtaCNV", s, "atacnv_result.rds")
    if (file.exists(f)) rownames(readRDS(f)$cnv$copy_ratio)
  },
  `Copy-scAT` = function(s) {
    f <- file.path("Result/benchmark/CopyscAT", s, paste0(s, "clean_cnv_cnv_scores.csv"))
    if (file.exists(f)) rownames(read.csv(f, row.names = 1, check.names = FALSE))
  },
  RIDDLER = function(s) {
    f <- file.path("Result/benchmark/RIDDLER", RIDDLER_DIR, s, "riddler_result.rds")
    if (file.exists(f)) sub("\\.(\\d+)$", "-\\1", colnames(readRDS(f)$riddler$CNV))
  },
  `CONGAS+` = function(s) {
    f <- file.path("Result/benchmark/CONGASplus", s, "congas_result.rds")
    if (file.exists(f))
      unique(sub("-(RNA|ATAC)$", "",
                 as.data.table(readRDS(f)$fit$best_fit$cluster_assignments)$cell))
  },
  inferCNV = function(s) {
    f <- file.path("Result/benchmark/inferCNV", INFERCNV_TAG, s, "infercnv_export.rds")
    if (!file.exists(f)) return(NULL)
    iv <- readRDS(f); intersect(iv$obs_cells, colnames(iv$expr))
  },
  MintCNA = function(s) {
    f <- file.path(sprintf("Result/MintCNA/%s_multi", s),
                   sprintf("%s_results_%s.RData", s, MINT_TAG))
    if (!file.exists(f)) return(NULL)
    e <- new.env(); load(f, e); colnames(get("cna_state_matrix", e))
  },
  `MintCNA(ATAC)` = function(s) {
    f <- sprintf("Result/MintCNA/%s_atac/%s_results_%s.RData", s, s, MINT_ATAC_TAG)
    if (!file.exists(f)) return(NULL)
    e <- new.env(); load(f, e); colnames(get("cna_state_matrix", e))
  })

arms <- as.data.table(readRDS("Result/benchmark/AtaCNV/P105T1/atacnv_result.rds")$bins)
arms <- arms[chr %in% CHRS & arm %in% c("p", "q"),
             .(from = min(start), to = max(end)), by = .(chr, arm)]
arms[, label := paste0(chr, arm)]

res <- rbindlist(lapply(SAMPLES, function(s) {
  tf <- file.path("Result/truth", paste0(s, "_truth_100kb.tsv.gz"))
  if (!file.exists(tf)) return(NULL)
  truth <- fread(tf)[chr %in% CHRS]
  cells <- fread(file.path("Data/cell_lists", paste0(s, "_cells.tsv")))
  keep  <- cells[cell_type == "Malignant", barcode]
  n_mal <- length(keep)
  if (CELLSET == "common") {
    cov <- lapply(cells_of, function(f) f(s))
    cov <- cov[!vapply(cov, is.null, logical(1))]
    keep <- Reduce(intersect, c(list(keep), cov))
    cat(sprintf("%s: %d of %d Malignant nuclei common to %d methods\n",
                s, length(keep), n_mal, length(cov)))
    if (length(keep) < 50) { cat("  too few; skipping\n"); return(NULL) }
  }

  P <- list(epiAneufinder = get_epi(s, truth, keep),
            AtaCNV        = get_atacnv(s, truth, keep),
            `Copy-scAT`   = get_copyscat(s, truth, keep, arms),
            RIDDLER       = get_riddler(s, truth, keep),
            `CONGAS+`     = get_congas(s, truth, keep),
            inferCNV      = get_infercnv(s, truth, keep),
            MintCNA       = get_mintcna(s, truth, keep),
            `MintCNA(ATAC)` = get_mintcna_atac(s, truth, keep))

  # Keep the projected pseudobulk tracks as well as the correlation: they are the
  # AtaCNA Fig 3C-F view (one mean-over-cells profile per method, stacked against
  # the bulk truth), and re-deriving them elsewhere would risk a second, subtly
  # different projection.
  prof <- rbindlist(lapply(names(P), function(nm) {
    v <- P[[nm]]
    if (is.null(v)) return(NULL)
    data.table(care_id = s, method = nm, chr = truth$chr,
               start = truth$start, end = truth$end, log2 = v)
  }))
  if (nrow(prof)) {
    prof <- rbind(prof, data.table(care_id = s, method = "bulk WGS", chr = truth$chr,
                                   start = truth$start, end = truth$end, log2 = truth$log2))
    dir.create("Result/qc/profiles", showWarnings = FALSE, recursive = TRUE)
    fwrite(prof, sprintf("Result/qc/profiles/%s_profiles_100kb.tsv.gz", s), sep = "\t")
  }

  # --- score on the GRID, three ways ------------------------------------------
  # own      : each method on the bins IT reports -- what the field does, but the
  #            methods are then graded on different exams (epiAneufinder covers
  #            80% of the genome, RIDDLER 100%), so it rewards aggressive filtering
  # common   : all scATAC methods on the identical bin set -- the fair comparison
  # coverage : the fraction each method reports on, which is itself a result and
  #            must not be allowed to hide inside `own`
  ATAC <- c("epiAneufinder", "AtaCNV", "Copy-scAT", "RIDDLER", "CONGAS+",
            "MintCNA", "MintCNA(ATAC)")
  agg <- rbindlist(c(
    lapply(names(P), function(nm) if (is.null(P[[nm]])) NULL else
      data.table(method = nm, chr = truth$chr, gb = truth$start %/% GRID, v = P[[nm]])),
    list(data.table(method = "truth", chr = truth$chr,
                    gb = truth$start %/% GRID, v = truth$log2))))
  agg <- agg[, .(v = mean(v, na.rm = TRUE)), by = .(method, chr, gb)]
  agg[is.nan(v), v := NA_real_]
  W <- dcast(agg, chr + gb ~ method, value.var = "v")
  tv <- W[["truth"]]
  mm <- setdiff(names(W), c("chr", "gb", "truth"))
  ok_atac <- intersect(ATAC, mm)
  kA <- if (length(ok_atac)) Reduce(`&`, c(lapply(ok_atac, function(m) is.finite(W[[m]])),
                                           list(is.finite(tv)))) else rep(FALSE, nrow(W))
  rbindlist(lapply(mm, function(nm) {
    v <- W[[nm]]; o <- is.finite(v) & is.finite(tv)
    data.table(care_id = s, method = nm, cellset = CELLSET,
               n_cell = length(keep), n_malignant = n_mal,
               n_bin = sum(o), n_common = sum(kA), grid_kb = GRID / 1e3,
               coverage = mean(is.finite(v)),
               r = if (sum(o) > 50 && sd(v[o]) > 0) cor(v[o], tv[o]) else NA_real_,
               r_common = if (nm %in% ok_atac && sum(kA) > 50 && sd(v[kA]) > 0)
                            cor(v[kA], tv[kA]) else NA_real_)
  }))
}))


# The MintCNA tag goes in the filename for any non-default configuration. A
# config sweep used to loop MINT_TAG over this script and silently overwrite the
# headline table with whichever tag ran last, so the benchmark rows and the
# MintCNA row came from different runs.
OUTF <- if (identical(MINT_TAG, CHOSEN_MINT))
          sprintf("Result/qc/pearson_vs_truth_%s.tsv", CELLSET) else
          sprintf("Result/qc/pearson_vs_truth_%s_mint-%s.tsv", CELLSET, MINT_TAG)
fwrite(res, OUTF, sep = "\t")
cat("wrote", OUTF, "\n")

w <- dcast(res, care_id ~ method, value.var = "r_common")
man <- fread("Data/meta/sample_manifest.csv"); setkey(man, care_id)
w[, subtype := substr(man[care_id, subtype], 1, 12)]
setcolorder(w, c("care_id", "subtype"))
cat(sprintf("=== Pearson vs bulk WGS, %g kb grid, COMMON scATAC bins, per sample ===\n\n", GRID/1e3))
print(w[order(subtype, care_id)], nrow = 30)

cat("\n=== per method ===\n")
print(res[!is.na(r), .(n = .N, median = round(median(r), 3), mean = round(mean(r), 3),
                       min = round(min(r), 3), max = round(max(r), 3),
                       `>0.5` = sum(r > 0.5)), by = method][order(-median)])
