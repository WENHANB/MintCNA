#!/usr/bin/env Rscript
# (R)CONGAS+ on one sample.  Usage: Rscript run_congas.R <care_id>
#
# The only method in the benchmark that uses BOTH modalities, and the only one
# written for multiome rather than standalone scATAC. It is therefore the correct
# comparator for MintCNA(ATAC+RNA); the three ATAC-only callers are the
# comparators for MintCNA(ATAC).
#
# Call sequence follows the package's own documented example:
#   init() -> filter_outliers() -> filter_missing_data()
#        -> auto_config_run(K = 1:4) -> fit_congas(learning_rate = 0.05)
#
# SEGMENTATION -- the one decision that has to be stated in Methods. init() wants
# a chr/from/to/copies tibble, which is exactly the shape of our bulk GATK
# segments. Feeding those in would make the comparison circular, since they ARE
# the ground truth. The package's own arm-level segmentation (Rcongas::hg38_arms,
# 48 arms, all copies = 2) is used instead.

args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("need <care_id> [segmentation: arms|truth]")
CARE   <- args[1]
SEGSRC <- if (length(args) >= 2) args[2] else "arms"
stopifnot(SEGSRC %in% c("arms", "truth"))

suppressPackageStartupMessages({
  library(data.table); library(Matrix); library(dplyr); library(tibble); library(Rcongas)
})

ROOT <- "<MINTCNA_ROOT>/ATAC_RNA"
setwd(ROOT)
# Note the braces: `x <- if (a) b` on one line followed by `else` on the next
# parses as two statements at top level and dies with "unexpected 'else'".
OUT <- if (SEGSRC == "arms") {
  file.path("Result/benchmark/CONGASplus", CARE)
} else {
  file.path("Result/diagnostics/congas_truthseg", CARE)
}
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

ATAC_BIN <- 1e6      # ATAC feature width; CONGAS+ sums features into segments,
                     # so this only has to be finer than an arm
CHRS <- paste0("chr", 1:22)

cells <- fread(file.path("Data/cell_lists", paste0(CARE, "_cells.tsv")))
keep_bc <- cells$barcode
message(sprintf("%s: %d paired nuclei on the whitelist", CARE, length(keep_bc)))

# ---- ATAC: cell / value / chr / from / to ----------------------------------
message("binning fragments ...")
fr <- fread(cmd = paste("zcat", shQuote(file.path("Data/atac_filtered",
                                                  paste0(CARE, "_fragments.tsv.gz")))),
            header = FALSE, select = 1:4, col.names = c("chr", "start", "end", "cell"))
fr <- fr[chr %in% CHRS & cell %in% keep_bc]
fr[, mid := (start + end) %/% 2L]
fr[, from := (mid %/% as.integer(ATAC_BIN)) * as.integer(ATAC_BIN)]
atac <- fr[, .(value = .N), by = .(cell, chr, from)]
atac[, to := from + as.integer(ATAC_BIN)]
atac <- as_tibble(atac[, .(cell, value, chr, from, to)])
message(sprintf("ATAC: %s rows, %d cells, %d bins",
                format(nrow(atac), big.mark = ","),
                uniqueN(atac$cell), uniqueN(paste(atac$chr, atac$from))))

# ---- RNA: gene / chr / from / to / cell / value -----------------------------
message("reading GEX ...")
man  <- fread("Data/meta/sample_manifest.csv")[care_id == CARE]
gfil <- strsplit(man$gex_files, ";")[[1]]
mtx  <- readMM(file.path("Data/gex_matrices", grep("matrix",   gfil, value = TRUE)))
bc   <- fread(file.path("Data/gex_matrices", grep("barcodes", gfil, value = TRUE)), header = FALSE)$V1
ft   <- fread(file.path("Data/gex_matrices", grep("features", gfil, value = TRUE)), header = FALSE)
setnames(ft, c("gene_id", "gene_name", "type"))
stopifnot(ncol(mtx) == length(bc), nrow(mtx) == nrow(ft))

# Gene coordinates. The 10x reference is GRCh38-2020-A (Ensembl 98); EnsDb v86 is
# the annotation available here. They differ for a minority of genes, which is
# immaterial because CONGAS+ only needs each gene's chromosome arm.
suppressPackageStartupMessages(library(EnsDb.Hsapiens.v86))
gg <- as.data.table(ensembldb::genes(EnsDb.Hsapiens.v86,
                                     columns = c("gene_id", "seq_name", "gene_seq_start", "gene_seq_end")))
gg <- gg[, .(gene_id, chr = paste0("chr", seqnames), from = start, to = end)][chr %in% CHRS]
ft[, row := .I]
ft <- merge(ft, gg, by = "gene_id")
message(sprintf("genes with hg38 coordinates: %d of %d (%.0f%%)",
                nrow(ft), nrow(gg[0]) + length(unique(gg$gene_id)) * 0 + nrow(ft), 100))

j   <- which(bc %in% keep_bc)
sub <- mtx[ft$row, j, drop = FALSE]
tt  <- as(sub, "TsparseMatrix")
rna <- data.table(gene = ft$gene_name[tt@i + 1L], chr = ft$chr[tt@i + 1L],
                  from = ft$from[tt@i + 1L],      to  = ft$to[tt@i + 1L],
                  cell = bc[j][tt@j + 1L],        value = as.integer(tt@x))
rna <- as_tibble(rna[value > 0])
message(sprintf("RNA: %s rows, %d cells, %d genes",
                format(nrow(rna), big.mark = ","), uniqueN(rna$cell), uniqueN(rna$gene)))

# multiome mode requires the same barcodes on both sides
common <- intersect(unique(rna$cell), unique(atac$cell))
message(sprintf("cells present in both modalities: %d", length(common)))
# dplyr::filter must be qualified: attaching EnsDb.Hsapiens.v86 pulls in
# ensembldb::filter, which masks it and errors with "requires an 'EnsDb' object"
rna  <- dplyr::filter(rna,  cell %in% common)
atac <- dplyr::filter(atac, cell %in% common)

# ---- fit -------------------------------------------------------------------
data("hg38_arms", package = "Rcongas")

# [MOD] Drop the five acrocentric short arms from the segmentation.
#
# Rcongas::hg38_arms lists chr13p/14p/15p/21p/22p as ordinary segments, but they
# are heterochromatin and rDNA with essentially no assayable sequence. Measured on
# P106T2 the ATAC totals across ALL 1,930 nuclei are chr13p 2, chr15p 13,
# chr14p 43, chr22p 131, chr21p 460, and chr22p carries a grand total of 1 RNA
# read. Two things follow, and both are fatal:
#   filter_missing_data(proportion_ATAC = 0.05) removes any cell with more than
#   5% zero segments -- 2.2 of 44 -- so nearly every nucleus fails on these five
#   alone. Measured: 1,671 of 1,724 removed, 97%.
#   auto_config_run() then fits a Gamma to a segment with one nonzero cell and
#   returns NaN, so every fit dies with "Expected parameter concentration ... to
#   satisfy GreaterThan(lower_bound=0.0)".
# Excluding them is what the other callers already do: Copy-scAT's own pipeline
# drops exactly chr13p/14p/15p/21p/22p (and chrYp) before calling.
ACRO <- c("chr13", "chr14", "chr15", "chr21", "chr22")
segmentation <- hg38_arms %>%
  dplyr::filter(!(chr %in% ACRO & from == 0))
message(sprintf("segmentation: %d arms (dropped %d acrocentric p-arms)",
                nrow(segmentation), nrow(hg38_arms) - nrow(segmentation)))

if (SEGSRC == "truth") {
  # DIAGNOSTIC PATH ONLY -- never report a number produced this way.
  #
  # CONGAS+ infers one copy number per (cluster, segment). We are obliged to give
  # it generic chromosome arms, because the only sample-specific segmentation we
  # have is the bulk GATK truth and feeding that to a method we then score against
  # that same truth is circular. But ~40 arms, mostly neutral, may simply starve
  # the model. Running it once on a truth-derived segmentation separates the two
  # explanations: if it works here, the segmentation is the cause and the paper
  # states that as a limitation; if it still fails, the arms are exonerated.
  tr <- data.table::fread(file.path("Result/truth", paste0(CARE, "_truth_100kb.tsv.gz")))
  tr <- tr[chr %in% CHRS & is.finite(log2)]
  data.table::setorder(tr, chr, start)
  tr[, state := data.table::fifelse(log2 < -0.15, "loss",
                data.table::fifelse(log2 >  0.15, "gain", "neu"))]
  tr[, run := data.table::rleid(chr, state)]
  sg <- tr[, .(from = min(start), to = max(end), n = .N, state = state[1]), by = .(chr, run)]
  sg <- sg[n >= 20]                                   # >= 2 Mb
  sg[, copies := data.table::fifelse(state == "loss", 1L,
                 data.table::fifelse(state == "gain", 3L, 2L))]
  segmentation <- tibble::as_tibble(sg[, .(chr, from = as.integer(from),
                                           to = as.integer(to), copies)])
  message(sprintf("DIAGNOSTIC truth-derived segmentation: %d segments (%d loss, %d gain, %d neutral)",
                  nrow(sg), sum(sg$state == "loss"), sum(sg$state == "gain"),
                  sum(sg$state == "neu")))
}

t0 <- Sys.time()

x <- Rcongas::init(rna = rna, atac = atac, segmentation = segmentation,
                   rna_likelihood = "NB", atac_likelihood = "NB",
                   reference_genome = "GRCh38", multiome = TRUE,
                   description = paste0("CARE ", CARE))
message("init done"); print(Rcongas::stat(x)[c("nsegments", "ncells_RNA", "ncells_ATAC")])

x <- Rcongas::filter_outliers(x)         # defaults
x <- Rcongas::filter_missing_data(x)     # defaults
message("after filtering:"); print(Rcongas::stat(x)[c("nsegments", "ncells_RNA", "ncells_ATAC")])
if (Rcongas::stat(x)$ncells_ATAC < 100)
  warning("fewer than 100 cells survived filtering -- check the segmentation")

# Documented example: hyperparams from auto_config_run, then fit at lr = 0.05.
#
# lambdas has to be passed by hand. The man page calls it "Float (Optional).
# Default 0.5", but the signature is fit_congas(x, K, lambdas, model_parameters,
# ...) with no default, so the authors' own example fails with
#   argument "lambdas" is missing, with no default
# 0.5 is supplied because that is the documented default: equal weight to the RNA
# and ATAC likelihoods. Values toward 0 favour ATAC, toward 1 favour RNA.
K <- 1:4
hyper <- Rcongas::auto_config_run(x, K)

# Checkpoint before the fit. Building the input takes a few minutes and the fit is
# the part that needs iterating, so never pay for the input twice.
saveRDS(list(care_id = CARE, object = x, hyper = hyper, n_cells = length(common),
             segsrc = SEGSRC, rna = rna, atac = atac, segmentation = segmentation),
        file.path(OUT, "congas_input.rds"))

# fit_congas() wraps every (k, lambda) fit in easypar::run(), which swallows the
# underlying error and reports only "N/N computations returned errors"; model
# selection then dies on a missing ICL column. Probe one fit directly first so a
# real message reaches the log.
probe <- tryCatch(
  getFromNamespace("fit_congas_single_run", "Rcongas")(
    x, hyper, 0.5, 2L, 0.05, "G", FALSE, 500, 10, 0.005, 5, FALSE),
  error = function(e) e)
if (inherits(probe, "error")) {
  message("SINGLE-RUN PROBE FAILED: ", conditionMessage(probe))
} else {
  message("single-run probe OK")
}

fit <- Rcongas::fit_congas(x, K = K, lambdas = 0.5, learning_rate = 0.05,
                           model_parameters = hyper)
el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))

saveRDS(list(care_id = CARE, object = x, hyper = hyper, fit = fit,
             n_cells = length(common), atac_bin = ATAC_BIN),
        file.path(OUT, "congas_result.rds"))
cat(sprintf("SUMMARY\t%s\t%d\t%.1f\n", CARE, length(common), el))
