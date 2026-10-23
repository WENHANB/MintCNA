#!/usr/bin/env Rscript
# Refit CONGAS+ from the saved input checkpoint under a different configuration.
#   Usage: Rscript refit_congas.R <care_id> <normal_cells TRUE|FALSE> <lambda> <tag>
#
# The input object in congas_input.rds is plain R tibbles and survives
# serialisation, so the expensive part (binning 40-100M fragments, joining gene
# coordinates) is paid once. The hyperparameters do NOT survive -- they are
# reticulate handles to Python tensors and come back as dead pointers -- so
# auto_config_run() is always re-run here, which is also where normal_cells acts.
#
# WHY normal_cells is the first thing to try. In LatentCategorical.py the flag
# only changes INITIALISATION:
#     X    = (data / sf) / nf
#     kmat = KMeans(K).fit(X.T).cluster_centers_
#     diploid_cluster = sum((round(kmat) == 2).T).argmax()  if normal_cells else -1
#     init[k,i,j] = 0.99 at CN=2      if k == diploid_cluster
#                   else high_prob at round(kmat[k,i])
# With the flag off, every cluster starts at whatever the KMeans centroid rounds
# to. The default run put essentially every segment at CN = 3 and stayed there,
# which is the signature of X being centred near 3 rather than 2. Turning the flag
# on pins one cluster at diploid and gives the fit a reference to move away from.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) stop("need <care_id> <normal_cells> <lambda> <tag> [select_segments]")
CARE   <- args[1]
NORMAL <- as.logical(args[2])
LAMBDA <- as.numeric(args[3])
TAG    <- args[4]
SELECT <- if (length(args) >= 5) as.logical(args[5]) else FALSE

suppressPackageStartupMessages({library(Rcongas); library(data.table)})
ROOT <- "<MINTCNA_ROOT>/ATAC_RNA"
setwd(ROOT)

IN  <- file.path("Result/benchmark/CONGASplus", CARE, "congas_input.rds")
stopifnot(file.exists(IN))
OUT <- file.path("Result/diagnostics/congas_variants", TAG, CARE)
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

x <- readRDS(IN)$object
cat(sprintf("[%s] %s  normal_cells=%s  lambda=%.2f\n", format(Sys.time()), CARE, NORMAL, LAMBDA))
print(Rcongas::stat(x)[c("nsegments", "ncells_RNA", "ncells_ATAC")])

t0 <- Sys.time()

# segments_selector_congas() -- "Filter segments where the optimal number of
# clusters is 1", i.e. drop the segments that are unimodal across cells. This is
# the package's own answer to the failure mode we diagnosed: with ~40 chromosome
# arms, most of them neutral, the uninformative segments dilute the few that
# carry signal. Skipping it was an omission on our side, not a defect in CONGAS+.
if (SELECT) {
  n_before <- Rcongas::stat(x)$nsegments
  x <- Rcongas::segments_selector_congas(x, K_max = 3, score = "BIC", lambda = LAMBDA)
  cat(sprintf("segments_selector_congas: %d -> %d segments\n",
              n_before, Rcongas::stat(x)$nsegments))
}

K     <- 1:4
hyper <- Rcongas::auto_config_run(x, K, normal_cells = NORMAL)
fit   <- Rcongas::fit_congas(x, K = K, lambdas = LAMBDA, learning_rate = 0.05,
                             model_parameters = hyper)
el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))

saveRDS(list(care_id = CARE, normal_cells = NORMAL, lambda = LAMBDA,
             object = x, fit = fit),
        file.path(OUT, "congas_result.rds"))

# immediate verdict: does the profile track the truth within this sample?
bf  <- fit$best_fit
cna <- as.data.table(if (!is.null(bf$CNA_real)) bf$CNA_real else bf$CNA)
ass <- as.data.table(bf$cluster_assignments)
ass[, cell := sub("-(RNA|ATAC)$", "", cell)]
w   <- unique(ass[, .(cell, cluster)])[, .(n = .N), by = cluster]
prof <- merge(cna, w, by = "cluster")[, .(cn = sum(value * n) / sum(n)), by = segment_id]
prof[, c("chr", "from", "to") := tstrsplit(segment_id, ":", fixed = TRUE)]
prof[, `:=`(from = as.numeric(from), to = as.numeric(to))]
tr <- fread(file.path("Result/truth", paste0(CARE, "_truth_100kb.tsv.gz")))
prof[, truth := vapply(seq_len(.N), function(i)
      tr[chr == prof$chr[i] & start >= prof$from[i] & end <= prof$to[i] &
         is.finite(log2), mean(log2)], numeric(1))]
ok  <- prof[is.finite(truth) & is.finite(cn)]
rho <- if (nrow(ok) > 5 && sd(ok$cn) > 0) cor(log2(ok$cn / 2), ok$truth) else NA_real_

cat(sprintf("SUMMARY\t%s\t%s\t%.2f\t%s\t%d\t%d\t%.2f\t%.3f\t%.1f\n",
            CARE, NORMAL, LAMBDA, SELECT, nrow(prof), uniqueN(cna$cluster),
            median(prof$cn, na.rm = TRUE), rho, el))
cat("  (care_id, normal_cells, lambda, select_segments, n_segments, n_clusters, median_CN, pearson, minutes)\n")
