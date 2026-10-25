#!/usr/bin/env Rscript
# epiAneufinder (1.1.5) on one multiome sample. The paper reports the published settings
# (windowSize 1e5, minFrags 20000, tag w100kb_mf20000); windowSize and minFrags are arguments because the
# package gives no single recommendation for 10x Multiome depth.
#   Rscript reproduce/analysis/multiome/methods/run_epiAneufinder.R <care_id> <windowSize> <minFrags> <tag> [k = 4]
#   e.g.    run_epiAneufinder.R P111T2 1e5 20000 w100kb_mf20000
# Input : ATAC_RNA/Data/atac_filtered/<care_id>_fragments.tsv.gz, ATAC_RNA/Data/ref/hg38-blacklist.v2.bed
# Output: ATAC_RNA/Result/benchmark/epiAneufinder/<tag>/<care_id>/epiAneufinder_results/results_table.tsv
# Run from the root of the MintCNA repository with MINTCNA_DATA set (see reproduce/README.md).
suppressPackageStartupMessages(library(epiAneufinder))
ROOT <- Sys.getenv("MINTCNA_DATA"); if (!nzchar(ROOT)) stop("set MINTCNA_DATA (see reproduce/README.md)")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) stop("need <care_id> <windowSize> <minFrags> <tag> [k]")
CARE <- args[1]; WINDOW <- as.numeric(args[2]); MINFRAGS <- as.integer(args[3]); TAG <- args[4]
KBP  <- if (length(args) >= 5) as.integer(args[5]) else 4L
ROOT <- file.path(ROOT, "ATAC_RNA")
IN  <- file.path(ROOT, "Data/atac_filtered", paste0(CARE, "_fragments.tsv.gz"))
BL  <- file.path(ROOT, "Data/ref/hg38-blacklist.v2.bed")
OUT <- file.path(ROOT, "Result/benchmark/epiAneufinder", TAG, CARE); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
stopifnot(file.exists(IN), file.exists(BL))

t0 <- Sys.time()
epiAneufinder(input = IN, outdir = OUT, blacklist = BL, windowSize = WINDOW,
              genome = "BSgenome.Hsapiens.UCSC.hg38", exclude = c("chrX", "chrY", "chrM"),
              reuse.existing = FALSE, title_karyo = CARE, ncores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4")),
              minFrags = MINFRAGS, minsizeCNV = 0, k = KBP, plotKaryo = TRUE)

res <- file.path(OUT, "epiAneufinder_results", "results_table.tsv")     # space separated
if (file.exists(res)) {
  h <- strsplit(readLines(res, n = 1), " ")[[1]]
  cat(sprintf("SUMMARY\t%s\t%s\t%d cells\t%d bins\t%.1f min\n", CARE, TAG, sum(grepl("^cell-", h)),
              length(readLines(res)) - 1L, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
} else cat(sprintf("SUMMARY\t%s\t%s\tno results table\n", CARE, TAG))
