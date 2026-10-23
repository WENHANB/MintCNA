#!/usr/bin/env Rscript
# epiAneufinder on one sample.
#   Usage: Rscript run_epiAneufinder.R <care_id> <windowSize> <minFrags> <tag>
#
# Two parameters are swept, because the package's own guidance does not resolve
# to a single setting on 10x Multiome depth.
#
#  windowSize  has NO DEFAULT in the function signature -- it is required. Every
#    piece of documentation (README, vignette, CLI wrapper) passes 1e5, but the
#    roxygen on the parameter itself reads "Reccomended for sparse data - 1e6".
#    This cohort is sparse by the authors' own standard: at 100 kb, lambda = 0.45
#    fragments per bin per nucleus and 64% of bins are empty. Running only 1e5
#    would test epiAneufinder in a regime its authors flag as unsuitable.
#
#  minFrags  defaults to 20000, tuned for standalone scATAC-seq. Multiome splits
#    the sequencing budget across two libraries from the same nuclei, so that gate
#    keeps just 23.4% of this cohort (16,257 of 69,384), and drops 4 of 20 samples
#    below 300 nuclei. 5000 is depth-appropriate.
#
# Everything else is README verbatim.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) stop("need <care_id> <windowSize> <minFrags> <tag> [k]")
CARE <- args[1]; WINDOW <- as.numeric(args[2]); MINFRAGS <- as.integer(args[3]); TAG <- args[4]
KBP  <- if (length(args) >= 5) as.integer(args[5]) else 4L

suppressPackageStartupMessages(library(epiAneufinder))

ROOT <- "<MINTCNA_ROOT>/ATAC_RNA"
IN   <- file.path(ROOT, "Data/atac_filtered", paste0(CARE, "_fragments.tsv.gz"))
OUT  <- file.path(ROOT, "Result/benchmark/epiAneufinder", TAG, CARE)
BL   <- file.path(ROOT, "Data/ref/hg38-blacklist.v2.bed")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

stopifnot(file.exists(IN), file.exists(BL))
NCORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))

cat(sprintf("[%s] %s  windowSize=%.0f  minFrags=%d  k=%d  ncores=%d\n",
            format(Sys.time()), CARE, WINDOW, MINFRAGS, KBP, NCORES))
t0 <- Sys.time()

epiAneufinder(
  input        = IN,
  outdir       = OUT,
  blacklist    = BL,
  windowSize   = WINDOW,
  genome       = "BSgenome.Hsapiens.UCSC.hg38",
  exclude      = c("chrX", "chrY", "chrM"), # README
  reuse.existing = FALSE,
  title_karyo  = CARE,
  ncores       = NCORES,
  minFrags     = MINFRAGS,
  minsizeCNV   = 0,                         # README default
  # k = breakpoints sought per chromosome per cell, so segments <= k+1. It matters
  # far more than it looks: assign_gainloss() z-scores the segment means WITHIN a
  # chromosome and hard-codes a |z| < 1 collapse to normal. A chromosome split into
  # exactly two segments gives z = +/-0.707 for both, so a clean whole-arm event is
  # forced to "normal" no matter how large it is. Raising k adds segments and can
  # break that degeneracy.
  k            = KBP,
  plotKaryo    = TRUE
)

el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))

# Record how many nuclei and bins actually survived -- these are the numbers the
# Methods section has to quote, and they are the whole reason both parameters are
# being swept.
#
# Two traps in reading this file, both of which cost a wasted pass earlier:
# it sits one directory BELOW the outdir passed to epiAneufinder(), and it is
# SPACE-separated, not tab-separated.
res <- file.path(OUT, "epiAneufinder_results", "results_table.tsv")
if (file.exists(res)) {
  h     <- strsplit(readLines(res, n = 1), " ")[[1]]
  ncell <- sum(grepl("^cell-", h))
  nbin  <- length(readLines(res)) - 1L
  cat(sprintf("SUMMARY\t%s\t%s\t%.0f\t%d\t%d\t%d\t%.1f\n",
              CARE, TAG, WINDOW, MINFRAGS, ncell, nbin, el))
} else {
  cat(sprintf("SUMMARY\t%s\t%s\t%.0f\t%d\tNA\tNA\t%.1f\tNO_RESULTS_TABLE\n",
              CARE, TAG, WINDOW, MINFRAGS, el))
}
