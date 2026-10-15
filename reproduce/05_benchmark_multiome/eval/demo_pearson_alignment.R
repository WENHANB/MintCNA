#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# How a cross-method Pearson correlation is actually built
#
# Reproduces the mechanism behind epiAneufinder Fig 3a, on our own GATK truth.
# Three mock methods stand in for the three native output formats you will meet:
#
#   A  per-bin, per-cell discrete states   (epiAneufinder, MintCNA)
#   B  chromosome-arm scalars              (Copy-scAT)
#   C  per-gene values                     (inferCNV, CaSpER, copyKat)
#
# The point of the script is the projection step, not the mock data.
# Run:  module load R/4.4 && Rscript Code/demo_pearson_alignment.R
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({library(GenomicRanges); library(data.table)})
set.seed(1)

ROOT <- "<MINTCNA_ROOT>/ATAC_RNA"
BIN  <- 1e5                      # 100 kb analysis grid
CHR  <- "1"                      # one chromosome keeps the demo readable
say  <- function(...) cat(sprintf(...), "\n")
rule <- function(t) cat("\n", strrep("=", 76), "\n", t, "\n", strrep("=", 76), "\n", sep="")

# ---------------------------------------------------------------------------
rule("STEP 1  truth: variable-length GATK segments -> GRanges")
# ---------------------------------------------------------------------------
seg <- fread(file.path(ROOT, "Data/synapse/variants_gatk_seg.csv"))
seg <- seg[aliquot_barcode == "OLIG-SJ-0006-R1-01D-WGS-606WVS" & chrom == CHR]
# pos is the string "[start,end)" -- parse it into two integers
seg[, `:=`(start = as.numeric(sub("\\[(\\d+),.*", "\\1", pos)),
           end   = as.numeric(sub(".*,(\\d+)\\)", "\\1", pos)))]
truth_gr <- GRanges(paste0("chr", CHR), IRanges(seg$start, seg$end - 1),
                    log2 = seg$log2_copy_ratio)
say("segments on chr%s: %d, spanning %.1f Mb, log2 range %.2f to %.2f",
    CHR, length(truth_gr), sum(width(truth_gr))/1e6,
    min(truth_gr$log2), max(truth_gr$log2))

# ---------------------------------------------------------------------------
rule("STEP 2  the anchor grid")
# ---------------------------------------------------------------------------
# Everything gets projected onto this. In the paper the anchor was
# epiAneufinder's own tileGenome output; here we build it explicitly.
chr_end <- max(end(truth_gr))
grid <- GRanges(paste0("chr", CHR),
                IRanges(seq(1, chr_end, by = BIN),
                        width = BIN))
say("grid: %d bins of %.0f kb", length(grid), BIN/1e3)

# --- project the truth onto the grid ---------------------------------------
# Overlap-length weighted mean: a bin straddling two segments gets the
# bp-weighted average, not whichever segment happens to be first.
project_weighted <- function(anchor, src, col) {
  hits <- findOverlaps(anchor, src)
  if (!length(hits)) return(rep(NA_real_, length(anchor)))
  w   <- width(pintersect(anchor[queryHits(hits)], src[subjectHits(hits)]))
  v   <- mcols(src)[[col]][subjectHits(hits)]
  num <- tapply(v * w, queryHits(hits), sum)
  den <- tapply(w,     queryHits(hits), sum)
  out <- rep(NA_real_, length(anchor))
  out[as.integer(names(num))] <- num / den
  out
}
grid$truth <- project_weighted(grid, truth_gr, "log2")
say("bins with truth coverage: %d / %d  (%d dropped -- centromere, telomere, low mappability)",
    sum(!is.na(grid$truth)), length(grid), sum(is.na(grid$truth)))

# ---------------------------------------------------------------------------
rule("STEP 3  three mock methods, each in its NATIVE output format")
# ---------------------------------------------------------------------------
NCELL <- 200
tr <- grid$truth

## --- A: per-bin per-cell discrete states, like epiAneufinder ---------------
# state 0 = loss, 1 = neutral, 2 = gain, one column per cell
pA <- ifelse(is.na(tr), 1, ifelse(tr < -0.3, 0, ifelse(tr > 0.3, 2, 1)))
A_states <- sapply(seq_len(NCELL), function(i) {
  s <- pA
  flip <- sample(length(s), round(0.15 * length(s)))   # 15% per-cell noise
  s[flip] <- sample(0:2, length(flip), replace = TRUE)
  s
})
say("A (epiAneufinder-like): %d bins x %d cells of {0,1,2}", nrow(A_states), ncol(A_states))

## --- B: chromosome-arm scalars, like Copy-scAT -----------------------------
# chr1 hg19 centromere ~125.0 Mb
arm_gr <- GRanges(paste0("chr", CHR),
                  IRanges(c(1, 125e6), c(125e6 - 1, chr_end)),
                  value = c(2 + 0.6 * (-1.4), 2 + 0.6 * (0.3)))   # 2 + sdCNV*z
say("B (Copy-scAT-like): %d arm-level values: %.3f, %.3f",
    length(arm_gr), arm_gr$value[1], arm_gr$value[2])

## --- C: per-gene values, like inferCNV -------------------------------------
# genes are unevenly spaced, and this is what shrinks the comparison set
NGENE <- 900
gs <- sort(sample(seq(1, chr_end - 4e4), NGENE))
gene_gr <- GRanges(paste0("chr", CHR), IRanges(gs, width = 3e4))
gi <- findOverlaps(gene_gr, grid, select = "first")
gv <- tr[gi]
gene_gr$value <- ifelse(is.na(gv), 2,
                        ifelse(gv < -0.3, 1, ifelse(gv > 0.3, 3, 2))) +
                 rnorm(NGENE, 0, 0.25)
say("C (inferCNV-like): %d genes, each with one value -- NOT on the grid", length(gene_gr))

# ---------------------------------------------------------------------------
rule("STEP 4  collapse per-cell output to pseudobulk on the 1/2/3 scale")
# ---------------------------------------------------------------------------
# The paper's formula:  mean over cells of (loss=1, neutral=2, gain=3)
# For a 0/1/2 state matrix that is simply rowMeans(states) + 1.
grid$A <- rowMeans(A_states) + 1
say("A pseudobulk: range %.3f - %.3f, %d distinct values",
    min(grid$A), max(grid$A), length(unique(round(grid$A, 6))))
say("  note: averaging 200 cells of {0,1,2} turned a 3-state caller into a continuous one")

# ---------------------------------------------------------------------------
rule("STEP 5  project B and C onto the anchor grid")
# ---------------------------------------------------------------------------
# This is combine_range_objects() from the paper: findOverlaps, then an
# UNWEIGHTED mean of every source feature touching the bin.
project_unweighted <- function(anchor, src, col) {
  hits <- findOverlaps(anchor, src, type = "any")
  out  <- rep(NA_real_, length(anchor))
  if (!length(hits)) return(out)
  m <- tapply(mcols(src)[[col]][subjectHits(hits)], queryHits(hits), mean)
  out[as.integer(names(m))] <- m
  out
}
grid$B <- project_unweighted(grid, arm_gr,  "value")   # broadcast: arm -> every bin
grid$C <- project_unweighted(grid, gene_gr, "value")   # average genes within bin

say("B after broadcast: %d bins covered, but only %d DISTINCT values (= number of arms)",
    sum(!is.na(grid$B)), length(unique(na.omit(grid$B))))
say("C after gene mapping: %d / %d bins covered -- %d bins contain NO gene and are dropped",
    sum(!is.na(grid$C)), length(grid), sum(is.na(grid$C)))

# ---------------------------------------------------------------------------
rule("STEP 6  intersect: keep only bins where EVERY method has a value")
# ---------------------------------------------------------------------------
df <- as.data.frame(mcols(grid))
keep <- complete.cases(df)
say("sequential attrition:")
say("  grid                        %5d bins", length(grid))
say("  + truth covered             %5d", sum(!is.na(df$truth)))
say("  + truth & A                 %5d", sum(!is.na(df$truth) & !is.na(df$A)))
say("  + truth & A & B             %5d", sum(!is.na(df$truth) & !is.na(df$A) & !is.na(df$B)))
say("  + truth & A & B & C         %5d   <-- final comparison set", sum(keep))
say("")
say("  %.1f%% of the grid survives. In the real Fig 3a it was 7,301 / 26,557 = 27.5%%,",
    100 * sum(keep) / length(grid))
say("  and the gene-based methods caused nearly all of the loss -- same mechanism as here.")
d <- df[keep, ]

# ---------------------------------------------------------------------------
rule("STEP 7  Pearson")
# ---------------------------------------------------------------------------
say("%-28s %8s", "method", "r vs truth")
for (m in c("A", "B", "C")) say("%-28s %8.3f", m, cor(d[[m]], d$truth))
say("")
say("Pearson is affine-invariant, so the 1/2/3 scale of A and B never has to be")
say("converted into log2. That is the whole trick -- and also the whole limitation:")
say("it measures the SHAPE of the profile, never its amplitude.")

# ---------------------------------------------------------------------------
rule("STEP 8  the diagnostic nobody publishes: centre within chromosome first")
# ---------------------------------------------------------------------------
# Much of a genome-wide correlation is just whole-chromosome ploidy offset.
# Removing the per-chromosome mean shows how much agreement is left at finer scale.
cen <- as.data.frame(lapply(d, function(x) x - mean(x)))
say("%-28s %10s %10s %8s", "method", "raw r", "centred r", "drop")
for (m in c("A", "B", "C")) {
  r0 <- cor(d[[m]], d$truth); r1 <- cor(cen[[m]], cen$truth)
  say("%-28s %10.3f %10.3f %8.3f", m, r0, r1, r0 - r1)
}
say("")
say("(single chromosome here, so the drop is small. Genome-wide in the real paper:")
say(" epiAneufinder 0.857 -> 0.796, Copy-scAT 0.734 -> 0.591, inferCNV 0.609 -> 0.461.)")
