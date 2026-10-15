#!/usr/bin/env Rscript
# The check that must pass before any benchmark number is believed.
#
# 1p/19q codeletion is the defining truncal event of oligodendroglioma and is
# absent from astrocytoma. It is present in essentially every malignant cell, so
# any method that cannot see it on this cohort has a pipeline fault -- wrong
# coordinates, wrong reference cells, a scrambled bin axis -- not a performance
# deficit. Running this first turns a whole class of silent errors into a
# visible one.
#
# The check runs on the TRUTH first. If the bulk WGS does not show the event in a
# given sample, a method that misses it there is not wrong, and grading it as a
# miss would be the actual error.

suppressPackageStartupMessages(library(data.table))
ROOT <- "<MINTCNA_ROOT>/ATAC_RNA"
setwd(ROOT)

# hg38 centromere midpoints
P1_END    <- 123.4e6
Q19_START <- 28.1e6

S <- c("P100T1","P100T3","P103T1","P103T2","P105T1","P105T3","P108T1","P108T2",
       "P109T1","P109T2","P110T1","P110T2","P111T2","P111T3",
       "P59T2","P59T3","P104T1","P104T2","P106T1","P106T2")
m <- fread("Data/meta/sample_manifest.csv"); setkey(m, care_id)

# ---- 1. the truth ----------------------------------------------------------
truth <- rbindlist(lapply(S, function(s) {
  f <- file.path("Result/truth", paste0(s, "_truth_100kb.tsv.gz"))
  if (!file.exists(f)) return(NULL)
  d <- fread(f)
  data.table(care_id = s, subtype = m[s, subtype],
             truth_1p  = d[chr == "chr1"  & start <  P1_END,    mean(log2, na.rm = TRUE)],
             truth_19q = d[chr == "chr19" & start >= Q19_START, mean(log2, na.rm = TRUE)])
}))
# a whole-arm one-copy loss reads about -0.3 to -1.0 depending on tumour purity
truth[, truth_codel := truth_1p < -0.15 & truth_19q < -0.15]
# Grade by the strength of the CODELETION, which requires both arms. Grading on
# min(1p, 19q) instead would label P104T1 "strong" -- it has a genuine 19q loss
# (-0.51) but no 1p loss (+0.05), so it is an astrocytoma with a single-arm event,
# not a weak codeletion, and Copy-scAT is right to call 19q and not 1p there.
truth[, truth_strength := fifelse(truth_codel & pmax(truth_1p, truth_19q) < -0.5, "strong",
                          fifelse(truth_codel, "weak", "absent"))]

# ---- 2. Copy-scAT (arm-level output, so it reads directly) -----------------
cs <- rbindlist(lapply(S, function(s) {
  f <- file.path("Result/benchmark/CopyscAT", s, paste0(s, "clean_cnv_cnv_scores.csv"))
  if (!file.exists(f)) return(NULL)
  d <- read.csv(f, row.names = 1, check.names = FALSE)
  # Copy-scAT reports only arms it called as altered; an arm absent from the
  # output is a non-call, which is NOT the same as a call of "normal".
  frac_lost <- function(a) if (a %in% colnames(d)) mean(d[[a]] < 1.5) else NA_real_
  data.table(care_id = s, cs_ncell = nrow(d), cs_narm = ncol(d),
             cs_1p = frac_lost("chr1p"), cs_19q = frac_lost("chr19q"))
}))

x <- merge(truth, cs, by = "care_id", all.x = TRUE)
x[, cs_codel := !is.na(cs_1p) & !is.na(cs_19q) & cs_1p > 0.5 & cs_19q > 0.5]
setorder(x, subtype, truth_1p)

fwrite(x, "Result/qc/sanity_1p19q.tsv", sep = "\t")

fmt <- function(v) fifelse(is.na(v), "  --", sprintf("%3.0f%%", 100 * v))
cat("=== 1p/19q sanity check ===\n\n")
print(x[, .(care_id, subtype = substr(subtype, 1, 12),
            t1p = round(truth_1p, 2), t19q = round(truth_19q, 2), truth = truth_strength,
            cs1p = fmt(cs_1p), cs19q = fmt(cs_19q), cs = cs_codel)], nrow = 25)

cat("\n--- graded against the truth, not against the diagnosis label ---\n")
for (lvl in c("strong", "weak", "absent")) {
  sub <- x[truth_strength == lvl]
  if (!nrow(sub)) next
  cat(sprintf("truth %-6s (n=%2d): Copy-scAT calls codeletion in %d\n",
              lvl, nrow(sub), sum(sub$cs_codel)))
}
cat("\ncodeletion false positives where the truth has none:",
    sum(x[truth_strength == "absent", cs_codel]), "of", nrow(x[truth_strength == "absent"]), "\n")

# Single-arm accuracy is the more informative view, since a codeletion call needs
# both arms and one missing arm hides an otherwise correct detection.
cat("\n--- per-arm, against the truth for that arm ---\n")
arm <- rbindlist(list(
  x[, .(care_id, arm = "1p",  t = truth_1p,  called = cs_1p)],
  x[, .(care_id, arm = "19q", t = truth_19q, called = cs_19q)]))
arm[, truth_loss := t < -0.15]
arm[, detected  := !is.na(called) & called > 0.5]
print(arm[, .(n = .N, detected = sum(detected)), by = .(arm, truth_loss)][order(arm, -truth_loss)])
