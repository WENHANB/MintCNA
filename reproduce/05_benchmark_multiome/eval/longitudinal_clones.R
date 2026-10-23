#!/usr/bin/env Rscript
# Longitudinal subclone analysis: what changes between the earlier and the later
# resection of the same tumour?
#
# WHY THIS EXISTS. Of the seven single-cell CNA callers we benchmark, none
# compares clones across timepoints -- and it is not for want of data: AtaCNV had
# pre/post anti-PD-1 BCC, Copy-scAT had five diagnosis-versus-relapse pGBM pairs,
# RIDDLER had a post-treatment liver metastasis. All of them stopped at the
# accuracy panel. This cohort is ten complete pairs, and the question needs
# per-cell copy number in the same cells where the phenotype is measured, which
# is exactly what a multiome caller provides and a bulk or scRNA method does not.
#
# WHAT IS MEASURED, and the confounds each measure has to survive:
#   burden      fraction of the genome called altered. The obvious confound is
#               sequencing depth -- a deeper library calls more events -- so the
#               paired depth difference is reported beside it and the burden
#               change is regressed on it.
#   clonality   number of subclones and the dominant clone's share. More cells
#               find more clones, so cell count is reported beside it too.
#   direction   gains and losses separately: a tumour that gains genome at
#               relapse is a different story from one that loses it.
#
# Nothing here is a per-cell test. The replicate unit is the PATIENT, n = 10, and
# every comparison is paired within patient.
.libPaths(c("<R_LIBS_USER>", .libPaths()))
suppressPackageStartupMessages(library(data.table))
setwd("<MINTCNA_ROOT>/ATAC_RNA")
TAG <- "h15a005e020"
OUT <- "Result/qc"

man <- fread("Data/meta/sample_manifest.csv")
have <- sub("_results.*", "", basename(Sys.glob(sprintf("Result/MintCNA/*/*_results_%s.RData", TAG))))
man <- man[care_id %in% have]
# Initial < Recurrence < 2nd Recurrence
man[, tp_rank := match(timepoint, c("Initial", "Recurrence", "2nd Recurrence"))]
man[, grade_num := as.integer(sub("^G", "", grade))]
pairs <- man[, .N, by = patient][N == 2, patient]
man <- man[patient %in% pairs][order(patient, tp_rank)]
cat(sprintf("%d patients, %d samples, all paired\n\n", length(pairs), nrow(man)))

## ---- per-sample summaries ---------------------------------------------------
one <- function(S) {
  e <- new.env()
  load(sprintf("Result/MintCNA/%s_multi/%s_results_%s.RData", S, S, TAG), e)
  cna <- base::get("cna_state_matrix", e); sub <- base::get("subclones", e)
  ref <- as.data.frame(base::get("ref", e))
  cl  <- fread(sprintf("Data/cell_lists/%s_cells.tsv", S))
  keep <- cl[cell_type == "Malignant", barcode]
  i <- colnames(cna) %in% keep
  cna <- cna[, i, drop = FALSE]; sub <- sub[i]
  bc  <- colnames(cna)
  autos <- sub("^chr", "", as.character(ref$chr)) %in% as.character(1:22)
  cna <- cna[autos, , drop = FALSE]

  # per-cell altered fraction, then the median cell -- not the mean, so one
  # pathological nucleus cannot move a patient-level number
  f_alt  <- colMeans(cna != 2L)
  f_gain <- colMeans(cna == 3L)
  f_loss <- colMeans(cna == 1L)
  # consensus: modal state per bin across malignant nuclei
  modal <- apply(cna, 1, function(v) as.integer(names(which.max(table(v)))))
  tb <- table(sub); frac <- sort(as.numeric(tb) / sum(tb), decreasing = TRUE)
  data.table(
    care_id = S, n_cell = ncol(cna), n_bin = nrow(cna),
    depth = median(cl$nFrags[match(bc, cl$barcode)], na.rm = TRUE),
    burden       = median(f_alt),
    burden_gain  = median(f_gain),
    burden_loss  = median(f_loss),
    burden_cons  = mean(modal != 2L),
    n_clone      = length(tb),
    dom_frac     = frac[1],
    shannon      = -sum(frac * log(frac)))
}
S <- rbindlist(lapply(man$care_id, one))
S <- merge(man[, .(care_id, patient, timepoint, tp_rank, grade_num, subtype)], S, by = "care_id")
setorder(S, patient, tp_rank)
cat("=== per-sample ===\n")
print(S[, .(care_id, patient, timepoint, grade_num, n_cell, depth,
            burden = round(burden, 3), gain = round(burden_gain, 3),
            loss = round(burden_loss, 3), n_clone, dom = round(dom_frac, 2))],
      row.names = FALSE)

## ---- paired deltas, later minus earlier ------------------------------------
W <- dcast(S, patient + subtype ~ tp_rank, value.var =
             c("care_id", "burden", "burden_gain", "burden_loss", "burden_cons",
               "n_clone", "dom_frac", "shannon", "depth", "n_cell", "grade_num"))
lo <- grep("_1$", names(W), value = TRUE); hi <- grep("_[23]$", names(W), value = TRUE)
# each patient has exactly two timepoints, but they are not always ranks 1 and 2
num <- function(pref) {
  cols <- grep(sprintf("^%s_[0-9]+$", pref), names(W), value = TRUE)
  a <- as.matrix(W[, ..cols])
  t(apply(a, 1, function(r) { k <- which(!is.na(r)); c(early = r[k[1]], late = r[k[2]]) }))
}
D <- data.table(patient = W$patient, subtype = W$subtype)
for (v in c("burden","burden_gain","burden_loss","burden_cons","n_clone",
            "dom_frac","shannon","depth","n_cell","grade_num")) {
  m <- num(v); D[[paste0(v, "_e")]] <- m[, 1]; D[[paste0(v, "_l")]] <- m[, 2]
  D[[paste0("d_", v)]] <- m[, 2] - m[, 1]
}
D[, progressed := d_grade_num > 0]

cat("\n=== paired change, later minus earlier ===\n")
print(D[, .(patient, subtype, grade = sprintf("G%d->G%d", grade_num_e, grade_num_l),
            d_burden = round(d_burden, 3), d_gain = round(d_burden_gain, 3),
            d_loss = round(d_burden_loss, 3), d_clone = d_n_clone,
            d_dom = round(d_dom_frac, 2), d_depth = round(d_depth))], row.names = FALSE)

cat("\n=== paired tests, n =", nrow(D), "patients ===\n")
tst <- function(lab, v) {
  x <- D[[v]]
  w <- wilcox.test(x, mu = 0, exact = FALSE)
  cat(sprintf("  %-28s median %+.4f   Wilcoxon P = %.3f   (%d of %d increase)\n",
              lab, median(x), w$p.value, sum(x > 0), length(x)))
}
tst("CNA burden",            "d_burden")
tst("  gained fraction",     "d_burden_gain")
tst("  lost fraction",       "d_burden_loss")
tst("consensus burden",      "d_burden_cons")
tst("number of subclones",   "d_n_clone")
tst("dominant clone share",  "d_dom_frac")
tst("clonal diversity (H)",  "d_shannon")

cat("\n=== confound check: is the burden change just depth or cell number? ===\n")
for (cv in c("d_depth", "d_n_cell")) {
  ct <- cor.test(D$d_burden, D[[cv]], method = "spearman", exact = FALSE)
  cat(sprintf("  d_burden vs %-9s Spearman rho = %+.3f  P = %.3f\n",
              cv, unname(ct$estimate), ct$p.value))
}
cat(sprintf("  depth itself: median change %+.0f fragments, %d of %d pairs deeper at the later timepoint\n",
            median(D$d_depth), sum(D$d_depth > 0), nrow(D)))

cat("\n=== does burden change track grade progression? ===\n")
print(D[, .(n = .N, median_d_burden = round(median(d_burden), 4)), by = progressed], row.names = FALSE)
if (D[, uniqueN(progressed)] > 1)
  cat(sprintf("  Mann-Whitney progressed vs stable: P = %.3f\n",
              wilcox.test(d_burden ~ progressed, data = D, exact = FALSE)$p.value))

fwrite(S, file.path(OUT, "longitudinal_per_sample.tsv"), sep = "\t")
fwrite(D, file.path(OUT, "longitudinal_paired.tsv"), sep = "\t")
cat(sprintf("\nwrote %s/longitudinal_{per_sample,paired}.tsv\n", OUT))
