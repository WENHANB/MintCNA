#!/usr/bin/env Rscript
# Which copy-number events are TRUNCAL and which are ACQUIRED at recurrence?
#
# The burden comparison in longitudinal_clones.R came back null, and what signal
# it had tracked the number of cells profiled (Spearman 0.709, P = 0.022). Burden
# is a count, and counts inherit every power difference between two libraries.
# This script asks a positional question instead -- WHICH bins are altered, not
# how many -- which is what the tumour-evolution literature actually reports and
# is far less sensitive to how deeply each timepoint was sampled.
#
# Per patient, each autosomal bin is classified from the two consensus profiles:
#   truncal   altered in the same direction at both timepoints
#   lost      altered early, neutral late
#   acquired  neutral early, altered late
#   switched  altered both, opposite directions
#
# The claim worth making is not about one patient. It is whether the SAME regions
# are acquired across patients more often than chance, which is tested against a
# within-patient circular permutation that preserves each tumour's own segment
# structure and event count -- shuffling bins independently would destroy the
# autocorrelation of CNA segments and manufacture significance.
.libPaths(c("<R_LIBS_USER>", .libPaths()))
suppressPackageStartupMessages(library(data.table))
setwd("<MINTCNA_ROOT>/ATAC_RNA")
set.seed(1)
TAG <- "h15a005e020"; OUT <- "Result/qc"; NPERM <- 2000

man <- fread("Data/meta/sample_manifest.csv")
have <- sub("_results.*", "", basename(Sys.glob(sprintf("Result/MintCNA/*/*_results_%s.RData", TAG))))
man <- man[care_id %in% have]
man[, tp_rank := match(timepoint, c("Initial", "Recurrence", "2nd Recurrence"))]
man[, grade_num := as.integer(sub("^G", "", grade))]
man <- man[patient %in% man[, .N, by = patient][N == 2, patient]][order(patient, tp_rank)]

consensus <- function(S) {
  e <- new.env()
  load(sprintf("Result/MintCNA/%s_multi/%s_results_%s.RData", S, S, TAG), e)
  cna <- base::get("cna_state_matrix", e); ref <- as.data.table(base::get("ref", e))
  keep <- fread(sprintf("Data/cell_lists/%s_cells.tsv", S))[cell_type == "Malignant", barcode]
  cna <- cna[, colnames(cna) %in% keep, drop = FALSE]
  ref[, chr := sub("^chr", "", as.character(chr))]
  k <- ref$chr %in% as.character(1:22)
  cna <- cna[k, , drop = FALSE]; ref <- ref[k]
  # modal state per bin; ties resolve to neutral, which is the conservative call
  md <- apply(cna, 1, function(v) { t <- tabulate(v, 3); if (sum(t == max(t)) > 1) 2L else which.max(t) })
  list(ref = ref[, .(chr, start, end)], state = as.integer(md))
}

cls <- function(a, b) fifelse(a == b & a != 2L, "truncal",
              fifelse(a != 2L & b == 2L, "lost",
              fifelse(a == 2L & b != 2L, "acquired",
              fifelse(a != 2L & b != 2L & a != b, "switched", "neutral"))))

# TWO PASSES. The shared grid is the intersection over ALL patients, not the
# first patient's bins: per-sample bin counts run 5,241-5,259 because each sample
# drops its own low-mappability bins, so keying off patient one left almost
# nobody on the grid and the recurrence test silently never ran.
key <- function(r) paste(r$chr, r$start, sep = ":")
CON <- list()
for (p in unique(man$patient)) {
  ids <- man[patient == p][order(tp_rank), care_id]
  CON[[p]] <- list(ids = ids, A = consensus(ids[1]), B = consensus(ids[2]))
}
GRID_KEY <- Reduce(intersect, lapply(CON, function(x) intersect(key(x$A$ref), key(x$B$ref))))
GRID <- CON[[1]]$A$ref[match(GRID_KEY, key(CON[[1]]$A$ref))]
cat(sprintf("shared grid across all %d patients: %d bins\n\n", length(CON), length(GRID_KEY)))

res <- list(); acq <- list()
for (p in names(CON)) {
  x <- CON[[p]]
  a <- x$A$state[match(GRID_KEY, key(x$A$ref))]
  b <- x$B$state[match(GRID_KEY, key(x$B$ref))]
  cc <- cls(a, b)
  res[[p]] <- data.table(patient = p, early = x$ids[1], late = x$ids[2], n_bin = length(GRID_KEY),
                         truncal = mean(cc == "truncal"), lost = mean(cc == "lost"),
                         acquired = mean(cc == "acquired"), switched = mean(cc == "switched"),
                         altered_early = mean(a != 2L), altered_late = mean(b != 2L))
  acq[[p]] <- as.integer(cc == "acquired")
}
R <- rbindlist(res)
cat("=== per patient, fraction of autosomal bins ===\n")
print(merge(R, man[, .(patient, grade_num, tp_rank)][, .(gr = paste0("G", grade_num, collapse = "->")), by = patient],
            by = "patient")[, .(patient, gr, n_bin,
        truncal = round(truncal, 3), lost = round(lost, 3),
        acquired = round(acquired, 3), switched = round(switched, 3))], row.names = FALSE)
cat(sprintf("\n  median truncal %.3f | lost %.3f | acquired %.3f | switched %.3f\n",
            median(R$truncal), median(R$lost), median(R$acquired), median(R$switched)))
cat(sprintf("  of everything altered at EITHER timepoint, %.1f%% is shared\n",
            100 * median(R$truncal / (R$truncal + R$lost + R$acquired + R$switched))))
w <- wilcox.test(R$acquired, R$lost, paired = TRUE, exact = FALSE)
cat(sprintf("  acquired vs lost, paired over %d patients: median %+.4f, Wilcoxon P = %.3f\n",
            nrow(R), median(R$acquired - R$lost), w$p.value))

## ---- recurrently acquired regions ------------------------------------------
usable <- names(acq)
if (length(usable) >= 5) {
  M <- do.call(cbind, acq[usable])          # bins x patients, 0/1
  obs <- rowSums(M)
  cat(sprintf("\n=== recurrent acquisition across %d patients on a shared %d-bin grid ===\n",
              ncol(M), nrow(M)))
  # circular shift per patient: preserves each patient's segment lengths and count
  perm <- replicate(NPERM, {
    rowSums(vapply(seq_len(ncol(M)), function(j) {
      v <- M[, j]; s <- sample.int(length(v), 1); c(v[-seq_len(s)], v[seq_len(s)])
    }, numeric(nrow(M))))
  })
  pv <- (rowSums(perm >= obs) + 1) / (NPERM + 1)
  G <- copy(GRID)[, `:=`(n_acq = obs, p = pv)][, fdr := p.adjust(p, "BH")]
  cat(sprintf("  max patients sharing an acquired bin: %d\n", max(obs)))
  cat(sprintf("  bins with FDR < 0.10: %d\n", G[fdr < 0.10, .N]))
  top <- G[order(-n_acq, p)][1:15]
  print(top[, .(chr, start = format(start, big.mark = ","), n_acq,
                p = signif(p, 3), fdr = signif(fdr, 3))], row.names = FALSE)
  if (G[fdr < 0.10, .N] > 0) {
    sig <- G[fdr < 0.10][order(chr, start)]
    sig[, grp := cumsum(c(1, diff(start) != 5e5 | diff(as.integer(factor(chr))) != 0))]
    cat("\n  contiguous significant regions:\n")
    print(sig[, .(chr = chr[1], from = format(min(start), big.mark = ","),
                  to = format(max(end), big.mark = ","), n_bin = .N,
                  n_patients = max(n_acq)), by = grp][, grp := NULL][], row.names = FALSE)
  }
  fwrite(G, file.path(OUT, "longitudinal_acquired_bins.tsv"), sep = "\t")
}
fwrite(R, file.path(OUT, "longitudinal_event_classes.tsv"), sep = "\t")
cat(sprintf("\nwrote %s/longitudinal_{event_classes,acquired_bins}.tsv\n", OUT))
