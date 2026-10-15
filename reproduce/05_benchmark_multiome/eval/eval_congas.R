#!/usr/bin/env Rscript
# Is CONGAS+ tracking the truth at all?
#
# CONGAS+ returns an integer copy number per (cluster, segment), not per cell, so
# the fair test is WITHIN a sample, ACROSS segments: does its per-segment call
# rank the genome the way the bulk WGS does? A pure ploidy offset (everything
# called 3 instead of 2) leaves that correlation intact; a scrambled ordering
# destroys it. That distinction is what decides whether the output is usable
# after recentring or not usable at all.
#
# Cluster-level calls are collapsed to one profile per sample by weighting each
# cluster by the number of cells assigned to it.

suppressPackageStartupMessages({library(Rcongas); library(data.table)})
ROOT <- "<MINTCNA_ROOT>/ATAC_RNA"
setwd(ROOT)

S <- c("P100T1","P100T3","P103T1","P103T2","P105T1","P105T3","P108T1","P108T2",
       "P109T1","P109T2","P110T1","P110T2","P111T2","P111T3",
       "P59T2","P59T3","P104T1","P104T2","P106T1","P106T2")
man <- fread("Data/meta/sample_manifest.csv"); setkey(man, care_id)

res <- rbindlist(lapply(S, function(s) {
  f <- file.path("Result/benchmark/CONGASplus", s, "congas_result.rds")
  if (!file.exists(f)) return(NULL)
  r <- readRDS(f)

  # best_fit$CNA is the ROUNDED integer call; CNA_real is the continuous estimate
  # behind it. Correlating the continuous version is the fairer test -- a uniform
  # ploidy offset leaves it intact, and rounding to integers can only destroy
  # information. Fall back to CNA if CNA_real is absent.
  cna <- if (!is.null(r$fit$best_fit$CNA_real))
           as.data.table(r$fit$best_fit$CNA_real)
         else as.data.table(Rcongas::get_fit(r$fit, what = "CNA"))
  ass <- tryCatch(as.data.table(Rcongas::get_fit(r$fit, what = "cluster_assignments")),
                  error = function(e) NULL)
  w <- if (!is.null(ass) && "cluster" %in% names(ass))
         ass[, .(n = .N), by = cluster] else data.table(cluster = unique(cna$cluster), n = 1)
  cna <- merge(cna, w, by = "cluster")
  prof <- cna[, .(cn = sum(value * n) / sum(n)), by = segment_id]   # cell-weighted

  prof[, c("chr", "from", "to") := tstrsplit(segment_id, ":", fixed = TRUE)]
  prof[, `:=`(from = as.numeric(from), to = as.numeric(to))]

  tr <- fread(file.path("Result/truth", paste0(s, "_truth_100kb.tsv.gz")))
  prof[, truth := vapply(seq_len(.N), function(i)
        tr[chr == prof$chr[i] & start >= prof$from[i] & end <= prof$to[i] &
           is.finite(log2), mean(log2)], numeric(1))]

  ok <- prof[is.finite(truth) & is.finite(cn)]
  # CONGAS+ integer CN vs truth log2: compare on a common scale
  rho <- if (nrow(ok) > 5 && sd(ok$cn) > 0) cor(log2(ok$cn / 2), ok$truth) else NA_real_

  a1p  <- prof[chr == "chr1"  & from < 1e8, cn][1]
  a19q <- prof[chr == "chr19" & from > 2e7, cn][1]
  med  <- median(prof$cn, na.rm = TRUE)

  data.table(care_id = s, subtype = substr(man[s, subtype], 1, 12),
             n_clust = uniqueN(cna$cluster), n_seg = nrow(prof),
             cn_median = round(med, 2),
             t1p  = round(tr[chr == "chr1"  & start <  123.4e6, mean(log2, na.rm = TRUE)], 2),
             t19q = round(tr[chr == "chr19" & start >= 28.1e6, mean(log2, na.rm = TRUE)], 2),
             cn_1p = round(a1p, 2), cn_19q = round(a19q, 2),
             pearson = round(rho, 3))
}))

res[, codel := t1p < -0.15 & t19q < -0.15]
# does CONGAS+ put 1p/19q BELOW that sample's own median segment?
res[, rel_1p  := cn_1p  - cn_median]
res[, rel_19q := cn_19q - cn_median]
setorder(res, subtype, t1p)

fwrite(res, "Result/qc/congas_vs_truth.tsv", sep = "\t")
cat("=== CONGAS+ per-sample, cell-weighted CN profile vs bulk WGS truth ===\n\n")
print(res[, .(care_id, subtype, n_clust, cn_median, t1p, cn_1p, rel_1p,
              t19q, cn_19q, rel_19q, pearson)], nrow = 25)

cat("\n--- Pearson (log2(CN/2) vs truth log2, across segments, within sample) ---\n")
cat(sprintf("median %.3f   mean %.3f   range %.3f to %.3f   n>0.3: %d/%d\n",
            median(res$pearson, na.rm = TRUE), mean(res$pearson, na.rm = TRUE),
            min(res$pearson, na.rm = TRUE), max(res$pearson, na.rm = TRUE),
            sum(res$pearson > 0.3, na.rm = TRUE), nrow(res)))

cat("\n--- on the 12 truth-positive samples, is 1p/19q called below the sample median? ---\n")
o <- res[codel == TRUE]
cat(sprintf("1p  below own median: %d/%d\n", sum(o$rel_1p  < 0, na.rm = TRUE), nrow(o)))
cat(sprintf("19q below own median: %d/%d\n", sum(o$rel_19q < 0, na.rm = TRUE), nrow(o)))
cat(sprintf("clusters found: %s\n", paste(sort(unique(res$n_clust)), collapse = ", ")))
