#!/usr/bin/env Rscript
# What does the COPY NUMBER itself explain in panels D and E?
#
# Panels D and E compare two subclones and find lineage differences. On their own
# that is a statement about cell state, and a reader is entitled to ask what the
# CNA contributed beyond splitting the cells into two groups the RNA annotation
# already separates. Two questions, and they pull in opposite directions:
#
#   1. DOSAGE. Where the two clones differ in copy number, do accessibility and
#      expression follow? If yes, the CNA is causal, not just a grouping label.
#
#      THE INDEPENDENCE HERE IS NARROWER THAN IT LOOKS, so AR_ARM exists.
#      The clone PARTITION is genuinely RNA-free: the AG-DCAE latent it is
#      clustered from was trained on the ATAC log2-ratio alone, and the proof is
#      that the subclone vectors of the multi arm and of the ATAC-only ablation
#      are byte-identical (ARI = 1.000, one shared AE training).
#      The differing REGIONS were not. In mode "multi" the change-point scan is a
#      multivariate CUSUM over an rbind(RNA, DNA) matrix and the del/nor/dup call
#      is a coverage-weighted fusion of both rows, so RNA -- which carries a value
#      in 77.6% of bins -- helped choose which bin a gene lands in. 13.3% of
#      (bin, cell) state calls change when it is removed. Reading a dosage
#      staircase off the multi arm and calling it independent corroboration is
#      therefore partly circular: RNA is on both axes.
#      AR_ARM=atac recomputes d_cn from the ATAC-only arm at the same tag, which
#      removes RNA from the x axis entirely. That is the version to report.
#
#   2. CONFOUNDING. Are the lineage genes of panel E sitting in those same
#      differing regions? If PDGFRA happens to be in a region gained in the
#      OPC-like clone, its expression shift is dosage, not lineage, and panel E's
#      reading would be wrong. This has to be checked, not assumed.
.libPaths(c("<R_LIBS_USER>", .libPaths()))
suppressPackageStartupMessages(library(data.table))
setwd("<MINTCNA_ROOT>/ATAC_RNA")
S <- Sys.getenv("AR_SAMPLE", "P111T2"); TAG <- "h15a005e020"
ARM <- Sys.getenv("AR_ARM", "multi")   # "atac" makes the x axis RNA-free; see above
SUF <- if (identical(ARM, "multi")) "" else paste0("_", ARM)

## ---- per-bin copy-number difference between the clones ----------------------
e <- new.env(); load(sprintf("Result/MintCNA/%s_%s/%s_results_%s.RData", S, ARM, S, TAG), e)
cat(sprintf("copy-number states from the %s arm, tag %s\n", toupper(ARM), TAG))
Sm <- base::get("cna_state_matrix", e); sub <- base::get("subclones", e)
rf <- as.data.table(base::get("ref", e))
keep <- fread(sprintf("Data/cell_lists/%s_cells.tsv", S))[cell_type == "Malignant", barcode]
i <- colnames(Sm) %in% keep; bcn <- colnames(Sm)[i]; Sm <- Sm[, i, drop = FALSE]; sub <- sub[i]
o  <- as.integer(names(sort(table(sub), decreasing = TRUE))[1:2])
cl <- ifelse(sub == o[1], "C1", ifelse(sub == o[2], "C2", NA_character_))

md <- as.data.table(readRDS("Data/synapse/snrna_malignant_meta_data.RDS"))[care_id == S]
md[, bare := sub("^.*-([ACGT]+-[0-9]+)$", "\\1", CellID)]
st <- md$cell_state[match(bcn, md$bare)]
dom <- sapply(c("C1","C2"), function(k) names(sort(table(st[cl == k]), decreasing = TRUE))[1])
OPC <- names(dom)[grepl("OPC", dom)]; AC <- names(dom)[grepl("AC", dom)]
cat(sprintf("%s: %s is OPC-like, %s is AC-like\n", S, OPC, AC))

# d_cn is on the same OPC-minus-AC orientation as d_atac and d_rna
rf[, d_cn := rowMeans(Sm[, cl == OPC & !is.na(cl), drop = FALSE]) -
             rowMeans(Sm[, cl == AC  & !is.na(cl), drop = FALSE])]
rf[, chr := sub("^chr", "", as.character(chr))]

## ---- attach d_cn to every gene already scored in both modalities ------------
ar <- fread(sprintf("Result/qc/atac_rna_clone_%s.tsv", S))
ga <- fread(sprintf("Data/mintcna_input/%s/gene_annotation.tsv", S))[, .(gene_id, chr, start)]
ga[, chr := sub("^chr", "", chr)]
ar <- merge(ar, ga, by = "gene_id")
setkey(rf, chr, start, end)
hit <- foverlaps(ar[, .(chr, s = start, e = start, gene_id)],
                 rf[, .(chr, start, end, d_cn)], by.x = c("chr","s","e"),
                 type = "within", nomatch = 0L)
ar <- merge(ar, hit[, .(gene_id, d_cn)], by = "gene_id")
cat(sprintf("genes with a copy-number difference attached: %d\n", nrow(ar)))

## ---- Q1: does copy number move accessibility and expression? ---------------
ar[, cn_bin := cut(d_cn, breaks = c(-Inf, -0.4, -0.2, 0.2, 0.4, Inf),
                   labels = c("<= -0.4", "-0.4..-0.2", "neutral", "0.2..0.4", ">= 0.4"))]
cat("\n=== Q1  DOSAGE: mean shift by copy-number difference (OPC minus AC) ===\n")
tab <- ar[, .(n_genes = .N,
              mean_dATAC = round(mean(d_atac), 3),
              mean_dRNA  = round(mean(d_rna), 3)), by = cn_bin][order(cn_bin)]
print(tab, row.names = FALSE)
neu <- ar[cn_bin == "neutral"]
for (m in c("d_atac", "d_rna")) {
  up <- ar[d_cn >= 0.4][[m]]; dn <- ar[d_cn <= -0.4][[m]]
  cat(sprintf("  %s : gained %+.3f vs neutral %+.3f (Wilcoxon P = %.2g) | lost %+.3f (P = %.2g)\n",
              m, mean(up), mean(neu[[m]]),
              wilcox.test(up, neu[[m]])$p.value, mean(dn),
              wilcox.test(dn, neu[[m]])$p.value))
}
cat(sprintf("\n  Spearman d_cn vs d_atac = %.3f ; d_cn vs d_rna = %.3f\n",
            cor(ar$d_cn, ar$d_atac, method = "spearman"),
            cor(ar$d_cn, ar$d_rna,  method = "spearman")))

## ---- Q2: are the lineage genes explained by dosage? ------------------------
LIN <- data.table(
  gene = c("OLIG1","OLIG2","SOX10","PDGFRA","CSPG4","NKX2-2","SOX8","ASCL1","DLL3","BCAN",
           "GFAP","AQP4","ALDH1L1","S100B","SLC1A2","SLC1A3","NFIA","NFIB","NFIX","SOX9","CLU","ID3"),
  fam  = c(rep("OPC", 10), rep("AC", 12)))
L <- merge(ar, LIN, by = "gene")
cat("\n=== Q2  CONFOUNDING: is each lineage gene in a region that differs in CN? ===\n")
print(L[order(fam, -d_rna), .(gene, fam, d_cn = round(d_cn, 3),
                              d_atac = round(d_atac, 2), d_rna = round(d_rna, 2))],
      row.names = FALSE)
cat(sprintf("\n  lineage genes with |d_cn| > 0.2 : %d of %d\n", L[abs(d_cn) > 0.2, .N], nrow(L)))
cat(sprintf("  median |d_cn|, lineage genes %.3f vs all genes %.3f\n",
            median(abs(L$d_cn)), median(abs(ar$d_cn))))

# how much of the RNA shift does dosage account for, genome-wide and for these genes
fit <- lm(d_rna ~ d_cn, data = ar)
cat(sprintf("\n  genome-wide  d_rna ~ d_cn :  slope %+.3f, R^2 = %.4f\n",
            coef(fit)[2], summary(fit)$r.squared))
L[, resid := d_rna - predict(fit, .SD)]
cat("  lineage-gene RNA shift remaining after subtracting the dosage fit:\n")
print(L[order(fam, -resid), .(gene, fam, d_rna = round(d_rna,2),
                              dosage_expected = round(predict(fit, .SD), 3),
                              residual = round(resid, 2))], row.names = FALSE)

fwrite(ar, sprintf("Result/qc/cna_dosage_link_%s%s.tsv", S, SUF), sep = "\t")
cat(sprintf("\nwrote Result/qc/cna_dosage_link_%s%s.tsv\n", S, SUF))
