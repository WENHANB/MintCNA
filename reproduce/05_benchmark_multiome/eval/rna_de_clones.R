#!/usr/bin/env Rscript
# Gene-level differential expression between the two MintCNA subclones, C1 vs C2.
# Per-gene Wilcoxon on log-normalised counts across nuclei, BH-adjusted. Written
# 2026-09-06 because panel F carried only clone MEANS per gene and no test; the
# question was whether a gene-level comparison is far more significant than the
# motif-level one (it is expected to be: it compares AC-like vs OPC-like nuclei).
.libPaths("<R_LIBS_USER>")
suppressPackageStartupMessages({ library(data.table); library(Matrix) })
setwd("<MINTCNA_ROOT>/ATAC_RNA")
S <- Sys.getenv("AR_SAMPLE", "P111T2"); GEX <- Sys.getenv("AR_GEX", "GSM9625483_SJ20-2"); TAG <- "h15a005e020"
e <- new.env(); load(sprintf("Result/MintCNA/%s_multi/%s_results_%s.RData", S, S, TAG), e)
Sm <- get("cna_state_matrix", e); sub <- get("subclones", e)
keep <- fread(sprintf("Data/cell_lists/%s_cells.tsv", S))[cell_type == "Malignant", barcode]
i <- colnames(Sm) %in% keep; bc <- colnames(Sm)[i]; sub <- sub[i]
o <- as.integer(names(sort(table(sub), decreasing = TRUE))[1:2])
cl <- ifelse(sub == o[1], "C1", ifelse(sub == o[2], "C2", NA)); use <- !is.na(cl); bc <- bc[use]; cl <- cl[use]
pre <- sprintf("Data/gex_matrices/%s", GEX)
X <- readMM(sprintf("%s_matrix.mtx.gz", pre))
feat <- fread(sprintf("%s_features.tsv.gz", pre), header = FALSE); bcs <- fread(sprintf("%s_barcodes.tsv.gz", pre), header = FALSE)[[1]]
rownames(X) <- make.unique(feat[[2]]); colnames(X) <- bcs
bare <- sub("^.*-([ACGT]+-[0-9]+)$", "\\1", colnames(X)); m <- match(bc, bare)
X <- X[, m[!is.na(m)]]; cl <- cl[!is.na(m)]
cat(sprintf("%s: %d nuclei with RNA (C1 %d, C2 %d), %d genes\n", S, ncol(X), sum(cl=="C1"), sum(cl=="C2"), nrow(X)))
# log-normalise (Seurat-style: per-cell 1e4 scaling, log1p)
X <- X[rowSums(X > 0) >= 0.05 * ncol(X), ]                       # expressed in >= 5 % of nuclei
Xn <- log1p(t(t(X) / colSums(X)) * 1e4)
a <- cl == "C1"; b <- cl == "C2"
res <- data.table(gene = rownames(Xn),
  mean_C1 = rowMeans(Xn[, a]), mean_C2 = rowMeans(Xn[, b]),
  pct_C1 = rowMeans(Xn[, a] > 0), pct_C2 = rowMeans(Xn[, b] > 0))
res[, log2FC_C2_vs_C1 := log2((expm1(mean_C2) + 1) / (expm1(mean_C1) + 1))]
res[, p := vapply(seq_len(nrow(Xn)), function(k) wilcox.test(Xn[k, a], Xn[k, b])$p.value, numeric(1))]
res[, fdr := p.adjust(p, "BH")]; setorder(res, fdr)
LIN <- c("OLIG1","OLIG2","SOX10","PDGFRA","CSPG4","NKX2-2","SOX8","ASCL1","DLL3","BCAN",
         "GFAP","AQP4","ALDH1L1","S100B","SLC1A2","SLC1A3","NFIA","NFIB","NFIX","SOX9","CLU","ID3")
res[, lineage := gene %in% LIN]; res[, rank := .I]
dir.create("Result/qc", showWarnings = FALSE)
fwrite(res, sprintf("Result/qc/rna_de_clones_%s.tsv", S), sep = "\t")
cat(sprintf("\ngenes tested %d; FDR<0.05: %d; FDR<0.05 & |log2FC|>1: %d; FDR<1e-10: %d\n",
  nrow(res), res[fdr < .05, .N], res[fdr < .05 & abs(log2FC_C2_vs_C1) > 1, .N], res[fdr < 1e-10, .N]))
cat("\ntop 12 by FDR:\n"); print(res[1:12, .(gene, log2FC = round(log2FC_C2_vs_C1, 2), fdr = signif(fdr, 2), lineage)], row.names = FALSE)
cat("\nthe 22 lineage genes:\n"); print(res[lineage == TRUE, .(gene, log2FC = round(log2FC_C2_vs_C1, 2), fdr = signif(fdr, 2), rank)][order(rank)], row.names = FALSE)
