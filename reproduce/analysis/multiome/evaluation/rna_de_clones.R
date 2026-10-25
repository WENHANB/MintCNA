#!/usr/bin/env Rscript
# Differential expression between the two MintCNA subclones of one sample (Fig. 5e): per-gene two-sided
# Wilcoxon rank-sum test on log-normalised counts (per-nucleus scaling to 10,000, log1p) across nuclei,
# genes expressed in at least 5 % of the nuclei, BH-adjusted.
#   Rscript reproduce/analysis/multiome/evaluation/rna_de_clones.R      (AR_SAMPLE, AR_GEX select the sample)
# Input : ATAC_RNA/Result/MintCNA/<S>_multi/<S>_results_h15a005e020.RData, Data/cell_lists/<S>_cells.tsv, Data/gex_matrices/<GEX>_*
# Output: ATAC_RNA/Result/qc/rna_de_clones_<S>.tsv
# Run from the root of the MintCNA repository with MINTCNA_DATA set (see reproduce/README.md).
suppressPackageStartupMessages({ library(data.table); library(Matrix) })
ROOT <- Sys.getenv("MINTCNA_DATA"); if (!nzchar(ROOT)) stop("set MINTCNA_DATA (see reproduce/README.md)")
setwd(file.path(ROOT, "ATAC_RNA"))
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
m <- match(bc, sub("^.*-([ACGT]+-[0-9]+)$", "\\1", colnames(X)))
X <- X[, m[!is.na(m)]]; cl <- cl[!is.na(m)]
cat(sprintf("%s: %d nuclei with RNA (C1 %d, C2 %d)\n", S, ncol(X), sum(cl == "C1"), sum(cl == "C2")))

X <- X[rowSums(X > 0) >= 0.05 * ncol(X), ]
Xn <- log1p(t(t(X) / colSums(X)) * 1e4)
a <- cl == "C1"; b <- cl == "C2"
res <- data.table(gene = rownames(Xn), mean_C1 = rowMeans(Xn[, a]), mean_C2 = rowMeans(Xn[, b]),
                  pct_C1 = rowMeans(Xn[, a] > 0), pct_C2 = rowMeans(Xn[, b] > 0))
res[, log2FC_C2_vs_C1 := log2((expm1(mean_C2) + 1) / (expm1(mean_C1) + 1))]
res[, p := vapply(seq_len(nrow(Xn)), function(k) wilcox.test(Xn[k, a], Xn[k, b])$p.value, numeric(1))]
res[, fdr := p.adjust(p, "BH")]; setorder(res, fdr)
LIN <- c("OLIG1","OLIG2","SOX10","PDGFRA","CSPG4","NKX2-2","SOX8","ASCL1","DLL3","BCAN",
         "GFAP","AQP4","ALDH1L1","S100B","SLC1A2","SLC1A3","NFIA","NFIB","NFIX","SOX9","CLU","ID3")
res[, lineage := gene %in% LIN]; res[, rank := .I]
dir.create("Result/qc", showWarnings = FALSE)
fwrite(res, sprintf("Result/qc/rna_de_clones_%s.tsv", S), sep = "\t")
cat(sprintf("genes tested %d; FDR < 0.05: %d\n", nrow(res), res[fdr < .05, .N]))
