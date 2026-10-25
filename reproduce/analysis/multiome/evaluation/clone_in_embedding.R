#!/usr/bin/env Rscript
# Do the two MintCNA subclones (defined from the ATAC-trained latent, no RNA involved) separate in the RNA
# expression embedding of the same nuclei? UMAP coordinates for Fig. 5d, plus the silhouette and the
# leave-one-out 15-nearest-neighbour accuracy of the clone labels in PCA space.
#   Rscript reproduce/analysis/multiome/evaluation/clone_in_embedding.R      (AR_SAMPLE, AR_GEX select the sample)
# Input : ATAC_RNA/Result/MintCNA/<S>_multi/<S>_results_h15a005e020.RData, Data/cell_lists/<S>_cells.tsv,
#         Data/synapse/snrna_malignant_meta_data.RDS (cell states), Data/gex_matrices/<GEX>_*
# Output: ATAC_RNA/Result/qc/clone_embedding_<S>.tsv, clone_embedding_<S>_summary.tsv
# Run from the root of the MintCNA repository with MINTCNA_DATA set (see reproduce/README.md).
suppressPackageStartupMessages({ library(data.table); library(Matrix); library(Seurat); library(cluster); library(FNN) })
ROOT <- Sys.getenv("MINTCNA_DATA"); if (!nzchar(ROOT)) stop("set MINTCNA_DATA (see reproduce/README.md)")
setwd(file.path(ROOT, "ATAC_RNA"))
set.seed(1)
S <- Sys.getenv("AR_SAMPLE", "P111T2"); GEX <- Sys.getenv("AR_GEX", "GSM9625483_SJ20-2")
TAG <- "h15a005e020"; OUT <- "Result/qc"

# clone labels (C1 = larger clone) of the malignant nuclei
e <- new.env(); load(sprintf("Result/MintCNA/%s_multi/%s_results_%s.RData", S, S, TAG), e)
Sm <- base::get("cna_state_matrix", e); sub <- base::get("subclones", e)
keep <- fread(sprintf("Data/cell_lists/%s_cells.tsv", S))[cell_type == "Malignant", barcode]
i <- colnames(Sm) %in% keep; bc <- colnames(Sm)[i]; sub <- sub[i]
o  <- as.integer(names(sort(table(sub), decreasing = TRUE))[1:2])
cl <- ifelse(sub == o[1], "C1", ifelse(sub == o[2], "C2", NA_character_))
use <- bc[!is.na(cl)]; grp <- cl[!is.na(cl)]
md <- as.data.table(readRDS("Data/synapse/snrna_malignant_meta_data.RDS"))[care_id == S]
md[, bare := sub("^.*-([ACGT]+-[0-9]+)$", "\\1", CellID)]
st <- md$cell_state[match(use, md$bare)]

# RNA of the same nuclei: standard Seurat embedding
pre <- sprintf("Data/gex_matrices/%s", GEX)
rb <- fread(cmd = sprintf("zcat %s_barcodes.tsv.gz", pre), header = FALSE)$V1
rf <- fread(cmd = sprintf("zcat %s_features.tsv.gz", pre), header = FALSE)
X  <- as(readMM(sprintf("%s_matrix.mtx.gz", pre)), "CsparseMatrix")
rownames(X) <- make.unique(rf$V2); colnames(X) <- rb
sel <- match(use, rb); ok <- !is.na(sel)
X <- X[, sel[ok]]; grp <- grp[ok]; st <- st[ok]
so <- CreateSeuratObject(X, min.cells = 3)
so$clone <- grp; so$state <- st
so <- NormalizeData(so, verbose = FALSE); so <- FindVariableFeatures(so, nfeatures = 2000, verbose = FALSE)
so <- ScaleData(so, verbose = FALSE); so <- RunPCA(so, npcs = 30, verbose = FALSE); so <- RunUMAP(so, dims = 1:30, verbose = FALSE)
pc <- Embeddings(so, "pca")[, 1:30]

# separation of the clones in PCA space, with the study's own cell states as the reference point
sil_clone <- mean(silhouette(as.integer(factor(grp)), dist(pc))[, "sil_width"])
acc <- mean(as.character(FNN::knn.cv(pc, factor(grp), k = 15)) == grp)
floor <- max(table(grp)) / length(grp)
ok2 <- !is.na(st) & st %in% c("OPC-like", "AC-like")
sil_st <- mean(silhouette(as.integer(factor(st[ok2])), dist(pc[ok2, ]))[, "sil_width"])
acc_st <- mean(as.character(FNN::knn.cv(pc[ok2, ], factor(st[ok2]), k = 15)) == st[ok2])
cat(sprintf("%s: %d nuclei; clone silhouette %.3f, 15-NN accuracy %.3f (majority floor %.3f); cell-state silhouette %.3f, 15-NN accuracy %.3f\n",
            S, ncol(X), sil_clone, acc, floor, sil_st, acc_st))

um <- as.data.table(Embeddings(so, "umap"), keep.rownames = "barcode"); setnames(um, c("barcode", "umap1", "umap2"))
um[, `:=`(clone = grp, state = st)]
fwrite(um, file.path(OUT, sprintf("clone_embedding_%s.tsv", S)), sep = "\t")
fwrite(data.table(sample = S, n = ncol(X), sil_clone = sil_clone, knn_acc_clone = acc, majority_floor = floor, sil_state = sil_st,
                  knn_acc_state = acc_st),
       file.path(OUT, sprintf("clone_embedding_%s_summary.tsv", S)), sep = "\t")
