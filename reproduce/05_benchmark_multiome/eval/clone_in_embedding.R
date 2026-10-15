#!/usr/bin/env Rscript
# Do the CNA-defined subclones separate in the RNA expression embedding?
#
# The field's cheapest downstream check (CopyKAT Fig 5d-f, CONGAS+ Fig 6D-E), and
# its inverse is also reported as a result (epiAneufinder, Copy-scAT: "clones are
# invisible in the embedding"). The clones here come from the ATAC/CNA side only,
# so this asks whether a partition made without any transcript is visible in the
# transcriptome. Because the clones are 95% / 87% lineage-pure, separation is
# expected; the point of running it is to put a number on it and to show the
# UMAP a reader expects to see.
#
# Separation is quantified two ways, both in PCA space (not on the UMAP, which is
# for looking at, not for measuring):
#   silhouette   how compact each clone is relative to the other
#   kNN accuracy leave-one-out 15-NN vote on clone label, against the majority-
#                class rate as the floor
.libPaths(c("<R_LIBS_USER>", .libPaths()))
suppressPackageStartupMessages({
  library(data.table); library(Matrix); library(Seurat); library(cluster); library(FNN)
})
setwd("<MINTCNA_ROOT>/ATAC_RNA")
set.seed(1)
S <- Sys.getenv("AR_SAMPLE", "P111T2"); GEX <- Sys.getenv("AR_GEX", "GSM9625483_SJ20-2")
TAG <- "h15a005e020"; OUT <- "Result/qc"

## ---- clone labels, from CNA only ------------------------------------------
e <- new.env()
load(sprintf("Result/MintCNA/%s_multi/%s_results_%s.RData", S, S, TAG), e)
Sm <- base::get("cna_state_matrix", e); sub <- base::get("subclones", e)
keep <- fread(sprintf("Data/cell_lists/%s_cells.tsv", S))[cell_type == "Malignant", barcode]
i <- colnames(Sm) %in% keep; bc <- colnames(Sm)[i]; sub <- sub[i]
o  <- as.integer(names(sort(table(sub), decreasing = TRUE))[1:2])
cl <- ifelse(sub == o[1], "C1", ifelse(sub == o[2], "C2", NA_character_))
use <- bc[!is.na(cl)]; grp <- cl[!is.na(cl)]

md <- as.data.table(readRDS("Data/synapse/snrna_malignant_meta_data.RDS"))[care_id == S]
md[, bare := sub("^.*-([ACGT]+-[0-9]+)$", "\\1", CellID)]
st <- md$cell_state[match(use, md$bare)]

## ---- RNA of the same nuclei -------------------------------------------------
pre <- sprintf("Data/gex_matrices/%s", GEX)
rb  <- fread(cmd = sprintf("zcat %s_barcodes.tsv.gz", pre), header = FALSE)$V1
rf  <- fread(cmd = sprintf("zcat %s_features.tsv.gz", pre), header = FALSE)
X   <- readMM(sprintf("%s_matrix.mtx.gz", pre)); X <- as(X, "CsparseMatrix")
rownames(X) <- make.unique(rf$V2); colnames(X) <- rb
sel <- match(use, rb); ok <- !is.na(sel)
X <- X[, sel[ok]]; grp <- grp[ok]; st <- st[ok]; use <- use[ok]
cat(sprintf("%s: %d clone-labelled nuclei with RNA (C1 %d, C2 %d)\n", S, ncol(X), sum(grp == "C1"), sum(grp == "C2")))

so <- CreateSeuratObject(X, min.cells = 3)
so$clone <- grp; so$state <- st
so <- NormalizeData(so, verbose = FALSE)
so <- FindVariableFeatures(so, nfeatures = 2000, verbose = FALSE)
so <- ScaleData(so, verbose = FALSE)
so <- RunPCA(so, npcs = 30, verbose = FALSE)
so <- RunUMAP(so, dims = 1:30, verbose = FALSE)
pc <- Embeddings(so, "pca")[, 1:30]

## ---- separation, measured in PCA space -------------------------------------
sil <- silhouette(as.integer(factor(grp)), dist(pc))
sil_clone <- mean(sil[, "sil_width"])
knn <- FNN::knn.cv(pc, factor(grp), k = 15)
acc <- mean(as.character(knn) == grp); floor <- max(table(grp)) / length(grp)
cat(sprintf("\nclone separation in RNA PCA space:\n"))
cat(sprintf("  mean silhouette (C1 vs C2)      = %.3f\n", sil_clone))
cat(sprintf("  15-NN leave-one-out accuracy    = %.3f   (majority-class floor %.3f)\n", acc, floor))

# the same measures for the study's own cell state, as the reference point: the
# states were DEFINED from this RNA, so they are the ceiling for what the RNA can separate
ok2 <- !is.na(st) & st %in% c("OPC-like", "AC-like")
sil_st <- mean(silhouette(as.integer(factor(st[ok2])), dist(pc[ok2, ]))[, "sil_width"])
knn2 <- FNN::knn.cv(pc[ok2, ], factor(st[ok2]), k = 15)
cat(sprintf("  for comparison, the study's own OPC-like vs AC-like labels on the same cells:\n"))
cat(sprintf("  mean silhouette                 = %.3f\n", sil_st))
cat(sprintf("  15-NN accuracy                  = %.3f\n", mean(as.character(knn2) == st[ok2])))

um <- as.data.table(Embeddings(so, "umap"), keep.rownames = "barcode")
setnames(um, c("barcode", "umap1", "umap2"))
um[, `:=`(clone = grp, state = st)]
fwrite(um, file.path(OUT, sprintf("clone_embedding_%s.tsv", S)), sep = "\t")
fwrite(data.table(sample = S, n = ncol(X), sil_clone = sil_clone, knn_acc_clone = acc,
                  majority_floor = floor, sil_state = sil_st),
       file.path(OUT, sprintf("clone_embedding_%s_summary.tsv", S)), sep = "\t")
cat(sprintf("\nwrote %s/clone_embedding_%s{,_summary}.tsv\n", OUT, S))
