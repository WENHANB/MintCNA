#!/usr/bin/env Rscript
# Is the subclone split copy number, or is it chromatin accessibility wearing a
# copy-number costume?
#
# The clusters come from the AG-DCAE latent, whose input is ATAC coverage, and
# coverage reflects BOTH copy number and accessibility. If the latent is really
# separating cell states, then "the CNA subclone predicts the cell state" is a
# tautology. The test: recluster using ONLY the discrete called states (1/2/3),
# which are the output of a change-point model and carry no continuous coverage
# signal, and ask whether the same two groups come back.
suppressPackageStartupMessages(library(data.table))
setwd("<MINTCNA_ROOT>/ATAC_RNA")
set.seed(1); TAG <- "h15a005e020"
md <- as.data.table(readRDS("Data/synapse/snrna_malignant_meta_data.RDS"))

for (S in c("P100T1","P111T2","P105T3")) {
  e <- new.env(); load(sprintf("Result/MintCNA/%s_multi/%s_results_%s.RData", S, S, TAG), e)
  Sm <- base::get("cna_state_matrix", e); sub <- base::get("subclones", e)
  keep <- fread(sprintf("Data/cell_lists/%s_cells.tsv", S))[cell_type == "Malignant", barcode]
  i <- colnames(Sm) %in% keep; bc <- colnames(Sm)[i]; Sm <- Sm[, i, drop = FALSE]; sub <- sub[i]
  tb <- sort(table(sub), decreasing = TRUE); o <- as.integer(names(tb)[1:2])
  latent_cl <- ifelse(sub == o[1], "C1", ifelse(sub == o[2], "C2", NA))
  ok <- !is.na(latent_cl)

  # cluster on the DISCRETE states alone: k-means, k = 2, on the bins that vary
  X <- t(Sm[, ok, drop = FALSE])
  v <- apply(X, 2, sd); X <- X[, which(v > 0), drop = FALSE]
  km <- kmeans(X, centers = 2, nstart = 50, iter.max = 1000)
  state_cl <- paste0("K", km$cluster)

  a <- table(latent = latent_cl[ok], state = state_cl)
  # agreement under the best of the two labellings
  agree <- max(sum(diag(a)), sum(a) - sum(diag(a))) / sum(a)
  ari <- {
    n <- sum(a); si <- sum(choose(a, 2)); sr <- sum(choose(rowSums(a), 2))
    sc <- sum(choose(colSums(a), 2)); ex <- sr*sc/choose(n, 2)
    (si - ex) / ((sr + sc)/2 - ex)
  }
  cat(sprintf("\n=== %s ===\n", S))
  print(a)
  cat(sprintf("  agreement %.1f%%   ARI %.3f\n", 100*agree, ari))

  # and does the state-only clustering ALSO predict the RNA cell state?
  mm <- md[care_id == S]; mm[, bare := sub("^.*-([ACGT]+-[0-9]+)$", "\\1", CellID)]
  st <- mm$cell_state[match(bc[ok], mm$bare)]
  if (sum(!is.na(st)) > 50) {
    t2 <- table(clone = state_cl, state = st)
    cat("  state-only clusters vs the paper's RNA annotation (% of row):\n")
    print(round(100*prop.table(t2, 1), 1))
    cat(sprintf("  chi-square p = %.3g\n", suppressWarnings(chisq.test(t2)$p.value)))
  }
}
