#!/usr/bin/env Rscript
# Panel-a sample choice, both criteria at once:
#   (1) is the subclone split just depth?  AUC of nFrags separating the clones
#   (2) do the clones correspond to the paper's own malignant cell states?
.libPaths(c("<R_LIBS_USER>", .libPaths()))
suppressPackageStartupMessages(library(data.table))
setwd("<MINTCNA_ROOT>/ATAC_RNA")
TAG <- "h15a005e020"
md <- as.data.table(readRDS("Data/synapse/snrna_malignant_meta_data.RDS"))
auc <- function(x, g) { r <- rank(x); n1 <- sum(g); n0 <- sum(!g)
  a <- (sum(r[g]) - n1*(n1+1)/2)/(n1*n0); max(a, 1-a) }

for (S in c("P111T3","P111T2","P100T1")) {
  e <- new.env(); load(sprintf("Result/MintCNA/%s_multi/%s_results_%s.RData", S, S, TAG), e)
  Sm <- base::get("cna_state_matrix", e); sub <- base::get("subclones", e)
  cl_tab <- fread(sprintf("Data/cell_lists/%s_cells.tsv", S))
  keep <- cl_tab[cell_type == "Malignant", barcode]
  i <- colnames(Sm) %in% keep; bc <- colnames(Sm)[i]; sub <- sub[i]
  g <- sub == as.integer(names(sort(table(sub), decreasing = TRUE))[1])
  nf <- cl_tab$nFrags[match(bc, cl_tab$barcode)]
  a  <- auc(nf, g)

  mm <- md[care_id == S]
  mm[, bare := sub("^.*-([ACGT]+-[0-9]+)$", "\\1", CellID)]
  st <- mm$cell_state[match(bc, mm$bare)]
  cat(sprintf("\n=== %s ===  n = %d  (C1 %d / C2 %d)\n", S, length(bc), sum(g), sum(!g)))
  cat(sprintf("  depth AUC (nFrags vs clone) = %.3f\n", a))
  cat(sprintf("  matched to paper annotation: %d / %d\n", sum(!is.na(st)), length(bc)))
  if (sum(!is.na(st)) > 50) {
    t <- table(clone = ifelse(g, "C1", "C2"), state = st)
    print(round(100*prop.table(t, 1), 1))
    cat(sprintf("  chi-square p = %.3g\n", suppressWarnings(chisq.test(t)$p.value)))
  }
}
