.libPaths(c("<R_LIBS_USER>", .libPaths()))
suppressPackageStartupMessages(library(data.table))
ROOT <- "<MINTCNA_ROOT>"; TAG <- "h15a005e020"
for (S in c("P111T2","P111T3")) {
  e <- new.env()
  load(sprintf("%s/ATAC_RNA/Result/MintCNA/%s_multi/%s_results_%s.RData",ROOT,S,S,TAG), e)
  cna <- get("cna_state_matrix", e); sub <- get("subclones", e)
  ref <- as.data.frame(get("ref", e))
  keep <- fread(sprintf("%s/ATAC_RNA/Data/cell_lists/%s_cells.tsv",ROOT,S))[cell_type=="Malignant", barcode]
  i <- colnames(cna) %in% keep; cna <- cna[,i,drop=FALSE]; sub <- sub[i]
  g <- sub == as.integer(names(sort(table(sub), decreasing=TRUE))[1])
  d <- rowMeans(cna[,g,drop=FALSE]) - rowMeans(cna[,!g,drop=FALSE])
  ref$chr <- sub("^chr","",as.character(ref$chr))
  dt <- data.table(chr=ref$chr, start=ref$start, d=d)[chr %in% as.character(1:22)]
  dt[, arm := ifelse(start < 1e8, "p/prox", "q/dist")]
  top <- dt[, .(n=.N, mean_d=mean(d), frac_big=mean(abs(d)>0.5)), by=.(chr,arm)][
            order(-abs(mean_d))][1:6]
  cat("\n==", S, "== largest clone differences (C1 - C2, state units):\n")
  print(top, row.names=FALSE)
  cat(sprintf("  longest run of |d|>0.4: %d consecutive bins\n",
      { r <- rle(abs(dt$d)>0.4); if(any(r$values)) max(r$lengths[r$values]) else 0L }))
}
