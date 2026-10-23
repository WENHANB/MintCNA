#!/usr/bin/env Rscript
# Do the two CNA-defined subclones differ in TF motif accessibility?
#
# The chain this closes:
#   ATAC coverage  -> CNA -> subclone          (AG-DCAE latent, no RNA involved)
#   subclone       -> cell state               (the paper's own RNA annotation)
#   subclone       -> TF motif accessibility   (ATAC again, at promoter resolution)
#
# The prediction is not open-ended. In P100T1 clone C1 is 90% OPC-like and C2 is
# 99% AC-like, so if the subclones are real lineages rather than a clustering
# artefact, C1 should open at oligodendrocyte-lineage motifs (OLIG1/2, SOX10) and
# C2 at astrocyte ones (the NFI family, SOX9). Naming the expected TFs before
# looking is what makes a hit meaningful.
#
# WINDOWS, not peaks. No peak calls exist for this cohort and MACS2 is not
# installed; promoters are where TF binding concentrates, so each gene gets
# TSS +/- 2 kb (strand-aware). This is less sensitive than a real peak set -- distal enhancers
# are missed -- and the result should be read as a floor, not a survey.
.libPaths(c("<R_LIBS_USER>", .libPaths()))
suppressPackageStartupMessages({
  library(data.table); library(GenomicRanges); library(SummarizedExperiment)
  library(chromVAR); library(motifmatchr); library(JASPAR2020); library(TFBSTools)
  library(BSgenome.Hsapiens.UCSC.hg38)
})
setwd("<MINTCNA_ROOT>/ATAC_RNA")
set.seed(1)
S <- Sys.getenv("TF_SAMPLE", "P100T1"); TAG <- "h15a005e020"
CHRS <- paste0("chr", 1:22)

## ---- subclone labels, from the CNA side only --------------------------------
e <- new.env(); load(sprintf("Result/MintCNA/%s_multi/%s_results_%s.RData", S, S, TAG), e)
Sm <- base::get("cna_state_matrix", e); sub <- base::get("subclones", e)
keep <- fread(sprintf("Data/cell_lists/%s_cells.tsv", S))[cell_type == "Malignant", barcode]
i <- colnames(Sm) %in% keep; bc <- colnames(Sm)[i]; sub <- sub[i]
tb <- sort(table(sub), decreasing = TRUE); o <- as.integer(names(tb)[1:2])
cl <- ifelse(sub == o[1], "C1", ifelse(sub == o[2], "C2", NA_character_))
use <- bc[!is.na(cl)]; grp <- cl[!is.na(cl)]
cat(sprintf("%s: C1 %d, C2 %d nuclei\n", S, sum(grp == "C1"), sum(grp == "C2")))

## ---- promoter windows -------------------------------------------------------
ga <- fread(sprintf("Data/mintcna_input/%s/gene_annotation.tsv", S))[chr %in% CHRS]
# TSS by STRAND. gene_annotation.tsv holds only the genomic left/right coordinates,
# so `start` is the 5' end for plus-strand genes but the 3' end for minus-strand
# ones (49.5% of genes). Until 2026-09-25 this line read `tss := start` under the
# comment "annotation is already 5' start", which was false. Same fix, same source
# (EnsDb.Hsapiens.v86), as atac_to_rna_clone.R on 2026-09-14, so Fig 5F and S7A
# now use the SAME promoter windows.
suppressPackageStartupMessages(library(EnsDb.Hsapiens.v86))
gs <- as.data.table(ensembldb::genes(EnsDb.Hsapiens.v86, columns = c("gene_id", "seq_strand")))[, .(gene_id, strand = as.character(strand))]
ga <- merge(ga, gs, by = "gene_id", all.x = TRUE)
stopifnot(all(ga$strand %in% c("+", "-")))
ga[, tss := fifelse(strand == "+", start, end)]
win <- GRanges(ga$chr, IRanges(pmax(ga$tss - 2000, 1), ga$tss + 2000), gene = ga$gene_name)
# genomic order: the strand join above re-orders rows by gene_id, and chromVAR's
# filterPeaks(non_overlapping = TRUE) refuses unsorted ranges
win <- sort(win)
win <- win[!duplicated(paste(seqnames(win), start(win)))]
cat(sprintf("promoter windows: %d\n", length(win)))

## ---- count fragments per window per cell ------------------------------------
cat("reading fragments ...\n")
fr <- fread(cmd = sprintf("zcat %s", shQuote(sprintf("Data/atac_filtered/%s_fragments.tsv.gz", S))),
            header = FALSE, select = 1:4, col.names = c("chr", "start", "end", "bc"),
            showProgress = FALSE)
fr <- fr[bc %in% use & chr %in% CHRS]
cat(sprintf("fragments in the two clones: %s\n", format(nrow(fr), big.mark = ",")))
fr[, `:=`(mid = (start + end) %/% 2L)]
wdt <- data.table(chr = as.character(seqnames(win)), s = start(win), e = end(win),
                  j = seq_along(win)); setkey(wdt, chr, s, e)
fr[, `:=`(s = mid, e = mid)]
ov <- foverlaps(fr[, .(chr, s, e, bc)], wdt, type = "within", nomatch = 0L)
cnt <- ov[, .N, by = .(j, bc)]
cells <- sort(unique(use))
M <- Matrix::sparseMatrix(i = cnt$j, j = match(cnt$bc, cells), x = cnt$N,
                          dims = c(length(win), length(cells)))
cat(sprintf("count matrix: %d windows x %d cells, %.1f%% non-zero\n",
            nrow(M), ncol(M), 100*Matrix::nnzero(M)/length(M)))

se <- SummarizedExperiment(assays = list(counts = M), rowRanges = win,
                           colData = DataFrame(clone = grp[match(cells, use)]))
se <- addGCBias(se, genome = BSgenome.Hsapiens.UCSC.hg38)
se <- filterPeaks(se, non_overlapping = TRUE)
cat(sprintf("after filterPeaks: %d windows\n", nrow(se)))

## ---- motifs and deviations --------------------------------------------------
mo <- getMatrixSet(JASPAR2020, list(species = 9606, collection = "CORE"))
cat(sprintf("JASPAR CORE motifs: %d\n", length(mo)))
ix <- matchMotifs(mo, se, genome = BSgenome.Hsapiens.UCSC.hg38)
dev <- computeDeviations(object = se, annotations = ix)
z <- assays(dev)$z
nm <- vapply(mo, TFBSTools::name, character(1))[rownames(z)]

g <- colData(dev)$clone
res <- rbindlist(lapply(seq_len(nrow(z)), function(k) {
  a <- z[k, g == "C1"]; b <- z[k, g == "C2"]
  a <- a[is.finite(a)]; b <- b[is.finite(b)]
  if (length(a) < 20 || length(b) < 20) return(NULL)
  data.table(motif = rownames(z)[k], TF = nm[k],
             mean_C1 = mean(a), mean_C2 = mean(b), diff = mean(a) - mean(b),
             p = wilcox.test(a, b)$p.value)
}))
res[, fdr := p.adjust(p, "BH")]
setorder(res, -diff)
fwrite(res, sprintf("Result/qc/tf_motif_%s.tsv", S), sep = "\t")

cat("\n=== most open in C1 (OPC-like) ===\n")
print(head(res[fdr < 0.05, .(TF, diff = round(diff,2), fdr = signif(fdr,2))], 12), row.names = FALSE)
cat("\n=== most open in C2 (AC-like) ===\n")
print(head(res[fdr < 0.05][order(diff), .(TF, diff = round(diff,2), fdr = signif(fdr,2))], 12),
      row.names = FALSE)
cat("\n=== the TFs named in advance ===\n")
want <- "OLIG1|OLIG2|SOX10|SOX8|NKX2-2|NFIA|NFIB|NFIC|NFIX|SOX9|STAT3|RFX4"
print(res[grepl(want, TF), .(TF, diff = round(diff,2), fdr = signif(fdr,2))][order(-diff)],
      row.names = FALSE)
