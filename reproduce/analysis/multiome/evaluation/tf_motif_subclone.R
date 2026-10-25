#!/usr/bin/env Rscript
# Transcription-factor motif accessibility of the two MintCNA subclones (Supplementary Fig. 8a): fragments in
# promoter windows (TSS +/- 2 kb, strand-aware, EnsDb.Hsapiens.v86) -> chromVAR deviation z-scores for the
# JASPAR 2020 CORE vertebrate motifs -> per-motif Wilcoxon test between the clones, BH-adjusted.
#   Rscript reproduce/analysis/multiome/evaluation/tf_motif_subclone.R      (TF_SAMPLE selects the sample, default P100T1)
# Input : ATAC_RNA/Result/MintCNA/<S>_multi/<S>_results_h15a005e020.RData, Data/cell_lists/<S>_cells.tsv,
#         Data/mintcna_input/<S>/gene_annotation.tsv, Data/atac_filtered/<S>_fragments.tsv.gz
# Output: ATAC_RNA/Result/qc/tf_motif_<S>.tsv
# Run from the root of the MintCNA repository with MINTCNA_DATA set (see reproduce/README.md).
suppressPackageStartupMessages({
  library(data.table); library(GenomicRanges); library(SummarizedExperiment); library(chromVAR); library(motifmatchr)
  library(JASPAR2020); library(TFBSTools); library(BSgenome.Hsapiens.UCSC.hg38); library(EnsDb.Hsapiens.v86)
})
ROOT <- Sys.getenv("MINTCNA_DATA"); if (!nzchar(ROOT)) stop("set MINTCNA_DATA (see reproduce/README.md)")
setwd(file.path(ROOT, "ATAC_RNA"))
set.seed(1)
S <- Sys.getenv("TF_SAMPLE", "P100T1"); TAG <- "h15a005e020"
CHRS <- paste0("chr", 1:22)

# subclone labels (C1 = larger clone) of the malignant nuclei
e <- new.env(); load(sprintf("Result/MintCNA/%s_multi/%s_results_%s.RData", S, S, TAG), e)
Sm <- base::get("cna_state_matrix", e); sub <- base::get("subclones", e)
keep <- fread(sprintf("Data/cell_lists/%s_cells.tsv", S))[cell_type == "Malignant", barcode]
i <- colnames(Sm) %in% keep; bc <- colnames(Sm)[i]; sub <- sub[i]
o <- as.integer(names(sort(table(sub), decreasing = TRUE))[1:2])
cl <- ifelse(sub == o[1], "C1", ifelse(sub == o[2], "C2", NA_character_))
use <- bc[!is.na(cl)]; grp <- cl[!is.na(cl)]
cat(sprintf("%s: C1 %d, C2 %d nuclei\n", S, sum(grp == "C1"), sum(grp == "C2")))

# promoter windows, TSS by strand
ga <- fread(sprintf("Data/mintcna_input/%s/gene_annotation.tsv", S))[chr %in% CHRS]
gs <- as.data.table(ensembldb::genes(EnsDb.Hsapiens.v86, columns = c("gene_id", "seq_strand")))[, .(gene_id, strand = as.character(strand))]
ga <- merge(ga, gs, by = "gene_id", all.x = TRUE)
stopifnot(all(ga$strand %in% c("+", "-")))
ga[, tss := fifelse(strand == "+", start, end)]
win <- sort(GRanges(ga$chr, IRanges(pmax(ga$tss - 2000, 1), ga$tss + 2000), gene = ga$gene_name))
win <- win[!duplicated(paste(seqnames(win), start(win)))]

# fragments (midpoints) per window and nucleus
fr <- fread(cmd = sprintf("zcat %s", shQuote(sprintf("Data/atac_filtered/%s_fragments.tsv.gz", S))),
            header = FALSE, select = 1:4, col.names = c("chr", "start", "end", "bc"), showProgress = FALSE)
fr <- fr[bc %in% use & chr %in% CHRS]
fr[, `:=`(s = (start + end) %/% 2L)][, e := s]
wdt <- data.table(chr = as.character(seqnames(win)), s = start(win), e = end(win), j = seq_along(win)); setkey(wdt, chr, s, e)
cnt <- foverlaps(fr[, .(chr, s, e, bc)], wdt, type = "within", nomatch = 0L)[, .N, by = .(j, bc)]
cells <- sort(unique(use))
M <- Matrix::sparseMatrix(i = cnt$j, j = match(cnt$bc, cells), x = cnt$N, dims = c(length(win), length(cells)))

se <- SummarizedExperiment(assays = list(counts = M), rowRanges = win, colData = DataFrame(clone = grp[match(cells, use)]))
se <- addGCBias(se, genome = BSgenome.Hsapiens.UCSC.hg38)
se <- filterPeaks(se, non_overlapping = TRUE)
mo <- getMatrixSet(JASPAR2020, list(species = 9606, collection = "CORE"))
dev <- computeDeviations(object = se, annotations = matchMotifs(mo, se, genome = BSgenome.Hsapiens.UCSC.hg38))
z <- assays(dev)$z
nm <- vapply(mo, TFBSTools::name, character(1))[rownames(z)]
g <- colData(dev)$clone
res <- rbindlist(lapply(seq_len(nrow(z)), function(k) {
  a <- z[k, g == "C1"]; b <- z[k, g == "C2"]; a <- a[is.finite(a)]; b <- b[is.finite(b)]
  if (length(a) < 20 || length(b) < 20) return(NULL)
  data.table(motif = rownames(z)[k], TF = nm[k], mean_C1 = mean(a), mean_C2 = mean(b), diff = mean(a) - mean(b), p = wilcox.test(a, b)$p.value)
}))
res[, fdr := p.adjust(p, "BH")]
setorder(res, -diff)
fwrite(res, sprintf("Result/qc/tf_motif_%s.tsv", S), sep = "\t")
cat("\nmost open in C1:\n"); print(head(res[fdr < 0.05, .(TF, diff = round(diff, 2), fdr = signif(fdr, 2))], 12), row.names = FALSE)
cat("\nmost open in C2:\n"); print(head(res[fdr < 0.05][order(diff), .(TF, diff = round(diff, 2), fdr = signif(fdr, 2))], 12), row.names = FALSE)
