#!/usr/bin/env Rscript
# Dosage: does the CNA difference between two subclones show up in their RNA?
#
# This is the test only a multiome assay can run. The subclones are defined from
# the AG-DCAE latent, which sees ATAC coverage and nothing else; the expression
# comes from the RNA of the SAME nuclei. Neither the clustering nor the CNA call
# has seen a single transcript, so an expression shift at exactly the regions
# where the clones differ in copy number is independent corroboration rather
# than a restatement.
#
# The comparison is within-gene and between-clones, so gene identity, length and
# GC are held constant; the control is the set of genes in regions where the two
# clones do NOT differ, which absorbs any global difference in library size or
# cell state that is not dosage.
suppressPackageStartupMessages(library(data.table))
setwd("<MINTCNA_ROOT>/ATAC_RNA")
S <- Sys.getenv("DOSE_SAMPLE", "P100T1"); TAG <- "h15a005e020"

e <- new.env(); load(sprintf("Result/MintCNA/%s_multi/%s_results_%s.RData", S, S, TAG), e)
Sm <- base::get("cna_state_matrix", e); rf <- as.data.table(base::get("ref", e))
sub <- base::get("subclones", e)
keep <- fread(sprintf("Data/cell_lists/%s_cells.tsv", S))[cell_type == "Malignant", barcode]
i <- colnames(Sm) %in% keep; bc <- colnames(Sm)[i]; Sm <- Sm[, i, drop = FALSE]; sub <- sub[i]
tb <- sort(table(sub), decreasing = TRUE); o <- as.integer(names(tb)[1:2])
cl <- ifelse(sub == o[1], "C1", ifelse(sub == o[2], "C2", NA_character_))
m1 <- rowMeans(Sm[, cl == "C1" & !is.na(cl), drop = FALSE])
m2 <- rowMeans(Sm[, cl == "C2" & !is.na(cl), drop = FALSE])
rf[, `:=`(d = m1 - m2, ad = abs(m1 - m2))]

ga <- fread(sprintf("Data/mintcna_input/%s/gene_annotation.tsv", S))
rna <- fread(sprintf("Data/mintcna_input/%s/rna_counts.tsv", S))
gid <- rna[[1]]; R <- as.matrix(rna[, -1]); rownames(R) <- gid
common <- intersect(colnames(R), bc[!is.na(cl)])
R <- R[, common, drop = FALSE]; g <- cl[match(common, bc)]
cpm <- log2(t(t(R) / pmax(colSums(R), 1)) * 1e6 + 1)      # per-cell library size out
cat(sprintf("%s: %d genes x %d nuclei;  C1 %d, C2 %d\n", S, nrow(cpm), ncol(cpm),
            sum(g == "C1"), sum(g == "C2")))

# assign every gene the copy-number difference of the bin it sits in
setkey(rf, chr, start, end)
ga[, `:=`(mid = (start + end)/2)]
gb <- rf[ga, on = .(chr, start <= mid, end >= mid), mult = "first",
         .(gene_id = i.gene_id, gene_name = i.gene_name, chr = x.chr, d = x.d, ad = x.ad)]
gb <- gb[!is.na(d) & gene_id %in% rownames(cpm)]
gb[, grp := fifelse(ad > 0.3, fifelse(d > 0, "higher in C1", "higher in C2"), "no CN difference")]
cat("\ngenes by the copy-number difference of their bin:\n"); print(table(gb$grp))

# per gene: mean expression in each clone, and the difference
mu1 <- rowMeans(cpm[gb$gene_id, g == "C1", drop = FALSE])
mu2 <- rowMeans(cpm[gb$gene_id, g == "C2", drop = FALSE])
gb[, `:=`(e1 = mu1, e2 = mu2, de = mu1 - mu2)]
gb <- gb[e1 + e2 > 0.2]                                    # drop genes off in both clones
cat(sprintf("%d genes expressed in at least one clone\n\n", nrow(gb)))

cat("=== RNA difference (C1 - C2, log2 CPM) by the DNA difference of the region ===\n\n")
print(gb[, .(n_genes = .N, median_dRNA = round(median(de), 4),
             mean_dRNA = round(mean(de), 4)), by = grp][order(-median_dRNA)], row.names = FALSE)
w <- wilcox.test(gb[grp == "higher in C1", de], gb[grp == "no CN difference", de])
w2 <- wilcox.test(gb[grp == "higher in C2", de], gb[grp == "no CN difference", de])
cat(sprintf("\nhigher-in-C1 regions vs no-difference regions: Wilcoxon p = %.3g\n", w$p.value))
cat(sprintf("higher-in-C2 regions vs no-difference regions: Wilcoxon p = %.3g\n", w2$p.value))
cat("\n=== by chromosome, where the clones differ most ===\n\n")
print(gb[ad > 0.3, .(n_genes = .N, median_dCN = round(median(d), 3),
                     median_dRNA = round(median(de), 4)), by = chr][order(-abs(median_dCN))][1:8],
      row.names = FALSE)
fwrite(gb, sprintf("Result/qc/dosage_%s.tsv", S), sep = "\t")
