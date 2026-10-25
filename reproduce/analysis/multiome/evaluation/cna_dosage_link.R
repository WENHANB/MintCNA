#!/usr/bin/env Rscript
# Copy-number dosage between the two MintCNA subclones of one sample (Supplementary Fig. 8b,c): for every
# gene scored by atac_to_rna_clone.R, the difference in mean copy-number state between the OPC-like and the
# AC-like clone (d_cn) is attached from the bin containing the gene start, so promoter accessibility and
# expression differences can be read against the copy-number difference. By default d_cn comes from the
# ATAC-only arm, whose states were called without the RNA (AR_ARM=multi uses the ATAC + RNA arm).
#   Rscript reproduce/analysis/multiome/evaluation/cna_dosage_link.R      (AR_SAMPLE selects the sample)
# Input : ATAC_RNA/Result/MintCNA/<S>_<arm>/<S>_results_h15a005e020.RData, Data/cell_lists/<S>_cells.tsv,
#         Data/synapse/snrna_malignant_meta_data.RDS, Result/qc/atac_rna_clone_<S>.tsv, Data/mintcna_input/<S>/gene_annotation.tsv
# Output: ATAC_RNA/Result/qc/cna_dosage_link_<S>_atac.tsv  (or cna_dosage_link_<S>.tsv for the multi arm)
# Run from the root of the MintCNA repository with MINTCNA_DATA set (see reproduce/README.md).
suppressPackageStartupMessages(library(data.table))
ROOT <- Sys.getenv("MINTCNA_DATA"); if (!nzchar(ROOT)) stop("set MINTCNA_DATA (see reproduce/README.md)")
setwd(file.path(ROOT, "ATAC_RNA"))
S <- Sys.getenv("AR_SAMPLE", "P111T2"); TAG <- "h15a005e020"
ARM <- Sys.getenv("AR_ARM", "atac")
SUF <- if (identical(ARM, "multi")) "" else paste0("_", ARM)

e <- new.env(); load(sprintf("Result/MintCNA/%s_%s/%s_results_%s.RData", S, ARM, S, TAG), e)
Sm <- base::get("cna_state_matrix", e); sub <- base::get("subclones", e); rf <- as.data.table(base::get("ref", e))
keep <- fread(sprintf("Data/cell_lists/%s_cells.tsv", S))[cell_type == "Malignant", barcode]
i <- colnames(Sm) %in% keep; bcn <- colnames(Sm)[i]; Sm <- Sm[, i, drop = FALSE]; sub <- sub[i]
o <- as.integer(names(sort(table(sub), decreasing = TRUE))[1:2])
cl <- ifelse(sub == o[1], "C1", ifelse(sub == o[2], "C2", NA_character_))
md <- as.data.table(readRDS("Data/synapse/snrna_malignant_meta_data.RDS"))[care_id == S]
md[, bare := sub("^.*-([ACGT]+-[0-9]+)$", "\\1", CellID)]
st <- md$cell_state[match(bcn, md$bare)]
dom <- sapply(c("C1", "C2"), function(k) names(sort(table(st[cl == k]), decreasing = TRUE))[1])
OPC <- names(dom)[grepl("OPC", dom)]; AC <- names(dom)[grepl("AC", dom)]
cat(sprintf("%s (%s arm): %s is OPC-like, %s is AC-like\n", S, ARM, OPC, AC))

# per-bin copy-number difference, OPC-like minus AC-like
rf[, d_cn := rowMeans(Sm[, cl == OPC & !is.na(cl), drop = FALSE]) - rowMeans(Sm[, cl == AC & !is.na(cl), drop = FALSE])]
rf[, chr := sub("^chr", "", as.character(chr))]

ar <- fread(sprintf("Result/qc/atac_rna_clone_%s.tsv", S))
ga <- fread(sprintf("Data/mintcna_input/%s/gene_annotation.tsv", S))[, .(gene_id, chr, start)]
ga[, chr := sub("^chr", "", chr)]
ar <- merge(ar, ga, by = "gene_id")
setkey(rf, chr, start, end)
hit <- foverlaps(ar[, .(chr, s = start, e = start, gene_id)], rf[, .(chr, start, end, d_cn)], by.x = c("chr", "s", "e"), type = "within", nomatch = 0L)
ar <- merge(ar, hit[, .(gene_id, d_cn)], by = "gene_id")

ar[, cn_bin := cut(d_cn, breaks = c(-Inf, -0.4, -0.2, 0.2, 0.4, Inf), labels = c("<= -0.4", "-0.4..-0.2", "neutral", "0.2..0.4", ">= 0.4"))]
print(ar[, .(n_genes = .N, mean_dATAC = round(mean(d_atac), 3), mean_dRNA = round(mean(d_rna), 3)), by = cn_bin][order(cn_bin)], row.names = FALSE)
cat(sprintf("Spearman d_cn vs d_atac = %.3f; d_cn vs d_rna = %.3f\n",
            cor(ar$d_cn, ar$d_atac, method = "spearman"), cor(ar$d_cn, ar$d_rna, method = "spearman")))
fwrite(ar, sprintf("Result/qc/cna_dosage_link_%s%s.tsv", S, SUF), sep = "\t")
