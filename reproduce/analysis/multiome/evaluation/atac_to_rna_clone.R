#!/usr/bin/env Rscript
# Promoter accessibility versus gene expression between the two MintCNA subclones of one sample (Fig. 5f):
# per gene, log2 fold change of pseudobulk promoter fragments (TSS +/- 2 kb, strand-aware) and of pseudobulk
# RNA counts, both as OPC-like clone over AC-like clone (clones are named by their dominant cell state).
#   Rscript reproduce/analysis/multiome/evaluation/atac_to_rna_clone.R      (AR_SAMPLE, AR_GEX select the sample)
# Input : ATAC_RNA/Result/MintCNA/<S>_multi/<S>_results_h15a005e020.RData, Data/cell_lists/<S>_cells.tsv,
#         Data/synapse/snrna_malignant_meta_data.RDS, Data/mintcna_input/<S>/gene_annotation.tsv,
#         Data/atac_filtered/<S>_fragments.tsv.gz, Data/gex_matrices/<GEX>_*
# Output: ATAC_RNA/Result/qc/atac_rna_clone_<S>.tsv  (gene_id, gene, pseudobulk CPMs, d_atac, d_rna)
# Run from the root of the MintCNA repository with MINTCNA_DATA set (see reproduce/README.md).
suppressPackageStartupMessages({ library(data.table); library(Matrix); library(EnsDb.Hsapiens.v86) })
ROOT <- Sys.getenv("MINTCNA_DATA"); if (!nzchar(ROOT)) stop("set MINTCNA_DATA (see reproduce/README.md)")
setwd(file.path(ROOT, "ATAC_RNA"))
S <- Sys.getenv("AR_SAMPLE", "P111T2"); GEX <- Sys.getenv("AR_GEX", "GSM9625483_SJ20-2")
TAG <- "h15a005e020"; CHRS <- paste0("chr", 1:22); OUT <- "Result/qc"

# subclones and their dominant cell state
e <- new.env(); load(sprintf("Result/MintCNA/%s_multi/%s_results_%s.RData", S, S, TAG), e)
Sm <- base::get("cna_state_matrix", e); sub <- base::get("subclones", e)
keep <- fread(sprintf("Data/cell_lists/%s_cells.tsv", S))[cell_type == "Malignant", barcode]
i <- colnames(Sm) %in% keep; bc <- colnames(Sm)[i]; sub <- sub[i]
o <- as.integer(names(sort(table(sub), decreasing = TRUE))[1:2])
cl <- ifelse(sub == o[1], "C1", ifelse(sub == o[2], "C2", NA_character_))
use <- bc[!is.na(cl)]; grp <- cl[!is.na(cl)]
md <- as.data.table(readRDS("Data/synapse/snrna_malignant_meta_data.RDS"))[care_id == S]
md[, bare := sub("^.*-([ACGT]+-[0-9]+)$", "\\1", CellID)]
st <- md$cell_state[match(use, md$bare)]
domin <- sapply(c("C1", "C2"), function(k) names(sort(table(st[grp == k]), decreasing = TRUE))[1])
OPC <- names(domin)[grepl("OPC", domin)]; AC <- names(domin)[grepl("AC", domin)]
stopifnot(length(OPC) == 1, length(AC) == 1)
cat(sprintf("%s: %s is OPC-like (n = %d), %s is AC-like (n = %d)\n", S, OPC, sum(grp == OPC), AC, sum(grp == AC)))

# ATAC: promoter fragments per clone
ga <- fread(sprintf("Data/mintcna_input/%s/gene_annotation.tsv", S))[chr %in% CHRS][!duplicated(gene_id)]
gs <- as.data.table(ensembldb::genes(EnsDb.Hsapiens.v86, columns = c("gene_id", "seq_strand")))[, .(gene_id, strand = as.character(strand))]
ga <- merge(ga, gs, by = "gene_id", all.x = TRUE)
stopifnot(all(ga$strand %in% c("+", "-")))
ga[, tss := fifelse(strand == "+", start, end)]
ga[, `:=`(s = pmax(tss - 2000, 1), e = tss + 2000, j = .I)]
setkey(ga, chr, s, e)
fr <- fread(cmd = sprintf("zcat %s", shQuote(sprintf("Data/atac_filtered/%s_fragments.tsv.gz", S))),
            header = FALSE, select = 1:4, col.names = c("chr", "start", "end", "bc"), showProgress = FALSE)
fr <- fr[bc %in% use & chr %in% CHRS]
fr[, mid := (start + end) %/% 2L][, `:=`(s = mid, e = mid)]
ov <- foverlaps(fr[, .(chr, s, e, bc)], ga[, .(chr, s, e, j)], type = "within", nomatch = 0L)
A <- dcast(ov[, .(n = .N), by = .(j, clone = grp[match(bc, use)])], j ~ clone, value.var = "n", fill = 0)
A <- merge(data.table(j = ga$j, gene_id = ga$gene_id, gene = ga$gene_name), A, by = "j", all.x = TRUE)
for (k in c("C1", "C2")) A[[k]][is.na(A[[k]])] <- 0

# RNA: counts per clone, same nuclei
pre <- sprintf("Data/gex_matrices/%s", GEX)
rb <- fread(sprintf("zcat %s_barcodes.tsv.gz", pre), header = FALSE)$V1
rf <- fread(sprintf("zcat %s_features.tsv.gz", pre), header = FALSE)
X  <- readMM(sprintf("%s_matrix.mtx.gz", pre))
sel <- match(use, rb); ok <- !is.na(sel)
Xs <- X[, sel[ok], drop = FALSE]; gs <- grp[ok]
R <- data.table(gene_id = rf$V1, C1 = as.numeric(Matrix::rowSums(Xs[, gs == "C1", drop = FALSE])),
                C2 = as.numeric(Matrix::rowSums(Xs[, gs == "C2", drop = FALSE])))

# CPM within modality and clone, log2 ratio OPC / AC with a pseudocount of 1 CPM; sparse genes dropped
cpm <- function(v) v / sum(v) * 1e6
A[, `:=`(a_opc = cpm(get(OPC)), a_ac = cpm(get(AC)), a_tot = C1 + C2)]
R[, `:=`(r_opc = cpm(get(OPC)), r_ac = cpm(get(AC)), r_tot = C1 + C2)]
D <- merge(A[, .(gene_id, gene, a_opc, a_ac, a_tot)], R[, .(gene_id, r_opc, r_ac, r_tot)], by = "gene_id")
D <- D[a_tot >= 20 & r_tot >= 20]
D[, `:=`(d_atac = log2((a_opc + 1) / (a_ac + 1)), d_rna = log2((r_opc + 1) / (r_ac + 1)))]
cat(sprintf("%d genes: Spearman rho = %.3f, Pearson r = %.3f\n", nrow(D),
            cor(D$d_atac, D$d_rna, method = "spearman"), cor(D$d_atac, D$d_rna)))
fwrite(D, file.path(OUT, sprintf("atac_rna_clone_%s.tsv", S)), sep = "\t")
