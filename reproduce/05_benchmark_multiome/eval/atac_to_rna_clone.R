#!/usr/bin/env Rscript
# Does the accessibility difference between the two CNA-defined subclones show up
# in gene EXPRESSION in the same nuclei?
#
# Panel D of Figure 4 shows the clones differ at lineage TF motifs. That is an
# ATAC-only statement. This script closes the loop on the other modality of the
# same assay:
#
#   subclone (from CNA)  ->  promoter accessibility (ATAC, TSS +/- 2 kb, strand-aware)
#                        ->  gene expression        (RNA, same nuclei, same cells)
#
# and asks whether the two move together gene by gene. Both axes are the SAME
# contrast (OPC-like clone minus AC-like clone), so a positive relation means an
# opening promoter is an expressed gene rather than an isolated ATAC artefact.
#
# Barcodes are shared between the two modalities in 10x Multiome, and genes are
# joined on Ensembl ID rather than symbol, so neither join is approximate.
.libPaths(c("<R_LIBS_USER>", .libPaths()))
suppressPackageStartupMessages({
  library(data.table); library(Matrix)
})
setwd("<MINTCNA_ROOT>/ATAC_RNA")
S    <- Sys.getenv("AR_SAMPLE", "P111T2")
GEX  <- Sys.getenv("AR_GEX",    "GSM9625483_SJ20-2")
TAG  <- "h15a005e020"
CHRS <- paste0("chr", 1:22)
OUT  <- "Result/qc"

## ---- subclones, and which one is which lineage ------------------------------
e <- new.env()
load(sprintf("Result/MintCNA/%s_multi/%s_results_%s.RData", S, S, TAG), e)
Sm  <- base::get("cna_state_matrix", e); sub <- base::get("subclones", e)
keep <- fread(sprintf("Data/cell_lists/%s_cells.tsv", S))[cell_type == "Malignant", barcode]
i   <- colnames(Sm) %in% keep; bc <- colnames(Sm)[i]; sub <- sub[i]
tb  <- sort(table(sub), decreasing = TRUE); o <- as.integer(names(tb)[1:2])
cl  <- ifelse(sub == o[1], "C1", ifelse(sub == o[2], "C2", NA_character_))
use <- bc[!is.na(cl)]; grp <- cl[!is.na(cl)]

md <- as.data.table(readRDS("Data/synapse/snrna_malignant_meta_data.RDS"))[care_id == S]
md[, bare := sub("^.*-([ACGT]+-[0-9]+)$", "\\1", CellID)]
st <- md$cell_state[match(use, md$bare)]
# name each clone by its dominant RNA state, so both axes below read in lineage
# terms rather than in arbitrary clone numbers
domin <- sapply(c("C1","C2"), function(k) names(sort(table(st[grp == k]), decreasing = TRUE))[1])
cat(sprintf("%s: C1 n=%d dominant %s | C2 n=%d dominant %s\n",
            S, sum(grp=="C1"), domin[["C1"]], sum(grp=="C2"), domin[["C2"]]))
OPC <- names(domin)[grepl("OPC", domin)]; AC <- names(domin)[grepl("AC", domin)]
stopifnot(length(OPC) == 1, length(AC) == 1)
cat(sprintf("  contrast direction: %s (OPC-like) minus %s (AC-like)\n", OPC, AC))

## ---- ATAC: promoter accessibility per clone ---------------------------------
ga <- fread(sprintf("Data/mintcna_input/%s/gene_annotation.tsv", S))[chr %in% CHRS]
ga <- ga[!duplicated(gene_id)]
# TSS +/- 2 kb, STRAND-AWARE. gene_annotation.tsv carries only the genomic left
# and right coordinates; until 2026-09-14 the window was start +/- 2 kb for every
# gene, which for the 49.5% of genes on the minus strand is the 3' end, not the
# promoter (Fig 5F: 2,873 of 8,869 genes). Strand comes from the same EnsDb the
# annotation was built from, joined on Ensembl ID.
suppressPackageStartupMessages(library(EnsDb.Hsapiens.v86))
gs <- as.data.table(ensembldb::genes(EnsDb.Hsapiens.v86, columns = c("gene_id", "seq_strand")))[, .(gene_id, strand = as.character(strand))]
ga <- merge(ga, gs, by = "gene_id", all.x = TRUE)
stopifnot(all(ga$strand %in% c("+", "-")))
ga[, tss := fifelse(strand == "+", start, end)]
ga[, `:=`(s = pmax(tss - 2000, 1), e = tss + 2000, j = .I)]
cat(sprintf("promoter windows: TSS +/- 2 kb; %d genes, %d on the minus strand\n", nrow(ga), ga[strand == "-", .N]))
setkey(ga, chr, s, e)
cat("reading fragments ...\n")
fr <- fread(cmd = sprintf("zcat %s", shQuote(sprintf("Data/atac_filtered/%s_fragments.tsv.gz", S))),
            header = FALSE, select = 1:4, col.names = c("chr","start","end","bc"),
            showProgress = FALSE)
fr <- fr[bc %in% use & chr %in% CHRS]
fr[, mid := (start + end) %/% 2L][, `:=`(s = mid, e = mid)]
ov <- foverlaps(fr[, .(chr, s, e, bc)], ga[, .(chr, s, e, j)], type = "within", nomatch = 0L)
acc <- ov[, .(n = .N), by = .(j, clone = grp[match(bc, use)])]
A <- dcast(acc, j ~ clone, value.var = "n", fill = 0)
A <- merge(data.table(j = ga$j, gene_id = ga$gene_id, gene = ga$gene_name), A, by = "j", all.x = TRUE)
for (k in c("C1","C2")) A[[k]][is.na(A[[k]])] <- 0
cat(sprintf("promoters: %d; ATAC fragments in promoters: %s\n",
            nrow(A), format(sum(A$C1 + A$C2), big.mark = ",")))

## ---- RNA: expression per clone, same nuclei ---------------------------------
pre <- sprintf("Data/gex_matrices/%s", GEX)
rb  <- fread(sprintf("zcat %s_barcodes.tsv.gz", pre), header = FALSE)$V1
rf  <- fread(sprintf("zcat %s_features.tsv.gz", pre), header = FALSE)
X   <- readMM(sprintf("%s_matrix.mtx.gz", pre))
cat(sprintf("RNA matrix: %d genes x %d cells\n", nrow(X), ncol(X)))
sel <- match(use, rb); ok <- !is.na(sel)
cat(sprintf("  %d of %d clone-labelled nuclei found in the RNA matrix\n", sum(ok), length(use)))
Xs  <- X[, sel[ok], drop = FALSE]; gs <- grp[ok]
R <- data.table(gene_id = rf$V1,
                C1 = as.numeric(Matrix::rowSums(Xs[, gs == "C1", drop = FALSE])),
                C2 = as.numeric(Matrix::rowSums(Xs[, gs == "C2", drop = FALSE])))

## ---- one table, one contrast ------------------------------------------------
# CPM inside each modality and each clone, then log2 of the OPC/AC ratio. A
# pseudocount of 1 CPM keeps zero-in-one-clone genes finite; genes too sparse in
# either modality are dropped rather than being carried as noise at the axes.
cpm <- function(v) v / sum(v) * 1e6
A[, `:=`(a_opc = cpm(get(OPC)), a_ac = cpm(get(AC)), a_tot = C1 + C2)]
R[, `:=`(r_opc = cpm(get(OPC)), r_ac = cpm(get(AC)), r_tot = C1 + C2)]
D <- merge(A[, .(gene_id, gene, a_opc, a_ac, a_tot)],
           R[, .(gene_id, r_opc, r_ac, r_tot)], by = "gene_id")
MIN_A <- 20; MIN_R <- 20
n0 <- nrow(D)
D <- D[a_tot >= MIN_A & r_tot >= MIN_R]
cat(sprintf("genes: %d joined, %d kept (ATAC >= %d and RNA >= %d counts)\n",
            n0, nrow(D), MIN_A, MIN_R))
D[, `:=`(d_atac = log2((a_opc + 1)/(a_ac + 1)),
         d_rna  = log2((r_opc + 1)/(r_ac + 1)))]

ct <- cor.test(D$d_atac, D$d_rna, method = "spearman")
cp <- cor.test(D$d_atac, D$d_rna, method = "pearson")
cat(sprintf("\nATAC promoter vs RNA, same contrast, %d genes:\n  Spearman rho = %.3f\n  Pearson  r   = %.3f\n",
            nrow(D), unname(ct$estimate), unname(cp$estimate)))
# concordance: of the genes that clearly move in RNA, how many move the same way in ATAC
mv <- D[abs(d_rna) > 1]
cat(sprintf("  of %d genes with |log2FC RNA| > 1, %.1f%% change promoter accessibility in the same direction\n",
            nrow(mv), 100*mean(sign(mv$d_atac) == sign(mv$d_rna))))

fwrite(D, file.path(OUT, sprintf("atac_rna_clone_%s.tsv", S)), sep = "\t")
cat(sprintf("\nwrote %s/atac_rna_clone_%s.tsv\n", OUT, S))

LIN <- list(
  OPC = c("OLIG1","OLIG2","SOX10","PDGFRA","CSPG4","NKX2-2","SOX8","ASCL1","DLL3","BCAN"),
  AC  = c("GFAP","AQP4","ALDH1L1","S100B","SLC1A2","SLC1A3","NFIA","NFIB","NFIX","SOX9","CLU","ID3"))
cat("\nlineage genes (log2FC OPC-like minus AC-like):\n")
for (k in names(LIN)) {
  x <- D[gene %in% LIN[[k]], .(gene, d_atac = round(d_atac,2), d_rna = round(d_rna,2))][order(-d_rna)]
  cat(sprintf("  -- %s lineage --\n", k)); print(x, row.names = FALSE)
}
