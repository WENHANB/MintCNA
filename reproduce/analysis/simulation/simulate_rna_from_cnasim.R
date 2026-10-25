#!/usr/bin/env Rscript
# Paired scRNA-seq counts for one CNAsim setting (Splatter 1.28.0, set.seed(1)).
#   Rscript reproduce/analysis/simulation/simulate_rna_from_cnasim.R <setting>
# Input : <setting>/profiles.tsv (CNAsim copy-number truth), Simulation/ref_hg38/gencode.hg38.annotation.gtf
# Output: <setting>/rna_counts_gene_by_cell.tsv  (copy-number effect in every gene, concordance p = 1)
#         <setting>/rna_counts_gene_by_cell_g1_p<k>.tsv, k = 1..4  (concordance 25 / 50 / 75 / 100 %)
#         <setting>/genes_hg38.tsv  (gene_id, chrom, start, end)
# Run from the root of the MintCNA repository with MINTCNA_DATA set (see reproduce/README.md).
suppressPackageStartupMessages({
  library(data.table); library(Matrix); library(rtracklayer); library(GenomicRanges); library(splatter)
})
ROOT <- Sys.getenv("MINTCNA_DATA"); if (!nzchar(ROOT)) stop("set MINTCNA_DATA (see reproduce/README.md)")
args <- commandArgs(trailingOnly = TRUE)
BASE <- if (length(args) >= 1) args[1] else "c3_N100_100kb_l1e4"
CNA_DIR <- file.path(ROOT, "Simulation/cnasim", BASE)
GTF <- file.path(ROOT, "Simulation/ref_hg38/gencode.hg38.annotation.gtf")
set.seed(1)

# ---- 1. bins x cells total copy number from profiles.tsv --------------------------------
profiles <- fread(file.path(CNA_DIR, "profiles.tsv"))
colnames(profiles)[5] <- "cn_states"
setnames(profiles, tolower(names(profiles)))
profiles[, cn_total := {                       # "a.b" = the two haplotype copy numbers
  x <- trimws(as.character(cn_states))
  ifelse(grepl("\\.", x),
         rowSums(matrix(as.numeric(unlist(tstrsplit(x, "\\.", perl = TRUE))), ncol = 2, byrow = FALSE)),
         round(as.numeric(x))) |> as.integer()
}]
profiles[!grepl("^chr", chrom), chrom := paste0("chr", chrom)]
profiles[, bin_id := paste0(chrom, ":", start, "-", end)]
profiles[, cell_num := as.integer(sub("^cell", "", cell))]
profiles[, chr_num := as.integer(sub("^chr", "", chrom))]
profiles <- profiles[grepl("^cell\\d+$", cell)]   # cells only, no ancestral nodes
setorder(profiles, chr_num, start, end, cell_num)

cells <- profiles[, unique(.SD[order(cell_num)]$cell), .SDcols = c("cell", "cell_num")]
CNdt <- dcast(profiles, bin_id ~ cell, value.var = "cn_total")
CNmat <- as.matrix(CNdt[, ..cells]); rownames(CNmat) <- CNdt$bin_id
ord_bin <- order(as.integer(sub("^chr", "", sub(":.*$", "", rownames(CNmat)))),
                 as.integer(sub("^.*:(\\d+)-.*$", "\\1", rownames(CNmat))))
CNmat <- CNmat[ord_bin, , drop = FALSE]

# ---- 2. genes (hg38 autosomes) and their overlap with the bins -------------------------
gtf <- import(GTF)
genes_gr <- gtf[gtf$type == "gene"]
genes_gr <- unique(genes_gr[seqnames(genes_gr) %in% paste0("chr", 1:22)])

bin_id <- rownames(CNmat)
bins_gr <- GRanges(seqnames = sub(":.*$", "", bin_id),
                   ranges = IRanges(start = as.integer(sub("^.*:(\\d+)-.*$", "\\1", bin_id)) + 1,
                                    end = as.integer(sub("^.*-(\\d+)$", "\\1", bin_id))))
hits <- findOverlaps(genes_gr, bins_gr, ignore.strand = TRUE)
A <- sparseMatrix(i = queryHits(hits), j = subjectHits(hits), x = 1,
                  dims = c(length(genes_gr), length(bins_gr)))       # genes x bins
den <- as.numeric(A %*% rep(1, length(bins_gr))); den[den == 0] <- 1

# ---- 3. baseline counts, then the copy-number dosage effect ------------------------------
set.seed(1)
sim <- splatSimulate(nGenes = length(genes_gr), batchCells = ncol(CNmat),
                     method = "groups", group.prob = 1, verbose = FALSE)
Y0 <- as.matrix(counts(sim))
rownames(Y0) <- mcols(genes_gr)$gene_id
colnames(Y0) <- colnames(CNmat)

# Per-gene dosage = mean (CN / 2) of the overlapping bins; a Bernoulli(p) mask decides which
# gene x cell entries carry the effect (p = concordance between the RNA and the DNA signal).
inject <- function(p, gamma = 1) {
  Dbin <- (CNmat / 2)^gamma
  Dbin[CNmat == 2L] <- 1.0
  Gdos <- (A %*% Dbin) / den
  mask <- matrix(rbinom(length(Gdos), 1, p), nrow = nrow(Gdos))
  Gdos_final <- mask * (1 + gamma * (Gdos - 1)) + (1 - mask)
  Y <- round(Y0 * Gdos_final); Y[Y < 0] <- 0
  Y
}
write_counts <- function(Y, name)
  fwrite(data.table(gene_id = rownames(Y), as.matrix(Y)), file.path(CNA_DIR, name), sep = "\t")

write_counts(inject(1), "rna_counts_gene_by_cell.tsv")
fwrite(data.table(gene_id = mcols(genes_gr)$gene_id, chrom = as.character(seqnames(genes_gr)),
                  start = start(genes_gr), end = end(genes_gr)),
       file.path(CNA_DIR, "genes_hg38.tsv"), sep = "\t")

# concordance sweep (Supplementary Fig. 5a): p = 1/4, 2/4, 3/4, 4/4
for (k in 1:4) write_counts(inject(k / 4), sprintf("rna_counts_gene_by_cell_g1_p%d.tsv", k))
