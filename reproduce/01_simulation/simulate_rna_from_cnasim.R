#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(rtracklayer)
  library(GenomicRanges)
  library(splatter)
})

# ---------------------------
# PATHS
# ---------------------------
args <- commandArgs(trailingOnly = TRUE)
BASE <- if (length(args) >= 1) args[1] else "c3_N100_50kb_l1e4"
print(BASE)
CNA_DIR <- file.path("<MINTCNA_ROOT>/Simulation/cnasim", BASE)  
GTF    <- file.path("<MINTCNA_ROOT>/Simulation", "ref_hg38", "hg38.annotation.gtf") 

# ---------------------------
# SETTINGS
# ---------------------------
set.seed(1)

# ---------------------------
# 1) Read CN bin truth (profiles.tsv)
# ---------------------------
profiles <- fread(file.path(CNA_DIR, "profiles.tsv"))
colnames(profiles)[5] <- "cn_states"
setnames(profiles, tolower(names(profiles)))

# cn_total from "1.1" (two haplotypes split by ".")
profiles[, cn_total := {
  x <- trimws(as.character(cn_states))
  ifelse(grepl("\\.", x),
         rowSums(matrix(as.numeric(unlist(tstrsplit(x, "\\.", perl=TRUE))), ncol=2, byrow=FALSE)),
         round(as.numeric(x))) |> as.integer()
}]

# standardize chrom and make ids
profiles[!grepl("^chr", chrom), chrom := paste0("chr", chrom)]
profiles[, bin_id := paste0(chrom, ":", start, "-", end)]

# numeric keys for proper ordering (cell1,cell2,... and chr1,chr2,...)
profiles[, cell_num := as.integer(sub("^cell", "", cell))]
profiles[, chr_num  := as.integer(sub("^chr",  "", chrom))]

# keep only real cells if ancestors are present (optional; comment out if not needed)
profiles <- profiles[grepl("^cell\\d+$", cell)]

# order rows (clean)
setorder(profiles, chr_num, start, end, cell_num)

# correct cell order
cells <- profiles[, unique(.SD[order(cell_num)]$cell), .SDcols = c("cell","cell_num")]

# bin x cell CN matrix with correct column order
CNdt  <- dcast(profiles, bin_id ~ cell, value.var = "cn_total")
bin_ids <- CNdt$bin_id
CNmat <- as.matrix(CNdt[, ..cells])
rownames(CNmat) <- bin_ids

# enforce row order by chr/start
ord_bin <- order(
  as.integer(sub("^chr", "", sub(":.*$", "", rownames(CNmat)))),
  as.integer(sub("^.*:(\\d+)-.*$", "\\1", rownames(CNmat)))
)
CNmat <- CNmat[ord_bin, , drop=FALSE]

# ---------------------------
# 2) Simulate scRNA counts (genes x cells)
# ---------------------------

## inputs from DNA step
# CNmat: bins x cells integer CN; rownames like "chr1:0-100000"; colnames like "cell1"
# Dmat : bins x cells dosage (CN/2)^gamma with CN=2 -> 1

## 1) genes with real coordinates (hg38 autosomes)
gtf <- import("<MINTCNA_ROOT>/Simulation/ref_hg38/gencode.hg38.annotation.gtf")
genes_gr <- gtf[gtf$type == "gene"]
#mcols(genes_gr)$gene_id <- ifelse(is.null(mcols(genes_gr)$gene_id), mcols(genes_gr)$ID, mcols(genes_gr)$gene_id)
genes_gr <- genes_gr[seqnames(genes_gr) %in% paste0("chr", 1:22)]
genes_gr <- unique(genes_gr)

## 2) bins GRanges from CNmat rownames
bin_id <- rownames(CNmat)
bin_chr <- sub(":.*$", "", bin_id)
bin_sta <- as.integer(sub("^.*:(\\d+)-.*$", "\\1", bin_id))
bin_end <- as.integer(sub("^.*-(\\d+)$", "\\1", bin_id))
bins_gr <- GRanges(seqnames=bin_chr, ranges=IRanges(start=bin_sta+1, end=bin_end))  # +1 for 1-based

## 3) gene–bin overlap sparse matrix A (genes x bins)
hits <- findOverlaps(genes_gr, bins_gr, ignore.strand=TRUE)
A <- sparseMatrix(i=queryHits(hits), j=subjectHits(hits), x=1,
                  dims=c(length(genes_gr), length(bins_gr)))
den <- as.numeric(A %*% rep(1, length(bins_gr)))
den[den==0] <- 1  # safety

## 4) baseline scRNA counts with Splatter (genes x cells, gene names from hg38)
set.seed(1)
ng <- length(genes_gr); nc <- ncol(CNmat)
sim <- splatSimulate(nGenes=ng, batchCells=nc, method="groups", group.prob=rep(1,1), verbose=FALSE)
Y0  <- as.matrix(counts(sim))  # dgCMatrix, genes x cells
rownames(Y0) <- mcols(genes_gr)$gene_id
colnames(Y0) <- colnames(CNmat)

## 5) inject CNA into RNA with partial correspondence p (0-1)
#Dbin <- Dmat  # bins x cells, rownames must match bins_gr order
# ensure bin order consistent
p <- 1        # presence/absence
gamma <- 1    # strength when present

# dosage per bin per cell
Dbin <- (CNmat / 2) ^ gamma
Dbin[CNmat == 2L] <- 1.0

Gdos <- (A %*% Dbin) / den        # genes x cells dosage from overlapping bins
mask <- matrix(rbinom(length(Gdos), 1, p), nrow=nrow(Gdos))
Gdos_eff <- 1 + gamma*(Gdos - 1)            # weak effect
Gdos_final <- mask * Gdos_eff + (1-mask) * 1
Y <- round(Y0 * Gdos_final)
Y[Y < 0] <- 0

## save
fwrite(data.table(gene_id=rownames(Y), as.matrix(Y)), file.path(CNA_DIR, "rna_counts_gene_by_cell.tsv"), sep="\t")
fwrite(data.table(gene_id=mcols(genes_gr)$gene_id,
                  chrom=as.character(seqnames(genes_gr)),
                  start=start(genes_gr), end=end(genes_gr)),
       file.path(CNA_DIR, "genes_hg38.tsv"), sep="\t")


######################### Grid for RNA ########################## 
## grid
p_grid     <- c(1/4,2/4,3/4,4/4)
gamma_grid <- c(1)

## (optional) save genes annotation once
# fwrite(data.table(
#   gene_id = mcols(genes_gr)$gene_id,
#   chrom   = as.character(seqnames(genes_gr)),
#   start   = start(genes_gr),
#   end     = end(genes_gr)
# ), file.path(CNA_DIR, "rna","genes_hg38.tsv"), sep = "\t")

## loop over all combinations
for (g_idx in seq_along(gamma_grid)) {
  for (p_idx in seq_along(p_grid)) {
    
    gamma <- gamma_grid[g_idx]
    p     <- p_grid[p_idx]
    
    ## --- your simulation block ---
    Dbin <- (CNmat / 2) ^ gamma
    Dbin[CNmat == 2L] <- 1.0
    
    Gdos <- (A %*% Dbin) / den
    mask <- matrix(rbinom(length(Gdos), 1, p), nrow = nrow(Gdos))
    
    Gdos_eff   <- 1 + gamma * (Gdos - 1)         # weak effect
    Gdos_final <- mask * Gdos_eff + (1 - mask) * 1
    
    Y <- round(Y0 * Gdos_final)
    Y[Y < 0] <- 0
    ## -----------------------------
    
    ## save with gX_pX
    out_fn <- sprintf("rna_counts_gene_by_cell_g%d_p%d.tsv", g_idx, p_idx)
    fwrite(data.table(gene_id = rownames(Y), as.matrix(Y)),
           file.path(CNA_DIR, out_fn), sep = "\t")
  }
}



