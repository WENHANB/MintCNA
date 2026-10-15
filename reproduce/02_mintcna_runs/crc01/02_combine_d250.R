# CRC01 re-run 2026-10-14: copy of Code/CRC/RNA_mapp_DNAmethy.R with the d=250 DNA denoising and RAW RNA
# (the published run combined an older per-chromosome AE output with an AE-denoised RNA matrix).
library(tidyverse)
library(GenomicRanges)
library(data.table)
setwd("<MINTCNA_ROOT>/")

DNAmethy_ref <- read.csv("./Data/CRC/DNAmeth/processed_l2r/ref_all_cells_100kb.csv",row.names = 1)
RNA_ref <- read.csv("./Data/CRC/RNA/processed_l2r/ref_all_cells_CopyKAT.csv",row.names = 1)

# DNAmethy_l2r <- read.csv("./Data/CRC/DNAmeth/processed_l2r/l2r_smoothed_all_cells_100kb.csv",row.names = 1)
# RNA_l2r <- read.csv("./Data/CRC/RNA/processed_l2r/l2r_all_cells_CopyKAT.csv",row.names = 1)
# 
# DNAmethy_l2r <- read.csv("./Data/CRC/multi_modal/DNAmeth_denoised_l2r_multimodal.csv",row.names = 1)
# RNA_l2r <- read.csv("./Data/CRC/multi_modal/RNA_denoised_l2r_multimodal.csv",row.names = 1)

DNAmethy_l2r <- read.csv("./Data/CRC/DNAmeth/processed_l2r/ae_dna_d250_197.csv", row.names = 1)   # AG-DCAE d=250, 197 cells
RNA_l2r <- read.csv("./Data/CRC/RNA/processed_l2r/l2r_all_cells_CopyKAT.csv",row.names = 1)   # raw RNA log2 ratios, as in the simulation pipeline
stopifnot(ncol(RNA_l2r) == 197, ncol(DNAmethy_l2r) == 197)

## Sample matching
# SraRunTable <- read.csv("./Data/CRC/SraRunTable.csv")
# SraRunTable <- SraRunTable[,c("Run","GEO_Accession..exp.")]
# colnames(SraRunTable)[2] <- "Accession"
# sample <- read.csv("./Data/CRC/sample.csv")
# sample <- sample[which(sample$Title %in% intersect(paste0("scTrioSeq2Rna_", colnames(DNAmethy_l2r)), sample$Title)), ]
# Sra_Sample_ref <- merge(SraRunTable, sample[,c("Accession","Title")], by = "Accession")
# Sra_Sample_ref$Title <- sub("^scTrioSeq2Rna_", "", Sra_Sample_ref$Title)
# name_map <- setNames(Sra_Sample_ref$Title, Sra_Sample_ref$Run)
# colnames(RNA_l2r) <- name_map[colnames(RNA_l2r)]

## Genomic Coordinate mapping (RNA to DNA methylation)
# Step 1: Rename columns to avoid duplication
RNA_ref <- dplyr::rename(RNA_ref,
                         chr_rna = chr,
                         start_rna = start,
                         end_rna = end)

DNAmethy_ref <- dplyr::rename(
  DNAmethy_ref,
  chr_dna   = chr,
  start_dna = start,
  end_dna   = end
)

# Step 2: Convert to GRanges
rna_gr <- GRanges(seqnames = RNA_ref$chr_rna,
                  ranges = IRanges(start = RNA_ref$start_rna, end = RNA_ref$end_rna))

dna_gr <- GRanges(seqnames = DNAmethy_ref$chr_dna,
                  ranges = IRanges(start = DNAmethy_ref$start_dna, end = DNAmethy_ref$end_dna))

# Step 3: Find overlaps
hits <- findOverlaps(rna_gr, dna_gr)

# Step 4: Create mapping table
rna_mapped <- cbind(RNA_ref[queryHits(hits), ], DNAmethy_ref[subjectHits(hits), ])

# Step 5: Summarize: how many RNA regions fall into each DNA methylation bin
bin_counts <- rna_mapped %>%
  group_by(chr_dna, start_dna, end_dna) %>%
  summarise(n_genes = n(), .groups = "drop")

# Join the counts back to the rna_mapped dataframe
rna_mapped <- left_join(rna_mapped, bin_counts,
                        by = c("chr_dna", "start_dna", "end_dna"))

library(GenomicRanges)
library(tidyverse)

# Step 1: Reshape RNA matrix (genes by sample) to long format
rna_long <- RNA_l2r %>%
  rownames_to_column("gene_id") %>%
  pivot_longer(-gene_id, names_to = "sample", values_to = "expr")

# Step 2: Merge with RNA_ref to get genomic info per gene
rna_annotated <- rna_long %>%
  left_join(RNA_ref, by = "gene_id")  # assumes gene_id column exists in RNA_ref

# Step 3: Convert to GRanges
rna_gr <- GRanges(seqnames = rna_annotated$chr_rna,
                  ranges = IRanges(start = rna_annotated$start_rna, end = rna_annotated$end_rna),
                  gene_id = rna_annotated$gene_id,
                  sample = rna_annotated$sample,
                  expr = rna_annotated$expr)

dna_gr <- GRanges(seqnames = DNAmethy_ref$chr,
                  ranges = IRanges(start = DNAmethy_ref$start, end = DNAmethy_ref$end))

# Step 4: Find overlaps
hits <- findOverlaps(rna_gr, dna_gr)

# Step 5: Map expression to DNA bins
mapped <- data.frame(
  chr = as.character(seqnames(dna_gr[subjectHits(hits)])),
  start = start(dna_gr[subjectHits(hits)]),
  end = end(dna_gr[subjectHits(hits)]),
  sample = mcols(rna_gr[queryHits(hits)])$sample,
  expr = mcols(rna_gr[queryHits(hits)])$expr
)

# Step 6: Aggregate expression per bin and sample (mean or median)
expr_by_bin_sample <- mapped %>%
  group_by(chr, start, end, sample) %>%
  summarise(mean_expr = median(expr, na.rm = TRUE), .groups = "drop")  # use median(expr, na.rm = TRUE) if preferred

# Create a unique bin ID for rows
expr_by_bin_sample <- expr_by_bin_sample %>%
  mutate(bin_id = paste(chr, start, end, sep = "_"))

# Convert to wide format
bin_by_sample_matrix <- expr_by_bin_sample %>%
  select(bin_id, sample, mean_expr) %>%
  pivot_wider(names_from = sample, values_from = mean_expr)

bin_ref <- expr_by_bin_sample %>%
  select(bin_id, chr, start, end) %>%
  distinct()

ord <- with(bin_ref, order(as.numeric(sub("chr", "", chr)), start))
bin_by_sample_matrix <- bin_by_sample_matrix[ord,]
bin_ref <- bin_ref[ord,]

## Make list
# Ensure bin_by_sample_matrix has bin_id as rownames
options(scipen = 999)
bin_ids <- bin_by_sample_matrix$bin_id
bin_by_sample_matrix <- bin_by_sample_matrix[, -1]  # remove bin_id column
rownames(bin_by_sample_matrix) <- bin_ids

DNAmethy_ref <- DNAmethy_ref %>%
  mutate(bin_id = paste(chr_dna, as.character(start_dna), as.character(end_dna), sep = "_"))
bin_ids <- DNAmethy_ref$bin_id
rownames(DNAmethy_l2r) <- bin_ids

## Combine
combined_multiomic_list <- list()
shared_samples <- intersect(names(RNA_l2r), names(DNAmethy_l2r))

ord <- with(DNAmethy_ref, order(as.numeric(sub("chr", "", chr_dna)), start_dna))
DNAmethy_l2r <- DNAmethy_l2r[ord,]
DNAmethy_ref <- DNAmethy_ref[ord,]

DNAmethy_l2r$bin_id <- rownames(DNAmethy_l2r)
bin_by_sample_matrix$bin_id <- rownames(bin_by_sample_matrix)

for (s in shared_samples) {
  rna_mat <- bin_by_sample_matrix[,c("bin_id",s)]
  colnames(rna_mat)[2] <- paste(s,"_RNA")
  dna_mat <- DNAmethy_l2r[,c("bin_id",s)]
  colnames(dna_mat)[2] <- paste(s,"_DNA")
  # Full join by bin_id to align rows
  merged <- merge(rna_mat,dna_mat,by = "bin_id", all = T) %>%
    as.matrix()
  combined_multiomic_list[[s]] <- merged[,]
}

parts <- do.call(rbind, strsplit(merged[,1], "_"))
ref <- data.frame(chr = parts[,1],
                  start = as.integer(parts[,2]),
                  end = as.integer(parts[,3]))
ord <- with(ref, order(as.numeric(sub("chr", "", chr)), start))
combined_multiomic_list <- lapply(combined_multiomic_list, function(data) data[ord,])
ref <- ref[ord,]

#save(ref,combined_multiomic_list, file = "Data/CRC/combined_multiomic_100kb_allcells.Rdata")
stopifnot(length(combined_multiomic_list) == 197)
save(ref,combined_multiomic_list, file = "Data/CRC/combined_multiomic_100kb_allcells_AE_d250.Rdata")
