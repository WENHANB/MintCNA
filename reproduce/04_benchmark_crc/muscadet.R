library(muscadet)

# Load example dataset inputs:
# Matrices of raw counts per features
data("mat_counts_atac_tumor", "mat_counts_rna_tumor")
# Table of raw counts per allele
data("allele_counts_atac_tumor", "allele_counts_rna_tumor")
# Table of feature coordinates
data("peaks", "genes")

# Create individual omic objects
atac <- CreateMuscomicObject(
  type = "ATAC",
  mat_counts = mat_counts_atac_tumor, 
  features = peaks)
rna <- CreateMuscomicObject(
  type = "RNA",
  mat_counts = mat_counts_rna_tumor,
  features = genes)
