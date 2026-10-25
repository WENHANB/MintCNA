#!/usr/bin/env Rscript
# Per-sample nucleus lists from the CARE authors' snATAC annotation: the QC-passing nuclei analysed in
# their paper, their cell type, and the reference (diploid) flag. Reference nuclei are the
# non-neuroectodermal lineages (Myeloid, Endo/Mural, Lymphocyte) whose RNA and ATAC labels agree;
# oligodendrocytes and astrocytes are not used as reference because they share lineage with the tumour.
#   Rscript reproduce/analysis/multiome/export_cell_lists.R
# Input : ATAC_RNA/Data/synapse/snatac_celltype_meta_data.RDS (download_synapse.py)
# Output: ATAC_RNA/Data/cell_lists/<care_id>_cells.tsv (barcode, cell_type, nFrags, TSSEnrichment, is_reference),
#         <care_id>_barcodes.txt, cell_list_summary.tsv
# Run from the root of the MintCNA repository with MINTCNA_DATA set (see reproduce/README.md).
suppressPackageStartupMessages(library(data.table))
ROOT <- Sys.getenv("MINTCNA_DATA"); if (!nzchar(ROOT)) stop("set MINTCNA_DATA (see reproduce/README.md)")
ROOT <- file.path(ROOT, "ATAC_RNA")
OUT  <- file.path(ROOT, "Data/cell_lists"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
SAMPLES <- c("P59T2","P59T3","P104T1","P104T2","P106T1","P106T2","P100T1","P100T3","P103T1","P103T2",
             "P105T1","P105T3","P108T1","P108T2","P109T1","P109T2","P110T1","P110T2","P111T2","P111T3")
REF_TYPES <- c("Myeloid", "Endo/Mural", "Lymphocyte")

meta <- readRDS(file.path(ROOT, "Data/synapse/snatac_celltype_meta_data.RDS")); setDT(meta)
meta[, barcode := sub(".*#", "", CellID)]
meta[, concordant := rna_cell_type == atac_cell_type_cluster |
       (rna_cell_type %in% c("ExcNeuron", "InhNeuron") & atac_cell_type_cluster == "Neuron") |
       (rna_cell_type %in% c("Endothelial", "Mural") & atac_cell_type_cluster == "Endo/Mural")]

summary_rows <- list()
for (s in SAMPLES) {
  d <- meta[care_id == s]
  if (!nrow(d)) { message("no metadata for ", s); next }
  d[, is_reference := atac_cell_type_cluster %in% REF_TYPES & concordant]
  fwrite(d[, .(barcode, cell_type = atac_cell_type_cluster, nFrags, TSSEnrichment, is_reference)],
         file.path(OUT, paste0(s, "_cells.tsv")), sep = "\t")
  writeLines(d$barcode, file.path(OUT, paste0(s, "_barcodes.txt")))
  summary_rows[[s]] <- data.table(care_id = s, n_cells = nrow(d), n_malig = sum(d$atac_cell_type_cluster == "Malignant"),
                                  n_ref = sum(d$is_reference), med_nFrags = as.integer(median(d$nFrags)))
}
sm <- rbindlist(summary_rows)
fwrite(sm, file.path(OUT, "cell_list_summary.tsv"), sep = "\t")
print(sm)
cat("samples with fewer than 100 reference nuclei:", paste(sm$care_id[sm$n_ref < 100], collapse = ", "), "\n")
