#!/usr/bin/env Rscript
# Export a per-sample cell whitelist from the authors' snATAC metadata.
#
# Why the RDS and not the GEX barcodes: the GEX barcodes.tsv.gz is the GEX
# library's called-cell list, which for P103T1 holds 5,264 barcodes. The authors'
# ArchR QC (doublet removal, TSS/nFrags thresholds) reduces that to 2,671, and
# those 2,671 are the nuclei the paper actually analysed -- the same set that
# sums to the published 118,180. Filtering fragments on the GEX list would leave
# doublets in, and their fragment counts run to 888k, which would dominate any
# per-bin normalisation.
#
# Emits Data/cell_lists/<care_id>_cells.tsv with columns:
#   barcode  cell_type  nFrags  TSSEnrichment  is_reference
# plus a bare barcode list for the awk filter.

suppressPackageStartupMessages(library(data.table))

ROOT <- "<MINTCNA_ROOT>/ATAC_RNA"
OUT  <- file.path(ROOT, "Data/cell_lists")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

SAMPLES <- c("P59T2","P59T3","P104T1","P104T2","P106T1","P106T2",
             "P100T1","P100T3","P103T1","P103T2","P105T1","P105T3",
             "P108T1","P108T2","P109T1","P109T2","P110T1","P110T2","P111T2","P111T3")

# Non-neuroectodermal lineages only. Oligodendrocyte and Astrocyte are excluded
# on purpose: glioma is glial, so a glial reference risks sharing lineage-specific
# open chromatin -- and for oligodendroglioma it risks carrying the very 1p/19q
# codeletion we are trying to detect.
REF_TYPES <- c("Myeloid", "Endo/Mural", "Lymphocyte")

meta <- readRDS(file.path(ROOT, "Data/synapse/snatac_celltype_meta_data.RDS"))
setDT(meta)
meta[, barcode := sub(".*#", "", CellID)]

# Nuclei whose RNA-side and ATAC-side labels disagree are dropped from the
# reference set (not from the analysis): 340 cohort-wide sit between Malignant
# and Astrocyte, and a single mis-annotated malignant nucleus shifts a mean.
meta[, concordant := rna_cell_type == atac_cell_type_cluster |
       (rna_cell_type %in% c("ExcNeuron","InhNeuron") & atac_cell_type_cluster == "Neuron") |
       (rna_cell_type %in% c("Endothelial","Mural")   & atac_cell_type_cluster == "Endo/Mural")]

summary_rows <- list()
for (s in SAMPLES) {
  d <- meta[care_id == s]
  if (!nrow(d)) { message("no metadata for ", s); next }
  d[, is_reference := atac_cell_type_cluster %in% REF_TYPES & concordant]

  fwrite(d[, .(barcode, cell_type = atac_cell_type_cluster, nFrags,
               TSSEnrichment, is_reference)],
         file.path(OUT, paste0(s, "_cells.tsv")), sep = "\t")
  writeLines(d$barcode, file.path(OUT, paste0(s, "_barcodes.txt")))

  summary_rows[[s]] <- data.table(
    care_id   = s,
    n_cells   = nrow(d),
    n_malig   = sum(d$atac_cell_type_cluster == "Malignant"),
    n_ref     = sum(d$is_reference),
    n_myeloid = sum(d$atac_cell_type_cluster == "Myeloid" & d$concordant),
    n_diploid_other = sum(!d$atac_cell_type_cluster %in% c("Malignant", REF_TYPES)),
    med_nFrags = as.integer(median(d$nFrags)))
}

sm <- rbindlist(summary_rows)
fwrite(sm, file.path(OUT, "cell_list_summary.tsv"), sep = "\t")
print(sm)
cat("\ntotals: cells", sum(sm$n_cells), " malignant", sum(sm$n_malig),
    " reference", sum(sm$n_ref), " other-diploid", sum(sm$n_diploid_other), "\n")
cat("samples with <100 reference nuclei: ",
    paste(sm$care_id[sm$n_ref < 100], collapse = ", "), "\n")
