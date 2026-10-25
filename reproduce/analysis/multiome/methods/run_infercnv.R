#!/usr/bin/env Rscript
# inferCNV (1.26.0) on the RNA arm of one multiome sample, package defaults (cutoff 1, denoise FALSE, HMM FALSE),
# reference group = the reference nuclei of the cell list, gene order from EnsDb.Hsapiens.v86.
#   Rscript reproduce/analysis/multiome/methods/run_infercnv.R <care_id>
#   INFERCNV_TAG (default "default"), INFERCNV_CUTOFF (1) and INFERCNV_DENOISE (FALSE) name and set the configuration.
# Input : ATAC_RNA/Data/gex_matrices/ (download_geo.sh), ATAC_RNA/Data/cell_lists/<care_id>_cells.tsv
# Output: ATAC_RNA/Result/benchmark/inferCNV/<tag>/<care_id>/  (inferCNV output; export_infercnv.R extracts the matrix)
# Run from the root of the MintCNA repository with MINTCNA_DATA set (see reproduce/README.md).
suppressPackageStartupMessages({ library(data.table); library(Matrix); library(infercnv); library(EnsDb.Hsapiens.v86) })
ROOT <- Sys.getenv("MINTCNA_DATA"); if (!nzchar(ROOT)) stop("set MINTCNA_DATA (see reproduce/README.md)")
args <- commandArgs(trailingOnly = TRUE); if (!length(args)) stop("need <care_id>")
CARE <- args[1]
setwd(file.path(ROOT, "ATAC_RNA"))
TAG <- Sys.getenv("INFERCNV_TAG", "default")
CUTOFF  <- as.numeric(Sys.getenv("INFERCNV_CUTOFF", "1"))
DENOISE <- as.logical(Sys.getenv("INFERCNV_DENOISE", "FALSE"))
OUT <- file.path("Result/benchmark/inferCNV", TAG, CARE); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
CHRS <- paste0("chr", 1:22)

# ---- counts of the whitelisted nuclei, two groups: reference and observed -------------------------
man  <- fread("Data/meta/sample_manifest.csv")[care_id == CARE]
gfil <- strsplit(man$gex_files, ";")[[1]]
mtx  <- readMM(file.path("Data/gex_matrices", grep("matrix",   gfil, value = TRUE)))
bc   <- fread(file.path("Data/gex_matrices", grep("barcodes", gfil, value = TRUE)), header = FALSE)$V1
ft   <- fread(file.path("Data/gex_matrices", grep("features", gfil, value = TRUE)), header = FALSE)
setnames(ft, c("gene_id", "gene_name", "type"))
stopifnot(ncol(mtx) == length(bc), nrow(mtx) == nrow(ft))
cells <- fread(file.path("Data/cell_lists", paste0(CARE, "_cells.tsv")))
keep  <- which(bc %in% cells$barcode)
mtx <- mtx[, keep, drop = FALSE]; bc <- bc[keep]
grp <- ifelse(cells[match(bc, barcode)]$is_reference, "reference", "observed")
message(sprintf("%s: %d cells (%d reference)", CARE, length(bc), sum(grp == "reference")))
if (sum(grp == "reference") < 20) stop("too few reference nuclei for inferCNV")

# ---- gene order (one row per gene symbol, the longest) --------------------------------------------
gg <- as.data.table(ensembldb::genes(EnsDb.Hsapiens.v86, columns = c("gene_id", "seq_name", "gene_seq_start", "gene_seq_end")))
gg <- gg[, .(gene_id, chr = paste0("chr", seqnames), start, end)][chr %in% CHRS]
ft[, row := .I]
ft <- merge(ft, gg, by = "gene_id")
ft[, len := end - start]; setorder(ft, gene_name, -len); ft <- ft[!duplicated(gene_name)]; setorder(ft, chr, start)
counts <- as.matrix(mtx[ft$row, , drop = FALSE]); rownames(counts) <- ft$gene_name; colnames(counts) <- bc
gene_file <- file.path(OUT, "gene_order.txt"); ann_file <- file.path(OUT, "cell_annot.txt")
fwrite(ft[, .(gene_name, chr, start, end)], gene_file, sep = "\t", col.names = FALSE)
fwrite(data.table(cell = bc, group = grp), ann_file, sep = "\t", col.names = FALSE)

# ---- run ----------------------------------------------------------------------------------------
t0 <- Sys.time()
obj <- CreateInfercnvObject(raw_counts_matrix = counts, annotations_file = ann_file, gene_order_file = gene_file,
                            ref_group_names = "reference", delim = "\t")
obj <- infercnv::run(obj, cutoff = CUTOFF, out_dir = OUT, cluster_by_groups = TRUE, denoise = DENOISE, HMM = FALSE,
                     num_threads = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4")), no_plot = FALSE, save_rds = TRUE)
saveRDS(list(care_id = CARE, n_cells = length(bc), n_ref = sum(grp == "reference"), n_genes = nrow(ft),
             genes = ft[, .(gene_name, chr, start, end)]), file.path(OUT, "infercnv_meta.rds"))
cat(sprintf("SUMMARY\t%s\t%d cells\t%d reference\t%d genes\t%.1f min\n", CARE, length(bc), sum(grp == "reference"), nrow(ft),
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
