#!/usr/bin/env Rscript
# inferCNV on one sample's RNA arm.  Usage: Rscript run_infercnv.R <care_id>
#
# The only RNA-side caller in the benchmark, and the one the CARE authors ran
# themselves -- so there is a published result on this exact cohort to check
# against. It is RNA-only: it belongs in the same tier as the ATAC-only callers
# (epiAneufinder, AtaCNV, Copy-scAT, RIDDLER), not with CONGAS+.
#
# Reference cells: our own annotation (Myeloid + Endo/Mural + Lymphocyte with
# concordant RNA and ATAC labels). This DIVERGES from the CARE paper, which used
# myeloid AND oligodendrocytes. Oligodendrocytes are excluded here because glioma
# is a glial tumour and, for oligodendroglioma, an oligodendrocyte reference risks
# carrying the 1p/19q codeletion that is the detection target. State this in Methods.
#
# Parameters are inferCNV's documented settings for 10x data:
#   cutoff = 0.1        (the package's own recommendation for 10x; 1 is for Smart-seq)
#   cluster_by_groups   TRUE
#   denoise             TRUE
#   HMM                 FALSE
# Everything else is left at the package default.
#
# HMM = FALSE matches the settings already used for this project's simulation
# benchmark (Simulation/cnasim/.../RNA_benchmark/InferCNV/run_infercnv.R), which
# ran inferCNV "non-HMM, continuous residual expression". Keeping the two
# sections on the same setting matters twice over: the results stay comparable
# across the paper, and the continuous residual is what Pearson and RMSE need --
# the discrete HMM states would put inferCNV in the same awkward position as
# epiAneufinder, where a magnitude metric requires an invented amplitude.

args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("need <care_id>")
CARE <- args[1]

suppressPackageStartupMessages({
  library(data.table); library(Matrix); library(infercnv)
})

ROOT <- "<MINTCNA_ROOT>/ATAC_RNA"
setwd(ROOT)
TAG <- Sys.getenv("INFERCNV_TAG", "default")
OUT <- file.path("Result/benchmark/inferCNV", TAG, CARE)
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
CHRS <- paste0("chr", 1:22)

# ---- counts -----------------------------------------------------------------
man  <- fread("Data/meta/sample_manifest.csv")[care_id == CARE]
gfil <- strsplit(man$gex_files, ";")[[1]]
mtx  <- readMM(file.path("Data/gex_matrices", grep("matrix",   gfil, value = TRUE)))
bc   <- fread(file.path("Data/gex_matrices", grep("barcodes", gfil, value = TRUE)), header = FALSE)$V1
ft   <- fread(file.path("Data/gex_matrices", grep("features", gfil, value = TRUE)), header = FALSE)
setnames(ft, c("gene_id", "gene_name", "type"))
stopifnot(ncol(mtx) == length(bc), nrow(mtx) == nrow(ft))

# ---- cells: our whitelist, with the reference flag ---------------------------
cells <- fread(file.path("Data/cell_lists", paste0(CARE, "_cells.tsv")))
keep  <- which(bc %in% cells$barcode)
mtx   <- mtx[, keep, drop = FALSE]
bc    <- bc[keep]

ann <- cells[match(bc, barcode)]
# Two groups only: reference nuclei versus everything else. inferCNV needs the
# reference group named explicitly; finer cell-type groups would split the
# reference and shrink each one below a usable size.
grp <- ifelse(ann$is_reference, "reference", "observed")
message(sprintf("%s: %d cells (%d reference, %d observed)",
                CARE, length(bc), sum(grp == "reference"), sum(grp == "observed")))
if (sum(grp == "reference") < 20) stop("too few reference nuclei for inferCNV")

# ---- gene order file --------------------------------------------------------
# inferCNV orders genes along the genome, so every gene needs coordinates. Same
# EnsDb source as the CONGAS+ runner, for consistency between the two.
suppressPackageStartupMessages(library(EnsDb.Hsapiens.v86))
gg <- as.data.table(ensembldb::genes(EnsDb.Hsapiens.v86,
                                     columns = c("gene_id", "seq_name", "gene_seq_start", "gene_seq_end")))
gg <- gg[, .(gene_id, chr = paste0("chr", seqnames), start, end)][chr %in% CHRS]
ft[, row := .I]
ft <- merge(ft, gg, by = "gene_id")
# Duplicate gene symbols break inferCNV's rowname-based join; keep the longest.
ft[, len := end - start]
setorder(ft, gene_name, -len)
ft <- ft[!duplicated(gene_name)]
setorder(ft, chr, start)
message(sprintf("genes with coordinates: %d", nrow(ft)))

counts <- as.matrix(mtx[ft$row, , drop = FALSE])
rownames(counts) <- ft$gene_name
colnames(counts) <- bc

gene_file <- file.path(OUT, "gene_order.txt")
ann_file  <- file.path(OUT, "cell_annot.txt")
fwrite(ft[, .(gene_name, chr, start, end)], gene_file, sep = "\t", col.names = FALSE)
fwrite(data.table(cell = bc, group = grp),  ann_file,  sep = "\t", col.names = FALSE)

# ---- run --------------------------------------------------------------------
t0 <- Sys.time()
obj <- CreateInfercnvObject(raw_counts_matrix   = counts,
                            annotations_file    = ann_file,
                            gene_order_file     = gene_file,
                            ref_group_names     = "reference",
                            delim               = "\t")

# PARAMETERS ARE THE PACKAGE DEFAULTS unless a flag says otherwise. Departures
# found in an audit on 2026-08-18 and reverted:
#   denoise  was TRUE, formals default is FALSE  -> now FALSE
#   cutoff   formals default is 1, but inferCNV's own documentation says
#            "cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for
#            10x Genomics". Both are run; CUTOFF selects which.
CUTOFF  <- as.numeric(Sys.getenv("INFERCNV_CUTOFF",  "1"))
DENOISE <- as.logical(Sys.getenv("INFERCNV_DENOISE", "FALSE"))
message(sprintf("[inferCNV] cutoff = %s   denoise = %s", CUTOFF, DENOISE))

obj <- infercnv::run(obj,
                     cutoff            = CUTOFF,
                     out_dir           = OUT,
                     cluster_by_groups = TRUE,
                     denoise           = DENOISE,
                     HMM               = FALSE,
                     num_threads       = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4")),
                     no_plot           = FALSE,
                     save_rds          = TRUE)
el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))

saveRDS(list(care_id = CARE, n_cells = length(bc), n_ref = sum(grp == "reference"),
             n_genes = nrow(ft), genes = ft[, .(gene_name, chr, start, end)]),
        file.path(OUT, "infercnv_meta.rds"))
cat(sprintf("SUMMARY\t%s\t%d\t%d\t%d\t%.1f\n",
            CARE, length(bc), sum(grp == "reference"), nrow(ft), el))
cat("  (care_id, n_cells, n_reference, n_genes, minutes)\n")
