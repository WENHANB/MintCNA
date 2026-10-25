#!/usr/bin/env Rscript
# CONGAS+ (Rcongas 0.2.0 / congas 0.0.77) on one multiome sample, the package's documented sequence:
# init() -> filter_outliers() -> filter_missing_data() -> auto_config_run(K = 1:4) -> fit_congas(learning_rate 0.05,
# lambdas 0.5). ATAC counts are summed in 1 Mb bins; the segmentation is the package's chromosome arms
# (Rcongas::hg38_arms) without the five acrocentric p arms, which carry almost no reads and would remove
# nearly every nucleus in filter_missing_data().
#   Rscript reproduce/analysis/multiome/methods/run_congas.R <care_id>
# Input : ATAC_RNA/Data/atac_filtered/<care_id>_fragments.tsv.gz, ATAC_RNA/Data/gex_matrices/,
#         ATAC_RNA/Data/cell_lists/<care_id>_cells.tsv
# Output: ATAC_RNA/Result/benchmark/CONGASplus/<care_id>/{congas_input.rds, congas_result.rds}
# Run from the root of the MintCNA repository with MINTCNA_DATA set (see reproduce/README.md).
suppressPackageStartupMessages({ library(data.table); library(Matrix); library(dplyr); library(tibble); library(Rcongas); library(EnsDb.Hsapiens.v86) })
ROOT <- Sys.getenv("MINTCNA_DATA"); if (!nzchar(ROOT)) stop("set MINTCNA_DATA (see reproduce/README.md)")
args <- commandArgs(trailingOnly = TRUE); if (!length(args)) stop("need <care_id>")
CARE <- args[1]
setwd(file.path(ROOT, "ATAC_RNA"))
OUT <- file.path("Result/benchmark/CONGASplus", CARE); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
ATAC_BIN <- 1e6
CHRS <- paste0("chr", 1:22)
keep_bc <- fread(file.path("Data/cell_lists", paste0(CARE, "_cells.tsv")))$barcode

# ---- ATAC: fragments per 1 Mb bin and nucleus ---------------------------------------------------
fr <- fread(cmd = paste("zcat", shQuote(file.path("Data/atac_filtered", paste0(CARE, "_fragments.tsv.gz")))),
            header = FALSE, select = 1:4, col.names = c("chr", "start", "end", "cell"))
fr <- fr[chr %in% CHRS & cell %in% keep_bc]
fr[, from := (((start + end) %/% 2L) %/% as.integer(ATAC_BIN)) * as.integer(ATAC_BIN)]
atac <- fr[, .(value = .N), by = .(cell, chr, from)][, to := from + as.integer(ATAC_BIN)]
atac <- as_tibble(atac[, .(cell, value, chr, from, to)])

# ---- RNA: counts per gene and nucleus, gene coordinates from EnsDb ------------------------------
man  <- fread("Data/meta/sample_manifest.csv")[care_id == CARE]
gfil <- strsplit(man$gex_files, ";")[[1]]
mtx  <- readMM(file.path("Data/gex_matrices", grep("matrix",   gfil, value = TRUE)))
bc   <- fread(file.path("Data/gex_matrices", grep("barcodes", gfil, value = TRUE)), header = FALSE)$V1
ft   <- fread(file.path("Data/gex_matrices", grep("features", gfil, value = TRUE)), header = FALSE)
setnames(ft, c("gene_id", "gene_name", "type"))
stopifnot(ncol(mtx) == length(bc), nrow(mtx) == nrow(ft))
gg <- as.data.table(ensembldb::genes(EnsDb.Hsapiens.v86, columns = c("gene_id", "seq_name", "gene_seq_start", "gene_seq_end")))
gg <- gg[, .(gene_id, chr = paste0("chr", seqnames), from = start, to = end)][chr %in% CHRS]
ft[, row := .I]; ft <- merge(ft, gg, by = "gene_id")
j  <- which(bc %in% keep_bc)
tt <- as(mtx[ft$row, j, drop = FALSE], "TsparseMatrix")
rna <- data.table(gene = ft$gene_name[tt@i + 1L], chr = ft$chr[tt@i + 1L], from = ft$from[tt@i + 1L],
                  to = ft$to[tt@i + 1L], cell = bc[j][tt@j + 1L], value = as.integer(tt@x))
rna <- as_tibble(rna[value > 0])

common <- intersect(unique(rna$cell), unique(atac$cell))               # multiome mode: same barcodes on both sides
rna  <- dplyr::filter(rna,  cell %in% common)
atac <- dplyr::filter(atac, cell %in% common)
message(sprintf("%s: %d nuclei in both modalities", CARE, length(common)))

# ---- segmentation: chromosome arms minus the acrocentric p arms ------------------------------------
data("hg38_arms", package = "Rcongas")
segmentation <- dplyr::filter(hg38_arms, !(chr %in% c("chr13", "chr14", "chr15", "chr21", "chr22") & from == 0))

# ---- fit ------------------------------------------------------------------------------------------
t0 <- Sys.time()
x <- Rcongas::init(rna = rna, atac = atac, segmentation = segmentation, rna_likelihood = "NB", atac_likelihood = "NB",
                   reference_genome = "GRCh38", multiome = TRUE, description = paste0("CARE ", CARE))
x <- Rcongas::filter_outliers(x)
x <- Rcongas::filter_missing_data(x)
print(Rcongas::stat(x)[c("nsegments", "ncells_RNA", "ncells_ATAC")])
K <- 1:4
hyper <- Rcongas::auto_config_run(x, K)
saveRDS(list(care_id = CARE, object = x, hyper = hyper, n_cells = length(common), segmentation = segmentation),
        file.path(OUT, "congas_input.rds"))
fit <- Rcongas::fit_congas(x, K = K, lambdas = 0.5, learning_rate = 0.05, model_parameters = hyper)
saveRDS(list(care_id = CARE, object = x, hyper = hyper, fit = fit, n_cells = length(common), atac_bin = ATAC_BIN),
        file.path(OUT, "congas_result.rds"))
cat(sprintf("SUMMARY\t%s\t%d cells\t%.1f min\n", CARE, length(common), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
