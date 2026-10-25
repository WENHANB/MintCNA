#!/usr/bin/env Rscript
# MintCNA - Build Demo Dataset
#
# Extracts a small chr1-22 subset from the c3_N100_100kb_l1e4_c0.2_n0.1
# CNAsim simulation directory, suitable for end-to-end pipeline testing.
#
# Usage:
#   Rscript build_demo.R [<source_dir>] [<output_dir>]
#
# Defaults:
#   source_dir: $MINTCNA_DATA/Simulation/cnasim/c3_N100_100kb_l1e4_c0.2_n0.1
#               (the Zenodo setting unpacked as described in reproduce/README.md)
#   output_dir: ../demo  (relative to this script)

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(data.table)
})

# -------------------- Args --------------------
args <- commandArgs(trailingOnly = TRUE)
SRC <- if (length(args) >= 1) args[1] else {
  if (!nzchar(Sys.getenv("MINTCNA_DATA"))) stop("give <source_dir> or set MINTCNA_DATA (see reproduce/README.md)")
  file.path(Sys.getenv("MINTCNA_DATA"), "Simulation/cnasim/c3_N100_100kb_l1e4_c0.2_n0.1")
}

# Auto-detect output dir from this script's location
script_args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", grep("--file=", script_args, value = TRUE))
if (length(script_path)) {
  OUT <- if (length(args) >= 2) args[2] else file.path(dirname(normalizePath(script_path)), "input")
} else {
  OUT <- if (length(args) >= 2) args[2] else "demo/input"
}

dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

message("[demo] Source: ", SRC)
message("[demo] Output: ", OUT)

KEEP_CHRS <- paste0("chr", 1:22)

# -------------------- 1. DNA: extract chr1-22 from DNA_RC.RData --------------------
message("\n[demo] Extracting DNA read counts (chr1-22)...")
e <- new.env()
load(file.path(SRC, "DNA_RC.RData"), envir = e)
Y_full  <- e$Y
gr_full <- e$gr

# Normalize chr names — gr may have "1","2",... or "chr1","chr2",...
gr_chrs <- as.character(seqnames(gr_full))
if (!any(grepl("^chr", gr_chrs))) gr_chrs <- paste0("chr", gr_chrs)

keep_idx <- which(gr_chrs %in% KEEP_CHRS)
message(sprintf("[demo] Keeping %d / %d DNA bins.", length(keep_idx), length(gr_full)))

Y_demo  <- Y_full[keep_idx, , drop = FALSE]
gr_demo <- gr_full[keep_idx]

# Force chr names to "chr<N>" form for downstream R code consistency
seqlevels(gr_demo) <- sub("^chr", "", seqlevels(gr_demo))   # strip if present
gr_demo <- renameSeqlevels(gr_demo, setNames(paste0("chr", seqlevels(gr_demo)), seqlevels(gr_demo)))

# Save same structure as original DNA_RC.RData (Y + gr)
Y  <- Y_demo
gr <- gr_demo
save(Y, gr, file = file.path(OUT, "DNA_RC.RData"))
message("[demo] Wrote ", file.path(OUT, "DNA_RC.RData"),
        " (", nrow(Y), " bins x ", ncol(Y), " cells)")

# -------------------- 2. RNA: extract chr1-22 genes --------------------
message("\n[demo] Extracting RNA gene-level counts (chr1-22)...")
gene_ref <- fread(file.path(SRC, "genes_hg38.tsv"))
# Expected columns: gene_id, chrom, start, end
chrom_col <- if ("chrom" %in% names(gene_ref)) "chrom" else "chr"
gene_ref_sub <- gene_ref[get(chrom_col) %in% KEEP_CHRS]

# Read RNA counts, subset to those genes
rna_counts <- fread(file.path(SRC, "rna_counts_gene_by_cell.tsv"))
gid_col <- names(rna_counts)[1]
rna_demo <- rna_counts[get(gid_col) %in% gene_ref_sub$gene_id]
gene_ref_demo <- gene_ref_sub[gene_id %in% rna_demo[[gid_col]]]

# Normalize the gene_ref column name to "chr" for downstream compatibility
if (chrom_col == "chrom") setnames(gene_ref_demo, "chrom", "chr")

message(sprintf("[demo] RNA: %d genes x %d cells.",
                nrow(rna_demo), ncol(rna_demo) - 1))

fwrite(rna_demo, file.path(OUT, "rna_counts.tsv"), sep = "\t")
fwrite(gene_ref_demo, file.path(OUT, "genes.tsv"), sep = "\t")

# -------------------- 3. cell_types.tsv (copy as-is) --------------------
message("\n[demo] Copying cell_types.tsv...")
file.copy(file.path(SRC, "cell_types.tsv"),
          file.path(OUT, "cell_types.tsv"), overwrite = TRUE)

# -------------------- 4. Summary --------------------
cell_types <- fread(file.path(OUT, "cell_types.tsv"))
n_normal <- sum(cell_types$normal == "normal")
n_tumor  <- sum(cell_types$normal != "normal")
message(sprintf("\n[demo] Done. Cells: %d normal + %d tumor = %d total.",
                n_normal, n_tumor, nrow(cell_types)))
message("[demo] Demo files in: ", OUT)
