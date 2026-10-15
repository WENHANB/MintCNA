#!/usr/bin/env Rscript
# Copy-scAT on one sample.  Usage: Rscript run_copyscat.R <care_id>
#
# The call sequence and every parameter below are copied verbatim from the
# authors' own copyscat_tutorial.R (fetched to Data/ref/copyscat/), OPTION 1 --
# the documented main path, which builds its 'normal' control from all cells via
# dummy cells rather than from a supplied normal list. OPTION 2 (NMF-detected or
# user-supplied normals, annotateCNV4B) is described by the authors as "still
# somewhat experimental - use if the copy number calls from A don't make sense",
# so it is not the default and is not what we run.
#
# The Python helper process_fragment_file.py is reimplemented in data.table for
# speed, keeping its two defining conventions: bins are assigned by fragment
# START coordinate, and a bin's value is the SUMMED FRAGMENT LENGTH in bp
# (end - start - 1), not a fragment count. Copying those is what makes this a
# faithful Copy-scAT run rather than a variant.

args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("need <care_id>")
CARE <- args[1]

suppressPackageStartupMessages({library(data.table); library(dplyr); library(CopyscAT)})

ROOT <- "<MINTCNA_ROOT>/ATAC_RNA"
OUT  <- file.path(ROOT, "Result/benchmark/CopyscAT", CARE)
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# The hg38 reference tables live in the repo's hg38_references/ folder, which is
# outside inst/ and therefore is NOT carried into the installed package.
REFDIR <- file.path(ROOT, "Data/ref/copyscat")
ref_chrom <- file.path(REFDIR, "hg38_chrom_sizes.tsv")
ref_cyto  <- file.path(REFDIR, "hg38_1e+06_cytoband_densities_granges.tsv")
ref_cpg   <- file.path(REFDIR, "hg38_1e+06_cpg_densities.tsv")
stopifnot(file.exists(ref_chrom), file.exists(ref_cyto), file.exists(ref_cpg))

BINSIZE  <- 1e6     # tutorial
MINFRAGS <- 1e4     # tutorial

# ---- bin fragments the Copy-scAT way --------------------------------------
# collapseChrom3N() attaches arm labels and CpG density by ROW INDEX,
#   mutate(chromArm = cytoband_data$V4, cpg = cpg_data$cpg_density)
# not by coordinate, so the bin axis must carry exactly the reference's 3,103
# entries in exactly its order -- including chrM, chrX, chrY and all-zero bins.
# Dropping empty bins or sex chromosomes makes dplyr error, or worse silently
# recycle and mislabel every bin. So the grid is taken FROM the reference file.
ref_grid <- fread(ref_cpg, header = FALSE, col.names = c("chr", "start", "end", "cpg"))
ref_grid[, bin := paste0(chr, "_", ((start - 1) %/% BINSIZE) + 1L)]

frag <- file.path(ROOT, "Data/atac_filtered", paste0(CARE, "_fragments.tsv.gz"))
message("reading fragments ...")
fr <- fread(cmd = paste("zcat", shQuote(frag)), header = FALSE, select = 1:4,
            col.names = c("chr", "start", "end", "barcode"))
fr <- fr[chr %in% unique(ref_grid$chr)]
fr[, len := end - start - 1L]                             # Copy-scAT's bin value
fr[, bin := paste0(chr, "_", (start %/% BINSIZE) + 1L)]   # assigned by START

tot  <- fr[, .(n = .N), by = barcode]
keep <- tot[n > MINFRAGS, barcode]                        # strict '>' as in the helper
message(sprintf("cells passing minFrags=%d: %d / %d", MINFRAGS, length(keep), nrow(tot)))
fr <- fr[barcode %in% keep]

agg <- fr[, .(v = sum(len)), by = .(bin, barcode)]

# ORIENTATION: the table must be CELLS AS ROWS, bins as columns, with the cell
# barcode in the first column. readInputTable() does
#   colnames(scData)[1] <- "Cell_id"; rownames(scData) <- scData$Cell_id
# and everything downstream assumes it: normalizeMatrixN starts with
# sc_t <- t(inputMatrix) and treats the COLUMNS of sc_t as cells. Writing
# bins-as-rows (the natural choice, and what an earlier version of this script
# did) is accepted silently and then computes library sizes over the bin axis --
# 3,103 pseudo-cells, one summing to zero, which is what made cpm() fail with
# "library sizes should be finite and non-negative".
wide <- dcast(agg, barcode ~ bin, value.var = "v", fill = 0)
miss <- setdiff(ref_grid$bin, names(wide))
if (length(miss)) wide[, (miss) := 0L]
setcolorder(wide, c("barcode", ref_grid$bin))
setnames(wide, "barcode", "chrom")        # first column header is discarded by readInputTable
stopifnot(ncol(wide) - 1L == nrow(ref_grid))
tsv <- file.path(OUT, "copyscat_input.tsv")
fwrite(wide, tsv, sep = "\t")
message(sprintf("input matrix: %d cells x %d bins", nrow(wide), ncol(wide) - 1L))

# ---- run Copy-scAT, tutorial OPTION 1 --------------------------------------
t0 <- Sys.time()
initialiseEnvironment(
  genomeFile   = ref_chrom,
  cytobandFile = ref_cyto,
  cpgFile      = ref_cpg,
  binSize      = BINSIZE,
  minFrags     = MINFRAGS,
  cellSuffix   = c("-1", "-2"),
  lowerTrim    = 0.5,
  upperTrim    = 0.8)

# Without this the package writes its PDFs and CSVs to the process working
# directory under an empty prefix -- that is where the stray _signal_distribution.pdf
# in the project root came from.
setOutputFile(OUT, CARE)

scData        <- readInputTable(tsv)
scData_k_norm <- normalizeMatrixN(scData, logNorm = FALSE, maxZero = 2000,
                                  imputeZeros = FALSE, blacklistProp = 0.8,
                                  blacklistCutoff = 125, dividingFactor = 1,
                                  upperFilterQuantile = 0.95)

summaryFunction <- cutAverage
scData_collapse <- collapseChrom3N(scData_k_norm, summaryFunction = summaryFunction,
                                   binExpand = 1, minimumChromValue = 100,
                                   logTrans = FALSE, tssEnrich = 1,
                                   logBase = 2, minCPG = 300, powVal = 0.73)
scData_collapse <- filterCells(scData_collapse, minimumSegments = 40, minDensity = 0.1)
message(sprintf("cells after filterCells: %d", ncol(scData_collapse) - 1L))

# identifyCNVClusters' SECOND POSITIONAL ARGUMENT is median_iqr, not
# summaryFunction -- and it comes from computeCenters, a step it is easy to skip.
median_iqr <- computeCenters(scData_collapse, summaryFunction = summaryFunction)

candidate_cnvs <- identifyCNVClusters(scData_collapse, median_iqr,
                                      useDummyCells = TRUE, propDummy = 0.25,
                                      minMix = 0.01, deltaMean = 0.03,
                                      deltaBIC2 = 0.25, bicMinimum = 0.1,
                                      subsetSize = 600, fakeCellSD = 0.08,
                                      uncertaintyCutoff = 0.55,
                                      summaryFunction = summaryFunction,
                                      maxClust = 4, mergeCutoff = 3,
                                      IQRCutoff = 0.2, medianQuantileCutoff = 0.4)

candidate_cnvs_clean <- clusterCNV(initialResultList = candidate_cnvs,
                                   medianIQR = candidate_cnvs[[3]], minDiff = 1.5)

final_cnv_list <- annotateCNV4(candidate_cnvs_clean, saveOutput = TRUE,
                               outputSuffix = "clean_cnv", sdCNV = 0.5,
                               filterResults = TRUE, filterRange = 0.8)
el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))

saveRDS(list(care_id = CARE, n_cells = nrow(wide), n_bins = ncol(wide) - 1L,
             collapse = scData_collapse, median_iqr = median_iqr,
             candidates = candidate_cnvs_clean, final = final_cnv_list),
        file.path(OUT, "copyscat_result.rds"))
cat(sprintf("SUMMARY\t%s\t%d\t%d\t%.1f\n", CARE, nrow(wide), ncol(wide) - 1L, el))
