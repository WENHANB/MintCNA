#!/usr/bin/env Rscript
# Copy-scAT (R package CopyscAT 0.40) on one multiome sample, following the authors' tutorial
# (copyscat_tutorial.R, option 1: control built from all cells with dummy cells) with its parameters.
# The Python fragment binning of the package is replaced by an equivalent data.table step: 1 Mb bins
# assigned by fragment start, bin value = summed fragment length; the bin axis is the package's hg38 grid.
#   Rscript reproduce/analysis/multiome/methods/run_copyscat.R <care_id>
# Input : ATAC_RNA/Data/atac_filtered/<care_id>_fragments.tsv.gz,
#         ATAC_RNA/Data/ref/copyscat/{hg38_chrom_sizes.tsv, hg38_1e+06_cytoband_densities_granges.tsv,
#         hg38_1e+06_cpg_densities.tsv}  (the package's hg38_references/ folder)
# Output: ATAC_RNA/Result/benchmark/CopyscAT/<care_id>/  (copyscat_result.rds and the package's CSV/PDF output)
# Run from the root of the MintCNA repository with MINTCNA_DATA set (see reproduce/README.md).
suppressPackageStartupMessages({ library(data.table); library(dplyr); library(CopyscAT) })
ROOT <- Sys.getenv("MINTCNA_DATA"); if (!nzchar(ROOT)) stop("set MINTCNA_DATA (see reproduce/README.md)")
args <- commandArgs(trailingOnly = TRUE); if (!length(args)) stop("need <care_id>")
CARE <- args[1]
ROOT <- file.path(ROOT, "ATAC_RNA")
OUT  <- file.path(ROOT, "Result/benchmark/CopyscAT", CARE); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
REFDIR <- file.path(ROOT, "Data/ref/copyscat")
ref_chrom <- file.path(REFDIR, "hg38_chrom_sizes.tsv")
ref_cyto  <- file.path(REFDIR, "hg38_1e+06_cytoband_densities_granges.tsv")
ref_cpg   <- file.path(REFDIR, "hg38_1e+06_cpg_densities.tsv")
stopifnot(file.exists(ref_chrom), file.exists(ref_cyto), file.exists(ref_cpg))
BINSIZE <- 1e6; MINFRAGS <- 1e4                                      # tutorial values

# ---- bin the fragments: cells x bins, bins exactly as in the reference tables -----------------------
ref_grid <- fread(ref_cpg, header = FALSE, col.names = c("chr", "start", "end", "cpg"))
ref_grid[, bin := paste0(chr, "_", ((start - 1) %/% BINSIZE) + 1L)]
fr <- fread(cmd = paste("zcat", shQuote(file.path(ROOT, "Data/atac_filtered", paste0(CARE, "_fragments.tsv.gz")))),
            header = FALSE, select = 1:4, col.names = c("chr", "start", "end", "barcode"))
fr <- fr[chr %in% unique(ref_grid$chr)]
fr[, len := end - start - 1L]
fr[, bin := paste0(chr, "_", (start %/% BINSIZE) + 1L)]
keep <- fr[, .N, by = barcode][N > MINFRAGS, barcode]
fr <- fr[barcode %in% keep]
wide <- dcast(fr[, .(v = sum(len)), by = .(bin, barcode)], barcode ~ bin, value.var = "v", fill = 0)
miss <- setdiff(ref_grid$bin, names(wide)); if (length(miss)) wide[, (miss) := 0L]
setcolorder(wide, c("barcode", ref_grid$bin))
setnames(wide, "barcode", "chrom")
stopifnot(ncol(wide) - 1L == nrow(ref_grid))
tsv <- file.path(OUT, "copyscat_input.tsv"); fwrite(wide, tsv, sep = "\t")
message(sprintf("%s: %d cells with > %d fragments x %d bins", CARE, nrow(wide), MINFRAGS, ncol(wide) - 1L))

# ---- Copy-scAT, tutorial option 1 ---------------------------------------------------------------
t0 <- Sys.time()
initialiseEnvironment(genomeFile = ref_chrom, cytobandFile = ref_cyto, cpgFile = ref_cpg, binSize = BINSIZE,
                      minFrags = MINFRAGS, cellSuffix = c("-1", "-2"), lowerTrim = 0.5, upperTrim = 0.8)
setOutputFile(OUT, CARE)
scData        <- readInputTable(tsv)
scData_k_norm <- normalizeMatrixN(scData, logNorm = FALSE, maxZero = 2000, imputeZeros = FALSE, blacklistProp = 0.8,
                                  blacklistCutoff = 125, dividingFactor = 1, upperFilterQuantile = 0.95)
summaryFunction <- cutAverage
scData_collapse <- collapseChrom3N(scData_k_norm, summaryFunction = summaryFunction, binExpand = 1, minimumChromValue = 100,
                                   logTrans = FALSE, tssEnrich = 1, logBase = 2, minCPG = 300, powVal = 0.73)
scData_collapse <- filterCells(scData_collapse, minimumSegments = 40, minDensity = 0.1)
median_iqr <- computeCenters(scData_collapse, summaryFunction = summaryFunction)
candidate_cnvs <- identifyCNVClusters(scData_collapse, median_iqr, useDummyCells = TRUE, propDummy = 0.25, minMix = 0.01,
                                      deltaMean = 0.03, deltaBIC2 = 0.25, bicMinimum = 0.1, subsetSize = 600, fakeCellSD = 0.08,
                                      uncertaintyCutoff = 0.55, summaryFunction = summaryFunction, maxClust = 4, mergeCutoff = 3,
                                      IQRCutoff = 0.2, medianQuantileCutoff = 0.4)
candidate_cnvs_clean <- clusterCNV(initialResultList = candidate_cnvs, medianIQR = candidate_cnvs[[3]], minDiff = 1.5)
final_cnv_list <- annotateCNV4(candidate_cnvs_clean, saveOutput = TRUE, outputSuffix = "clean_cnv", sdCNV = 0.5,
                               filterResults = TRUE, filterRange = 0.8)

saveRDS(list(care_id = CARE, n_cells = nrow(wide), n_bins = ncol(wide) - 1L, collapse = scData_collapse,
             median_iqr = median_iqr, candidates = candidate_cnvs_clean, final = final_cnv_list),
        file.path(OUT, "copyscat_result.rds"))
cat(sprintf("SUMMARY\t%s\t%d\t%d\t%.1f min\n", CARE, nrow(wide), ncol(wide) - 1L, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
