#!/usr/bin/env Rscript
# rcCAE input for one simulation setting: bins x cells read counts with GC and mappability columns,
# comma separated, first line "#bin size = <bp>" (the format of the rcCAE example data).
#   Rscript reproduce/analysis/simulation/methods/rcCAE_input.R <setting>      (called by rcCAE_run.sh)
# Input : <setting>/DNA_RC.RData (Y, gr)
# Output: <setting>/Calling/rcCAE/rcCAE_input.txt
# Run from the root of the MintCNA repository with MINTCNA_DATA set (see reproduce/README.md).
suppressPackageStartupMessages({ library(data.table); library(GenomicRanges) })
ROOT <- Sys.getenv("MINTCNA_DATA"); if (!nzchar(ROOT)) stop("set MINTCNA_DATA (see reproduce/README.md)")
args <- commandArgs(trailingOnly = TRUE)
BASE <- if (length(args) >= 1) args[1] else "c3_N100_100kb_l1e4"
CNA_DIR <- file.path(ROOT, "Simulation/cnasim", BASE)
load(file.path(CNA_DIR, "DNA_RC.RData"))
BIN_SIZE <- 100000L

dt <- data.table(chr = as.integer(sub("^chr", "", as.character(seqnames(gr)), ignore.case = TRUE)),
                 `gc-content` = as.numeric(gr$gc), mappability = as.numeric(gr$mapp))
Ydt <- as.data.table(Y); setnames(Ydt, paste0("cell", seq_len(ncol(Ydt))))
dt <- cbind(dt, Ydt)

out_dir <- file.path(CNA_DIR, "Calling/rcCAE"); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
tmp <- tempfile(); fwrite(dt, tmp, sep = ",", col.names = TRUE)
writeLines(c(paste0("#bin size = ", BIN_SIZE), readLines(tmp)), file.path(out_dir, "rcCAE_input.txt"))
unlink(tmp)
cat("wrote", file.path(out_dir, "rcCAE_input.txt"), "\n")
