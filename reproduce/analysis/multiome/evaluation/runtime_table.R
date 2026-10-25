#!/usr/bin/env Rscript
# Runtime and peak memory of every method's benchmark jobs from the SLURM accounting database
# (Supplementary Fig. 9b). CPU time (TotalCPU) is the comparable quantity: the jobs ran with different
# core counts, so wall time is not. Job names are matched to methods by the patterns in MAP.
#   Rscript reproduce/analysis/multiome/evaluation/runtime_table.R
# Output: ATAC_RNA/Result/qc/runtime_by_array.tsv (one row per SLURM array), runtime_by_method.tsv
# Run from the root of the MintCNA repository with MINTCNA_DATA set (see reproduce/README.md).
suppressPackageStartupMessages(library(data.table))
ROOT <- Sys.getenv("MINTCNA_DATA"); if (!nzchar(ROOT)) stop("set MINTCNA_DATA (see reproduce/README.md)")
setwd(file.path(ROOT, "ATAC_RNA"))

raw <- fread(cmd = paste("sacct -S 2026-08-01 -u $USER --parsable2 --noheader",
                         "--format=JobName,JobID,ReqCPUS,ReqMem,Elapsed,TotalCPU,MaxRSS,State"),
             header = FALSE, sep = "|", col.names = c("job", "id", "cpus", "mem", "elapsed", "cputime", "maxrss", "state"))
raw[, base := sub("\\.(batch|extern|[0-9]+)$", "", id)]            # the .batch step carries MaxRSS, the parent the name
parent <- raw[id == base, .(base, job, cpus, mem, elapsed, cputime, state)]
rss    <- raw[grepl("\\.batch$", id), .(base, maxrss)]
d <- merge(parent, rss, by = "base", all.x = TRUE)[state == "COMPLETED"]

tosec <- function(x) {                                              # [DD-]HH:MM:SS[.ss] -> seconds
  x <- sub("\\..*$", "", x)
  dd <- ifelse(grepl("-", x), as.numeric(sub("-.*", "", x)), 0)
  p <- do.call(rbind, lapply(strsplit(sub(".*-", "", x), ":"), function(v) as.numeric(rev(v))[1:3])); p[is.na(p)] <- 0
  dd * 86400 + p[, 3] * 3600 + p[, 2] * 60 + p[, 1]
}
tokb <- function(x) { n <- as.numeric(sub("[KMGT]$", "", x)); u <- sub(".*[0-9]", "", x); n * c(K = 1, M = 1024, G = 1024^2, T = 1024^3)[u] }
d[, `:=`(wall_s = tosec(elapsed), cpu_s = tosec(cputime), rss_gb = tokb(maxrss) / 1024^2)]

MAP <- c(epiAneu = "epiAneufinder", atacnv = "AtaCNV", copyscat = "Copy-scAT", riddler = "RIDDLER", congas = "CONGAS+",
         infercnv = "inferCNV", agdcae = "MintCNA (AG-DCAE)", mint_scan = "MintCNA (scan)", mint_call = "MintCNA (call)")
d[, method := NA_character_]
for (k in names(MAP)) d[is.na(method) & grepl(k, job, ignore.case = TRUE), method := MAP[k]]
d <- d[!is.na(method)]

d[, array_id := sub("_.*", "", base)]
byarr <- d[, .(n_task = .N, wall_med_min = round(median(wall_s) / 60, 1), cpu_med_min = round(median(cpu_s) / 60, 1),
               cpus = as.integer(median(cpus)), peak_gb = round(max(rss_gb, na.rm = TRUE), 1)), by = .(method, array_id)]
fwrite(byarr[order(method, -n_task)], "Result/qc/runtime_by_array.tsv", sep = "\t")
t <- d[, .(n_jobs = .N, wall_median_min = round(median(wall_s) / 60, 1), cpu_median_min = round(median(cpu_s) / 60, 1),
           cpu_total_h = round(sum(cpu_s) / 3600, 1), cpus = as.integer(median(cpus)), peak_gb = round(max(rss_gb, na.rm = TRUE), 1)),
       by = method][order(cpu_median_min)]
fwrite(t, "Result/qc/runtime_by_method.tsv", sep = "\t")
print(t)
