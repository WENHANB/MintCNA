#!/usr/bin/env Rscript
# Per-method runtime and memory, pulled from the SLURM accounting database.
#
# WHAT THIS CAN AND CANNOT SAY. The methods were not run under a common
# allocation -- each was sized from its own measured peak, which was the right
# call for throughput and the wrong one for a runtime comparison. So:
#   - WALL time is not comparable across methods (different --cpus-per-task).
#   - CPU time (TotalCPU) is the closest thing to comparable, since it counts
#     core-seconds regardless of how many cores were given. It still rewards a
#     method that parallelises badly, because idle cores cost nothing in it.
#   - Two methods do NOT run their own published input path: RIDDLER's
#     window_matrix and Copy-scAT's binning were replaced with data.table
#     equivalents (verified identical) because the originals were unusable on a
#     10x Multiome fragments file. Their times are therefore FASTER than the
#     published code would give.
# All three caveats belong in any table that reaches the manuscript.

suppressPackageStartupMessages(library(data.table))
setwd("<MINTCNA_ROOT>/ATAC_RNA")

raw <- fread(cmd = paste("sacct -S 2026-08-01 -u $USER --parsable2 --noheader",
  "--format=JobName,JobID,ReqCPUS,ReqMem,Elapsed,TotalCPU,MaxRSS,State"),
  header = FALSE, sep = "|",
  col.names = c("job","id","cpus","mem","elapsed","cputime","maxrss","state"))

# the .batch step carries MaxRSS; the parent carries the name -- join them
raw[, base := sub("\\.(batch|extern|[0-9]+)$", "", id)]
parent <- raw[id == base, .(base, job, cpus, mem, elapsed, cputime, state)]
rss    <- raw[grepl("\\.batch$", id), .(base, maxrss)]
d <- merge(parent, rss, by = "base", all.x = TRUE)[state == "COMPLETED"]

tosec <- function(x) {                       # [DD-]HH:MM:SS[.ss]
  x <- sub("\\..*$", "", x)
  dd <- ifelse(grepl("-", x), as.numeric(sub("-.*", "", x)), 0)
  hms <- sub(".*-", "", x)
  p <- do.call(rbind, lapply(strsplit(hms, ":"), function(v) as.numeric(rev(v))[1:3]))
  p[is.na(p)] <- 0
  dd*86400 + p[,3]*3600 + p[,2]*60 + p[,1]
}
tokb <- function(x) { n <- as.numeric(sub("[KMGT]$", "", x))
  u <- sub(".*[0-9]", "", x); n * c(K=1, M=1024, G=1024^2, T=1024^3)[u] }

d[, `:=`(wall_s = tosec(elapsed), cpu_s = tosec(cputime), rss_gb = tokb(maxrss)/1024^2)]

MAP <- c(epiAneu="epiAneufinder", atacnv="AtaCNV", copyscat="Copy-scAT",
         riddler="RIDDLER", congas="CONGAS+", infercnv="inferCNV",
         agdcae="MintCNA (AG-DCAE)", mint_scan="MintCNA (scan)", mint_call="MintCNA (call)")
d[, method := NA_character_]
for (k in names(MAP)) d[is.na(method) & grepl(k, job, ignore.case = TRUE), method := MAP[k]]
d <- d[!is.na(method)]

# One SLURM ARRAY is one run of a method over a set of samples; a task inside it
# is one sample -- EXCEPT for the parameter sweeps, which hold the sample fixed
# and vary a setting. Group by array so the two can be told apart, and report the
# per-sample median only from arrays that actually span samples.
d[, array_id := sub("_.*", "", base)]
byarr <- d[, .(n_task = .N, wall_med_min = round(median(wall_s)/60,1),
               cpu_med_min = round(median(cpu_s)/60,1), cpus = as.integer(median(cpus)),
               peak_gb = round(max(rss_gb, na.rm=TRUE),1)), by = .(method, array_id)]
fwrite(byarr[order(method, -n_task)], "Result/qc/runtime_by_array.tsv", sep = "\t")
cat("=== every array job, so sweeps can be excluded by hand ===\n")
print(byarr[order(method, -n_task)])
cat("\n")

t <- d[, .(n_jobs = .N,
           wall_median_min = round(median(wall_s)/60, 1),
           cpu_median_min  = round(median(cpu_s)/60, 1),
           cpu_total_h     = round(sum(cpu_s)/3600, 1),
           cpus            = as.integer(median(cpus)),
           peak_gb         = round(max(rss_gb, na.rm = TRUE), 1)), by = method][order(cpu_median_min)]
fwrite(t, "Result/qc/runtime_by_method.tsv", sep = "\t")
print(t)
cat("\nwall_median_min is NOT comparable across rows (different --cpus-per-task).\n")
cat("cpu_median_min is the comparable column. RIDDLER and Copy-scAT ran on\n")
cat("re-implemented input steps and are faster here than their published code.\n")
