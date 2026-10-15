library(tidyverse)
library(GenomicRanges)
library(data.table)
library(HMMcopy)
setwd("<MINTCNA_ROOT>/")

args <- commandArgs(trailingOnly = TRUE)
BASE <- if (length(args) >= 1) args[1] else "100kb_c2_N50_2"
print(BASE)
CNA_DIR <- file.path("<MINTCNA_ROOT>/Simulation/cnasim", BASE)  

ref <- read.csv(file.path(CNA_DIR,"ref_all_cells_100kb.csv"),row.names = 1)
l2r <- read.csv(file.path(CNA_DIR,"l2r_smoothed_all_cells_100kb.csv"), row.names = 1)

ref$chr <- paste0("chr", ref$chr)

ord <- with(ref, order(as.numeric(sub("chr", "", chr)), start))
l2r <- l2r[ord,]
ref <- ref[ord,]


# ref: data.frame with columns chr,start,end (bins in order)
# l2r: data.frame/matrix with rows=bins, cols=cells (log2 ratios already)
build_hmm_inputs_raw <- function(ref, l2r, cell_id,
                                 chr_levels = c(paste0("chr",1:22))) {
  stopifnot(nrow(ref) == nrow(l2r))
  stopifnot(all(c("chr","start","end") %in% names(ref)))
  x <- as.numeric(l2r[,cell_id])
  if (is.null(x)) stop("`cell_id` not found in `l2r`")
  
  # chromosomes (no relabeling of values)
  chr <- factor(ref$chr, levels = unique(c(chr_levels, as.character(ref$chr))))
  
  # keep rows required for HMMsegment to run (no value changes)
  keep <- !is.na(x) & !is.na(chr)
  copy <- x[keep]
  chr  <- droplevels(chr[keep])
  autosomes <- chr %in% paste0("chr", 1:22)
  
  somatic_dt <- data.table(
    chr   = as.factor(chr),
    start = ref$start[keep],
    end   = ref$end[keep],
    copy  = copy
  )
  
  somatic_dt
}

# Map one cell’s segments -> per-bin state vector (length = nrow(ref))
seg_to_state_vector <- function(ref, segs) {
  bins <- as.data.table(ref)[, .(chr, start, end)][, bin_id := .I]
  segs <- as.data.table(segs)[, .(chr, start, end, state)]
  # 5/6 -> 5, others unchanged
  segs[, state5 := ifelse(state >= 5, 5L, as.integer(state))]
  
  setkey(bins, chr, start, end)
  setkey(segs, chr, start, end)
  hits <- foverlaps(bins, segs, type = "within", nomatch = NA)
  
  ans <- integer(nrow(bins)); ans[] <- NA_integer_
  ans[hits$bin_id] <- hits$state5
  ans
}

# Build a bin × cell CNA status matrix from a named list of segment tables
# seg_list: e.g., list(cellA = segA, cellB = segB, ...)
segments_to_matrix <- function(ref, seg_list) {
  mat <- sapply(seg_list, function(segs) seg_to_state_vector(ref, segs))
  storage.mode(mat) <- "integer"
  rownames(mat) <- NULL
  mat
}

# 0) list of all cell IDs (columns in your l2r)
cell_ids <- colnames(l2r)

# 1) run HMMsegment for every cell and collect segment tables
seg_list <- vector("list", length(cell_ids))
names(seg_list) <- cell_ids

for (id in cell_ids) {
  print(id)
  inp <- build_hmm_inputs_raw(ref, l2r, cell_id = id)  # your function
  fit <- HMMsegment(inp, verbose = FALSE)
  seg_list[[id]] <- fit$segs
}

# 2) convert segment tables -> bin×cell CNA matrix (5/6 collapsed to 5)
cna_hmmcopy <- segments_to_matrix(ref, seg_list)
ref_hmmcopy <- ref

save(cna_hmmcopy, ref_hmmcopy, file = file.path(CNA_DIR,"Calling/HMMcopy.RData"))
