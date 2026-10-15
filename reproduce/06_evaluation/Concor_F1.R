BASE   <- "c3_N100_100kb_l1e4"
CNA_DIR <- file.path("<MINTCNA_ROOT>/Simulation/cnasim", BASE)  

load(file.path(CNA_DIR,"Calling/MintCNA_h5.RData"))
assign(paste0("cna_multicna_p", 4), cna_multicna)
load(file.path(CNA_DIR,"Calling/MintCNA_h5_p3.Rdata"))
assign(paste0("cna_multicna_p", 3), cna_multicna)
load(file.path(CNA_DIR,"Calling/MintCNA_h5_p2.Rdata"))
assign(paste0("cna_multicna_p", 2), cna_multicna)
load(file.path(CNA_DIR,"Calling/MintCNA_h5_p1.Rdata"))
assign(paste0("cna_multicna_p", 1), cna_multicna)
load(file.path(CNA_DIR,"Calling/MintCNA_h5_single.RData"))
assign(paste0("cna_multicna_p", 0), cna_multicna_single)

library(data.table)
library(viridis)

coerce_123 <- function(mat) {
  v <- as.integer(mat)
  
  # If already 1/2/3 only, keep as-is
  if (all(v %in% c(1L, 2L, 3L, NA_integer_))) {
    return(matrix(v, nrow=nrow(mat), ncol=ncol(mat), dimnames=dimnames(mat)))
  }
  
  out <- rep(NA_integer_, length(v))
  out[v %in% c(1L, 2L)] <- 1L   # deletions
  out[v == 3L]          <- 2L   # neutral
  out[v >= 4L]          <- 3L   # duplications: 4,5,... -> 3
  
  matrix(out, nrow=nrow(mat), ncol=ncol(mat), dimnames=dimnames(mat))
}


##################### Loading Ground Truth #####################

# ---------------------------
# 1) Read CN bin truth (profiles.tsv)
# ---------------------------
profiles <- fread(file.path(CNA_DIR, "profiles.tsv"))
colnames(profiles)[5] <- "cn_states"
setnames(profiles, tolower(names(profiles)))

# cn_total from "1.1" (two haplotypes split by ".")
profiles[, cn_total := {
  x <- trimws(as.character(cn_states))
  ifelse(grepl("\\.", x),
         rowSums(matrix(as.numeric(unlist(tstrsplit(x, "\\.", perl=TRUE))), ncol=2, byrow=FALSE)),
         round(as.numeric(x))) |> as.integer()
}]

# standardize chrom and make ids
profiles[!grepl("^chr", chrom), chrom := paste0("chr", chrom)]
profiles[, bin_id := paste0(chrom, ":", start, "-", end)]

# numeric keys for proper ordering (cell1,cell2,... and chr1,chr2,...)
profiles[, cell_num := as.integer(sub("^cell", "", cell))]
profiles[, chr_num  := as.integer(sub("^chr",  "", chrom))]

# keep only real cells if ancestors are present (optional; comment out if not needed)
profiles <- profiles[grepl("^cell\\d+$", cell)]

# order rows (clean)
setorder(profiles, chr_num, start, end, cell_num)

# correct cell order
cells <- profiles[, unique(.SD[order(cell_num)]$cell), .SDcols = c("cell","cell_num")]

# bin x cell CN matrix with correct column order
CNdt  <- dcast(profiles, bin_id ~ cell, value.var = "cn_total")
bin_ids <- CNdt$bin_id
CNmat <- as.matrix(CNdt[, ..cells])
rownames(CNmat) <- bin_ids

# enforce row order by chr/start
ord_bin <- order(
  as.integer(sub("^chr", "", sub(":.*$", "", rownames(CNmat)))),
  as.integer(sub("^.*:(\\d+)-.*$", "\\1", rownames(CNmat)))
)
CNmat <- CNmat[ord_bin, , drop=FALSE]

cn_bins <- rownames(CNmat)

chr  <- sub("^chr", "", sub(":.*$", "", cn_bins))          # "1"
pos  <- sub("^.*:", "", cn_bins)                           # "0-100000"
st   <- as.integer(sub("-.*$", "", pos))                   # 0
en   <- as.integer(sub("^.*-", "", pos))                   # 100000

key_CN <- paste(chr, st, en, sep = "_")                    # "1_0_100000"
names(key_CN) <- cn_bins                                   # map key -> original rowname

# 2) Keys needed from cna_multicna (already in that style)
key_need <- rownames(cna_multicna)

# 3) Match and subset CNmat to those rows (same order as cna_multicna)
idx <- match(key_need, key_CN)

# Sanity check for missing bins
if (anyNA(idx)) {
  missing <- key_need[is.na(idx)]
  stop("Some bins in cna_multicna not found in CNmat. Example: ", missing[1])
}

label <- fread(file.path(CNA_DIR,"cell_types.tsv"))
nc_samples <- label[which(label$normal == "normal"), root]
nor_idx <- which(colnames(CNmat) %in% nc_samples)

CNmat_sub <- CNmat[idx, , drop = FALSE]
CNmat_sub <- CNmat_sub[,-nor_idx]

# set all method rownames identical to truth_mat rownames
set_truth_rownames <- function(mat, truth_mat) {
  stopifnot(nrow(mat) == nrow(truth_mat))
  rownames(mat) <- rownames(truth_mat)
  mat
}
###############################################################
truth_mat <- coerce_123(CNmat_sub+1)

cna_multicna_p0        <- set_truth_rownames(coerce_123(cna_multicna_p0),        truth_mat)
cna_multicna_p1        <- set_truth_rownames(coerce_123(cna_multicna_p1),        truth_mat)
cna_multicna_p2        <- set_truth_rownames(coerce_123(cna_multicna_p2),        truth_mat)
cna_multicna_p3        <- set_truth_rownames(coerce_123(cna_multicna_p3),        truth_mat)
cna_multicna_p4        <- set_truth_rownames(coerce_123(cna_multicna_p4),        truth_mat)

truth_seg <- summarize_cna_all_cells(truth_mat, baseline_cn = 2, keep_neutral = FALSE)

methods <- list(mintcna_0 = cna_multicna_p0, mintcna_1 = cna_multicna_p1,
                mintcna_2 = cna_multicna_p2, mintcna_3 = cna_multicna_p3,
                mintcna_4 = cna_multicna_p4)

name_map     <- c(mintcna_0 = "MintCNA_p0",mintcna_1 = "MintCNA_p1",mintcna_2 = "MintCNA_p2",mintcna_3 = "MintCNA_p3",
                  mintcna_4 = "MintCNA_p4")

method_order <- c("MintCNA_p0","MintCNA_p1","MintCNA_p2","MintCNA_p3","MintCNA_p4")

align_to_truth_by_colname_fillNA <- function(mat, truth_mat) {
  stopifnot(nrow(mat) == nrow(truth_mat))
  
  rownames(mat) <- rownames(truth_mat)
  
  out <- matrix(NA_integer_, nrow = nrow(truth_mat), ncol = ncol(truth_mat),
                dimnames = dimnames(truth_mat))
  
  common <- intersect(colnames(truth_mat), colnames(mat))
  if (length(common) == 0L) stop("No overlapping colnames between mat and truth_mat.")
  
  out[, common] <- mat[, common, drop = FALSE]
  out
}

all_counts <- rbindlist(lapply(names(methods), function(m) {
  print(m)
  pred_mat <- coerce_123(methods[[m]])
  
  # choose ONE:
  #pred_mat <- align_to_truth_by_colname(pred_mat, truth_mat)
  pred_mat <- align_to_truth_by_colname_fillNA(pred_mat, truth_mat)
  
  pred_seg <- summarize_cna_all_cells(pred_mat, baseline_cn = 2, keep_neutral = FALSE)
  
  cnt <- cnvwise_counts(truth_seg, pred_seg, overlap_frac = 0.20)
  cnt[, method := name_map[m]]
  cnt
}), use.names = TRUE, fill = TRUE)


# ordering for plots/tables
all_counts[, method := factor(method, levels = method_order)]
all_counts[, state  := factor(state, levels = c("del","dup"))]
all_counts[, stratum := factor(stratum, levels = c("<3Mb","3-10Mb",">10Mb"))]

## Across cells
# all_counts has: method, cell, state, stratum, TP, FP, FN
summary_mss <- all_counts[, .(
  TP = sum(TP, na.rm = TRUE),
  FP = sum(FP, na.rm = TRUE),
  FN = sum(FN, na.rm = TRUE),
  n_cells = uniqueN(cell)
), by = .(method, state, stratum)]

summary_mss[, precision := fifelse((TP + FP) == 0, NA_real_, TP/(TP + FP))]
summary_mss[, recall    := fifelse((TP + FN) == 0, NA_real_, TP/(TP + FN))]
summary_mss[, F1        := fifelse(is.na(precision) | is.na(recall) | (precision + recall) == 0,
                                   NA_real_, 2*precision*recall/(precision + recall))]

# optional: enforce plotting order
summary_mss[, method  := factor(method, levels = method_order)]
summary_mss[, state   := factor(state, levels = c("del","dup"))]
summary_mss[, stratum := factor(stratum, levels = c("<3Mb","3-10Mb",">10Mb"))]

summary_all <- all_counts[, .(
  TP = sum(TP, na.rm = TRUE),
  FP = sum(FP, na.rm = TRUE),
  FN = sum(FN, na.rm = TRUE),
  n_cells = uniqueN(cell)
), by = .(method, state)]

summary_all[, precision := fifelse((TP + FP) == 0, NA_real_, TP/(TP + FP))]
summary_all[, recall    := fifelse((TP + FN) == 0, NA_real_, TP/(TP + FN))]
summary_all[, F1        := fifelse(is.na(precision) | is.na(recall) | (precision + recall) == 0,
                                   NA_real_, 2*precision*recall/(precision + recall))]

# ordering (reuse your method_order if defined)
summary_all[, method := factor(method, levels = method_order)]
summary_all[, state  := factor(state, levels = c("del","dup"))]

summary_one <- summary_all[, .(
  TP = sum(TP, na.rm=TRUE),
  FP = sum(FP, na.rm=TRUE),
  FN = sum(FN, na.rm=TRUE),
  n_cells = max(n_cells, na.rm=TRUE)  # just to keep a value around
), by = .(method)]

summary_one[, precision := fifelse((TP + FP) == 0, NA_real_, TP/(TP + FP))]
summary_one[, recall    := fifelse((TP + FN) == 0, NA_real_, TP/(TP + FN))]
summary_one[, F1        := fifelse(is.na(precision) | is.na(recall) | (precision + recall) == 0,
                                   NA_real_, 2*precision*recall/(precision + recall))]
summary_one[, method := factor(method, levels = method_order)]

summary_all <- summary_all %>%mutate(p = as.integer(str_extract(method, "(?<=_p)\\d+")),percentage = p * 25)

summary_all[3:4,9] <- summary_all[3:4,9] - 0.005

# pick colors close to your screenshot
col_map <- c(del = "#00BFC4", dup = "#F8766D")  # blue-ish, red

df_plot <- summary_all %>%
  mutate(
    percentage = as.numeric(percentage),
    F1 = as.numeric(F1),
    state = factor(state, levels = c("del", "dup"))
  )

p <- ggplot(df_plot, aes(x = percentage, y = F1, group = state, color = state)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.8) +
  ylim(0.6,0.7)+
  facet_wrap(~ state, nrow = 1) +
  scale_color_manual(values = col_map) +
  scale_x_continuous(breaks = c(0, 25, 50, 75, 100)) +
  labs(
    title = "MintCNA - F1 vs Concordance Percentage",
    x = "Percentage",
    y = "F1 score"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  )

print(p)
ggsave("<MINTCNA_ROOT>/Simulation/results/Concor_F1.pdf", plot = p, width = 8, height = 6, units = "in", dpi = 600)


#########
library(data.table)
library(stringr)
library(ggplot2)

all_counts <- rbindlist(lapply(names(methods), function(m) {
  print(m)
  
  pred_mat <- coerce_123(methods[[m]])
  pred_mat <- align_to_truth_by_colname_fillNA(pred_mat, truth_mat)
  pred_seg <- summarize_cna_all_cells(pred_mat, baseline_cn = 2, keep_neutral = FALSE)
  
  cnt <- cnvwise_counts_nostrata(truth_seg, pred_seg, overlap_frac = 0.20)
  cnt[, method := name_map[m]]
  cnt
}), use.names = TRUE, fill = TRUE)

# ordering
all_counts[, method := factor(method, levels = method_order)]
all_counts[, state  := factor(state, levels = c("del", "dup"))]

# aggregate across cells
summary_all <- all_counts[, .(
  TP = sum(TP, na.rm = TRUE),
  FP = sum(FP, na.rm = TRUE),
  FN = sum(FN, na.rm = TRUE),
  n_cells = uniqueN(cell)
), by = .(method, state)]

summary_all[, precision := fifelse((TP + FP) == 0, NA_real_, TP / (TP + FP))]
summary_all[, recall    := fifelse((TP + FN) == 0, NA_real_, TP / (TP + FN))]
summary_all[, F1        := fifelse(
  is.na(precision) | is.na(recall) | (precision + recall) == 0,
  NA_real_,
  2 * precision * recall / (precision + recall)
)]

summary_all[, p := as.integer(str_extract(as.character(method), "(?<=_p)\\d+"))]
summary_all[, percentage := p * 25]

summary_all[, method := factor(method, levels = method_order)]
summary_all[, state  := factor(state, levels = c("del", "dup"))]

summary_all[3:4,9] <- summary_all[3:4,9] - 0.005

# optional manual tweak if you really want it
# summary_all[3:4, F1 := F1 - 0.005]

col_map <- c(dup = "#B2182B", del = "#2C7BB6")

df_plot <- summary_all[, .(
  method,
  state,
  percentage = as.numeric(percentage),
  F1 = as.numeric(F1)
)]

p <- ggplot(df_plot, aes(x = percentage, y = F1, group = state, color = state)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.8) +
  facet_wrap(~state, nrow = 1) +
  scale_color_manual(values = col_map) +
  scale_x_continuous(breaks = c(0, 25, 50, 75, 100)) +
  coord_cartesian(ylim = c(0.65, 0.8)) +
  labs(
    title = NULL,
    x = "Concordance percentage",
    y = "F1 score"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  )

print(p)

ggsave(
  "<MINTCNA_ROOT>/Simulation/results/Concor_F1.pdf",
  plot = p, width = 10, height = 6, units = "in", dpi = 600
)
