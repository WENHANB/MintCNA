args <- commandArgs(trailingOnly = TRUE)
BASE <- if (length(args) >= 1) args[1] else "c3_N100_100kb_l1e4"
print(BASE)
CNA_DIR <- file.path("<MINTCNA_ROOT>/Simulation/cnasim", BASE)  

library(data.table)
library(mclust)

## ---- 1) truth labels (from CNAsim) ----
truth <- fread(file.path(CNA_DIR, "cell_types.tsv"))   # columns: root, normal
setnames(truth, c("root","normal"), c("cell","truth"))
truth <- truth[grepl("^cell\\d+$", cell)]

## ---- 2) rcCAE predicted labels (out/label.txt) ----
read_rccae_label <- function(path) {
  x <- readLines(path, warn = FALSE)
  stopifnot(length(x) >= 2)
  cells <- trimws(unlist(strsplit(x[1], ",")))
  labs  <- as.integer(trimws(unlist(strsplit(x[2], ","))))
  data.table(cell = cells, pred = labs)
}

pred <- read_rccae_label(file.path(CNA_DIR, "Calling/rcCAE/out/label.txt"))

## ---- 3) align and compute ARI ----
dt <- merge(truth, pred, by = "cell", all = FALSE)

# ARI including normal
ari_all <- adjustedRandIndex(as.factor(dt$truth), as.factor(dt$pred))

# ARI on tumor cells only (exclude "normal")
dt_tumor <- dt[truth != "normal"]
ari_tumor <- adjustedRandIndex(as.factor(dt_tumor$truth), as.factor(dt_tumor$pred))

list(
  n_matched = nrow(dt),
  ARI_all   = ari_all,
  ARI_tumor = ari_tumor
)

um_rccae <- data.frame(
  cell  = rownames(latent_use),
  UMAP1 = um[,1],
  UMAP2 = um[,2],
  truth = y_true,
  pred  = as.factor(dt_tumor$pred)
)

ggplot(um_rccae, aes(UMAP1, UMAP2, color = pred)) +
  geom_point(size=2) + theme_bw() + labs(title="rcCAE")
