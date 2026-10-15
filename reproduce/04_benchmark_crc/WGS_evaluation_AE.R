load("<MINTCNA_ROOT>/Result/CRC/Benchmark/SCOPE.RData")
load("<MINTCNA_ROOT>/Result/CRC/Benchmark/MultiCNA_AE_197_d250.RData")
load("<MINTCNA_ROOT>/Result/CRC/Benchmark/HMMcopy.RData")
load("<MINTCNA_ROOT>/Result/CRC/Benchmark/DNAcopy.RData")
load("<MINTCNA_ROOT>/Result/CRC/Benchmark/SCCNV.RData")
load("<MINTCNA_ROOT>/Result/CRC/Benchmark/SCONCE.RData")
load("<MINTCNA_ROOT>/Result/CRC/Benchmark/rcCAE.RData")

colnames(ref_scope)[1] <- "chr";
cna_scope <- cna_scope[,colnames(cna_scope) %in% colnames(cna_multicna_ae)]
cna_scope <- cna_scope + 1L; cna_scope[cna_scope >= 5] <- 5L 
cna_hmmcopy <- cna_hmmcopy[,colnames(cna_hmmcopy) %in% colnames(cna_multicna_ae)]
cna_sconce <- cna_sconce[,colnames(cna_sconce) %in% colnames(cna_multicna_ae)]
cna_sccnv <- cna_sccnv[,colnames(cna_sccnv) %in% colnames(cna_multicna_ae)]
cna_rccae <- cna_rccae[,colnames(cna_rccae) %in% colnames(cna_multicna_ae)]

library(data.table)
library(viridis)
library(ggplot2)

## --- helpers from before (keep as-is) ---
coerce_123 <- function(mat) {
  v <- as.integer(mat)
  # If already 1/2/3 only, keep as-is
  if (all(v %in% c(1L,2L,3L, NA_integer_))) {
    return(matrix(v, nrow = nrow(mat), ncol = ncol(mat),
                  dimnames = dimnames(mat)))
  }
  # Otherwise map 1..5 -> 1/2/3
  out <- rep(NA_integer_, length(v))
  out[v %in% c(1L,2L)] <- 1L  # deletions
  out[v %in% 3L]        <- 2L  # neutral
  out[v %in% c(4L,5L)] <- 3L  # duplications
  matrix(out, nrow = nrow(mat), ncol = ncol(mat),
         dimnames = dimnames(mat))
}

cna_hmmcopy  = coerce_123(cna_hmmcopy)
cna_multicna_ae = coerce_123(cna_multicna_ae)
cna_scope    = coerce_123(cna_scope)
cna_sconce    = coerce_123(cna_sconce)
cna_sccnv    = coerce_123(cna_sccnv)
cna_rccae    = coerce_123(cna_rccae)

############################## Macro-F1 ###############################################

# Map cell column name -> truth site
site_from_col <- function(x) {
  if (grepl("_PT", x, fixed=TRUE)) return("PT")
  if (grepl("_LN", x, fixed=TRUE)) return("LN1")
  if (grepl("_ML", x, fixed=TRUE)) return("ML1")
  if (grepl("_MP", x, fixed=TRUE)) return("MP1")
  NA_character_
}

# 3-class per-cell F1s and Macro-F1
macro_f1_3class <- function(pred, truth) {
  ok <- is.finite(pred) & is.finite(truth)
  classes <- c(1L,2L,3L)  # DEL, NOR, DUP
  f1s <- sapply(classes, function(k) {
    TP <- sum(pred[ok]==k & truth[ok]==k)
    FP <- sum(pred[ok]==k & truth[ok]!=k)
    FN <- sum(pred[ok]!=k & truth[ok]==k)
    pr <- ifelse(TP+FP>0, TP/(TP+FP), NA_real_)
    rc <- ifelse(TP+FN>0, TP/(TP+FN), NA_real_)
    ifelse(pr+rc>0, 2*pr*rc/(pr+rc), NA_real_)
  })
  data.table(F1_DEL=f1s[1], F1_NOR=f1s[2], F1_DUP=f1s[3],
             MacroF1=mean(f1s, na.rm=TRUE),
             n_eval=sum(ok))
}

# Truth (bins already correspond across matrices)
truth_mat <- coerce_123(cna_dnacopy)

# Methods
methods <- list(
  hmmcopy  = cna_hmmcopy,
  multicna_ae = cna_multicna_ae,
  scope    = cna_scope,
  sconce    = cna_sconce,
  sccnv    = cna_sccnv,
  rccae = cna_rccae
)

# display labels & plotting order
name_map     <- c(hmmcopy="HMMcopy", multicna_ae="MintCNA",scope="SCOPE", sconce ="SCONCE",sccnv="SCCNV",rccae="rcCAE")
method_order <- c("MintCNA","SCOPE","HMMcopy","SCONCE","SCCNV","rcCAE")
site_order   <- c("PT","LN1","ML1","MP1")
pal          <- setNames(viridis(length(method_order), option="D"), method_order)

# sanity: same nrow across all matrices
stopifnot(all(sapply(methods, nrow) == nrow(truth_mat)))

# ===================== evaluate (no alignment) =====================

# For each method column (cell), compare to the corresponding truth site column
res_macro <- rbindlist(lapply(names(methods), function(mn) {
  pred_mat <- methods[[mn]]
  pred_cols <- colnames(pred_mat)
  pred_site <- vapply(pred_cols, site_from_col, character(1))
  keep <- !is.na(pred_site) & pred_site %in% colnames(truth_mat)
  if (!all(keep)) message(sprintf("[site map] dropped %d/%d columns for %s.",
                                  sum(!keep), length(keep), mn))
  pred_cols <- pred_cols[keep]; pred_site <- pred_site[keep]
  
  rbindlist(lapply(seq_along(pred_cols), function(i) {
    cell <- pred_cols[i]; site <- pred_site[i]
    out <- macro_f1_3class(pred_mat[, cell], truth_mat[, site])
    out[, `:=`(method = name_map[mn], cell_id = cell, site = site, n_bins = nrow(truth_mat))]
  }))
}), use.names = TRUE)

## ---------- prettify site names + add ALL panel ----------
# map LN1→LN, ML1→ML, MP1→MP; keep PT
site_lab <- c("PT"="PT","LN1"="LN","ML1"="ML","MP1"="MP")
res_plot <- copy(res_macro)
res_plot[, facet := unname(site_lab[as.character(site)])]
res_plot[, facet := factor(facet, levels = c("PT","LN","ML","MP"))]

# add an "ALL" panel by duplicating rows with facet == "ALL"
res_all <- copy(res_plot)[, facet := factor("ALL", levels = c("PT","LN","ML","MP","ALL"))]
res_plot <- rbind(res_plot, res_all, fill = TRUE)
res_plot[, facet := factor(facet, levels = c("PT","LN","ML","MP","ALL"))]

# method palette & order (you already had these)
res_plot[, method := factor(method, levels = method_order)]
pal <- setNames(viridis(length(methods), option = "D"), levels(res_plot$method))

# count sample size per facet
facet_n <- res_plot[, .N, by = facet]
lab_names <- setNames(paste0(facet_n$facet, " (N = ", facet_n$N/length(methods), ")"),facet_n$facet)

p_macro <- ggplot(res_plot, aes(x = method, y = MacroF1, fill = method)) +
  geom_boxplot(width = 0.8, color = "black", outlier.size = 0.6) +
  scale_fill_manual(values = pal, breaks = levels(res_plot$method)) +
  facet_wrap(~ facet, nrow = 1,labeller = labeller(facet = lab_names)) +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = NULL, y = "Macro-F1") +
  theme_bw(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92"),
        axis.text.x = element_blank(),            # hide x tick labels
        axis.ticks.x = element_blank(),
        legend.position = "right",
        legend.title = element_blank())

print(p_macro)
ggsave("<MINTCNA_ROOT>/Result/CRC/Segmentation/Evaluation/wgs/macro_f1_ae.png",
       width=10, height=5, dpi=600)

############################## F1 for DUP and DEL ###############################################

suppressPackageStartupMessages({library(data.table); library(ggplot2); library(viridisLite) })

## --------- Binary F1 for one positive class (DEL or DUP) -------------------
# y_pred, y_true: integer vectors in {1,2,3}
# pos = 1 for DEL, 3 for DUP
f1_binary <- function(y_pred, y_true, pos) {
  ok <- is.finite(y_pred) & is.finite(y_true)
  y_pred <- y_pred[ok]; y_true <- y_true[ok]
  if (!length(y_pred)) return(NA_real_)
  tp <- sum(y_pred == pos & y_true == pos)
  fp <- sum(y_pred == pos & y_true != pos)
  fn <- sum(y_pred != pos & y_true == pos)
  denom <- (2*tp + fp + fn)
  if (denom == 0) return(NA_real_)  # undefined if no pos in pred and truth
  2*tp / denom
}

## --------- map column name -> site (matches your earlier code) -------------
site_from_col <- function(x){
  if (grepl("_PT", x, fixed=TRUE)) return("PT")
  if (grepl("_LN", x, fixed=TRUE)) return("LN1")
  if (grepl("_ML", x, fixed=TRUE)) return("ML1")
  if (grepl("_MP", x, fixed=TRUE)) return("MP1")
  NA_character_
}
site_lab <- c("PT"="PT","LN1"="LN","ML1"="ML","MP1"="MP")

## --------- inputs: bin-aligned CNA matrices (values 1/2/3) -----------------
truth_mat <- cna_dnacopy
methods <- list(
  hmmcopy  = cna_hmmcopy,
  multicna_ae = cna_multicna_ae,
  scope    = cna_scope,
  sconce    = cna_sconce,
  sccnv    = cna_sccnv,
  rccae = cna_rccae
)
method_order <- c("MintCNA","SCOPE","HMMcopy","SCONCE","SCCNV","rcCAE")
## --------- compute F1 for DEL & DUP per cell -------------------------------
res_f1 <- rbindlist(lapply(names(methods), function(mn){
  pred <- methods[[mn]]
  stopifnot(nrow(pred) == nrow(truth_mat))
  cols  <- colnames(pred)
  sites <- vapply(cols, site_from_col, character(1))
  keep  <- !is.na(sites) & sites %in% colnames(truth_mat)
  cols  <- cols[keep]; sites <- sites[keep]
  rbindlist(lapply(seq_along(cols), function(i){
    cell <- cols[i]; st <- sites[i]
    data.table(
      method = name_map[mn], site = st, cell = cell,
      DEL = f1_binary(pred[, cell], truth_mat[, st], pos = 1L),
      DUP = f1_binary(pred[, cell], truth_mat[, st], pos = 3L)
    )
  }))
}), use.names = TRUE)

## --------- reshape for plotting & add ALL column ---------------------------
res_plot <- melt(res_f1,
                 id.vars = c("method","site","cell"),
                 variable.name = "facet", value.name = "F1")

# site labels, add ALL pooled column
res_plot[, site := factor(site_lab[site], levels = c("PT","LN","ML","MP"))]
res_all <- copy(res_plot)[, site := factor("ALL", levels = c("PT","LN","ML","MP","ALL"))]
res_plot <- rbind(res_plot, res_all, fill = TRUE)
res_plot[, site := factor(site, levels = c("PT","LN","ML","MP","ALL"))]

# facet rows = DEL/DUP
res_plot[, facet := factor(facet, levels = c("DEL","DUP"))]

# method order + palette (MultiCNA first)
res_plot[, method := factor(method, levels = method_order)]
pal <- setNames(viridis(length(levels(res_plot$method)), option = "D"),
                levels(res_plot$method))


# ---- RECOMPUTE N for this plot ONLY ----
n_methods <- length(levels(res_plot$method))
n_types   <- length(unique(res_plot$facet))   # DEL / DUP => 2

# count regions per site
site_n <- res_plot[, .N, by = site][, N := N / (n_methods * n_types)]

lab_names_site <- setNames(
  paste0(site_n$site, " (N = ", site_n$N, ")"),
  site_n$site
)

## --------- boxplot (rows: DEL/DUP; cols: PT/LN/ML/MP/ALL) ------------------
p_f1 <- ggplot(res_plot, aes(x = method, y = F1, fill = method)) +
  geom_boxplot(width = 0.8, color = "black", outlier.size = 0.6) +
  scale_fill_manual(values = pal, breaks = levels(res_plot$method)) +
  facet_grid(facet ~ site, labeller = labeller(site = lab_names_site)) +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = NULL, y = "F1 score") +
  theme_bw(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank())

print(p_f1)

ggsave("<MINTCNA_ROOT>/Result/CRC/Segmentation/Evaluation/wgs/F1_DEL_DUP_ae.png",
       width=10, height=5, dpi=600)

############################## Kappa ###############################################

suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(viridisLite) })

## --- quadratic weighted kappa (works for K=3) ---
kappa_quadratic <- function(pred, truth, classes = 1:3) {
  ok <- is.finite(pred) & is.finite(truth)
  pred  <- factor(pred[ok],  levels = classes)
  truth <- factor(truth[ok], levels = classes)
  M <- table(truth, pred)
  n <- sum(M); if (n == 0) return(NA_real_)
  K <- nrow(M)
  O <- M / n
  r <- rowSums(O); c <- colSums(O); E <- outer(r, c)
  W <- outer(seq_len(K), seq_len(K), function(i,j) ((i-j)^2)/((K-1)^2))
  denom <- sum(W * E); if (denom == 0) return(NA_real_)
  1 - sum(W * O) / denom
}

## --- helpers to get site from cell id and pretty facet label ---
site_from_col <- function(x){
  if (grepl("_PT", x, fixed=TRUE)) return("PT")
  if (grepl("_LN", x, fixed=TRUE)) return("LN1")
  if (grepl("_ML", x, fixed=TRUE)) return("ML1")
  if (grepl("_MP", x, fixed=TRUE)) return("MP1")
  NA_character_
}
site_lab <- c("PT"="PT","LN1"="LN","ML1"="ML","MP1"="MP")

## --- inputs you already have (bin-wise aligned 1/2/3 matrices) ---
# truth
truth_mat <- cna_dnacopy

methods <- list(
  hmmcopy  = cna_hmmcopy,
  multicna_ae = cna_multicna_ae,
  scope    = cna_scope,
  sconce    = cna_sconce,
  sccnv    = cna_sccnv,
  rccae = cna_rccae
)
method_order <- c("MintCNA","SCOPE","HMMcopy","SCONCE","SCCNV","rcCAE")

## --- compute κ_w per cell for each method ---
res_kappa <- rbindlist(lapply(names(methods), function(mn){
  pred <- methods[[mn]]
  stopifnot(nrow(pred) == nrow(truth_mat))
  # map each column to its site in the truth matrix
  cols  <- colnames(pred)
  sites <- vapply(cols, site_from_col, character(1))
  keep  <- !is.na(sites) & sites %in% colnames(truth_mat)
  cols  <- cols[keep]; sites <- sites[keep]
  rbindlist(lapply(seq_along(cols), function(i){
    cell <- cols[i]; st <- sites[i]
    data.table(
      method = name_map[mn],
      site   = st,
      cell   = cell,
      KappaW = kappa_quadratic(pred[, cell], truth_mat[, st])
    )
  }))
}), use.names = TRUE)

## --- prep for plotting (facets & legend) ---
res_plot <- copy(res_kappa)
res_plot[, facet := factor(site_lab[site], levels = c("PT","LN","ML","MP"))]
# add ALL panel
res_all <- copy(res_plot)[, facet := factor("ALL", levels = c("PT","LN","ML","MP","ALL"))]
res_plot <- rbind(res_plot, res_all, fill = TRUE)
res_plot[, facet := factor(facet, levels = c("PT","LN","ML","MP","ALL"))]

# method order + palette (MultiCNA first)
res_plot[, method := factor(method, levels = method_order)]
pal <- setNames(viridis(length(levels(res_plot$method)), option = "D"),
                levels(res_plot$method))

# count sample size per facet
facet_n <- res_plot[, .N, by = facet]
lab_names <- setNames(paste0(facet_n$facet, " (N = ", facet_n$N/length(methods), ")"),facet_n$facet)

## --- boxplot (same style as your Macro-F1 plot) ---
p_kappa <- ggplot(res_plot, aes(x = method, y = KappaW, fill = method)) +
  geom_boxplot(width = 0.8, color = "black", outlier.size = 0.6) +
  scale_fill_manual(values = pal, breaks = levels(res_plot$method)) +
  facet_wrap(~ facet, nrow = 1, labeller = labeller(facet = lab_names)) +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = NULL, y = "Quadratic weighted Cohen's \u03BA") +
  theme_bw(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank())

print(p_kappa)
ggsave("<MINTCNA_ROOT>/Result/CRC/Segmentation/Evaluation/wgs/kappa.png",
       width=10, height=5, dpi=600)

p_kappa <- ggplot(res_plot[which(res_plot$facet =="ALL"),], aes(x = method, y = KappaW, fill = method)) +
  geom_boxplot(width = 0.8, color = "black", outlier.size = 0.6) +
  scale_fill_manual(values = pal, breaks = levels(res_plot$method)) +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = NULL, y = "Quadratic weighted Cohen's \u03BA") +
  theme_bw(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank(),
        axis.title.y = element_text(size = 14, face = "bold"))

print(p_kappa)
ggsave("<MINTCNA_ROOT>/Result/CRC/Segmentation/Evaluation/wgs/kappa_all.png",
       width= 6, height=6, dpi=600)

 
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(viridisLite) })

## ---------- helpers ----------
# map column name -> site in truth
site_from_col <- function(x){
  if (grepl("_PT", x, fixed=TRUE)) return("PT")
  if (grepl("_LN", x, fixed=TRUE)) return("LN1")
  if (grepl("_ML", x, fixed=TRUE)) return("ML1")
  if (grepl("_MP", x, fixed=TRUE)) return("MP1")
  NA_character_
}
site_lab <- c("PT"="PT","LN1"="LN","ML1"="ML","MP1"="MP")

# signed encoding: 1=DEL -> -1; 2=NOR -> 0; 3=DUP -> +1
signed_encode <- function(x) ifelse(x==1L, -1, ifelse(x==3L, +1, 0))

# cosine similarity on signed vectors
cosine_signed <- function(pred_states, truth_states) {
  ok <- is.finite(pred_states) & is.finite(truth_states)
  if (!any(ok)) return(NA_real_)
  x <- signed_encode(pred_states[ok])
  y <- signed_encode(truth_states[ok])
  den <- sqrt(sum(x^2)) * sqrt(sum(y^2))
  if (den == 0) return(NA_real_)   # no variation -> undefined
  sum(x * y) / den
}

## ---------- inputs (bin-wise aligned, values ∈ {1,2,3}) ----------
truth_mat <- cna_dnacopy
methods <- list(
  hmmcopy  = cna_hmmcopy,
  multicna_ae = cna_multicna_ae,
  scope    = cna_scope,
  sconce    = cna_sconce,
  sccnv    = cna_sccnv,
  rccae = cna_rccae
)
method_order <- c("MintCNA","SCOPE","HMMcopy","SCONCE","SCCNV","rcCAE")


## ---------- compute cosine per cell ----------
res_cos <- rbindlist(lapply(names(methods), function(mn){
  pred <- methods[[mn]]
  stopifnot(nrow(pred) == nrow(truth_mat))
  cols  <- colnames(pred)
  sites <- vapply(cols, site_from_col, character(1))
  keep  <- !is.na(sites) & sites %in% colnames(truth_mat)
  cols  <- cols[keep]; sites <- sites[keep]
  rbindlist(lapply(seq_along(cols), function(i){
    cell <- cols[i]; st <- sites[i]
    data.table(
      method = name_map[mn],
      site   = st,
      cell   = cell,
      Cosine = cosine_signed(pred[, cell], truth_mat[, st])
    )
  }))
}), use.names = TRUE)

## ---------- prep for plotting (PT/LN/ML/MP + ALL) ----------
res_plot <- copy(res_cos)
res_plot[, facet := factor(site_lab[site], levels = c("PT","LN","ML","MP"))]
res_all <- copy(res_plot)[, facet := factor("ALL", levels = c("PT","LN","ML","MP","ALL"))]
res_plot <- rbind(res_plot, res_all, fill = TRUE)
res_plot[, facet := factor(facet, levels = c("PT","LN","ML","MP","ALL"))]

# method order + palette (MultiCNA first), legend instead of x labels
res_plot[, method := factor(method, levels = method_order)]
pal <- setNames(viridis(length(levels(res_plot$method)), option = "D"),
                levels(res_plot$method))

## ---------- boxplot ----------
p_cos <- ggplot(res_plot, aes(x = method, y = Cosine, fill = method)) +
  geom_boxplot(width = 0.8, color = "black", outlier.size = 0.6) +
  scale_fill_manual(values = pal, breaks = levels(res_plot$method)) +
  facet_wrap(~ facet, nrow = 1) +
  coord_cartesian(ylim = c(-.5, 1)) +
  labs(title = "Cosine similarity on signed states (−1/0/+1)",
       x = NULL, y = "Cosine similarity") +
  theme_bw(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank())

print(p_cos)
ggsave("<MINTCNA_ROOT>/Result/CRC/Segmentation/Evaluation/wgs/cosine.png",
       width=10, height=5, dpi=320)

p_cos <- ggplot(res_plot[which(res_plot$facet =="ALL"),], aes(x = method, y = Cosine, fill = method)) +
  geom_boxplot(width = 0.8, color = "black", outlier.size = 0.6) +
  scale_fill_manual(values = pal, breaks = levels(res_plot$method)) +
  coord_cartesian(ylim = c(-.5, 1)) +
  labs(x = NULL, y = "Cosine similarity") +
  theme_bw(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank(),
        axis.title.y = element_text(size = 14, face = "bold"))

print(p_cos)
ggsave("<MINTCNA_ROOT>/Result/CRC/Segmentation/Evaluation/wgs/cosine_all.png",
       width=10, height=5, dpi=320)
