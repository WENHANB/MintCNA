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

# Truth (bins already correspond across matrices)
truth_mat <- coerce_123(cna_dnacopy)
cna_hmmcopy  = coerce_123(cna_hmmcopy)
cna_multicna_ae = coerce_123(cna_multicna_ae)
cna_scope    = coerce_123(cna_scope)
cna_sconce    = coerce_123(cna_sconce)
cna_sccnv    = coerce_123(cna_sccnv)
cna_rccae    = coerce_123(cna_rccae)

############################ F1-score ############################ 
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
site_lab <- c("PT"="PT","LN1"="LN","ML1"="ML","MP1"="MP")
pal          <- setNames(viridis(length(method_order), option="D"), method_order)

# F1
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


