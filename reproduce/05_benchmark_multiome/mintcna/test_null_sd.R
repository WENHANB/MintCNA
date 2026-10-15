suppressPackageStartupMessages(library(data.table))
RR <- "<MINTCNA_ROOT>/MintCNA_git"
for (f in list.files(file.path(RR,"R"), "\\.[rR]$", recursive=TRUE, full.names=TRUE))
  try(suppressWarnings(sys.source(f, globalenv())), silent=TRUE)
setwd("<MINTCNA_ROOT>/ATAC_RNA")
e <- new.env(); load("Result/MintCNA/P111T3_multi/P111T3_combined.RData", e)
cl <- base::get("combined_list", e); ref <- base::get("ref", e)
qn <- function(m, s){ mu<-apply(m,1,median,na.rm=TRUE); s0<-apply(m,1,sd,na.rm=TRUE)
  s0[!is.finite(s0)|s0==0]<-1; (m-mu)/s0*s+mu }
sdt <- sapply(seq_len(nrow(cl[[1]])), function(r) median(sapply(cl, function(m) sd(m[r,],na.rm=TRUE))))
x <- qn(cl[[length(cl)]], sdt)
cat("ROWSD", paste(round(apply(x,1,sd,na.rm=TRUE),3), collapse=" "), "\n")

# deparse() rewraps the call, so the literal "rep(sd_row, each = n)" does not
# survive it. Patch the SOURCE FILE text instead, which is what actually runs.
raw <- readLines(file.path(RR, "R/segmentation/null.R"))
stopifnot(sum(grepl("each = n", raw)) == 1)
raw2 <- sub("each = n", "times = n", raw, fixed = TRUE)
tmp <- tempfile(fileext = ".R"); writeLines(raw2, tmp)
env2 <- new.env(parent = globalenv()); sys.source(tmp, env2)
mi_fix <- base::get("multiInverse_chr", env2)

for (lab in c("shipped","fixed")) {
  set.seed(1)
  q <- if (lab=="shipped") multiInverse_chr(x, ref, h=15, simT=10, penalty="l2")
       else                 mi_fix(x, ref, h=15, simT=10, penalty="l2")
  t10 <- q$q[which.min(abs(q$probs-0.90))]; t05 <- q$q[which.min(abs(q$probs-0.95))]
  cat(sprintf("TAU %-8s alpha=0.10 -> %.4f   alpha=0.05 -> %.4f\n", lab, t10, t05))
}
