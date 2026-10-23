#!/usr/bin/env Rscript
# Does a genome-direction smoother add anything on top of AG-DCAE?
#
# AG-DCAE shares information ACROSS CELLS (a latent representation of the
# bins x cells matrix). AtaCNV's DLM shares information ALONG THE GENOME within
# one cell. Measured on P111T3 they are not the same operation at all: AG-DCAE
# leaves roughness at 0.946 (raw 0.948) while lifting per-cell r 0.208 -> 0.305.
# This script asks whether adding the second one helps, and at what cost to the
# breakpoints -- a smoother that blurs edges would show up as a LOWER amplitude
# at truly altered bins.
suppressPackageStartupMessages({library(data.table); library(matrixStats); library(dlm)})
setwd("<MINTCNA_ROOT>/ATAC_RNA")
args <- commandArgs(trailingOnly = TRUE); S <- if (length(args)) args[1] else "P111T3"
CHRS <- paste0("chr",1:22); G <- 1e6

tr <- fread(sprintf("Result/truth/%s_truth_100kb.tsv.gz", S))[chr %in% CHRS]
tr <- tr[, .(t = mean(log2, na.rm=TRUE)), by = .(chr, gb = start %/% G)]; setkey(tr, chr, gb)
ae  <- as.matrix(fread(sprintf("Result/MintCNA/%s_multi/%s_ae_dna.csv", S, S)), rownames=1)
ref <- fread(sprintf("Result/MintCNA/%s_multi/%s_ref_dna.csv", S, S))

# DLM local-level smoother, applied WITHIN each chromosome so it cannot smooth
# across a centromere or a chromosome boundary. AtaCNV's own dV/dW are the
# starting point; the ratio dV/dW sets how hard it smooths.
dlm_chr <- function(m, chr, dV, dW) {
  mod <- dlmModPoly(order = 1, dV = dV, dW = dW)
  idx <- split(seq_len(nrow(m)), chr)
  out <- m
  for (i in idx) for (j in seq_len(ncol(m))) {
    y <- m[i, j]; if (all(is.na(y))) next
    out[i, j] <- dlmSmooth(y, mod)$s[-1]
  }
  out
}
qn_modal <- function(m, sd_target) {
  mu <- rowMedians(m, na.rm=TRUE); s <- rowSds(m, na.rm=TRUE); s[!is.finite(s)|s==0] <- 1
  (m - mu)/s * sd_target + mu
}
agg <- function(M) {
  k <- paste(ref$chr, ref$start %/% G); i <- split(seq_len(nrow(M)), k)
  list(m = t(vapply(i, function(r) colMeans2(M[r,,drop=FALSE], na.rm=TRUE), numeric(ncol(M)))),
       key = do.call(rbind, strsplit(names(i), " ")))
}
stat <- function(M, nm) {
  a <- agg(M); d <- data.table(chr=a$key[,1], gb=as.numeric(a$key[,2]))
  ok <- d$chr %in% CHRS; Mm <- a$m[ok,,drop=FALSE]; d <- d[ok]
  j <- tr[d, on=.(chr,gb)]$t; keep <- is.finite(j); Mm <- Mm[keep,,drop=FALSE]; j <- j[keep]
  pb <- rowMeans2(Mm, na.rm=TRUE); alt <- abs(j) > 0.1
  pc <- suppressWarnings(apply(Mm[, sample(ncol(Mm), min(300,ncol(Mm)))], 2,
                               function(v) cor(v, j, use="complete.obs")))
  cat(sprintf("%-26s pb r %6.3f | per-cell r %6.3f | amp@alt %6.3f  sd@neut %6.3f  SNR %5.2f | rough %5.3f\n",
      nm, cor(pb,j,use="complete.obs"), mean(pc,na.rm=TRUE),
      mean(abs(pb[alt]),na.rm=TRUE), sd(pb[!alt],na.rm=TRUE),
      mean(abs(pb[alt]),na.rm=TRUE)/sd(pb[!alt],na.rm=TRUE),
      mean(abs(diff(pb)),na.rm=TRUE)/sd(pb,na.rm=TRUE)))
}
cat("===", S, "— does a genome-direction smoother add to AG-DCAE? ===\n\n")
sdt <- median(rowSds(ae, na.rm=TRUE), na.rm=TRUE)
stat(qn_modal(ae, sdt), "AG-DCAE + QN (current)")
for (p in list(c(0.3,0.01), c(0.3,0.03), c(0.3,0.003), c(1.0,0.01))) {
  t0 <- Sys.time()
  sm <- dlm_chr(ae, ref$chr, p[1], p[2])
  stat(qn_modal(sm, median(rowSds(sm, na.rm=TRUE), na.rm=TRUE)),
       sprintf("+ DLM dV=%.1f dW=%.3f", p[1], p[2]))
  cat(sprintf("%26s (%.1f min)\n", "", as.numeric(difftime(Sys.time(), t0, units="mins"))))
}
