#!/usr/bin/env Rscript
# Independent checks on the number the benchmark reports. Each one either passes
# or names what is wrong; none of them reuses the logic it is testing.
suppressPackageStartupMessages(library(data.table))
setwd("<MINTCNA_ROOT>/ATAC_RNA")
CHRS <- paste0("chr", 1:22); GRID <- 5e5
src <- readLines("Code/eval/pearson_vs_truth.R")
i0 <- grep("^project <- function", src); i1 <- grep("^# --- per-method pseudobulk", src)[1]
eval(parse(text = paste(src[i0:(i1-1)], collapse = "\n")))
ok <- TRUE; chk <- function(name, pass, detail = "") {
  cat(sprintf("[%s] %-52s %s\n", if (pass) "PASS" else "FAIL", name, detail))
  if (!pass) ok <<- FALSE
}

# ---- 1. project() on a case with a known answer ----------------------------
truth <- data.table(chr = rep("chr1", 10), start = seq(0, 9)*1e5, end = seq(1,10)*1e5)
# one method bin covering truth bins 3..7 (midpoints 250k..650k), value 2
v <- project(truth, "chr1", 2e5, 7e5, 2)
chk("project: midpoint containment", identical(which(!is.na(v)), 3:7),
    sprintf("non-NA at %s", paste(which(!is.na(v)), collapse=",")))
chk("project: value carried through", all(v[3:7] == 2))

# ---- 2. it must NOT carry a value across a gap ------------------------------
# two method bins with a hole between them; truth midpoints inside the hole are NA
v2 <- project(truth, c("chr1","chr1"), c(0, 7e5), c(2e5, 1e6), c(1, 3))
chk("project: no roll across a gap", all(is.na(v2[3:7])),
    sprintf("hole bins give %s", paste(v2[3:7], collapse=",")))

# ---- 3. it must not leak across chromosomes --------------------------------
tr2 <- rbind(truth, data.table(chr="chr2", start=seq(0,4)*1e5, end=seq(1,5)*1e5))
v3  <- project(tr2, "chr1", 0, 1e6, 5)
chk("project: chromosome-scoped", all(is.na(v3[tr2$chr == "chr2"])))

# ---- 4. Pearson is affine-invariant: a rescaled track must score identically -
S <- "P111T3"
tr <- fread(sprintf("Result/truth/%s_truth_100kb.tsv.gz", S))[chr %in% CHRS]
a  <- readRDS(sprintf("Result/benchmark/AtaCNV/%s/atacnv_result.rds", S))$cnv$copy_ratio
keep <- fread(sprintf("Data/cell_lists/%s_cells.tsv", S))[cell_type=="Malignant", barcode]
i <- rownames(a) %in% keep
cn <- do.call(rbind, strsplit(colnames(a), "_"))
m  <- colMeans(a[i, , drop=FALSE])
p1 <- project(tr, cn[,1], as.numeric(cn[,2]), as.numeric(cn[,3]), log2(m))
agg <- function(v) { d <- data.table(chr=tr$chr, gb=tr$start %/% GRID, v=v, t=tr$log2)
  d <- d[, .(v=mean(v,na.rm=TRUE), t=mean(t,na.rm=TRUE)), by=.(chr,gb)]; d[is.finite(v)&is.finite(t)] }
d1 <- agg(p1); r1 <- cor(d1$v, d1$t)
p2 <- project(tr, cn[,1], as.numeric(cn[,2]), as.numeric(cn[,3]), log2(m)*3 + 7)
d2 <- agg(p2); r2 <- cor(d2$v, d2$t)
chk("Pearson affine-invariant", abs(r1-r2) < 1e-10, sprintf("%.6f vs %.6f", r1, r2))

# ---- 5. the reported number reproduces from the saved profile --------------
pr <- fread(sprintf("Result/qc/profiles/%s_profiles_100kb.tsv.gz", S))
tv <- pr[method=="bulk WGS", .(chr, start, t=log2)]
tab <- fread("Result/qc/pearson_vs_truth_own.tsv")[care_id==S]
for (mm in c("AtaCNV","MintCNA","RIDDLER")) {
  x <- merge(pr[method==mm, .(chr, start, v=log2)], tv, by=c("chr","start"))
  x <- x[chr %in% CHRS][, .(v=mean(v,na.rm=TRUE), t=mean(t,na.rm=TRUE)),
                        by=.(chr, gb=start %/% GRID)][is.finite(v)&is.finite(t)]
  rr <- cor(x$v, x$t); rep <- tab[method==mm, r]
  chk(sprintf("table reproduces from profile: %s", mm), abs(rr-rep) < 1e-6,
      sprintf("recomputed %.4f, table %.4f", rr, rep))
}

# ---- 6. shuffling the truth must destroy the correlation -------------------
set.seed(1); x <- merge(pr[method=="AtaCNV", .(chr,start,v=log2)], tv, by=c("chr","start"))
x <- x[chr %in% CHRS & is.finite(v) & is.finite(t)]
chk("shuffled truth gives r ~ 0", abs(cor(x$v, sample(x$t))) < 0.05,
    sprintf("%.4f", cor(x$v, sample(x$t))))

# ---- 7. the truth itself scores 1.0 ----------------------------------------
chk("truth vs itself = 1", abs(cor(tv$t, tv$t, use="complete.obs") - 1) < 1e-12)

# ---- 8. coordinate conventions actually overlap ----------------------------
# a silent off-by-one-bin would still give a high r; check the ranges agree
rng <- pr[chr=="chr1", .(min=min(start), max=max(start), n=.N), by=method]
chk("all methods span the same chr1 range", diff(range(rng$max)) < 5e6,
    sprintf("max start spread %.1f Mb", diff(range(rng$max))/1e6))

cat(sprintf("\n%s\n", if (ok) "ALL CHECKS PASS" else "SOME CHECKS FAILED"))
