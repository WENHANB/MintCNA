#!/usr/bin/env Rscript
# Numerical equivalence: run the ORIGINAL (Code/MultiCNA, what the simulation and
# CRC01 sections used) and the PACKAGE (MintCNA_git) on identical inputs and
# compare outputs.
#
# Source-level diffing is unreliable in both directions -- it flags `return(x)`
# -> `x` as a difference, and it would miss a change hidden behind identical
# formatting. Only the numbers settle it.
#
# Written after `multiMissInspect_scores` was found to have `n <- 2*h` instead of
# `n <- ncol(x)`, which silently reduced the scan to two positions per chromosome.

set.seed(7)
O <- new.env(); P <- new.env()
for (f in list.files("<MINTCNA_ROOT>/Code/MultiCNA",
                     "\\.[rR]$", full.names = TRUE))
  try(suppressWarnings(suppressMessages(sys.source(f, envir = O))), silent = TRUE)
for (f in list.files("<MINTCNA_ROOT>/MintCNA_git/R",
                     "\\.[rR]$", recursive = TRUE, full.names = TRUE))
  try(suppressWarnings(suppressMessages(sys.source(f, envir = P))), silent = TRUE)

pass <- fail <- skip <- character(0)
chk <- function(name, expr_o, expr_p, tol = 1e-8) {
  a <- tryCatch(eval(expr_o, envir = list()), error = function(e) e)
  b <- tryCatch(eval(expr_p, envir = list()), error = function(e) e)
  if (inherits(a, "error") || inherits(b, "error")) {
    skip <<- c(skip, name)
    cat(sprintf("  %-26s SKIP  (%s)\n", name,
                conditionMessage(if (inherits(a, "error")) a else b)))
    return(invisible())
  }
  ok <- isTRUE(all.equal(a, b, tolerance = tol))
  if (ok) pass <<- c(pass, name) else fail <<- c(fail, name)
  cat(sprintf("  %-26s %s\n", name, if (ok) "identical" else "*** DIFFERS ***"))
}

# ---- fixtures ---------------------------------------------------------------
p <- 2; nb <- 300
X  <- rbind(RNA = c(rnorm(150), rnorm(150, -1)),
            DNA = c(rnorm(150), rnorm(150, -1)))
cus <- O$cusum(X[, 1:11])
cps <- c(60, 150, 220)
ref <- data.frame(chr = rep(paste0("chr", 1:3), each = 100),
                  start = rep(seq(0, 99) * 1e5, 3),
                  end   = rep(seq(1, 100) * 1e5, 3), stringsAsFactors = FALSE)

cat("=== scan / projection ===\n")
chk("cusum",              quote(O$cusum(X[, 1:11])),          quote(P$cusum(X[, 1:11])))
chk("cusum.vector",       quote(O$cusum.vector(X[1, 1:11])),  quote(P$cusum.vector(X[1, 1:11])))
chk("get.v.proj_l2",      quote(O$get.v.proj_l2(cus, 5)),     quote(P$get.v.proj_l2(cus, 5)))
chk("get.v.proj_l1",      quote(O$get.v.proj_l1(cus, 5)),     quote(P$get.v.proj_l1(cus, 5)))
for (h in c(5, 10, 15))
  chk(sprintf("multiMissInspect h=%d", h),
      bquote(O$multiMissInspect_scores(X, h = .(h), penalty = "l2")),
      bquote(P$multiMissInspect_scores(X, h = .(h), penalty = "l2")))

cat("\n=== null distribution ===\n")
for (h in c(5, 10)) {
  chk(sprintf("multiInverse_chr h=%d", h),
      bquote({set.seed(1); O$multiInverse_chr(X, ref = ref, h = .(h), simT = 1)}),
      bquote({set.seed(1); P$multiInverse_chr(X, ref = ref, h = .(h), simT = 1)}))
}

cat("\n=== fusion / calling ===\n")
chk("fusion_omics $z",
    quote(O$fusion_omics(X, cps, method = "count")$z),
    quote(P$fusion_omics(X, cps)$z))
chk("fusion_omics $sum.sq",
    quote(O$fusion_omics(X, cps, method = "count")$sum.sq),
    quote(P$fusion_omics(X, cps)$sum.sq))
for (st in c(3, 5)) {
  chk(sprintf("CNVcluster st=%d $state", st),
      bquote({set.seed(3); O$CNVcluster(Y = X, cp = cps, L = 30, method = "count", st = .(st))$cnv.state}),
      bquote({set.seed(3); P$CNVcluster(Y = X, cp = cps, L = 30, st = .(st))$cnv.state}))
  chk(sprintf("CNVcluster st=%d $iter", st),
      bquote({set.seed(3); O$CNVcluster(Y = X, cp = cps, L = 30, method = "count", st = .(st))$iter}),
      bquote({set.seed(3); P$CNVcluster(Y = X, cp = cps, L = 30, st = .(st))$iter}))
}

cat("\n=== transforms / output ===\n")
chk("qn_modal",      quote(O$qn_modal(X, sd_target = c(1, 1))), quote(P$qn_modal(X, sd_target = c(1, 1))))
chk("reindex_cps",   quote(O$reindex_cps(list(chr1 = c(20, 60), chr2 = c(30), chr3 = integer(0)), c(100, 100, 100))),
                     quote(P$reindex_cps(list(chr1 = c(20, 60), chr2 = c(30), chr3 = integer(0)), c(100, 100, 100))))

cat(sprintf("\n\nPASS %d   DIFFER %d   SKIP %d\n", length(pass), length(fail), length(skip)))
if (length(fail)) cat("DIFFERING:", paste(fail, collapse = ", "), "\n")
if (length(skip)) cat("SKIPPED:  ", paste(skip, collapse = ", "), "\n")
