#!/usr/bin/env Rscript
# Compare every function MintCNA_git ships against the original in Code/MultiCNA,
# which is what the simulation and CRC01 sections actually ran.
#
# The refactor into a package was not verified function-by-function, and
# multiMissInspect_scores already turned out to differ in a way that changes the
# result, not just the speed. This finds the rest before anything is run.

ORIG <- "<MINTCNA_ROOT>/Code/MultiCNA"
PKG  <- "<MINTCNA_ROOT>/MintCNA_git/R"

grab <- function(dir) {
  e <- new.env()
  for (f in list.files(dir, "\\.[rR]$", recursive = TRUE, full.names = TRUE)) {
    try(suppressWarnings(suppressMessages(sys.source(f, envir = e))), silent = TRUE)
  }
  e
}
eo <- grab(ORIG); ep <- grab(PKG)
fo <- ls(eo)[vapply(ls(eo), function(n) is.function(get(n, eo)), logical(1))]
fp <- ls(ep)[vapply(ls(ep), function(n) is.function(get(n, ep)), logical(1))]

norm <- function(f) paste(trimws(deparse(f)), collapse = "\n")
both <- intersect(fo, fp)
same <- diff <- character(0)
for (n in both) {
  if (identical(norm(get(n, eo)), norm(get(n, ep)))) same <- c(same, n) else diff <- c(diff, n)
}

cat(sprintf("functions in original: %d   in package: %d   shared: %d\n\n",
            length(fo), length(fp), length(both)))
cat(sprintf("IDENTICAL (%d): %s\n\n", length(same), paste(same, collapse = ", ")))
cat(sprintf("DIFFERENT (%d):\n", length(diff)))
for (n in diff) {
  a <- trimws(deparse(get(n, eo))); b <- trimws(deparse(get(n, ep)))
  # crude line-level diff: lines present in one and not the other
  only_o <- setdiff(a, b); only_p <- setdiff(b, a)
  cat(sprintf("\n=== %s ===  (original %d lines, package %d)\n", n, length(a), length(b)))
  if (length(only_o)) { cat("  ONLY IN ORIGINAL:\n"); cat(paste0("    ", head(only_o, 8)), sep = "\n") }
  if (length(only_p)) { cat("  ONLY IN PACKAGE:\n");  cat(paste0("    ", head(only_p, 8)), sep = "\n") }
}
cat(sprintf("\n\nIN ORIGINAL BUT NOT PACKAGED (%d): %s\n",
            length(setdiff(fo, fp)), paste(setdiff(fo, fp), collapse = ", ")))
cat(sprintf("\nNEW IN PACKAGE (%d): %s\n",
            length(setdiff(fp, fo)), paste(setdiff(fp, fo), collapse = ", ")))
