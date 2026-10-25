# Loaders shared by f1_states.R: every method's PER-CELL profile of one sample, projected onto the
# truth's 100 kb grid (each truth bin takes the value of the method bin containing its midpoint) and
# expressed as a log2 ratio to the method's own neutral value.
# Expects the working directory ATAC_RNA (Result/benchmark/..., Result/MintCNA/...).

# index of the method bin containing each truth bin's midpoint (NA when none does)
bin_index <- function(truth, chr, from, to) {
  d <- data.table(chr = as.character(chr), from = as.numeric(from), bin_from = as.numeric(from), bin_to = as.numeric(to), j = seq_along(chr))
  d <- d[is.finite(from) & is.finite(to)]
  setkey(d, chr, from)
  q <- data.table(chr = as.character(truth$chr), from = (as.numeric(truth$start) + as.numeric(truth$end)) / 2, i = seq_len(nrow(truth)))
  setkey(q, chr, from)
  j <- d[q, roll = TRUE]
  j[is.na(bin_to) | from < bin_from | from > bin_to, j := NA_integer_]      # a rolling join must not cross gaps
  out <- rep(NA_integer_, nrow(truth)); out[j$i] <- j$j
  out
}
expand_to_grid <- function(M, idx) {
  out <- matrix(NA_real_, nrow(M), length(idx), dimnames = list(rownames(M), NULL))
  ok <- !is.na(idx); out[, ok] <- M[, idx[ok], drop = FALSE]
  out
}

# chromosome arms (from the AtaCNV bin table), used to place Copy-scAT's arm-level calls
ARMS <- as.data.table(readRDS("Result/benchmark/AtaCNV/P105T1/atacnv_result.rds")$bins)
ARMS <- ARMS[chr %in% paste0("chr", 1:22) & arm %in% c("p", "q"), .(from = min(start), to = max(end)), by = .(chr, arm)]
ARMS[, label := paste0(chr, arm)]

# list of cells x truth-bins matrices (log2 ratio), one per method, for the cells in `keep`
percell_matrices <- function(s, truth, keep) {
  L <- list()

  # inferCNV: relative expression per gene, neutral = 1; genes mapped to the grid by midpoint
  f <- file.path("Result/benchmark/inferCNV", Sys.getenv("INFERCNV_TAG", "default"), s, "infercnv_export.rds")
  if (file.exists(f)) {
    iv <- readRDS(f)
    cc <- intersect(iv$obs_cells, colnames(iv$expr)); cc <- cc[cc %in% keep]
    if (length(cc)) {
      go <- iv$gene_order
      L$inferCNV <- expand_to_grid(log2(pmax(t(iv$expr[, cc, drop = FALSE]), 1e-3)), bin_index(truth, as.character(go$chr), go$start, go$stop))
    }
  }
  # MintCNA, both arms: states 1/2/3 -> log2(state / 2)
  for (a in list(c("MintCNA", "multi", Sys.getenv("MINT_TAG", "h15a005e020")),
                 c("MintCNA(ATAC)", "atac", Sys.getenv("MINT_ATAC_TAG", "h15a005e025")))) {
    f <- sprintf("Result/MintCNA/%s_%s/%s_results_%s.RData", s, a[2], s, a[3])
    if (!file.exists(f)) next
    e <- new.env(); if (inherits(try(load(f, e), silent = TRUE), "try-error")) next
    Sm <- base::get("cna_state_matrix", e); rf <- base::get("ref", e)
    i <- colnames(Sm) %in% keep; if (!any(i)) next
    L[[a[1]]] <- expand_to_grid(log2(t(Sm[, i, drop = FALSE]) / 2), bin_index(truth, rf$chr, rf$start, rf$end))
  }
  # epiAneufinder: states 0/1/2 -> copy number 1/2/3
  f <- file.path("Result/benchmark/epiAneufinder", Sys.getenv("EPI_TAG", "w100kb_mf20000"), s, "epiAneufinder_results/results_table.tsv")
  if (file.exists(f)) {
    ef <- fread(f); cc <- grep("^cell-", names(ef), value = TRUE); use <- cc[sub("^cell-", "", cc) %in% keep]
    if (length(use)) {
      M <- t(as.matrix(ef[, ..use])) + 1; rownames(M) <- sub("^cell-", "", use)
      L$epiAneufinder <- expand_to_grid(log2(M / 2), bin_index(truth, ef$seq, ef$start, ef$end))
    }
  }
  # AtaCNV: copy ratio, neutral = 1
  f <- file.path("Result/benchmark/AtaCNV", s, "atacnv_result.rds")
  if (file.exists(f)) {
    A <- readRDS(f)$cnv$copy_ratio; i <- rownames(A) %in% keep
    if (any(i)) {
      cn <- do.call(rbind, strsplit(colnames(A), "_"))
      L$AtaCNV <- expand_to_grid(log2(pmax(A[i, , drop = FALSE], 1e-3)), bin_index(truth, cn[, 1], as.numeric(cn[, 2]), as.numeric(cn[, 3])))
    }
  }
  # Copy-scAT: arm-level pseudo copy number, neutral = 2; arms absent from the output were called neutral
  f <- file.path("Result/benchmark/CopyscAT", s, paste0(s, "clean_cnv_cnv_scores.csv"))
  if (file.exists(f)) {
    cs <- read.csv(f, row.names = 1, check.names = FALSE); i <- rownames(cs) %in% keep
    if (any(i)) {
      tested <- as.character(readRDS(file.path("Result/benchmark/CopyscAT", s, "copyscat_result.rds"))$collapse[[1]])
      full <- matrix(2, nrow = sum(i), ncol = length(tested), dimnames = list(rownames(cs)[i], tested))
      shared <- intersect(colnames(cs), tested); full[, shared] <- as.matrix(cs[i, shared, drop = FALSE])
      a <- ARMS[label %in% tested]; full <- full[, a$label, drop = FALSE]
      L$`Copy-scAT` <- expand_to_grid(log2(pmax(full, 1e-3) / 2), bin_index(truth, a$chr, a$from, a$to))
    }
  }
  # RIDDLER: multiplier of the cell's baseline, neutral = 1 (barcodes come back with "." for "-")
  f <- file.path("Result/benchmark/RIDDLER/1MB", s, "riddler_result.rds")
  if (file.exists(f)) {
    rd <- readRDS(f)$riddler
    bc <- sub("\\.(\\d+)$", "-\\1", colnames(rd$CNV)); i <- bc %in% keep
    if (any(i)) {
      w <- rd$window; mi <- match(rownames(rd$CNV), as.character(w[[4]]))
      M <- t(rd$CNV[, i, drop = FALSE]); rownames(M) <- bc[i]
      L$RIDDLER <- expand_to_grid(log2(pmax(M, 1e-3)), bin_index(truth, as.character(w[[1]])[mi], w[[2]][mi], w[[3]][mi]))
    }
  }
  # CONGAS+: one profile per clone, assigned to its cells, neutral = 2
  f <- file.path("Result/benchmark/CONGASplus", s, "congas_result.rds")
  if (file.exists(f)) {
    bf <- readRDS(f)$fit$best_fit
    cna <- as.data.table(if (!is.null(bf$CNA_real)) bf$CNA_real else bf$CNA)
    ass <- as.data.table(bf$cluster_assignments); ass[, cell := sub("-(RNA|ATAC)$", "", cell)]
    ass <- unique(ass[, .(cell, cluster)])[cell %in% keep]
    if (nrow(ass)) {
      cna[, c("c", "fr", "to") := tstrsplit(segment_id, ":", fixed = TRUE)]
      segs <- unique(cna[, .(segment_id, c, fr = as.numeric(fr), to = as.numeric(to))])
      wide <- dcast(cna, cluster ~ segment_id, value.var = "value")
      cl <- as.matrix(wide[, -1]); rownames(cl) <- wide$cluster; cl <- cl[, segs$segment_id, drop = FALSE]
      M <- cl[match(ass$cluster, rownames(cl)), , drop = FALSE]; rownames(M) <- ass$cell
      L$`CONGAS+` <- expand_to_grid(log2(pmax(M, 1e-3) / 2), bin_index(truth, segs$c, segs$fr, segs$to))
    }
  }
  L
}
