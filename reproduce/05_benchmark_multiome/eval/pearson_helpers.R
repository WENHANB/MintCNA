# Shared loaders: return each method's PER-CELL profile projected onto the
# truth's 100 kb grid, in log2 ratio, for the cells requested.
#
# Everything is projected onto the common grid rather than compared at each
# method's native resolution. That choice matters for the arm-level callers: at
# arm resolution every arm counts once, so chr21q weighs as much as chr1q, while
# on the 100 kb grid each method bin is weighted by the genome it actually
# covers. The second is the fairer weighting and is what the field does.
#
# The mapping truth-bin -> method-bin is computed once as an index vector and
# then applied by column indexing; building it per cell would be ~1000x slower.

# index of the method bin containing each truth bin's midpoint (NA if none)
bin_index <- function(truth, chr, from, to) {
  d <- data.table(chr = as.character(chr), from = as.numeric(from),
                  bin_from = as.numeric(from), bin_to = as.numeric(to),
                  j = seq_along(chr))
  d <- d[is.finite(from) & is.finite(to)]
  setkey(d, chr, from)
  q <- data.table(chr = as.character(truth$chr),
                  from = (as.numeric(truth$start) + as.numeric(truth$end)) / 2,
                  i = seq_len(nrow(truth)))
  setkey(q, chr, from)
  j <- d[q, roll = TRUE]
  # a plain roll carries the last bin of a chromosome across centromere gaps and
  # past the chromosome end, so require real containment
  j[is.na(bin_to) | from < bin_from | from > bin_to, j := NA_integer_]
  out <- rep(NA_integer_, nrow(truth))
  out[j$i] <- j$j
  out
}

# expand a cells x methodbins matrix onto the truth grid using that index
expand_to_grid <- function(M, idx) {
  out <- matrix(NA_real_, nrow(M), length(idx),
                dimnames = list(rownames(M), NULL))
  ok <- !is.na(idx)
  out[, ok] <- M[, idx[ok], drop = FALSE]
  out
}

percell_matrices <- function(s, truth, keep) {
  L <- list()

  # --- inferCNV: gene-level relative expression, neutral = 1 ------------------
  # It was missing from this helper entirely, so every per-cell table built from
  # it silently compared six methods while the pseudobulk table compared seven.
  # Genes are mapped to the truth grid by midpoint, exactly as the coarser
  # methods' bins are; several genes can fall in one bin, and the last one wins
  # rather than being averaged, because averaging a ratio across genes is not the
  # same operation as inferCNV's own smoothing and would invent a value the
  # method never produced.
  f <- file.path("Result/benchmark/inferCNV", Sys.getenv("INFERCNV_TAG", "default"), s,
                 "infercnv_export.rds")
  if (file.exists(f)) {
    iv <- readRDS(f)
    cc <- intersect(iv$obs_cells, colnames(iv$expr))
    cc <- cc[cc %in% keep]
    if (length(cc)) {
      go  <- iv$gene_order
      idx <- bin_index(truth, as.character(go$chr), go$start, go$stop)
      L$inferCNV <- expand_to_grid(log2(pmax(t(iv$expr[, cc, drop = FALSE]), 1e-3)), idx)
    }
  }

  # --- MintCNA, both arms: discrete states 1/2/3 -> log2(state / 2) ----------
  # Same convention the pseudobulk table uses. The two arms differ only in mode
  # and in eps (0.20 multi, 0.25 ATAC-only, the value at which the two call the
  # same share of the genome); everything else, including the AG-DCAE output, is
  # shared.
  for (a in list(c("MintCNA", "multi", Sys.getenv("MINT_TAG", "h15a005e020")),
                 c("MintCNA(ATAC)", "atac", Sys.getenv("MINT_ATAC_TAG", "h15a005e025")))) {
    f <- sprintf("Result/MintCNA/%s_%s/%s_results_%s.RData", s, a[2], s, a[3])
    if (!file.exists(f)) next
    e <- new.env(); if (inherits(try(load(f, e), silent = TRUE), "try-error")) next
    Sm <- base::get("cna_state_matrix", e); rf <- base::get("ref", e)
    i <- colnames(Sm) %in% keep
    if (!any(i)) next
    M <- t(Sm[, i, drop = FALSE])                       # cells x bins
    idx <- bin_index(truth, rf$chr, rf$start, rf$end)
    L[[a[1]]] <- expand_to_grid(log2(M / 2), idx)
  }

  # --- epiAneufinder: states 0/1/2 -> copy number 1/2/3 ---------------------
  # As published (windowSize 1e5, minFrags 20000), matching
  # Code/eval/pearson_vs_truth.R. This file still pointed at the tuned
  # w1mb_mf5000 run after the 2026-08-18 revert to package defaults, so the
  # per-cell table was scoring a configuration the pseudobulk table had dropped.
  f <- file.path("Result/benchmark/epiAneufinder",
                 Sys.getenv("EPI_TAG", "w100kb_mf20000"), s,
                 "epiAneufinder_results/results_table.tsv")
  if (file.exists(f)) {
    ef <- fread(f); cc <- grep("^cell-", names(ef), value = TRUE)
    use <- cc[sub("^cell-", "", cc) %in% keep]
    if (length(use)) {
      M <- t(as.matrix(ef[, ..use])) + 1                 # cells x bins, CN 1/2/3
      rownames(M) <- sub("^cell-", "", use)              # align barcodes across methods
      idx <- bin_index(truth, ef$seq, ef$start, ef$end)
      L$epiAneufinder <- expand_to_grid(log2(M / 2), idx)
    }
  }

  # --- AtaCNV: linear copy ratio, neutral = 1 -------------------------------
  f <- file.path("Result/benchmark/AtaCNV", s, "atacnv_result.rds")
  if (file.exists(f)) {
    A <- readRDS(f)$cnv$copy_ratio
    i <- rownames(A) %in% keep
    if (any(i)) {
      cn <- do.call(rbind, strsplit(colnames(A), "_"))
      idx <- bin_index(truth, cn[, 1], as.numeric(cn[, 2]), as.numeric(cn[, 3]))
      L$AtaCNV <- expand_to_grid(log2(pmax(A[i, , drop = FALSE], 1e-3)), idx)
    }
  }

  # --- Copy-scAT: arm-level pseudo copy number, neutral = 2 -----------------
  f <- file.path("Result/benchmark/CopyscAT", s, paste0(s, "clean_cnv_cnv_scores.csv"))
  if (file.exists(f)) {
    cs <- read.csv(f, row.names = 1, check.names = FALSE)
    i <- rownames(cs) %in% keep
    if (any(i)) {
      # arms absent from the output were tested and called neutral -> 2
      tested <- as.character(readRDS(file.path("Result/benchmark/CopyscAT", s,
                                               "copyscat_result.rds"))$collapse[[1]])
      full <- matrix(2, nrow = sum(i), ncol = length(tested),
                     dimnames = list(rownames(cs)[i], tested))
      shared <- intersect(colnames(cs), tested)
      full[, shared] <- as.matrix(cs[i, shared, drop = FALSE])
      a <- ARMS[label %in% tested]
      full <- full[, a$label, drop = FALSE]
      idx <- bin_index(truth, a$chr, a$from, a$to)
      L$`Copy-scAT` <- expand_to_grid(log2(pmax(full, 1e-3) / 2), idx)
    }
  }

  # --- RIDDLER: multiplier of the cell's own baseline, neutral = 1 ----------
  f <- file.path("Result/benchmark/RIDDLER/1MB", s, "riddler_result.rds")
  if (file.exists(f)) {
    rd <- readRDS(f)$riddler
    # riddler_fit reads with read.table(header = TRUE), so "AAACAGCC-1" comes
    # back as "AAACAGCC.1"; matching without undoing that selects nothing
    bc <- sub("\\.(\\d+)$", "-\\1", colnames(rd$CNV))
    i <- bc %in% keep
    if (any(i)) {
      w  <- rd$window
      mi <- match(rownames(rd$CNV), as.character(w[[4]]))   # join by name, not position
      idx <- bin_index(truth, as.character(w[[1]])[mi], w[[2]][mi], w[[3]][mi])
      M <- t(rd$CNV[, i, drop = FALSE]); rownames(M) <- bc[i]
      L$RIDDLER <- expand_to_grid(log2(pmax(M, 1e-3)), idx)
    }
  }

  # --- CONGAS+: one profile per CLONE, painted onto its cells ---------------
  f <- file.path("Result/benchmark/CONGASplus", s, "congas_result.rds")
  if (file.exists(f)) {
    bf  <- readRDS(f)$fit$best_fit
    cna <- as.data.table(if (!is.null(bf$CNA_real)) bf$CNA_real else bf$CNA)
    ass <- as.data.table(bf$cluster_assignments)
    ass[, cell := sub("-(RNA|ATAC)$", "", cell)]
    ass <- unique(ass[, .(cell, cluster)])[cell %in% keep]
    if (nrow(ass)) {
      cna[, c("c", "fr", "to") := tstrsplit(segment_id, ":", fixed = TRUE)]
      segs <- unique(cna[, .(segment_id, c, fr = as.numeric(fr), to = as.numeric(to))])
      wide <- dcast(cna, cluster ~ segment_id, value.var = "value")
      cl <- as.matrix(wide[, -1]); rownames(cl) <- wide$cluster
      cl <- cl[, segs$segment_id, drop = FALSE]
      idx <- bin_index(truth, segs$c, segs$fr, segs$to)
      M <- cl[match(ass$cluster, rownames(cl)), , drop = FALSE]
      rownames(M) <- ass$cell
      L$`CONGAS+` <- expand_to_grid(log2(pmax(M, 1e-3) / 2), idx)
    }
  }

  # --- inferCNV: residual expression, neutral = 1 ---------------------------
  f <- file.path("Result/benchmark/inferCNV", s, "infercnv_export.rds")
  if (file.exists(f)) {
    iv <- readRDS(f)
    i <- colnames(iv$expr) %in% keep
    if (any(i)) {
      go <- iv$gene_order
      idx <- bin_index(truth, as.character(go$chr), go$start, go$stop)
      L$inferCNV <- expand_to_grid(log2(pmax(t(iv$expr[, i, drop = FALSE]), 1e-3)), idx)
    }
  }
  L
}

ARMS <- as.data.table(readRDS("Result/benchmark/AtaCNV/P105T1/atacnv_result.rds")$bins)
ARMS <- ARMS[chr %in% paste0("chr", 1:22) & arm %in% c("p", "q"),
             .(from = min(start), to = max(end)), by = .(chr, arm)]
ARMS[, label := paste0(chr, arm)]
