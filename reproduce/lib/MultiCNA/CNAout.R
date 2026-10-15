build_cna_matrix <- function(cna_df_all, ref,
                             code_map = c(del = 1, nor = 2, dup = 3),
                             default_state = "nor") {
  n_bins <- nrow(ref)
  cells  <- unique(cna_df_all$cell)
  M <- matrix(code_map[default_state], n_bins, length(cells),
              dimnames = list(sprintf("%s_%d_%d", ref$chr, ref$start, ref$end), cells))
  
  for (cell in cells) {
    segs <- cna_df_all[cna_df_all$cell == cell, ]
    for (i in seq_len(nrow(segs))) {
      M[segs$start_idx[i]:segs$end_idx[i], cell] <- code_map[segs$state[i]]
    }
  }
  M
}

library(dplyr); library(tidyr); library(ggplot2)
library(ggnewscale); library(viridisLite)

make_genome_map <- function(ref){
  ref2 <- ref %>% 
    mutate(chr = as.character(chr),
           chr_num = suppressWarnings(as.integer(sub("^chr","", chr)))) %>%
    arrange(is.na(chr_num), chr_num, chr, start) %>%
    mutate(bin_idx = row_number())
  chr_last  <- tapply(ref2$bin_idx, ref2$chr, max)
  chr_first <- tapply(ref2$bin_idx, ref2$chr, min)
  chr_mids  <- floor((chr_first + chr_last)/2)
  list(ref = ref2,
       chr_last = chr_last,
       chr_mids = chr_mids,
       chr_labels = sub("^chr","", names(chr_last)))
}

# S: bins×cells integer codes (from your build_cna_matrix), ref with chr/start/end
# code_map: named integer vector, e.g.
#   3-state: c(del=1, nor=2, dup=3)
#   5-state: c(del=1, loss=2, nor=3, gain=4, amp=5)
plot_S_heatmap <- function(S, ref, code_map,
                           file = "CNA_heatmap_from_S.pdf",
                           row_order = NULL, width = 14, height = 4,
                           strip_bins = NULL) {
  stopifnot(is.matrix(S), nrow(S) == nrow(ref))
  gm <- make_genome_map(ref)
  
  # reorder cells if provided
  cells <- colnames(S)
  if (!is.null(row_order)) cells <- intersect(row_order, cells)
  S <- S[, cells, drop = FALSE]
  
  # integer -> state label, using code_map (no collapsing)
  lab_levels <- names(code_map)
  lab_mat <- matrix(lab_levels[match(as.integer(S), code_map)],
                    nrow = nrow(S), dimnames = dimnames(S))
  # fill NAs as "nor" if present, else first level
  fallback <- if ("nor" %in% lab_levels) "nor" else lab_levels[1]
  lab_mat[is.na(lab_mat)] <- fallback
  
  # compress runs per cell to rectangles
  rle_rects <- function(v){
    r <- rle(v); ends <- cumsum(r$lengths); starts <- ends - r$lengths + 1
    data.frame(state = vapply(r$values, toupper, ""), start_idx = starts, end_idx = ends)
  }
  dd <- lapply(seq_len(ncol(lab_mat)), function(j){
    x <- rle_rects(lab_mat[, j]); x$cell <- colnames(lab_mat)[j]; x
  }) %>% bind_rows()
  
  # map bin indices → x
  n_bins <- nrow(gm$ref); x_map <- gm$ref$bin_idx
  dd <- dd %>%
    mutate(xmin = x_map[pmax(1L, start_idx)] - 0.5,
           xmax = x_map[pmin(n_bins, end_idx)] + 0.5,
           cell_f = factor(cell, levels = rev(cells)))
  
  # left strip: site from cell ID (LN/ML/MP/PT)
  if (is.null(strip_bins)) strip_bins <- max(50, round(0.02 * n_bins))
  site_df <- data.frame(cell = cells) %>%
    mutate(site = sub("^.*?_(.*?)_.*$", "\\1", cell),
           site = gsub("[0-9]+$", "", site),
           cell_f = factor(cell, levels = rev(cells)),
           ymin = as.numeric(cell_f) - 0.5, ymax = as.numeric(cell_f) + 0.5,
           xmin = 0.5, xmax = strip_bins + 0.5)
  
  # shift genome by strip width
  dd <- dd %>% mutate(xmin = xmin + strip_bins, xmax = xmax + strip_bins)
  bg <- dd %>% distinct(cell_f) %>%
    mutate(xmin = 0.5 + strip_bins, xmax = n_bins + 0.5 + strip_bins, state = toupper(fallback))
  chr_last <- as.numeric(gm$chr_last) + strip_bins
  chr_mids <- as.numeric(gm$chr_mids) + strip_bins
  
  # color palette and label order
  pal_full <- c(
    DUP  = "brown",
    AMP  = "brown",
    GAIN = "lightsalmon1",
    NOR  = "oldlace",
    LOSS = "deepskyblue",
    DEL  = "dodgerblue3"
  )
  pal_order <- c("DUP","AMP", "GAIN", "NOR", "LOSS", "DEL")
  
  states_present <- intersect(pal_order, unique(c(toupper(lab_levels), unique(dd$state))))
  pal_cna <- pal_full[states_present]
  
  site_levels <- unique(site_df$site)
  site_pal <- setNames(viridis(length(site_levels), option = "D"), site_levels)
  
  p <- ggplot() +
    geom_rect(data = site_df,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = site)) +
    scale_fill_manual(values = site_pal, drop = FALSE, name = "Tumor site") +
    ggnewscale::new_scale_fill() +
    geom_rect(data = bg,
              aes(xmin = xmin, xmax = xmax,
                  ymin = as.numeric(cell_f) - 0.5, ymax = as.numeric(cell_f) + 0.5,
                  fill = state)) +
    geom_rect(data = dd,
              aes(xmin = xmin, xmax = xmax,
                  ymin = as.numeric(cell_f) - 0.5, ymax = as.numeric(cell_f) + 0.5,
                  fill = state)) +
    scale_fill_manual(values = pal_cna, drop = FALSE, name = "CNA",
                      breaks = pal_order, labels = pal_order) +
    scale_y_continuous(breaks = seq_along(cells), labels = rev(cells), expand = c(0, 0)) +
    scale_x_continuous(limits = c(0.5, n_bins + strip_bins + 0.5),
                       breaks = chr_mids, labels = gm$chr_labels, expand = c(0, 0)) +
    geom_vline(xintercept = strip_bins + 0.5, color = "grey60", linewidth = 0.3) +
    geom_vline(xintercept = chr_last + 0.1, color = "grey70", linewidth = 0.2) +
    labs(x = "Chromosome", y = "Cell") +
    theme_minimal(base_size = 11) +
    theme(panel.grid = element_blank(), axis.ticks = element_blank())
  
  ggsave(file, plot = p, width = width, height = height, units = "in")
  invisible(p)
}

plot_S_heatmap <- function(S, ref, code_map,
                           file = "CNA_heatmap_from_S.pdf",
                           row_order = NULL, width = 14, height = 4,
                           strip_bins = NULL,
                           cluster = NULL) {
  
  stopifnot(is.matrix(S), nrow(S) == nrow(ref))
  gm <- make_genome_map(ref)
  
  cells <- colnames(S)
  if (!is.null(row_order)) cells <- intersect(row_order, cells)
  S <- S[, cells, drop = FALSE]
  
  # --- map codes to labels ---
  labs <- names(code_map)
  L <- matrix(labs[match(as.integer(S), code_map)], nrow = nrow(S), dimnames = dimnames(S))
  fb <- if ("nor" %in% labs) "nor" else labs[1]
  L[is.na(L)] <- fb
  L <- apply(L, 2, function(x){x <- toupper(x); x[x=="DUP"]<-"AMP"; x[x %in% c("NEU","NORMAL")]<-"NOR"; x})
  
  # --- rectangles per run ---
  rle_rects <- function(v){ r <- rle(v); e <- cumsum(r$lengths); s <- e - r$lengths + 1
  data.frame(state = r$values, start_idx = s, end_idx = e) }
  dd <- lapply(seq_len(ncol(L)), function(j){ x <- rle_rects(L[,j]); x$cell <- colnames(L)[j]; x }) %>% dplyr::bind_rows()
  
  # x map
  n_bins <- nrow(gm$ref); x_map <- gm$ref$bin_idx
  dd <- dd %>% dplyr::mutate(xmin = x_map[pmax(1L, start_idx)] - .5,
                             xmax = x_map[pmin(n_bins, end_idx)] + .5,
                             cell_f = factor(cell, levels = rev(cells)))
  
  # left strips
  if (is.null(strip_bins)) strip_bins <- max(50, round(0.02 * n_bins))
  
  site_df <- data.frame(cell = cells) %>%
    dplyr::mutate(site = gsub("[0-9]+$", "", sub("^.*?_(.*?)_.*$", "\\1", cell)),
                  cell_f = factor(cell, levels = rev(cells)),
                  ymin = as.numeric(cell_f)-.5, ymax = as.numeric(cell_f)+.5,
                  xmin = .5, xmax = strip_bins + .5)
  
  # align cluster to S’s columns + make the strip as wide as site strip
  cluster_strip_bins <- if (is.null(cluster)) 0 else strip_bins
  if (!is.null(cluster)) {
    cluster <- as.factor(cluster[cells])        # <- alignment fix
    cluster_df <- data.frame(cell = cells, cluster = cluster) %>%
      dplyr::mutate(cell_f = factor(cell, levels = rev(cells)),
                    ymin = as.numeric(cell_f)-.5, ymax = as.numeric(cell_f)+.5,
                    xmin = strip_bins + .5, xmax = strip_bins + cluster_strip_bins + .5)
  }
  shift_bins <- strip_bins + cluster_strip_bins
  
  dd <- dd %>% dplyr::mutate(xmin = xmin + shift_bins, xmax = xmax + shift_bins)
  bg <- dd %>% dplyr::distinct(cell_f) %>%
    dplyr::mutate(xmin = .5 + shift_bins, xmax = n_bins + .5 + shift_bins, state = "NOR")
  chr_last <- as.numeric(gm$chr_last) + shift_bins
  chr_mids <- as.numeric(gm$chr_mids) + shift_bins
  
  # palettes + factor order
  pal_order <- c("AMP","GAIN","NOR","LOSS","DEL")
  pal_cna <- c(AMP="brown", GAIN="lightsalmon1", NOR="oldlace", LOSS="deepskyblue", DEL="dodgerblue3")
  dd$state <- factor(dd$state, levels = pal_order)
  bg$state <- factor(bg$state, levels = pal_order)
  
  site_pal <- setNames(viridisLite::viridis(length(unique(site_df$site)), option = "D"),
                       unique(site_df$site))
  if (!is.null(cluster)) {
    cl_lv <- levels(cluster)
    cl_pal <- setNames(
      scales::hue_pal(l = 70, c = 100)(length(cl_lv)),  # evenly spaced around color wheel
      cl_lv
    )
  }
  
  
  p <- ggplot2::ggplot() +
    geom_rect(data = site_df, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=site)) +
    scale_fill_manual(values = site_pal, drop = FALSE, name = "Tumor site") +
    ggnewscale::new_scale_fill() +
    { if (!is.null(cluster)) geom_rect(data = cluster_df, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=cluster)) } +
    { if (!is.null(cluster)) scale_fill_manual(values = cl_pal, drop = FALSE, name = "Cluster") } +
    ggnewscale::new_scale_fill() +
    geom_rect(data = bg, aes(xmin=xmin, xmax=xmax, ymin=as.numeric(cell_f)-.5, ymax=as.numeric(cell_f)+.5, fill=state)) +
    geom_rect(data = dd, aes(xmin=xmin, xmax=xmax, ymin=as.numeric(cell_f)-.5, ymax=as.numeric(cell_f)+.5, fill=state)) +
    scale_fill_manual(values = pal_cna, drop = TRUE, name = "CNA",
                      breaks = pal_order, labels = pal_order) +
    scale_y_continuous(breaks = seq_along(cells), labels = rev(cells), expand = c(0,0)) +
    scale_x_continuous(limits = c(.5, n_bins + shift_bins + .5),
                       breaks = chr_mids, labels = gm$chr_labels, expand = c(0,0)) +
    geom_vline(xintercept = strip_bins + .5, color = "grey60", linewidth = .3) +
    { if (!is.null(cluster) && cluster_strip_bins > 0)
      geom_vline(xintercept = strip_bins + cluster_strip_bins + .5, color = "grey60", linewidth = .3) } +
    geom_vline(xintercept = chr_last + .1, color = "grey70", linewidth = .2) +
    labs(x = "Chromosome", y = "Cell") +
    theme_minimal(base_size = 11) +
    theme(panel.grid = element_blank(), axis.ticks = element_blank())
  
  ggsave(file, plot = p, width = width, height = height, units = "in")
  invisible(p)
}
