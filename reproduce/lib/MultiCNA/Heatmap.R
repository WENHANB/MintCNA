#' Fast genome-wide CNA heatmap (rectangles per segment)
#'
#' Draws a compact genome-wide heatmap where each row is a cell and each CNA
#' segment is a colored rectangle (dup/del). Chromosome boundaries and labels
#' are shown along the x-axis for orientation.
#' @param cna_df CNA table with columns `chr`, `start_idx`, `end_idx`, `state`
#'   (factor or character: "dup"/"del"), and `cell`.
#' @param ref Bin reference (see \link{make_genome_map}) with columns `chr`, `start`.
#' @param file Output filename (PDF/PNG supported by \code{ggsave}).
#' @param row_order Optional character vector to order rows (cells). Defaults
#'   to the order of appearance in `cna_df`.
#' @param width,height Plot size passed to \code{ggsave} (in inches).
#' @return Invisibly returns the \code{ggplot} object; also writes the figure to \code{file}.
#' @seealso \link{make_genome_map}, \link{map_cna_to_genome}, \link[ggplot2]{ggsave}
#' @examples
#' \dontrun{
#' gm <- make_genome_map(ref)
#' dd <- map_cna_to_genome(cna_df, gm)
#' plot_cna_gg(cna_df, ref, file = "CNA_heatmap.pdf")
#' }
#' @export
library(dplyr)
library(ggplot2)
# plot_cna_gg <- function(cna_df, ref, file = "CNA_heatmap.pdf",
#                         row_order = NULL, width = 14, height = 4) {
#   gm <- make_genome_map(ref)
#   dd <- map_cna_to_genome(cna_df, gm)
#   
#   # order rows (cells)
#   cells <- if (is.null(row_order)) unique(dd$cell) else row_order
#   cells <- intersect(cells, unique(dd$cell))
#   dd$cell_f <- factor(dd$cell, levels = rev(cells))  # first cell at top
#   
#   # draw a neutral background per row, then overlay dup/del
#   bg <- distinct(dd, cell_f) %>%
#     mutate(xmin = 0.5, xmax = nrow(gm$ref) + 0.5,
#            state = "neutral")
#   
#   pal <- c(dup = "brown", del = "dodgerblue3", neutral = "oldlace")
#   
#   p <- ggplot() +
#     geom_rect(data = bg,
#               aes(xmin = xmin, xmax = xmax,
#                   ymin = as.numeric(cell_f) - 0.5,
#                   ymax = as.numeric(cell_f) + 0.5,
#                   fill = state)) +
#     geom_rect(data = dd,
#               aes(xmin = xmin, xmax = xmax,
#                   ymin = as.numeric(cell_f) - 0.5,
#                   ymax = as.numeric(cell_f) + 0.5,
#                   fill = state)) +
#     scale_fill_manual(values = pal, drop = FALSE) +
#     scale_y_continuous(breaks = seq_along(cells),
#                        labels = rev(cells), expand = c(0, 0)) +
#     scale_x_continuous(breaks = as.numeric(gm$chr_mids),
#                        labels = gm$chr_labels, expand = c(0, 0)) +
#     geom_vline(xintercept = as.numeric(gm$chr_last) + 0.1,
#                color = "grey70", linewidth = 0.2) +
#     labs(x = "Chromosome", y = "Cell", fill = "CNA") +
#     theme_minimal(base_size = 11) +
#     theme(panel.grid = element_blank(),
#           axis.ticks = element_blank())
#   
#   ggsave(file, plot = p, width = width, height = height, units = "in")
# }

library(dplyr)
library(ggplot2)
library(ggnewscale)
library(RColorBrewer)
library(viridisLite)

plot_cna_gg <- function(cna_df, ref, file = "CNA_heatmap.pdf",
                        row_order = NULL, width = 14, height = 4,
                        strip_bins = NULL) {
  
  gm <- make_genome_map(ref)
  dd <- map_cna_to_genome(cna_df, gm)
  
  # order rows
  cells <- if (is.null(row_order)) unique(dd$cell) else row_order
  cells <- intersect(cells, unique(dd$cell))
  dd$cell_f <- factor(dd$cell, levels = rev(cells))
  
  # choose a visible strip width (default ~2% of genome bins, min 50)
  n_bins <- nrow(gm$ref)
  if (is.null(strip_bins)) strip_bins <- max(50, round(0.02 * n_bins))
  
  # site from cell (middle token), collapse LN1->LN, PT1->PT, etc.
  site_df <- cna_df %>%
    distinct(cell) %>%
    mutate(site = sub("^.*?_(.*?)_.*$", "\\1", cell),
           site = gsub("[0-9]+$", "", site),
           cell_f = factor(cell, levels = rev(cells)),
           ymin = as.numeric(cell_f) - 0.5,
           ymax = as.numeric(cell_f) + 0.5,
           xmin = 0.5,
           xmax = strip_bins + 0.5)
  

  site_levels <- unique(site_df$site)
  n_sites <- length(site_levels)
  # use viridis discrete palette
  site_cols <- viridis(n_sites, option = "D")   # options: "A","B","C","D","E","F","G","H"
  site_pal <- setNames(site_cols, site_levels)
  
  # shift genome by strip_bins to the right
  dd <- dd %>% mutate(xmin = xmin + strip_bins, xmax = xmax + strip_bins)
  bg <- distinct(dd, cell_f) %>%
    mutate(xmin = 0.5 + strip_bins, xmax = n_bins + 0.5 + strip_bins, state = "neutral")
  chr_last <- as.numeric(gm$chr_last) + strip_bins
  chr_mids <- as.numeric(gm$chr_mids) + strip_bins
  
  pal_cna <- c(dup = "brown", del = "dodgerblue3", neutral = "oldlace")
  
  p <- ggplot() +
    # left color block
    geom_rect(data = site_df,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = site)) +
    scale_fill_manual(values = site_pal, drop = FALSE, name = "Tumor site") +
    ggnewscale::new_scale_fill() +
    # heatmap
    geom_rect(data = bg,
              aes(xmin = xmin, xmax = xmax,
                  ymin = as.numeric(cell_f)-0.5, ymax = as.numeric(cell_f)+0.5,
                  fill = state)) +
    geom_rect(data = dd,
              aes(xmin = xmin, xmax = xmax,
                  ymin = as.numeric(cell_f)-0.5, ymax = as.numeric(cell_f)+0.5,
                  fill = state)) +
    scale_fill_manual(values = pal_cna, drop = FALSE, name = "CNA") +
    # axes
    scale_y_continuous(breaks = seq_along(cells), labels = rev(cells), expand = c(0, 0)) +
    scale_x_continuous(limits = c(0.5, n_bins + strip_bins + 0.5),
                       breaks = chr_mids, labels = gm$chr_labels, expand = c(0, 0)) +
    # divider between strip and genome
    geom_vline(xintercept = strip_bins + 0.5, color = "grey60", linewidth = 0.3) +
    geom_vline(xintercept = chr_last + 0.1, color = "grey70", linewidth = 0.2) +
    labs(x = "Chromosome", y = "Cell") +
    theme_minimal(base_size = 11) +
    theme(panel.grid = element_blank(), axis.ticks = element_blank())
  
  ggsave(file, plot = p, width = width, height = height, units = "in")
}


make_genome_map <- function(ref) {
  ref2 <- ref %>%
    mutate(chr = as.character(chr),
           chr_num = suppressWarnings(as.integer(sub("^chr", "", chr)))) %>%
    arrange(is.na(chr_num), chr_num, chr, start) %>%
    mutate(bin_idx = row_number())
  
  chr_last <- tapply(ref2$bin_idx, ref2$chr, max)
  chr_first <- tapply(ref2$bin_idx, ref2$chr, min)
  chr_levels <- names(chr_last)
  chr_mids <- floor((chr_first + chr_last) / 2)
  list(ref = ref2,
       chr_last = chr_last,
       chr_mids = chr_mids,
       chr_labels = sub("^chr", "", chr_levels))
}


map_cna_to_genome <- function(cna_df, gm) {
  # lookup vectors: global x for each bin index
  x_map <- gm$ref$bin_idx
  names(x_map) <- seq_len(nrow(gm$ref))
  
  cna_df %>%
    mutate(xmin = x_map[pmax(1L, start_idx)] - 0.5,
           xmax = x_map[pmin(nrow(gm$ref), end_idx)] + 0.5)
}

# cna_df must have: cell, state, start_idx, end_idx (aligned to rows in ref)
# plot_cna_gg(cna_df_all[which(cna_df_all$cell %in% converge_cell),], 
#             ref, file = "Result/CRC/Segmentation/CNV/CNA_all_heatmap_100kb_h10_TGMM_QNALL.pdf", width = 12, height = 8)