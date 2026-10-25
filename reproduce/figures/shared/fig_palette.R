# Shared palette, type settings and label helpers for the MintCNA manuscript figures.
# Sourced by every figure script under reproduce/figures/ (main figures, fig4_crc/, fig5_multiome/
# and supp/) as source("reproduce/figures/shared/fig_palette.R"), from the root of the repository.
#
# One colour channel per meaning, never reused: blue/red for copy-number direction (loss/gain),
# green for MintCNA and its single-modality arms, muted hues or grey for the competing methods,
# violet/teal/ochre for sample metadata (site, subclone, cell state). Nothing outside the
# direction channel is blue or red, so a colour never has to be read as either "loss" or "method".

## ---- copy-number direction: the only blue/red in the figures -----------------------------
CN5_COLORS <- c("1" = "#2C7BB6", "2" = "#8CC5E3", "3" = "#F2F2F2",
                "4" = "#F4A582", "5" = "#B2182B")
CN5_LABELS <- c("Del", "Loss", "Neu", "Gain", "Amp")   # "Neu" matches the state labels of Fig. 1
# Fig 5 calls three states rather than five: the ATAC arm has neither the depth nor the bin width
# to separate one-copy from homozygous events. Its hues are the two poles and the midpoint of the
# five-state ramp, so a deletion is the same blue in Fig 4 and Fig 5.
CN3_COLORS <- c("1" = "#2C7BB6", "2" = "#F2F2F2", "3" = "#B2182B")
CN3_LABELS <- c("Loss", "Neutral", "Gain")
GAIN <- "#B2182B"; LOSS <- "#2166AC"
TRUTH_GREY <- "#B8BCC2"          # bulk WGS silhouette

## ---- method colours ------------------------------------------------------------------------
# Figs 2 and 3 identify the methods by the eight hues of METHOD_HUES. In the accuracy panels the
# methods sit at fixed x with direct labels, so colour identifies nothing that position does not;
# method_fill() therefore offers three modes for the competitors: "full" keeps their hues, "hero"
# paints them one flat grey, and "muted" desaturates them so the MintCNA family leads while the
# cross-figure identity is kept. Fig 4b uses the default below; Fig 5b asks for "full".
B_COLOR_MODE <- "muted"

# MintCNA(DNA) and MintCNA(ATAC) share one olive, the single-modality arm in every figure, and
# MintCNA(RNA) (Supp Fig 5) is a lighter green, so the family reads as one method in three
# configurations and no competitor hue is spent on it. The three greens stay apart in greyscale
# and under deuteranopia. The multiome competitors (Fig 5) reuse the competitor hues of Figs 2-4.
METHOD_HUES <- c("MintCNA" = "#1B9E77", "MintCNA(DNA)" = "#7A8B1E",
                 "MintCNA(ATAC)" = "#7A8B1E", "MintCNA(RNA)" = "#A6D96A",
                 "SCOPE" = "#0072B2", "HMMcopy" = "#E69F00", "SCONCE" = "#984EA3",
                 "SCCNV" = "#56B4E9", "rcCAE" = "#D55E00", "SCICoNE" = "#CC79A7",
                 "AtaCNV" = "#0072B2", "Copy-scAT" = "#E69F00",
                 "RIDDLER" = "#56B4E9", "epiAneufinder" = "#984EA3",
                 "inferCNV" = "#D55E00", "CONGAS+" = "#CC79A7")
HERO       <- METHOD_HUES[c("MintCNA", "MintCNA(DNA)", "MintCNA(ATAC)", "MintCNA(RNA)")]
COMPETITOR <- "#B4B8BE"

method_fill <- function(ms, mode = B_COLOR_MODE) {
  hero <- intersect(ms, names(HERO))
  comp <- setdiff(ms, names(HERO))
  out <- setNames(character(length(ms)), ms)
  out[hero] <- HERO[hero]
  out[comp] <- switch(mode,
    full  = METHOD_HUES[comp],
    hero  = rep(COMPETITOR, length(comp)),
    muted = colorspace::lighten(colorspace::desaturate(METHOD_HUES[comp], 0.55), 0.30),
    stop("unknown colour mode: ", mode))
  out
}
# Label ink: the MintCNA family is named in its own green, the competitors in dark grey, in every mode.
label_ink <- function(ms) ifelse(ms %in% names(HERO), "#046B44", "grey20")

## ---- sample metadata: neither blue nor red -------------------------------------------------
# Clone identity is an arbitrary label and shares no palette with cell state (violet / ochre
# below) or with copy-number direction. The pair is a near-black neutral and a teal accent: two
# greys differing only in lightness do not separate at the small point sizes of the UMAP panel,
# and teal is the hue furthest from everything else the figures already mean.
CLONE_COLORS2 <- c(C1 = "#3A3A3A", C2 = "#42949E")        # neutral / teal accent
SITE_COLORS2  <- c(PT = "#2D2352", LN = "#6250A6",         # one violet ramp, PT (the primary) darkest,
                   ML = "#9C8FD1", MP = "#D3CCEC")         # spaced so ML and MP stay distinct in the legend

## ---- malignant cell state: Fig 5 biology panels --------------------------------------------
# Three categories: OPC-like and AC-like are the contrast the panels show; NPC-like, MES-like and
# Undifferentiated are each <= 5 % of nuclei in the sample shown and are merged into one grey
# (the five-way breakdown is in the source data).
STATE_COLORS <- c("OPC-like" = "#5B2D8E", "AC-like" = "#E8961E", "other" = "#BEBAB2")

# The five malignant cell states as two families rather than five unrelated hues: violet for the
# progenitor side (OPC, NPC), ochre for the astrocytic/mesenchymal side (AC, MES), grey for
# undifferentiated. Within a family the second member differs in lightness so the pair stays
# separable in greyscale.
CELLSTATE_COLORS <- c(
  "OPC-like"         = "#5B2D8E",
  "NPC-like"         = "#A88AD0",
  "AC-like"          = "#E8961E",
  "MES-like"         = "#9C5A16",
  "Undifferentiated" = "#BEBAB2")

## ---- type ----------------------------------------------------------------------------------
# The figures use Arial (FIG_FAMILY) and fall back to whatever fontconfig substitutes for it. Every
# script writes through cairo_pdf; the default device is pointed at cairo_pdf as well because
# cowplot measures text on the default device while assembling, and R's base pdf() device does
# not know Arial.
FIG_FAMILY <- "Arial"
options(device = function(file = NULL, ...) grDevices::cairo_pdf(if (is.null(file)) tempfile(fileext = ".pdf") else file, ...))
# ComplexHeatmap lays a heatmap out on a null base pdf() device, which knows Helvetica but not
# Arial. The two are metric-compatible, so "Arial" is registered there as an alias of the
# Helvetica metrics; it is only used for measuring, and the output is written by cairo_pdf.
if (!"Arial" %in% names(grDevices::pdfFonts()))
  grDevices::pdfFonts(Arial = grDevices::pdfFonts()[["Helvetica"]])
if (!"Arial" %in% names(grDevices::postscriptFonts()))
  grDevices::postscriptFonts(Arial = grDevices::postscriptFonts()[["Helvetica"]])

# Sizes are printed points: the panels are drawn at final physical size and nothing is rescaled
# at assembly. plotmath renders a sub/superscript at about 0.7x its parent, so text carrying one
# (MintCNA^DNA) uses PT_SUP as its parent to keep the superscript above 5 pt.
BASE_PT   <- 7
PT_SUP    <- 7.4        # 7.4 * 0.7 = 5.2 pt
FS_LGT <- 6.5           # legend titles and secondary labels
PT_SMALL  <- 5.5        # smallest plain text (no sub/superscript)

## ---- deterministic gene-label placement (Fig 5, panel_bio.R and panel_de.R) ----------------
# ggrepel avoids only the points in its own data and settles in local optima when the whole cloud
# is handed to it as obstacles. place_labels() instead scores every candidate position on rings
# around each labelled point against all points and every label already placed, and takes the
# cheapest. It is greedy in the order given, deterministic, and works in mm (box sizes follow from
# the panel size and the axis ranges), so the result is what the exported figure shows.
#   x, y      all points (data units);  col  logical, TRUE for a coloured point
#   idx       row indices to label, most important first;  txt  their labels
#   xr, yr    axis ranges including expansion;  panel_mm  c(width, height) of the plot area;
#   pt_mm     point radius in mm per row;  text_pt  label size
# Cost per candidate: grey points under the padded box (cheap, the cloud reads through at 45 %
# alpha), coloured points under it, overlap with labels already placed, coloured points a leader
# would cross, leader length (expensive, so leaders appear only where a point is boxed in), and a
# large penalty for leaving the panel or entering a caller-supplied keep-out box.
place_labels <- function(x, y, col, idx, txt, xr, yr, panel_mm, pt_mm, text_pt = 6.5,
                         pad_mm = 0.35, radii_mm = c(1.8, 2.4, 3.2, 4.2, 5.5, 7, 9, 12),
                         angles = seq(0, 345, by = 15), keepout = NULL,
                         w_grey = 1, w_col = 200, w_lab = 400, w_cross = 150, w_len = 6,
                         seg_min_mm = 1.6) {
  sx <- panel_mm[1] / diff(xr); sy <- panel_mm[2] / diff(yr)        # mm per unit
  X <- (x - xr[1]) * sx; Y <- (y - yr[1]) * sy                        # everything in mm
  ch_w <- text_pt * 0.3528 * 0.58; ch_h <- text_pt * 0.3528 * 0.95
  placed <- list(); out <- vector("list", length(idx))
  seg_hit <- function(ax, ay, bx, by, px, py, r) {                    # point within r of segment
    vx <- bx - ax; vy <- by - ay; L2 <- vx^2 + vy^2
    t <- pmin(1, pmax(0, ((px - ax) * vx + (py - ay) * vy) / L2))
    sqrt((ax + t * vx - px)^2 + (ay + t * vy - py)^2) < r
  }
  for (k in seq_along(idx)) {
    i <- idx[k]; hw <- nchar(txt[k]) * ch_w / 2 + pad_mm; hh <- ch_h / 2 + pad_mm
    best <- NULL
    for (r in radii_mm) for (a in angles) {
      cx <- X[i] + r * cos(a * pi / 180); cy <- Y[i] + r * sin(a * pi / 180)
      x0 <- cx - hw; x1 <- cx + hw; y0 <- cy - hh; y1 <- cy + hh
      cost <- 0
      if (x0 < 0.3 || x1 > panel_mm[1] - 0.3 || y0 < 0.3 || y1 > panel_mm[2] - 0.3) cost <- cost + 1e4
      if (!is.null(keepout) && x0 < keepout[2] && x1 > keepout[1] && y0 < keepout[4] && y1 > keepout[3]) cost <- cost + 1e4
      inb <- X > x0 - pt_mm & X < x1 + pt_mm & Y > y0 - pt_mm & Y < y1 + pt_mm
      inb[i] <- FALSE
      cost <- cost + w_grey * sum(inb & !col) + w_col * sum(inb & col)
      # the box may not cover or touch the point it names
      if (X[i] > x0 - pt_mm[i] - 0.3 && X[i] < x1 + pt_mm[i] + 0.3 &&
          Y[i] > y0 - pt_mm[i] - 0.3 && Y[i] < y1 + pt_mm[i] + 0.3) cost <- cost + w_col
      for (b in placed) if (x0 < b[2] && x1 > b[1] && y0 < b[4] && y1 > b[3]) cost <- cost + w_lab
      # the box must be closer to its own point than to any other coloured point
      dbox <- sqrt(pmax(0, x0 - X, X - x1)^2 + pmax(0, y0 - Y, Y - y1)^2)
      if (any(col & seq_along(X) != i & dbox < dbox[i] + 0.6)) cost <- cost + w_col
      # leader: from the box edge nearest the point, to the point
      ax <- pmin(pmax(X[i], x0), x1); ay <- pmin(pmax(Y[i], y0), y1)
      len <- sqrt((ax - X[i])^2 + (ay - Y[i])^2)
      if (len > 0.5) {
        hit <- seg_hit(ax, ay, X[i], Y[i], X, Y, pt_mm + 0.25); hit[i] <- FALSE
        cost <- cost + w_cross * sum(hit & col) + 0.5 * w_grey * sum(hit & !col)
      }
      cost <- cost + w_len * len
      if (is.null(best) || cost < best$cost)
        best <- list(cost = cost, cx = cx, cy = cy, box = c(x0, x1, y0, y1), ax = ax, ay = ay, len = len)
    }
    placed[[k]] <- best$box
    out[[k]] <- data.frame(gene = txt[k], lx = best$cx / sx + xr[1], ly = best$cy / sy + yr[1],
                           ax = best$ax / sx + xr[1], ay = best$ay / sy + yr[1],
                           px = x[i], py = y[i], seg = best$len > seg_min_mm, cost = best$cost)
  }
  do.call(rbind, out)
}

## ---- plotmath labels -----------------------------------------------------------------------
# Display labels keyed by method name, never by position, so a reordered axis cannot mislabel a
# method. "AtaCNV" is the key of the result tables and of the R package; the method was published
# as AtaCNA (Wang, Jin, Shi & Xi, Cell Rep Methods 2025;5(1):100939), so only the displayed label
# carries the published name and no join is touched.
METHOD_PLOTMATH <- c("MintCNA" = "MintCNA", "MintCNA(DNA)" = "MintCNA^DNA",
                     "MintCNA(ATAC)" = "MintCNA^ATAC", "MintCNA(RNA)" = "MintCNA^RNA", "Copy-scAT" = "'Copy-scAT'",
                     "AtaCNV" = "'AtaCNA'",
                     "CONGAS+" = "'CONGAS+'")
method_plotmath_chr <- function(ms)
  unname(ifelse(ms %in% names(METHOD_PLOTMATH), METHOD_PLOTMATH[ms],
                sprintf("'%s'", ms)))
method_labels <- function(ms) parse(text = method_plotmath_chr(ms))

# Line weight and point size for series plots: the MintCNA family heavier than the competitors.
LINE_W_HERO <- 0.85
LINE_W_COMP <- 0.42
PT_HERO     <- 1.5
PT_COMP     <- 1.0
method_line_w  <- function(ms) ifelse(ms %in% names(HERO), LINE_W_HERO, LINE_W_COMP)
method_pt_size <- function(ms) ifelse(ms %in% names(HERO), PT_HERO,     PT_COMP)
