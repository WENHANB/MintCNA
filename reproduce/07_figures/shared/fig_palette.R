# Shared palette for the manuscript figures (Fig 2 simulation, Fig 3 CRC).
#
# The problem this file fixes: across panels a/b/c the same hue currently carries
# three different meanings. #1B9E77 is both "MintCNA the method" (b) and
# "subclone C1" (a); #D95F02 (subclone C2) is one JND from #D55E00 (rcCAE);
# and the Site strip in a is blue/green/orange/red, sitting 2 mm from a CN ramp
# where blue already means loss and red already means amplification.
#
# The rule here is one channel per family, never reused:
#   DIRECTION  blue <-> red   CN state, gain/loss.        panels a, c   (signal)
#   HERO       green          MintCNA only.               panel  b      (accent)
#   NEUTRAL    greys          every competitor method.    panel  b
#   METADATA   purple/teal    subclone, site.             panel  a
# Nothing outside DIRECTION is allowed to be blue or red, so a reader never has
# to ask whether a colour means "loss" or "which method".

## ---- DIRECTION: the only blue/red in the figure -------------------------------
CN5_COLORS <- c("1" = "#2C7BB6", "2" = "#8CC5E3", "3" = "#F2F2F2",
                "4" = "#F4A582", "5" = "#B2182B")
CN5_LABELS <- c("Del", "Loss", "Neu", "Gain", "Amp")   # "Nor" until 2026-10-09, "Neutral" until 2026-10-14; "Neu" matches Fig. 1 (author). CN3 below keeps the full word for Fig. 5c
# Figure 4 calls three states, not five: the ATAC arm has neither the depth nor
# the bin width to separate one-copy from homozygous events. The three hues are
# the outer poles and the midpoint of the five-state ramp, so a deletion is the
# same blue in both figures and a reader crossing from Fig 3 to Fig 4 does not
# have to re-learn the direction channel.
CN3_COLORS <- c("1" = "#2C7BB6", "2" = "#F2F2F2", "3" = "#B2182B")
CN3_LABELS <- c("Loss", "Neutral", "Gain")
GAIN <- "#B2182B"; LOSS <- "#2166AC"
TRUTH_GREY <- "#B8BCC2"          # bulk WGS silhouette in panel c

## ---- panel b method colours ---------------------------------------------------
# Three modes, because there is a real trade-off here and it is not obvious which
# side wins. Fig 2 (radar) and Fig 3 (raincloud) identify methods by these eight
# hues, so greying the competitors in Fig 3b means rcCAE is orange in one figure
# and grey in another. Against that: in Fig 3b the methods sit at fixed x with
# direct labels, so colour identifies nothing that position does not already.
#   "full"  keep all eight hues            -> cross-figure identity intact
#   "hero"  competitors one flat grey      -> maximum hero salience
#   "muted" competitor hues desaturated    -> identity kept, hero still leads
B_COLOR_MODE <- Sys.getenv("FIG4_B_COLOR", "muted")

# Figure 4 adds the multiome competitors and the ATAC-only ablation of our own
# method. MintCNA(ATAC) takes the same dark green as MintCNA(DNA) and sits in the
# same position relative to MintCNA in every figure, so a reader who has learned
# "dark green = the single-modality version of the hero" in Fig 2/3 reads Fig 4
# without being told. The six multiome competitors reuse the Fig 2/3 competitor
# hues; they are different methods, but they are muted to the same grey family in
# panel b, so no hue means two things on the same page.
# MintCNA(RNA) (Supp Fig 5 only) is the third member of the hero family: the
# single-modality arms are the hero's own green, darker for DNA and lighter for
# RNA, so the family reads as one method in three configurations and no competitor
# hue is spent on it. L* 81 / 57 / 39 keeps the three apart in greyscale and under
# deuteranopia; in Supp Fig 5 the RNA point also sits far from the other two.
# 2026-10-13: the single-modality arms moved from dark green #046B44 to OLIVE #7A8B1E
# at Dr. Xiao's request -- the two greens could not be told apart (dE 20 in CIELAB).
# Olive is dE 47 from the hero green and >= 46 from every competitor hue (nearest:
# HMMcopy), still reads as the green family, and survives the 35 % alpha fills.
# Candidates compared in shared/qa/MintCNA_DNA_colour_candidates.png.
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
    stop("unknown FIG4_B_COLOR: ", mode))
  out
}
# Label ink follows the same logic: the hero pair is named in its own green, the
# competitors in plain dark grey, in every mode.
label_ink <- function(ms) ifelse(ms %in% names(HERO), "#046B44", "grey20")

## ---- METADATA: panel a, deliberately not blue and not red ---------------------
# Clone identity is an arbitrary LABEL; cell state is biology. They must not
# share a palette. The violet/ochre family belongs to cell state, and in P111T2
# the clone numbering happens to run against it -- C1 is the AC-like (ochre)
# clone -- so a violet C1 block beside an ochre state strip read backwards.
# Two neutrals keep clone as a separate channel that competes with nothing.
# DEEPENED 2026-09-03, and the reason is rendering, not hue. At full tone the old
# pair (#3F444B / #9BA3AD) was fine -- dE00 37.4, twice this project's floor. But
# panel D draws it at size 0.22 mm and alpha 0.7, and compositing onto white took
# C2 to #B9BFC6: 1.85:1 against the panel surface, below the 3:1 floor
# validate_palette.R applies, and the pair collapsed to dE00 20.9. The clouds read
# as "dark dots vs faint dots", not as two clone colours.
# It also fixes a quieter defect nobody reported: the old C2 was only dE00 12.0
# from the "Undifferentiated" cell-state key #BEBAB2 sitting in the SAME panel C
# legend, and 6.5 from COMPETITOR grey. Now 17.8 and 15.5.
# The pair stays NEUTRAL. A hue scan at L* 52 / C* 22 found every hue window
# within 15 dE00 of something already reserved except olive (which flirts with the
# hero green at 19.0) and mauve (with the gain red at 20.4) -- the palette has
# spent every family, so relaxing the neutral rule buys nothing. The contrast is
# bought with alpha in panel_umap.R instead, where it is free.
# NEUTRAL RULE RELAXED 2026-09-03, after the author twice reported C1 and C2 do
# not separate. The earlier note claimed the palette had spent every hue; a full
# scan of CIELAB says otherwise -- 5,213 pairs clear a 15 dE00 floor against every
# reserved colour AND a 20 dE00 floor against the ones that actually sit next to
# panel D (the copy-number blue/red and the cell-state violet/ochre in panel C's
# legend, a centimetre above). Two colours differing only in LIGHTNESS give the eye
# no categorical cue at 0.28 mm; that is why grey-on-grey failed twice.
# Teal is the pick because it is the hue furthest from everything this figure
# already means: violet/ochre are cell state, blue/red are copy-number direction,
# green is the method. Rendered at alpha 0.9 the pair is dE00 40.8, both members
# clear 3:1 against white, and deuteranopia/protanopia cost it under 4 dE00.
# NEUTRAL + ONE COOL ACCENT. Checked against skill/nature-figure design-theory.md,
# which was NOT consulted for the first four attempts and should have been:
#   "For dense multi-panel figures, family consistency beats maximal hue separation"
#   "Reserve green/red for arrows, gains, drops, thresholds, or signed biological
#    direction"
# The earlier petrol/coral pair separated best of anything measured (dE00 50.0) but
# spent the WARM RED family on an arbitrary label, on a page where red already
# means Gain -- a rule-3 violation. And beside the cell-state UMAP its coral read
# as a cousin of the AC-like ochre.
# The skill's own "one cool family, two steps" cannot be used here either: its
# blue_main #0F4D92 is 8.7 dE00 from our copy-number Loss blue, and a two-step cool
# pair (petrol + teal) separates by only 27.8 -- the same failure as grey-on-grey.
# That rule assumes a palette not already spent on biology; ours has blue/red on
# direction, violet/ochre on cell state and green on the method.
# What survives is the palette file's ORIGINAL intent -- clone as a channel that
# competes with nothing -- with the second member promoted from a grey to the
# skill's sanctioned teal accent (#42949E). C1 is #3A3A3A, not #272727: the neutral
# must sit BELOW the structural-grey band (chromosome fills #5A6069, COMPETITOR
# #B4B8BE, TRUTH_GREY #B8BCC2), and the skill's own neutral_dark #4D4D4D is 8.4
# dE00 from the chromosome fill, so mid-grey is not available.
# Measured: dE00 33.2 as rendered at alpha 0.9, 27.3 under dichromacy, 26.2 from
# the nearest cell state, 18.3 from the nearest copy-number colour, 13.8 from the
# nearest structural grey. Block-label ink is unchanged (white on C1, grey10 on C2
# at 4.94:1).
CLONE_COLORS2 <- c(C1 = "#3A3A3A", C2 = "#42949E")        # neutral / teal accent
SITE_COLORS2  <- c(PT = "#2D2352", LN = "#6250A6",         # one violet ramp, widened:
                   ML = "#9C8FD1", MP = "#D3CCEC")         # at the first spacing ML and MP
                                                           # were one step apart and the
                                                           # legend keys read as one colour.
                                                           # PT darkest = the primary.

## ---- malignant cell state: Figure 4 biology panels -----------------------------
# Three categories, not five. OPC-like and AC-like are the contrast the panels
# exist to show; NPC-like, MES-like and Undifferentiated are each <=5% of nuclei
# in the sample shown and are merged into one grey. The full five-way breakdown
# is in the panel's source data. Five categorical hues also cannot all clear the
# normal-vision separation floor once blue, red and green are reserved.
# Checked with Code/figures/validate_palette.R: pairwise normal-vision dE 42.9 /
# 41.3 / 15.1, every CVD type above target.
STATE_COLORS <- c("OPC-like" = "#5B2D8E", "AC-like" = "#E8961E", "other" = "#BEBAB2")

# The paper's own five malignant cell states, as TWO FAMILIES rather than five
# unrelated hues: violet for the progenitor/oligodendroglial side (OPC, NPC),
# ochre for the astrocytic/mesenchymal side (AC, MES), grey for undifferentiated.
# The two family hues are deliberately the clone hues, because the claim panel a
# makes is that the clone blocks and this strip line up -- shared colour is the
# argument, not a collision. Within a family the second member is the lighter
# (OPC/NPC) or darker (AC/MES) member so the pair stays separable in greyscale.
CELLSTATE_COLORS <- c(
  "OPC-like"         = "#5B2D8E",
  "NPC-like"         = "#A88AD0",
  "AC-like"          = "#E8961E",
  "MES-like"         = "#9C5A16",
  "Undifferentiated" = "#BEBAB2")

## ---- type ---------------------------------------------------------------------
# Nature Communications: in-figure text 5-7 pt, Arial/Helvetica. No Arial on this
# host, so ask for "Helvetica": fontconfig resolves it to Nimbus Sans, the URW
# Helvetica clone with identical metrics, and -- unlike the literal name
# "Nimbus Sans" -- R's base pdf() device knows it as a standard font. That
# matters because cowplot measures text on the default device while assembling;
# an unknown family sends it down the CID path and the build dies with
# "failed to find or load PDF CID font". Panels are written at FINAL physical
# size, so these are true printed points and nothing is rescaled at assembly.
# Verified on this host: "Liberation Sans" (what Fig 2 ask for), "Arial" and
# "Helvetica" ALL resolve to NimbusSans-Regular.otf, so the typeface is already
# consistent with the other main-text figures. Bare "sans" does NOT -- it lands on
# DejaVu Sans, which is how DejaVu ended up embedded alongside Nimbus in the first
# build: ComplexHeatmap's gpar() calls carry no family and inherit the device
# default. Passing family = FIG_FAMILY to cairo_pdf() fixes that at the device.
# 2026-10-13: the real Arial (Microsoft core fonts, arial32.exe via the corefonts
# SourceForge package, EULA-licensed) is installed in ~/.fonts, so the figures now embed
# genuine Arial (TrueType) instead of the Nimbus Sans clone -- NC's "Arial or Helvetica"
# read literally. R's base pdf() device does not know "Arial", and cowplot measures text
# on the DEFAULT device while assembling, so the default device is pointed at cairo_pdf
# (fontconfig-backed) below; every script writes through cairo_pdf anyway.
FIG_FAMILY <- "Arial"
options(device = function(file = NULL, ...) grDevices::cairo_pdf(if (is.null(file)) tempfile(fileext = ".pdf") else file, ...))
# ComplexHeatmap's grid.grabExpr() opens a null base pdf() device to lay the heatmap out,
# and that device knows Helvetica but not Arial ("failed to find or load PDF CID font",
# Fig 4 build 2026-10-13). Arial and Helvetica are metric-compatible, so "Arial" is
# registered there as an alias of the Helvetica metrics; it is only ever used for
# measuring -- the figure itself is written by cairo_pdf, which embeds the real Arial.
if (!"Arial" %in% names(grDevices::pdfFonts()))
  grDevices::pdfFonts(Arial = grDevices::pdfFonts()[["Helvetica"]])
if (!"Arial" %in% names(grDevices::postscriptFonts()))
  grDevices::postscriptFonts(Arial = grDevices::postscriptFonts()[["Helvetica"]])

# Font-size floor. Nature Communications asks 5-7 pt for in-figure text. plotmath
# renders a sub/superscript at ~0.7x its parent, so a 6 pt axis label carrying
# MintCNA^DNA puts "DNA" on the page at 4.2 pt -- under the floor even though the
# parent is compliant. Any text that carries a sub/superscript therefore uses
# PT_SUP as its parent, not BASE_PT.
BASE_PT   <- 7
PT_SUP    <- 7.4        # 7.4 * 0.7 = 5.2 pt, clear of the 5 pt floor
# Legend-title / secondary-label size, 6.5 pt. It lived only inside
# panel_a_heatmap.R (sourced LAST by the assembly), so a panel sourced earlier that
# used it -- E and F's gene labels from 2026-09-06 -- saw an undefined name. The
# heatmap keeps its own identical assignment; this one makes it global.
FS_LGT <- 6.5
PT_SMALL  <- 5.5        # smallest plain text (no sub/superscript) in the figure

## ---- deterministic gene-label placement (Fig 4 E/F) ----------------------------
# ggrepel avoids only the points in ITS OWN data, and when the whole cloud is handed
# to it as obstacles (2026-09-06) it settles in local optima: OLIG1 on the PDGFRA
# point, CLU across BCAN. Nine labels do not need a physics simulation. place_labels()
# scores every candidate position on rings around each labelled point against ALL
# points and every label already placed, and takes the cheapest. It is greedy in the
# order given, deterministic, and unit-aware: box sizes are computed in mm from the
# panel size and the axis ranges, so the result is the exported figure, not a guess.
#   x, y      all points (data units);  col  logical, TRUE for a coloured point
#   idx       row indices to label, most important first;  txt  their labels
#   xr, yr    axis ranges INCLUDING expansion;  panel_mm  c(width, height) of the
#             plot area;  pt_mm  point radius in mm per row;  text_pt  label size
#   Cost per candidate: grey points under the (padded) box, coloured points under it,
#   labels already placed, coloured points a leader would cross, leader length,
#   and a large penalty for leaving the panel (or a caller-supplied keep-out box).
# Weights (2026-09-06, second pass): a label wants to sit BESIDE its point. Grey
# points under a name cost little (the cloud is 45% alpha and reads through), a
# millimetre of leader costs a lot, so leaders appear only where the point is
# boxed in by coloured points or other names. The first pass had the weights the
# other way round and drew nine leaders across the panel.
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
      # the box may not cover, or touch, the point it names (LHFPL3 sat on its own
      # point at the top edge until this line existed)
      if (X[i] > x0 - pt_mm[i] - 0.3 && X[i] < x1 + pt_mm[i] + 0.3 &&
          Y[i] > y0 - pt_mm[i] - 0.3 && Y[i] < y1 + pt_mm[i] + 0.3) cost <- cost + w_col
      for (b in placed) if (x0 < b[2] && x1 > b[1] && y0 < b[4] && y1 > b[3]) cost <- cost + w_lab
      # no ambiguity: the box must be closer to its own point than to any other
      # coloured point (PDGFRA's name sat beside the OLIG1 point until this line)
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

## ---- line-chart identity (Fig 2 A/B) ------------------------------------------
# Fig 3B can mute the competitors because the methods sit at fixed x with direct
# labels, so colour there carries no information. In a line chart the series
# overlap and colour is the ONLY identity channel, so Fig 2 uses the full hues.
# Hero salience is bought with line WEIGHT instead of saturation, which keeps the
# hue mapping identical across the two figures.
## ---- plotmath labels ----------------------------------------------------------
# Superscript rendering for the single-modality arms. Keyed BY METHOD, never by
# position: a vector aligned to a performance ordering mislabels a violin the
# moment the ordering changes.
# "AtaCNV" is the name on the 2023 preprint and on the R package; the method was
# PUBLISHED as AtaCNA (Wang, Jin, Shi & Xi, Cell Rep Methods 2025;5(1):100939).
# The figure must carry the name the citation carries. The internal key stays
# "AtaCNV" because that is what the result tables and the package are keyed on --
# only the DISPLAYED label changes, so no join is touched.
METHOD_PLOTMATH <- c("MintCNA" = "MintCNA", "MintCNA(DNA)" = "MintCNA^DNA",
                     "MintCNA(ATAC)" = "MintCNA^ATAC", "MintCNA(RNA)" = "MintCNA^RNA", "Copy-scAT" = "'Copy-scAT'",
                     "AtaCNV" = "'AtaCNA'",
                     "CONGAS+" = "'CONGAS+'")
method_plotmath_chr <- function(ms)
  unname(ifelse(ms %in% names(METHOD_PLOTMATH), METHOD_PLOTMATH[ms],
                sprintf("'%s'", ms)))
method_labels <- function(ms) parse(text = method_plotmath_chr(ms))

LINE_W_HERO <- 0.85
LINE_W_COMP <- 0.42
PT_HERO     <- 1.5
PT_COMP     <- 1.0
method_line_w  <- function(ms) ifelse(ms %in% names(HERO), LINE_W_HERO, LINE_W_COMP)
method_pt_size <- function(ms) ifelse(ms %in% names(HERO), PT_HERO,     PT_COMP)
