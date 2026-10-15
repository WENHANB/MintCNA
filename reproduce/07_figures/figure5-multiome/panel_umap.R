#!/usr/bin/env Rscript
# Panel D — is the AG-DCAE subclone partition accurate?
#
# The same RNA UMAP drawn twice: coloured by the subclone MintCNA infers from the
# AG-DCAE latent space, and coloured
# by the source study's own RNA-derived cell state. The two colourings are the
# same picture, which is the point.
#
# WHAT "WITHOUT RNA" IS ALLOWED TO MEAN HERE -- checked against the pipeline,
# because the earlier comment said "copy number only, no transcript ever read"
# and that is not literally true.
#   TRUE  the AG-DCAE latent is trained on the ATAC log2-ratio matrix alone
#         (run_all.slurm passes only l2r_dna.csv; denoise.py takes one matrix),
#         and cluster_latent_kmeans sees only latent columns.
#   TRUE  and stronger, because it is MEASURED rather than argued: the subclone
#         vectors of the full run and of the ATAC-only ablation are identical,
#         ARI = 1.000, 690/305 both ways. The partition does not depend on RNA.
#         (Caveat for honesty: the two arms share one AE training by symlink, so
#         this is invariance confirmed, not an independent re-derivation.)
#   FALSE as an absolute. export_cell_lists.R keeps a nucleus in the DIPLOID
#         REFERENCE only if its RNA and ATAC lineage labels agree, and the
#         reference is the denominator of the log2 ratio -- so an RNA-derived
#         variable is evaluated upstream of the latent. For P111T2 it excluded
#         4 nuclei, all 4 of which fail the >3,000-fragment ATAC depth QC
#         anyway; the reference is the same 815 nuclei either way.
# So: no expression feature enters the partition, and the partition is provably
# invariant to RNA -- but do not write "no transcript was ever read".
#
# The titles NAME THE TOOL rather than the quantity. "Subclone (copy number)"
# described what the labels are made of; the panel's job is to say whose labels
# they are, and MintCNA is what the reader is being asked to trust. The right-hand
# title matches panel C's legend wording exactly, so the same thing is called the
# same thing twice on one page.
#
# NO ACCURACY NUMBERS ON THE PANEL. The titles used to carry the 15-NN agreements
# (0.983 and 0.994) and they raised more questions than they answered: 0.994 looks
# like a competitor beating us until you know the cell-state labels were themselves
# defined from this same RNA and are therefore a CEILING, and that takes a sentence,
# not a number floating over a plot. Both numbers, the 0.693 majority floor, and
# the ceiling argument now live in the caption, where there is room to say what
# they mean. The panel makes the qualitative claim; the caption quantifies it.
# (They are still computed and printed to the console below, and written into
# Result/qc/clone_embedding_*_summary.tsv.)
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(cowplot) })

UMAP_SAMPLE <- Sys.getenv("FIG5_UMAP_SAMPLE", "P111T2")
QC <- file.path(ROOT, "ATAC_RNA/Result/qc")

.ud <- fread(file.path(QC, sprintf("clone_embedding_%s.tsv", UMAP_SAMPLE)))
.us <- fread(file.path(QC, sprintf("clone_embedding_%s_summary.tsv", UMAP_SAMPLE)))

# NO LEGEND IN EITHER UMAP.
#
# "A legend never sits on top of data" is unanimous across the six benchmarked
# papers, and the five-entry cell-state key printed straight through the orange
# cluster. It is also unnecessary: panel C's legend strip is directly above, in
# the same column, and already defines both the cell-state colours and the clone
# blocks. Repeating a key millimetres from the original is the thing RIDDLER Fig 3a
# is criticised for.
#
# The two subclones are DIRECT-LABELLED in whitespace beside their clouds instead
# (CopyKAT Fig 3). Direct annotation of two well-separated groups needs no key at
# all, and it is the one form of on-data text the literature does use -- for
# biological groups, never for methods.
#
# Clone colours are C's own CLONE_COLORS2, so the same two clones are the same two
# colours in both panels. With direct labels the pair does not have to carry the
# identification on its own.

# max(umap2) is the wrong anchor: a handful of C2 nuclei sit up inside the C1
# cloud, so the C2 label was pinned to THEIR height and floated ~8 mm above its own
# cluster while "OPC-like", which has no such strays, hugged its points. The label
# should track the top of the cluster's BULK, so the y anchor is the 98th
# percentile and the x anchor is the median of only the points in the top decile --
# i.e. where the cluster actually is at the height the label sits.
.lab_pos <- function(d, grp_col, grp) {
  x <- d[get(grp_col) == grp]
  hi <- as.numeric(stats::quantile(x$umap2, 0.90, names = FALSE))
  data.table(umap1 = median(x[umap2 >= hi, umap1]),
             umap2 = as.numeric(stats::quantile(x$umap2, 0.98, names = FALSE)),
             lab   = grp)
}
.cl <- rbind(.lab_pos(.ud, "clone", "C1"), .lab_pos(.ud, "clone", "C2"))
# The state panel is direct-labelled too, and only for the two states that
# actually dominate (95.0% and 87.1% of the two clones). Labelling the left panel
# and leaving the right one bare made the reader hunt for a key that is not there;
# labelling both puts "C1" and "AC-like" on the same cluster, which is the whole
# comparison. The three minor states stay decoded by panel C's legend above.
.dom <- .ud[!is.na(state) & nzchar(state), .N, by = state][order(-N)][1:2, state]
.st <- rbindlist(lapply(.dom, function(k) .lab_pos(.ud[!is.na(state)], "state", k)))

.umap_base <- function(p, ttl)
  # alpha 0.9 / size 0.28, not 0.7 / 0.22. A deeper palette alone cannot fix panel
  # D: at alpha 0.7 the light clone composites to 1.85:1 against white and stops
  # being a colour. This is the half of the fix that lives in the rendering, and it
  # goes in .umap_base so BOTH UMAPs change together -- they must stay comparable.
  p + geom_point(size = 0.28, alpha = 0.9, show.legend = FALSE) +
    labs(title = ttl, x = "UMAP 1", y = NULL) +
    theme_cowplot(BASE_PT, font_family = FIG_FAMILY) +
    # PT_SUP, not PT_SMALL. Figure 3's only title-tier text (its panel C facet
    # strips) renders at 7.4 pt; these two were at 5.5, the size this manuscript
    # reserves for dense chromosome ticks and direct data labels. A subpanel title
    # must never be the smallest text on the page.
    theme(plot.title = element_text(size = PT_SUP, face = "plain", colour = "grey15",
                                    hjust = 0.5, margin = margin(b = 1)),
          axis.text = element_blank(), axis.ticks = element_blank(),
          axis.line = element_line(linewidth = 0.3, colour = "grey30"),
          # BASE_PT, not PT_SMALL. EVERY axis title in Figure 3 renders at 7.0;
          # "UMAP 1"/"UMAP 2" were the only axis titles in either figure at 5.5.
          axis.title = element_text(size = BASE_PT, colour = "grey40"),
          panel.grid = element_blank(), legend.position = "none",
          plot.margin = margin(1, 1, 1, 1), aspect.ratio = 1) +
    # Room for the direct labels. "OPC-like" is centred on a cluster that sits
    # against the right edge and was clipped to "OPC-lik"; ggplot clips at the
    # panel, so the loss is silent. The SAME expansion goes on both UMAPs because
    # they are the same embedding and their axes must stay comparable.
    scale_x_continuous(expand = expansion(mult = 0.16))

p_umap_clone <- .umap_base(
  ggplot(.ud[!is.na(clone)], aes(umap1, umap2,
                                 colour = factor(clone, levels = c("C1", "C2")))),
  "MintCNA subclone") +
  scale_colour_manual(values = CLONE_COLORS2) +
  geom_text(data = .cl, aes(umap1, umap2, label = lab), inherit.aes = FALSE,
            size = PT_SUP/.pt, fontface = "bold", colour = "grey15", vjust = -0.3) +
  scale_y_continuous(expand = expansion(mult = c(0.04, 0.13))) +
  labs(y = "UMAP 2")

p_umap_state <- .umap_base(
  ggplot(.ud[!is.na(state) & nzchar(state)],
         aes(umap1, umap2, colour = factor(state, levels = names(CELLSTATE_COLORS)))),
  "Cell state (RNA)") +
  scale_colour_manual(values = CELLSTATE_COLORS, drop = TRUE) +
  geom_text(data = .st, aes(umap1, umap2, label = lab), inherit.aes = FALSE,
            size = PT_SUP/.pt, fontface = "bold", colour = "grey15", vjust = -0.3) +
  scale_y_continuous(expand = expansion(mult = c(0.04, 0.13)))

cat(sprintf("panel D: %s, %d nuclei, 15-NN clone %.3f vs state %.3f (floor %.3f)\n",
            UMAP_SAMPLE, nrow(.ud), .us$knn_acc_clone, .us$knn_acc_state, .us$majority_floor))
