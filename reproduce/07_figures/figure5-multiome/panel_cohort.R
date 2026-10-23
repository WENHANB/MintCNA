#!/usr/bin/env Rscript
# Panel E -- does the subclone / cell-state coincidence hold across the cohort?
#
# Panels C and D show one tumour. Without this panel the figure makes a claim about
# IDH-mutant glioma from a sample that was chosen BECAUSE the claim holds in it
# (selection rule 2 in the handoff: "the clones must be lineage-pure"). This panel
# is the cohort: every tumour in which MintCNA returned two or more subclones,
# placed by (x) whether the split is sequencing depth and (y) how strongly it
# follows the RNA cell-state boundary. Mark size is the copy-number difference
# between the two clones, so a reader can see which splits carry a real CNA block.
#
# The classification is a RULE, applied to all 20, not a curated list:
#   genuine split      depth AUC < 0.75  AND  >= 100 bins differ by > 0.3 state
#     ..follows state  dominant cell state differs between the two clones   (6)
#     ..within a state dominant cell state is the same                       (1)
#   not established    everything else: depth-driven (AUC > 0.96 in seven), or a
#                      split with no contiguous CNA difference behind it    (13)
#
# P59T2 is the panel's most important point. A 187-bin (93 Mb) subclonal block
# partitioning nuclei WITHIN a single AC-like state is what rules out the reading
# that the AG-DCAE clustering merely tracks chromatin accessibility -- if it did,
# a within-state split could not exist. It is labelled, not hidden.
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(cowplot); library(ggrepel) })

ch <- fread(file.path(DATA, "fig5_clone_state_all_samples.tsv"))
ch[, genuine := depth_AUC < 0.75 & bins_differ >= 100]
ch[, class := fifelse(!genuine, "not established",
              fifelse(states_differ, "follows cell state", "within one state"))]
ch[, class := factor(class, levels = c("follows cell state", "within one state", "not established"))]
cat(sprintf("panel E (cohort): %d tumours; genuine splits %d, of which follow state %d, within one state %d; not established %d\n",
            nrow(ch), sum(ch$genuine), ch[class == "follows cell state", .N],
            ch[class == "within one state", .N], ch[class == "not established", .N]))

# named: the figure's exemplar, the strongest tumour by every number, the counterexample
NAME <- c("P111T2", "P105T1", "P59T2")
lab <- ch[sample %in% NAME]

# clone C2 teal is deliberately reused for "follows cell state": that is the class
# the whole left column illustrates. The counterexample takes the clone C1 neutral.
# "Not established" is hollow so it reads as absence, not as a third biology.
CLS <- c(`follows cell state` = CLONE_COLORS2[["C2"]],
         `within one state`   = CLONE_COLORS2[["C1"]],
         `not established`    = "grey55")

# ONE point layer, fill AND outline both mapped, so every class reaches the legend
# with a glyph. Drawing the hollow class as a separate layer with a literal
# fill = "white" is what left "not established" as a key with no mark -- the same
# bug the motif volcano had. Hollow is encoded in the fill scale as white.
OUTL <- c(`follows cell state` = "grey25", `within one state` = "grey25",
          `not established`    = "grey55")
FILL <- c(CLS[1:2], `not established` = "white")
p_cohort <- ggplot(ch, aes(depth_AUC, cramerV)) +
  geom_vline(xintercept = 0.75, linewidth = 0.22, linetype = "22", colour = "grey70") +
  # ROTATED along the rule and placed on its LEFT, at mid-height: the only band
  # near x = 0.75 with no points is y 0.35-0.75. Written horizontally at the top
  # it clipped at the panel edge and ran into the P105T1 label.
  # centred higher and shorter than the first draft: at y = 0.55 the 24-character
  # run reached down into the stacked legend and overlapped its top row by 2 mm.
  annotate("text", x = 0.745, y = 0.62, angle = 90, hjust = 0.5, vjust = 0,
           size = PT_SMALL/.pt, colour = "grey45", label = "depth-driven split") +
  geom_point(aes(size = bins_differ, fill = class, colour = class),
             shape = 21, stroke = 0.35) +
  # P59T2 is labelled to the LEFT of its point: it sits at the bottom-left corner
  # right above the stacked legend, and the repel's default up-right placement put
  # "P59T2" 1.6 mm onto the size-key title. The two top labels keep the default.
  geom_text_repel(data = lab[sample != "P59T2"], aes(label = sample), size = PT_SMALL/.pt,
                  colour = "grey15", segment.size = 0.16, segment.colour = "grey65",
                  min.segment.length = 0, box.padding = 0.35, point.padding = 0.25,
                  seed = 4, show.legend = FALSE) +
  geom_text_repel(data = lab[sample == "P59T2"], aes(label = sample), size = PT_SMALL/.pt,
                  colour = "grey15", segment.size = 0.16, segment.colour = "grey65",
                  min.segment.length = 0, nudge_x = -0.005, nudge_y = 0.09,
                  box.padding = 0.25, point.padding = 0.25, seed = 4, show.legend = FALSE) +
  scale_fill_manual(values = FILL, breaks = names(FILL), name = NULL, drop = FALSE) +
  scale_colour_manual(values = OUTL, breaks = names(OUTL), name = NULL, drop = FALSE) +
  guides(fill = guide_legend(order = 1, override.aes = list(size = 2.2, shape = 21)),
         colour = guide_legend(order = 1),
         # the size key runs HORIZONTALLY under the class key; two stacked vertical
         # keys were 43 % of the panel height, and at the right edge the size key
         # clipped by 4.7 mm and sat on the depth-artefact points
         size = guide_legend(order = 2, direction = "horizontal", nrow = 1,
                             title.position = "left", label.position = "bottom",
                             override.aes = list(shape = 21, fill = "grey85", colour = "grey40"))) +
  scale_size_area(max_size = 4.2, breaks = c(100, 400, 800), name = "bins differing") +
  scale_x_continuous(limits = c(0.48, 1.0), breaks = c(0.5, 0.75, 1.0),
                     expand = expansion(mult = c(0.02, 0.04))) +
  scale_y_continuous(limits = c(0, 1.02), breaks = seq(0, 1, 0.25), expand = expansion(mult = c(0.02, 0.02))) +
  labs(x = "Depth AUC  (fragments vs subclone)",
       y = "Subclone × cell state  (Cramér's V)") +
  theme_cowplot(BASE_PT, font_family = FIG_FAMILY) +
  theme(axis.text = element_text(size = BASE_PT - 1, colour = "grey20"),
        axis.title = element_text(size = BASE_PT),
        axis.line = element_line(linewidth = 0.3, colour = "grey30"),
        axis.ticks = element_line(linewidth = 0.3, colour = "grey30"),
        panel.grid = element_blank(),
        # SQUARE: both axes are unit-scale association measures, like B and F.
        aspect.ratio = 1,
        # lower-middle: the coincident tumours sit top-left, the depth artefacts
        # hug the right edge, the counterexample is bottom-left. The band between
        # x 0.55-0.9 and y 0-0.3 is the one empty region.
        legend.position = "inside", legend.position.inside = c(0.24, 0.02),
        legend.justification = c(0, 0), legend.direction = "vertical",
        legend.box = "vertical", legend.box.just = "left",
        legend.background = element_blank(),
        legend.margin = margin(0, 0, 0, 0), legend.spacing.y = unit(0.3, "mm"),
        legend.spacing.x = unit(1.2, "mm"), legend.box.spacing = unit(0, "mm"),
        legend.title = element_text(size = BASE_PT - 1, colour = "grey30"),
        legend.text = element_text(size = BASE_PT - 1), legend.key.size = unit(2.6, "mm"),
        plot.margin = margin(2, 2, 1, 2))

# Supplementary since 2026-09-06: the main figure follows find -> compare -> interpret.
if (!exists("FIG5_ASSEMBLING")) {
  cairo_pdf(file.path(OUT, "figS_cohort.pdf"), width = 70/25.4, height = 70/25.4, family = FIG_FAMILY)
  print(p_cohort); dev.off(); cat("wrote figS_cohort.pdf\n")
}
