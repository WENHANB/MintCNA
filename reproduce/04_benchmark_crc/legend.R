library(cowplot)
library(ggpubr)

method_order <- c(
  "MintCNA",
  "MintCNA(DNA)",
  "SCOPE",
  "HMMcopy",
  "SCONCE",
  "SCCNV",
  "rcCAE",
  "SCICoNE"
)


leg_df <- data.frame(
  method_disp = factor(method_order, levels = method_order),
  y = 1
)

p_leg <- ggplot(leg_df, aes(x = method_disp, y = y, fill = method_disp)) +
  geom_col() +
  scale_fill_manual(
    name = "Method:",
    values = color_map,
    breaks = method_order
  ) +
  theme_void(base_size = 11) +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.key.size = grid::unit(0.45, "cm")
  ) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

leg <- cowplot::get_legend(p_leg)

ggsave(
  "<MINTCNA_ROOT>/Result/CRC/Segmentation/Evaluation/wgs/legend_method.pdf",
  plot = ggpubr::as_ggplot(leg),
  width = 6.5,
  height = 1.4,
  units = "in",
  dpi = 600
)
