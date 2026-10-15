library(readxl)
library(tidyr)
library(dplyr)
library(ggplot2)

# Read data
df <- read_excel("<MINTCNA_ROOT>/Simulation/cnasim/Computational Time.xlsx")

# Reshape to long format
df_long <- df %>%
  pivot_longer(cols = -Dataset, names_to = "Method", values_to = "Time") %>%
  mutate(Time = as.numeric(Time)) %>%
  filter(!is.na(Time)) %>%
  filter(Dataset != "Breast cancer (T10)")

# Set method order consistent with other figures
df_long$Method <- factor(df_long$Method, levels = c(
  "MintCNA", "MintCNA(DNA)", "SCOPE", "HMMcopy", 
  "SCONCE", "SCCNV", "rcCAE", "SCICoNE"
))

# Set dataset order
df_long$Dataset <- factor(df_long$Dataset, levels = c(
  "Simulation(N=100)", "Colorectal Cancer (scTrio-seq2)"
))

# x-axis tick labels: keep "MintCNA(DNA)" as the data key but render it as MintCNA with a
# superscript DNA (plotmath), aligned to the Method factor levels above.
method_axis_labels <- parse(text = c("MintCNA", "MintCNA^DNA", "SCOPE", "HMMcopy",
                                     "SCONCE", "SCCNV", "rcCAE", "SCICoNE"))

# Color palette matching your other figures
method_colors <- c(
  "MintCNA"      = "#1B9E77",
  "MintCNA(DNA)" = "#D95F02",
  "SCOPE"        = "#7570B3",
  "HMMcopy"      = "#E7298A",
  "SCONCE"       = "#66C2A5",
  "SCCNV"        = "#E6AB02",
  "rcCAE"        = "#A6761D",
  "SCICoNE"      = "#666666"
)

p <- ggplot(df_long, aes(x = Method, y = Time, fill = Dataset)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  scale_y_log10(
    breaks = c(1, 2, 5, 10, 20, 50, 100, 200, 500, 1000),
    labels = c("1", "2", "5", "10", "20", "50", "100", "200", "500", "1000")
  ) +
  scale_x_discrete(labels = method_axis_labels) +
  scale_fill_manual(values = c(
    "Simulation(N=100)" = "#3C5488",
    "Colorectal Cancer (scTrio-seq2)" = "#E64B35"
  )) +
  labs(
    x = NULL,
    y = "CPU time (minutes, log scale)",
    fill = "Dataset"
  ) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 13),
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(
  "<MINTCNA_ROOT>/Simulation/results/computational_time_barplot.pdf",
  p, width = 10, height = 5
)

