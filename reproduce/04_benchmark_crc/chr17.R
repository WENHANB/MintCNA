load("<MINTCNA_ROOT>/Result/CRC/Benchmark/MultiCNA_AE_197_d250.RData")

################### 17p #####################
## define 17p region in hg19
p17_end <- 25248165

## 1) find bins in 17p
idx_17p <- which(
  ref_multicna_ae$chr == "chr17" &
    ref_multicna_ae$start <= p17_end
)

length(idx_17p)
head(ref_multicna_ae[idx_17p, ])

## 2) subset CNA matrix for 17p
cna_17p <- cna_multicna_ae[idx_17p, , drop = FALSE]

## 3) summarize 17p status for each cell
status_17p <- data.frame(
  cell = colnames(cna_multicna_ae),
  
  prop_loss = apply(cna_17p, 2, function(x) mean(x == 1, na.rm = TRUE)),
  prop_neutral = apply(cna_17p, 2, function(x) mean(x == 2, na.rm = TRUE)),
  prop_gain = apply(cna_17p, 2, function(x) mean(x == 3, na.rm = TRUE)),
  
  status_17p = apply(cna_17p, 2, function(x) {
    x <- x[!is.na(x)]
    if (length(x) == 0) return(NA)
    
    if (all(x == 2)) {
      "neutral"
    } else if (mean(x == 1) > 0.5) {
      "loss"
    } else if (mean(x == 3) > 0.5) {
      "gain"
    } else {
      "mixed"
    }
  }),
  
  row.names = NULL
)

status_17p_loss <- data.frame(
  cell = colnames(cna_multicna_ae),
  loss_17p = apply(cna_17p, 2, function(x) mean(x == 1, na.rm = TRUE) > 0.5),
  prop_loss_17p = apply(cna_17p, 2, function(x) mean(x == 1, na.rm = TRUE)),
  row.names = NULL
)

library(dplyr)
library(ggplot2)
library(stringr)
library(scales)

## ---- 1) extract site from cell name ----
plot_df <- status_17p_loss %>%
  mutate(
    site = str_extract(cell, "LN|PT|ML|MP")
  ) %>%
  filter(!is.na(site))

## ---- 2) summarize percentage of 17p loss by site ----
site_17p <- plot_df %>%
  group_by(site) %>%
  summarise(
    n_cell = n(),
    n_loss = sum(loss_17p, na.rm = TRUE),
    pct_loss = 100 * n_loss / n_cell,
    .groups = "drop"
  ) %>%
  mutate(
    site = factor(site, levels = c("PT", "LN", "ML", "MP"))
  )

p <- ggplot(site_17p, aes(x = site, y = pct_loss)) +
  geom_col(width = 0.7, fill = "#FCAE91") +
  geom_text(
    aes(label = sprintf("%.1f%%", pct_loss)),
    vjust = -0.5,
    size = 5
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 20),
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0.08))
  ) +
  labs(
    x = "Tumor site",
    y = "Cells with 17p loss (%)"
  ) +
  theme_classic(base_size = 16) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14, color = "black"),
    axis.line = element_line(linewidth = 0.8),
    axis.ticks = element_line(linewidth = 0.8),
    plot.margin = margin(10, 15, 10, 10)
  )

p

ggsave("<MINTCNA_ROOT>/Result/CRC/Segmentation/Evaluation/Fig_17p_loss_by_site.pdf", 
       p, width = 5.5, height = 4.5)

############# 17 q #############
library(dplyr)
library(ggplot2)
library(stringr)
library(scales)

## -------------------------------
## 1) define 17q bins (hg19)
## -------------------------------
q17_start <- 25248166

idx_17q <- which(
  ref_multicna_ae$chr == "chr17" &
    ref_multicna_ae$end >= q17_start
)

cna_17q <- cna_multicna_ae[idx_17q, , drop = FALSE]

## -------------------------------
## 2) define 17q loss per cell
##    here: >50% of 17q bins are loss
## -------------------------------
status_17q_loss <- data.frame(
  cell = colnames(cna_multicna_ae),
  prop_loss_17q = apply(cna_17q, 2, function(x) mean(x == 1, na.rm = TRUE)),
  loss_17q = apply(cna_17q, 2, function(x) mean(x == 1, na.rm = TRUE) > 0.5),
  row.names = NULL
)

## -------------------------------
## 3) extract site from cell name
## -------------------------------
plot_df_17q <- status_17q_loss %>%
  mutate(site = str_extract(cell, "LN|PT|ML|MP")) %>%
  filter(!is.na(site))

## -------------------------------
## 4) summarize percentage by site
## -------------------------------
site_17q <- plot_df_17q %>%
  group_by(site) %>%
  summarise(
    n_cell = n(),
    n_loss = sum(loss_17q, na.rm = TRUE),
    pct_loss = 100 * n_loss / n_cell,
    .groups = "drop"
  ) %>%
  mutate(
    site = factor(site, levels = c("PT", "LN", "ML", "MP")),
    site_label = paste0(site, "\n(n=", n_cell, ")")
  )

p_17q <- ggplot(site_17q, aes(x = site, y = pct_loss)) +
  geom_col(width = 0.7, fill = "#FCAE91") +
  geom_text(
    aes(label = sprintf("%.1f%%", pct_loss)),
    vjust = -0.5,
    size = 5
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 20),
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0.08))
  ) +
  labs(
    x = "Tumor site",
    y = "Cells with 17q loss (%)"
  ) +
  theme_classic(base_size = 16) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14, color = "black"),
    axis.line = element_line(linewidth = 0.8),
    axis.ticks = element_line(linewidth = 0.8),
    plot.margin = margin(10, 15, 10, 10)
  )


ggsave("<MINTCNA_ROOT>/Result/CRC/Segmentation/Evaluation/Fig_17q_loss_by_site.pdf", 
       p_17q, width = 5.5, height = 4.5)
