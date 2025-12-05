library(ggplot2)
library(dplyr)
library(ggrepel)
library(readr)

de_results <- read_csv("/Users/cmdb/QB Project/global_DE_d11_vs_d1_scanpy.csv")

# filter
MAX_LOGFC <- 10  # maximum absolute log2 fold change
de_results_filtered <- de_results %>%
  filter(abs(logfoldchanges) <= MAX_LOGFC)

# png saving volcano plot 
png(
  filename = "/Users/cmdb/QB Project/global_volcano_d11_vs_d1_2.png",
  width = 2500,
  height = 2000,
  res = 300
)

# define significance categories for coloring
de_results_filtered <- de_results_filtered %>%
  mutate(
    negLogP = -log10(pvals_adj),
    significance = case_when(
      pvals_adj < 0.05 & abs(logfoldchanges) > 1 ~ "Significant",
      pvals_adj < 0.05 ~ "P only",
      abs(logfoldchanges) > 1 ~ "FC only",
      TRUE ~ "Not significant"
    )
  )

# select genes to label (highly significant)
label_genes <- de_results_filtered %>%
  filter(pvals_adj < 0.05 & abs(logfoldchanges) > 1.5)

# Volcano Plot
ggplot(de_results_filtered, aes(x = logfoldchanges, y = negLogP)) +
  geom_point(aes(color = significance), alpha = 0.8, size = 3) +  # points colored by significance
  scale_color_manual(values = c(
    "Not significant" = "grey30",
    "FC only" = "grey30",
    "P only" = "grey30",
    "Significant" = "red2"
  )) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +  # fold change cutoff lines
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +  # p val cutoff line
  geom_text_repel(
    data = label_genes,  # only label highly significant genes
    aes(label = names),
    size = 3.5,
    box.padding = 0.5,
    max.overlaps = 20,
    segment.color = "black",  # connectors
    segment.size = 0.5
  ) +
  labs(
    title = "Global Stress Genes DE: d11 vs d1",
    x = bquote(~Log[2]~ 'fold change (d11 vs d1)'),
    y = bquote(~-Log[10]~'adjusted P'),
    color = "Significance"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

dev.off()
