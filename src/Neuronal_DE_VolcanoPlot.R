library(ggplot2)
library(ggrepel)
library(dplyr)
library(readr)

de_results <- read_csv("/Users/cmdb/QB Project/neuronal_DE_stress_d11_vs_d1_scanpy.csv")

# Ensure numeric columns
de_results <- de_results %>%
  mutate(
    logfoldchanges = as.numeric(logfoldchanges),
    pvals_adj = as.numeric(pvals_adj)
  )

# Filter extreme outliers
MAX_LOGFC <- 10  # maximum absolute log2 fold change
de_results_filtered <- de_results %>%
  filter(abs(logfoldchanges) <= MAX_LOGFC)

# png saving volcano plot
png(
  filename = "/Users/cmdb/QB Project/neuronal_volcano_d11_vs_d1_2.png",
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
  geom_point(aes(color = significance), alpha = 0.8, size = 3) +
  scale_color_manual(values = c(
    "Not significant" = "grey30",
    "FC only" = "grey30",
    "P only" = "grey30",
    "Significant" = "red2"
  )) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_text_repel(
    data = label_genes,
    aes(label = names),
    size = 3.5,
    box.padding = 0.5,
    max.overlaps = 20,
    segment.color = "black",
    segment.size = 0.5
  ) +
  labs(
    title = "Neuronal Stress Genes DE: d11 vs d1",
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
