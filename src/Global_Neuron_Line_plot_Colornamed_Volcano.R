library(tidyverse)
install.packages("kernlab")
library(kernlab)

library(dplyr)
library(palmerpenguins)
library(matrixStats)
library(ggplot2)

age1 <- read.csv("/Users/cmdb/qb_project/global_timepoint_single_cell.csv")
age2 <- read.csv("/Users/cmdb/qb_project/neuronal_timepoint_single_cell.csv")
neuron_data <- read.csv("/Users/cmdb/qb_project/neuron_stress_gene_expression_by_timepoint.csv")
colnames(age1)
colnames(age2)
colnames(neuron_data)

#To convert data from wide format so that the gene expression can be easily plotted.
graph_age <- age1 %>%
  pivot_longer(
    cols = -timepoint,
    names_to = "gene",
    values_to = "expression"
  )
#group by gene and expression and calcuate for the mean expression.
sum <- graph_age %>%
  group_by(timepoint, gene)%>%
  summarize(mean_expression = mean(expression, na.rm = TRUE))
#To organize the timepoint from day 1 to day 15 in order.
sum$timepoint <- factor(
  sum$timepoint,
  levels = c("d1", "d3", "d5", "d8", "d11", "d15"),
  ordered = TRUE
)

ggplot(sum, aes(x = timepoint, y= mean_expression, color = gene, group = gene))+
  geom_line(size = 1) +
  geom_point(size =1)+
  facet_wrap(~ gene, scales = "free_y") +
  labs(
    title ="Changes in stress gene expression during C. elegans aging",
    x = "Timepoint",
    y= "Mean expression"
  )

#Now the y-axis is no longer free 
ggplot(sum, aes(x = timepoint, y= mean_expression, color = gene, group = gene))+
  geom_line(size = 1) +
  geom_point(size =1)+
  facet_wrap(~ gene) +
  labs(
    title ="Changes in stress gene expression during C. elegans aging",
    x = "Timepoint",
    y= "Mean expression"
  )


ggplot(sum, aes(x = timepoint, y= mean_expression, color = gene, group = gene))+
  geom_line(size = 1) +
  geom_point(size =1)+
  labs(
    color = "Genes",
    title ="Changes in stress gene expression during C. elegans aging",
    x = "Timepoint",
    y= "Mean expression"
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 16, face = "bold"),
    legend.key.size = unit(0.3, "lines"),
    legend.spacing = unit(0.1, "lines"),
    legend.text = element_text(size = 8)
  )

  




#group by gene and expression and calculate for the mean expression.
graph_age$timepoint <- factor(
  graph_age$timepoint,
  levels = c("d1", "d3", "d5", "d8", "d11", "d15"),
  ordered = TRUE
)
#group the gene vs timepoint to find specific gene expression at specific timepoints
aging_k <- graph_age %>%
  group_by(gene, timepoint) %>%
  summarize(expression = mean(expression, na.rm= TRUE))

#Was running into errors and used google to learn that the data needs to be in wide format for one row per gene, and one column per timepoint
wide <- aging_k %>%
  pivot_wider(names_from = timepoint, values_from = expression) %>%
  column_to_rownames("gene")
#Using google, incorporated the chronological order of timepoint from before
wide <- wide[, levels(graph_age$timepoint)]
#Was running into error, and troubleshooted with ChatGPT to learn that a numeric matrix need to be generated for following functions.
wide_matrix <- as.matrix(wide)
mode(wide_matrix) <- "numeric"
#Learned from week7, on finding the standard deviation across timepoints for each genes
standarddevi <- rowSds(wide_matrix)
filtered <-wide_matrix[standarddevi >0.01,]
num_clusters <- min(3, nrow(filtered))

#Learned from week7, on generating the kmeans heatplot
set.seed(42)
result <- kmeans(scale(as.matrix(filtered)),centers =num_clusters, nstart=100)

labels_kmeans <-result$cluster
ordered_labels <-labels_kmeans[order(labels_kmeans)]
fil <- filtered[order(labels_kmeans),]

#The heatmap generated shows the changes of gene expression changing overtime, now with the genes of similar expressions grouped together
heatmap(fil, Rowv=NA, Colv=NA, RowSideColors=RColorBrewer::brewer.pal(n= num_clusters,name= "Paired")[ordered_labels], ylab="Gene")




#I wanted to create a section of code just for the Volcano plot 
# Filter for the two timepoints you want to compare
volcano_input <- graph_age %>%
  #To compare the two oppsosite timepoints of d1 and d11
  filter(timepoint %in% c("d1", "d11"))

# Calculate fold change and p-value per gene
i <- volcano_input %>%
  #To calculate the fold change and p-value needed for the plot for each gene individually.
  group_by(gene) %>%
  summarize(
    #Calculating the average expression of all replicates at day 11, doing +1 to avoid taking the log of zero. (Searched help using google and chatGPT)
    log2_fold_change = log2(mean(expression[timepoint == "d11"]) + 1) -
      log2(mean(expression[timepoint == "d1"]) + 1),
    p_value = t.test(expression[timepoint == "d11"], 
                     expression[timepoint == "d1"])$p.value
  ) %>%
  # Creating the y-axis for volcano plot
  mutate(
    adjusted_p = p.adjust(p_value, method = "BH"),
    neg_log10_adjp = -log10(adjusted_p),                     
    significance = ifelse(p_value < 0.05 & abs(log2_fold_change) > 1,
                          "Significant", "Not Significant")
  )

#Creating the volcano plot 
ggplot(i, aes(x = log2_fold_change, y = neg_log10_adjp, color = gene)) +
  geom_point(size = 3) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  labs(
    
    title = "Volcano Plot: Gene Expression d11 vs d1",
    x = "log2 Fold Change",
    y = "-log10(adjusted p-value)",
    color = "Gene"
  ) +
  coord_cartesian(ylim = c(0, 250)) +

  theme_bw(base_size = 14) +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 16, face = "bold"),
    legend.key.size = unit(0.3, "lines"),
    legend.spacing = unit(0.1, "lines"),
    legend.text = element_text(size = 8)
  )




# Pivot to long format
neuron_long <- neuron_data %>%
  pivot_longer(
    cols = -timepoint,
    names_to = "gene",
    values_to = "expression"
  )%>%
  group_by(timepoint, gene) %>%          
  summarize(mean_expression = mean(expression, na.rm = TRUE)) %>% 
  ungroup()

# Order timepoints
neuron_long$timepoint <- factor(
  neuron_long$timepoint,
  levels = c("d1", "d3", "d5", "d8", "d11", "d15"),
  ordered = TRUE
)

# Line plot
ggplot(neuron_long, aes(x = timepoint, y = mean_expression, color = gene, group = gene)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(
    title = "Neuron-specific Stress Gene Expression Over Time",
    x = "Timepoint",
    y = "Mean Expression",
    color = "Gene"
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 16, face = "bold"),
    legend.key.size = unit(0.3, "lines"),
    legend.spacing = unit(0.1, "lines"),
    legend.text = element_text(size = 8)
  )





high_genes <- sum %>%
  filter(mean_expression > 2
         ) %>%
  pull(gene)  

high_genes


high_genesdd <- neuron_long  %>%
  filter(mean_expression > 3
  ) %>%
  pull(gene)  

high_genesdd
