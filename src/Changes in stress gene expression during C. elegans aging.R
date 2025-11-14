library(tidyverse)
install.packages("kernlab")
library(kernlab)

library(dplyr)
library(palmerpenguins)
library(matrixStats)
library(ggplot2)

age1 <- read.csv("/Users/cmdb/qb_project/expression_df_old.csv")

colnames(age1)



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

ggplot(sum, aes(x = timepoint, y= mean_expression, color = gene, group = gene))+
  geom_line(size = 1) +
  geom_point(size =1)+
  labs(
    title ="Changes in stress gene expression during C. elegans aging",
    x = "Timepoint",
    y= "Mean expression"
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


