library(tidyverse)

age <- read.csv("/Users/cmdb/qb_project/expression_df.csv")


#To convert data from wide format so that the gene expression can be easily plotted.

graph_age <- age %>%
  
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

