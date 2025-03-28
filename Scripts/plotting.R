#script for plotting the worst fit genes
#loading necessary libraries

library(ggplot2)

#take the ten worst genes, aka the genes with the lowest coefficients
worst_genes <- head(gene_effects, 10)

#create a data frame from the worst_genes object
worst_df <- data.frame(
  Gene = names(worst_genes),
  Coefficient = as.numeric(worst_genes)
)

#plot the worst_df
ggplot(worst_df, aes(x = reorder(Gene, Coefficient), y = Coefficient)) +
  geom_col(fill = "tomato") +
  coord_flip() +  #flips plot so we can read gene names
  labs(
    title = "10 Worst-Fit Genes (Lowest Log-Likelihood)",
    x = "Gene",
    y = "Log-Likelihood Deviation (from reference)"
  ) +
  theme_minimal()
