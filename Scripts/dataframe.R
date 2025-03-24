#creating a dataframe for mutation and selection traces 



#creating a data frame for easy access using the objects for the loglikelihood functions 

access.df <- data.frame(
  dM_trace = dM_trace,
  dE_trace = dE_trace,
  gene_names = gene_names,
  gene_index = gene_index,
  samples = 100
)

