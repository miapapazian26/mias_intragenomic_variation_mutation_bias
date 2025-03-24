#creating a data frame for easy access using the objects for the loglikelihood function 
access.df <- data.frame(
  dM_trace = dM_trace,
  dE_trace = dE_trace,
  gene_names = gene_names,
  gene_index = gene_index,
  samples = 100
)