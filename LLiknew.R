#getting data 
mutationTrace <- trace$getCodonSpecificParameterTrace(0) #mu
delta_M <- mutationTrace[[1]][[1]]
selectionTrace <- trace$getCodonSpecificParameterTrace(1) #eta
delta_eta <- selectionTrace[[1]][[1]]

dM_trace <- delta_M
dE_trace <- delta_eta
gene_index <- match(gene_names, gene_names)
gene_names <- getNames(genome, FALSE)

phi_trace <- synth_trace_mix


access.df <- data.frame(
  dM_trace = dM_trace,
  dE_trace = dE_trace,
  gene_names = gene_names,
  gene_index = gene_index,
  samples = 100
)


for(i in gene_index) {
  phi <- phi_trace[i, ]
  cc <- genome$getCodonCountsPerGene(i)
  
  for(s in samples) {
    dM <- dM_trace[mix][[s]]
    dE <- dE_trace[mix][[s]]
    p <- phi[[s]]
llik[[index,s]]=LLik(cc, dM, dE, p)
  }
}


