#genome is: candida_tenuis.max.cds

#getting mutation trace 
mutationTrace <- trace$getCodonSpecificParameterTrace(0) #mu
delta_M <- mutationTrace[[1]]
#renaming mutation trace 
dM_trace <- delta_M

#getting selection trace
selectionTrace <- trace$getCodonSpecificParameterTrace(1) #eta
delta_eta <- selectionTrace[[1]]
#renaming selection trace 
dE_trace <- delta_eta

#getting codon names 
codon_names <- colnames(codon_counts)

#creating gene indexes 
gene_index <- match(gene_names, gene_names)

#creating gene names 
gene_names <- getNames(genome, FALSE)

#creating phi values 
phi_trace <- synth_trace_mix

#getting codon counts and gene ids
codon_counts <- getCodonCounts(genome)

