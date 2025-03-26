#genome is: candida_tenuis.max.cds

#getting mutation trace 
mutationTrace <- trace$getCodonSpecificParameterTrace(0) #full trace
delta_M <- mutationTrace[[1]] #40 codon vectors 
#renaming mutation trace 
dM_trace <- delta_M 

#getting selection trace
selectionTrace <- trace$getCodonSpecificParameterTrace(1) #full trace
delta_eta <- selectionTrace[[1]] #40 codon vectors
#renaming selection trace 
dE_trace <- delta_eta 

#getting codon names 
codon_names <- c()
for (aa in aminoAcids()) {
  if (aa %in% c("M", "W", "X")) next
  codons <- AAToCodon(aa, TRUE)  # TRUE = excludes reference codon
  codon_names <- c(codon_names, codons)
}

#creating gene names 
gene_names <- getNames(genome, FALSE)

#creating gene indexes 
gene_index <- match(gene_names, gene_names)


#creating phi values 
phi_trace <- synth_trace_mix

#getting codon counts and gene ids
codon_counts <- codon_counts[, codon_names] #filtered 40

#full 64-codon count
full_codon_counts <- getCodonCounts(genome)

#matching codon names with vectors for delta m and delta eta 
names(delta_M) <- codon_names
names(delta_eta) <- codon_names

