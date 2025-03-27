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

#creating gene names 
gene_names <- getNames(genome, FALSE)

#creating gene indexes 
gene_index <- 1:length(gene_names)

#creating phi values 
phi_trace <- synth_trace #changing from synth_trace_mix to synth_trace
#reshaping into a tidy format
phi_long <- phi_trace %>%
  pivot_longer(
    cols = -id_names,
    names_to = "Sample",
    values_to = "Phi"
  ) %>%
  mutate(
    Sample = as.integer(Sample) + 1  # convert from "0" to 1-based integer
  )

#getting codon names 
codon_names <- csp_codons
#getting rid of the loop and using csp_codons

#getting codon counts and gene ids
codon_names_vec <- codon_names$Codon
codon_counts <- full_codon_counts[, codon_names_vec] #filtered 40

#full 64-codon count
full_codon_counts <- getCodonCounts(genome)

#matching codon names with vectors for delta m and delta eta 
names(delta_M) <- codon_names
names(delta_eta) <- codon_names

