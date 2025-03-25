#checking and usage 
head(delta_M) #looking delta M **fix name
head(delta_eta) #looking delta eta **fix name
head(codon_counts) #looking at codon counts 
head(phi_trace) #looking at phi values 
head(codon_names) #looking at codon names
head(gene_index) #looking at gene indexes
head(gene_names) #looking at gene names

LLikGene(gene_name, delta_M, delta_eta, phi, codon_counts)

LLikAA(delta_M, delta_eta, phi, codon_counts) 


