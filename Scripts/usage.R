#checking and usage 
head(dM_trace) #looking delta M 
head(dE_trace) #looking delta eta 
head(codon_counts) #looking at codon counts 
head(phi_trace) #looking at phi values 
head(codon_names) #looking at codon names
head(gene_index) #looking at gene indexes **can be dropped using seq_along(gene_names)
head(gene_names) #looking at gene names

str(dM_trace) #looking delta M 
str(dE_trace) #looking delta eta 
str(codon_counts) #looking at codon counts 
str(phi_trace) #looking at phi values 
str(codon_names) #looking at codon names
str(gene_index) #looking at gene indexes
str(gene_names) #looking at gene names 

LLikGene(gene_name, dM_trace, dE_trace, phi, codon_counts)
##do i have to call the function like this instead?
LLikGene(gene, dM_trace[[s]][names(counts)], dE_trace[[s]][names(counts)], phi, counts)

LLikAA(dM_trace, dE_trace, phi, codon_counts) 


