#tests to ensure that the number of genes and samples match across objects 
stopifnot(length(gene_names) == nrow(codon_counts))         # 4945
stopifnot(length(gene_names) == nrow(phi_trace))            # 4945
stopifnot(ncol(phi_trace) == length(dM_trace[[1]]))         # 101 samples
stopifnot(all(colnames(codon_counts) %in% names(dM_trace[[1]])))  # codon match

#how to manually test the gene likelihood function
#extract one sample, then call the function
#dM and dE are named for SINGULAR SAMPLES using dM trace and dE trace
test_gene <- gene_names[1]
test_sample <- 1
counts <- codon_counts[test_gene, ]
phi <- phi_trace[1, test_sample]
dM <- dM_trace[[test_sample]][names(counts)]
dE <- dE_trace[[test_sample]][names(counts)]

# call the log-likelihood function using the sample-specific values
LLikGene(test_gene, dM, dE, phi, counts)
