#tests to ensure that the number of genes and samples match across objects 
stopifnot(length(gene_names) == nrow(codon_counts))         # 4945
stopifnot(length(gene_names) == nrow(phi_trace))            # 4945
stopifnot(ncol(phi_trace) == length(dM_trace[[1]]))         # 101 samples
stopifnot(all(colnames(codon_counts) %in% names(dM_trace[[1]])))  # codon match

#how to manually test the gene likelihood function
#extract one sample, then call the function
#dM and dE are named for SINGULAR SAMPLES using dM trace and dE trace


# call the log-likelihood function using the sample-specific values
test_gene <- gene_names[1]
test_counts <- codon_counts[test_gene, ]
test_phi <- phi_trace[1, 10]
test_dM <- sapply(dM_trace, `[`, 10)[names(test_counts)]
test_dE <- sapply(dE_trace, `[`, 10)[names(test_counts)]
#try running the function
LLikGene(test_gene, test_dM, test_dE, test_phi, test_counts)

#check the first few gene names and rownames 
head(gene_names)
head(rownames(codon_counts))
#are they the same? 
all(gene_names == rownames(codon_counts)) #should return TRUE


#check names of dM_trace and length of each entry
names(dM_trace)  # should be codons like "GCA", "GCC", ...
length(dM_trace[[1]])  # should be 101 (number of samples)
#confirm all entries are length 101
all(sapply(dM_trace, length) == 101)  # should be TRUE

#check names of dE_trace and length of each entry
names(dE_trace)  # should be codons like "GCA", "GCC", ...
length(dE_trace[[1]])  # should be 101 (number of samples)
#confirm all entries are length 101
all(sapply(dE_trace, length) == 101)  # should be TRUE

#checking phi_trace is a 4945 x 101 matrix 
dim(phi_trace)  # should return c(4945, 101)
