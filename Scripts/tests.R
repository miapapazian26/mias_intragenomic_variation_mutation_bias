#tests to ensure that the number of genes and samples match across objects 
stopifnot(length(gene_names) == nrow(codon_counts))         # 4945
stopifnot(length(gene_names) == nrow(phi_trace))            # 4945
stopifnot(ncol(phi_trace) == length(dM_trace[[1]]))         # 101 samples
stopifnot(all(colnames(codon_counts) %in% names(dM_trace[[1]])))  # codon match

