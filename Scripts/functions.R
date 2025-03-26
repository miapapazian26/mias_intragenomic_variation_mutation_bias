##
#FUNCTION: compute log-likelihood for one gene
##
LLikGene <- function(gene_name, dM_trace, dE_trace, phi, codon_counts) {
  # compute log probabilities using the correct sign and exponent
  log_P <- -dM_trace - dE_trace * phi
  
  # convert to probabilities
  P <- exp(log_P)
  P <- P / sum(P)  # normalize (ensure probabilities sum to 1)
  
  # compute log-likelihood
  LL <- dmultinom(x = codon_counts, prob = P, log = TRUE)
  
  # return gene name and loglikelihood
  return(list(Gene = gene_name, LogLikelihood = LL))
}

##
#FUNCTION: logliklihood amino acid function 
##
LLikAA <- function(delta_M, delta_eta, phi, codon_counts) {
  log_P = delta_M - delta_eta * phi
  log_P = log_P - max(log_P)  
  p = exp(log_P)
  
  
  LL = dmultinom(x = codon_counts, prob = p, log = TRUE)
  return(LL)
}

##
#LOOP: log-likelihood loop analysis (original)
##
for(i in gene_index) {
  phi <- phi_trace[i, ]
  cc <- genome$getCodonCountsPerGene(i)
  
  for(s in samples) {
    dM <- dM_trace[mix][[s]]
    dE <- dE_trace[mix][[s]]
    p <- phi[[s]]
    llik[[index,s]] = LLik(cc, dM, dE, p)
  }
}

##
#LOOP: editted log-likelihood loop analysis
##
for (i in seq_along(gene_names)) {
  gene <- gene_names[i]
  counts <- codon_counts[gene, ]
  
  for (s in 1:ncol(phi_trace)) {
    phi <- phi_trace[i, s]
    dM <- dM_trace[[s]][names(counts)]
    dE <- dE_trace[[s]][names(counts)]
    
    result <- LLikGene(gene, dM, dE, phi, counts)
    log_likelihood_matrix[i, s] <- result$LogLikelihood
  }
}

##
#LOOP: codon usage analysis
##
aa_list <- as.character(aa.bar$data$AA)  # convert amino acids to characters
codon_usage_per_gene <- list()

for (gene in gene_names) {
  c_counts_gene <- codon_counts[gene, , drop = FALSE]  # keep as dataframe
  codon_usage_per_gene[[gene]] <- list()
  
  for (AA in aa_list) {
    codons <- AAToCodon(AA)  # get codons for this amino acid
    valid_codons <- intersect(codons, colnames(codon_counts))  # ensure valid codons
    
    if (length(valid_codons) > 0) {
      codon_counts_subset <- c_counts_gene[, valid_codons, drop = FALSE]
      codon_usage_per_gene[[gene]][[AA]] <- codon_counts_subset
    } else {
      codon_usage_per_gene[[gene]][[AA]] <- NULL
    }
  }
}


##
##CALLING
##
LLikGene(gene_name, dM_trace, dE_trace, phi, codon_counts)
##do i have to call the function like this instead?
dM <- dM_trace[[s]][names(counts)]     ##don't assign to name until right before calling loop. 
dE <- dE_trace[[s]][names(counts)]     ##these names are necessary for the loop, not the function
LLikGene(gene, dM, dE, phi, counts)

LLikAA(dM_trace, dE_trace, phi, codon_counts) 

