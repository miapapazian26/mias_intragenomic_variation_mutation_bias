##
#FUNCTION: linear model function
##
lm_fit <- lm(LogLikelihood ~ Gene + Sample, data = loglik_df)



##
# FUNCTION: compute log-likelihood for one gene
# takes one sample's worth of parameters and one gene's codon counts
# returns the log-likelihood for that gene under that sample
## 
LLikGene <- function(gene_name, dM, dE, phi, codon_counts) {
  # compute log probabilities using the correct sign and exponent
  log_P <- -dM - dE * phi
  
  # convert to probabilities
  P <- exp(log_P)
  P <- P / sum(P)  # normalize (ensure probabilities sum to 1)
  
  # compute log-likelihood
  LL <- dmultinom(x = codon_counts, prob = P, log = TRUE)
  
  # return gene name and loglikelihood
  return(list(Gene = gene_name, LogLikelihood = LL))
}

##
# FUNCTION: loglikelihood amino acid function 
##
LLikAA <- function(dM, dE, phi, codon_counts) {
  # compute log probabilities for one amino acid group
  log_P <- -dM - dE * phi
  log_P <- log_P - max(log_P)  # numerical stability
  
  # convert to probabilities
  p <- exp(log_P)
  p <- p / sum(p)
  
  # compute log-likelihood
  LL <- dmultinom(x = codon_counts, prob = p, log = TRUE)
  return(LL)
}

##
# LOOP: log-likelihood loop analysis (updated)
##
#create an empty list to store results
loglik_list <- list()

#start counter for rows
row_index <- 1

for (i in seq_along(gene_names)) {
  gene <- gene_names[i]
  counts <- codon_counts[gene, ]
  
  for (s in 1:ncol(phi_trace)) {
    phi <- phi_trace[i, s]
    
    dM <- sapply(dM_trace, `[`, s)[names(counts)]
    dE <- sapply(dE_trace, `[`, s)[names(counts)]
    
    ll_result <- LLikGene(gene, dM, dE, phi, counts)
    
    loglik_list[[row_index]] <- data.frame(
      Gene = gene,
      Sample = s,
      LogLikelihood = ll_result$LogLikelihood
    )
    
    row_index <- row_index + 1
  }
}
#combine all rows into one data frame
loglik_df <- do.call(rbind, loglik_list)


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
# LOOP: codon usage analysis
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
#CALLING
#need to specify a single sample for dM and dE
##

LLikGene(test_gene, dM, dE, phi, counts)

LLikAA(dM, dE, phi, codon_counts)

