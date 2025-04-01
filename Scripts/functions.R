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
#FUNCTION: compute one LL value per gene per sample based on amino acid group level
##
LLikGene <- function(gene_name, dM, dE, phi, codon_counts, aa_by_codon) {
  # ensure all inputs are properly named
  names(dM) <- names(dE) <- names(codon_counts)
  
  # group codons by their associated amino acid
  codon_groups <- split(names(codon_counts), aa_by_codon[names(codon_counts)])
  
  LL_total <- 0
  
  for (aa in names(codon_groups)) {
    codons <- codon_groups[[aa]]
    counts <- codon_counts[codons]
    if (sum(counts) == 0) next
    
    dM_aa <- dM[codons]
    dE_aa <- dE[codons]
    
    # add amino acid–level log-likelihood
    LL_total <- LL_total + LLikAA(dM_aa, dE_aa, phi, counts)
  }
  
  return(list(Gene = gene_name, LogLikelihood = LL_total))
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
loglik_list <- list()
row_index <- 1

for (gene in gene_names) {
  # get codon counts for this gene
  counts <- as.numeric(codon_counts[gene, , drop = TRUE])
  gene_codons <- colnames(codon_counts)
  names(counts) <- gene_codons
  
  for (s in 1:100) {
    # get φ (synthesis rate) for this gene and sample
    phi <- synth_trace %>%
      filter(id_names == gene) %>%
      pull(as.character(s - 1)) %>%
      as.numeric()
    
    # get mutation and selection params for this gene/sample
    dM <- csp_long %>%
      filter(Codon %in% gene_codons, Sample == s) %>%
      arrange(match(Codon, gene_codons)) %>%
      pull(Mutation)
    
    dE <- csp_long %>%
      filter(Codon %in% gene_codons, Sample == s) %>%
      arrange(match(Codon, gene_codons)) %>%
      pull(Selection)
    
    # compute log-likelihood
    ll_result <- LLikGene(gene, dM, dE, phi, counts)
    
    # store result
    loglik_list[[row_index]] <- data.frame(
      Gene = gene,
      Sample = s,
      LogLikelihood = ll_result$LogLikelihood
    )
    
    row_index <- row_index + 1
  }
}

# Combine results into a tibble
loglik_df <- dplyr::bind_rows(loglik_list)

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

