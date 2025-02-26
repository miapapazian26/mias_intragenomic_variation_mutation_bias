LLikGenome <- function(gene_name, delta_mu, delta_eta, phi, codon_counts) {
  # compute log probabilities
  log_P <- delta_mu + delta_eta * log(phi)
  
  # convert to probabilities
  P <- exp(log_P)
  P <- P / sum(P)  # Normalize so probabilities sum to 1
  
  # compute log-likelihood
  LL <- sum(codon_counts * log(P))
  
  # return gene name and log-likelihood
  return(list(Gene = gene_name, LogLikelihood = LL))
}

# example usage
gene_name <- "GeneX"
delta_mu <- c(-0.5, 0.1, -0.3)  # Δμ for codons
delta_eta <- c(0.2, -0.1, 0.5)  # Δη for codons
phi <- 50                       # expression level
codon_counts <- c(10, 5, 20)     # codon counts

LLikGenome(gene_name, delta_mu, delta_eta, phi, codon_counts)

