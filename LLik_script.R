LLikGenome <- function(gene_name, delta_M, delta_eta, phi, codon_counts) {
  # compute log probabilities using the correct sign and exponent
  log_P <- -delta_M - delta_eta * phi  # direct multiplication with Φ
  
  # convert to probabilities
  P <- exp(log_P)
  P <- P / sum(P)  # normalize (ensure probabilities sum to 1)
  
  # compute log-likelihood
  LL <- sum(codon_counts * log(P))
  
  # return gene name and log-likelihood
  return(list(Gene = gene_name, LogLikelihood = LL))
}

# example usage
gene_name <- "GeneX"
delta_M <- c(-0.5, 0.1, -0.3)  # ΔM for codons
delta_eta <- c(0.2, -0.1, 0.5)  # Δη for codons
phi <- 50                       # protein synthesis rate (Φ)
codon_counts <- c(10, 5, 20)     # codon counts

LLikGenome(gene_name, delta_M, delta_eta, phi, codon_counts)

