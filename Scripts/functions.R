#loglikelihood gene function
LLikGene <- function(gene_name, delta_M, delta_eta, phi, codon_counts) {
  # compute log probabilities using the correct sign and exponent
  log_P <- -delta_M - delta_eta * phi  # direct multiplication with Φ
  
  # convert to probabilities
  P <- exp(log_P)
  P <- P / sum(P)  # normalize (ensure probabilities sum to 1)
  
  # compute log-likelihood
  LL = dmultinom(x = codon_counts, prob = P, log = TRUE)
  
  # return gene name and loglikelihood
  return(list(Gene = gene_name, LogLikelihood = LL))
}


#logliklihood amino acid function 
LLikAA <- function(delta_M, delta_eta, phi, codon_counts) {
  log_P = delta_M - delta_eta * phi
  log_P = log_P - max(log_P)  
  p = exp(log_P)
  
  
  LL = dmultinom(x = codon_counts, prob = p, log = TRUE)
  return(LL)
}
