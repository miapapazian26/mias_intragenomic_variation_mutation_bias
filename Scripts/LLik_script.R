#usage
gene_names <- getNames(genome = genome)
delta_M <- mutationTrace[[1]][[1]]  # ΔM for codons
delta_eta <- selectionTrace[[1]][[1]]  # Δη for codons
phi <- synth_trace_mix  # protein synthesis rate (Φ)
codon_counts <- #working (photo from 3/5)   # codon counts and gene ids

#delta eta and delta mu
#need to run loop first 

mutationTrace <- trace$getCodonSpecificParameterTrace(0) #mu

selectionTrace <- trace$getCodonSpecificParameterTrace(1) #eta

#checking 
head(mutationTrace) #looking delta M
head(selectionTrace) #looking delta eta
head(codon_counts) #looking at codon counts 
head(phi) #looking at phi values 

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


LLikGene(gene_name, delta_M, delta_eta, phi, codon_counts)

#logliklihood amino acid function 
LLikAA <- function(delta_M, delta_eta, phi, codon_counts) {
  log_P = delta_M - delta_eta * phi
  log_P = log_P - max(log_P)  
  p = exp(log_P)

  
  LL = dmultinom(x = codon_counts, prob = p, log = TRUE)
  return(LL)
}

LLikAA(delta_M, delta_eta, phi, codon_counts) 


