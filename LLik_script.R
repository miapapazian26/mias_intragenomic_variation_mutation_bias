#found this in the code
# calculate the log-marginal likelihood
parameter <- loadParameterObject("parameter.Rda")
mcmc <- loadMCMCObject("mcmc.Rda")
calculate_marginal_likelihood(parameter, mcmc, mixture = 1, samples = 500, scaling = 1.5)

calculateMarginalLogLikelihood <- function(parameter, mcmc, mixture, n.samples, divisor,warnings=TRUE)

#found this in the code
# calculate the bayes factor for two models
parameter1 <- loadParameterObject("parameter1.Rda")
parameter2 <- loadParameterObject("parameter2.Rda")
mcmc1 <- loadMCMCObject("mcmc1.Rda")
mcmc2 <- loadMCMCObject("mcmc2.Rda")
mll1 <- calculate_marginal_likelihood(parameter1, mcmc1, mixture = 1, samples = 500, scaling = 1.5)
mll2 <- calculate_marginal_likelihood(parameter2, mcmc2, mixture = 1, samples = 500, scaling = 1.5)
cat("Bayes factor: ", mll1 - mll2, "\n")


#ai generated function
LLikGenome <- function(gene_name, delta_M, delta_eta, phi, codon_counts) {
  # compute log probabilities using the correct sign and exponent
  log_P <- -delta_M - delta_eta * phi  # direct multiplication with Φ
  
  # convert to probabilities
  P <- exp(log_P)
  P <- P / sum(P)  # normalize (ensure probabilities sum to 1)
  
  # compute log-likelihood
  LL <- sum(codon_counts * log(P))
  
  # return gene name and loglikelihood
  return(list(Gene = gene_name, LogLikelihood = LL))
}

# example usage
gene_name <- "GeneX"
delta_M <- c(-0.5, 0.1, -0.3)  # ΔM for codons
delta_eta <- c(0.2, -0.1, 0.5)  # Δη for codons
phi <- 50                       # protein synthesis rate (Φ)
codon_counts <- c(10, 5, 20)     # codon counts

LLikGenome(gene_name, delta_M, delta_eta, phi, codon_counts)

