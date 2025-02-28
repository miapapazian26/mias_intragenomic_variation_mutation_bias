#found this in the code
# calculate the log-marginal likelihood
parameter <- loadParameterObject("parameter.Rda")
mcmc <- loadMCMCObject("mcmc.Rda")
calculate_marginal_likelihood(parameter, mcmc, mixture = 1, samples = 100, scaling = 1.5)

#code function 
calculateMarginalLogLikelihood <- function(parameter, mcmc, mixture, n.samples, divisor,warnings=TRUE)
{  
  if(divisor < 1) stop("Generalized Harmonic Mean Estimation of Marginal Likelihood requires importance sampling distribution variance divisor be greater than 1")
  
  ## Collect information from AnaCoDa objects
  trace <- parameter$getTraceObject()
  ## This should be the posterior instead of the log_posterior but this causes an overflow, find fix!!!
  log_posterior <- mcmc$getLogPosteriorTrace()
  log_posterior <- log_posterior[(length(log_posterior) - n.samples+1):(length(log_posterior))]
  
  
  ### HANDLE CODON SPECIFIC PARAMETERS
  log_imp_dens_sample <- rep(0, n.samples)
  for (k in 1:mixture)
  {
    for(ptype in 0:1) # for all parameter types (mutation/selection parameters)
    {
      for(aa in AnaCoDa::aminoAcids()) # for all amino acids
      {
        if(aa == "M" || aa == "W" || aa == "X") next # ignore amino acids with only one codon or stop codons 
        codons <- AnaCoDa::AAToCodon(aa, focal = T)
        ## get covariance matrix and mean of importance distribution
        sample_mat <- matrix(NA, ncol = length(codons), nrow = n.samples)
        mean_vals <- rep(NA, length(codons))
        for(i in 1:length(codons)) # for all codons
        {
          vec <- trace$getCodonSpecificParameterTraceByMixtureElementForCodon(k, codons[i], ptype, TRUE)
          vec <- vec[(length(vec) - n.samples+1):(length(vec))]
          sample_mat[,i] <- vec
          mean_vals[i] <- mean(vec)
        }
        ## scale/shrinked covariance matrix
        cov_mat <- cov(sample_mat) / divisor
        if(all(cov_mat == 0))
        {
          if(warnings) print(paste("Covariance matrix for codons in amino acid",aa,"has 0 for all values. Skipping."))
          next
        }
        for(i in 1:n.samples)
        {
          ## calculate importance density for collected samples
          ## mikeg: We should double check with Russ that it is okay to use the full covariance matrix of the sample (even though we have no prior on the cov structure) when constructing the importance density function.
          ## It seems logical to do so
          log_imp_dens_aa <- dmvnorm(x = sample_mat[i,], mean = mean_vals, sigma = cov_mat, log = TRUE)
          log_imp_dens_sample[i] = log_imp_dens_sample[i] + log_imp_dens_aa
        }
      }
    }
  }
  ## HANDLE GENE SPECIFIC PARAMETERS
  
  # phi values are stored on natural scale.
  
  synt_trace <- trace$getSynthesisRateTrace()[[mixture]]
  n_genes <- length(synt_trace);
  
  sd_vals <- rep(NA, n_genes)
  mean_vals <- rep(NA, n_genes)
  for(i in 1:n_genes) ## i is indexing across genes
  {
    vec <- synt_trace[[i]]
    vec <- vec[(length(vec) - n.samples+1):(length(vec))]
    sd_vals[i] <- sd(vec)
    if (all(sd_vals[i] == 0))
    {
      if(warnings) print(paste("Variance of gene",i,"is 0. Skipping."))
      next
    }
    mean_vals[i] <- mean(vec)
    log_mean_vals <- log(mean_vals) - 0.5 * log(1+(sd_vals^2/mean_vals^2))
    log_sd_vals <- sqrt(log(1+(sd_vals^2/mean_vals^2)))
    
    ## Calculate vector of importance density for entire \phi trace of gene.
    log_imp_dens_phi <- dlnorm(x = vec, meanlog = log_mean_vals[i], sdlog = log_sd_vals[i]/divisor, log = TRUE)
    
    
    ## update importance density function vector of sample with current gene;
    log_imp_dens_sample <- log_imp_dens_sample + log_imp_dens_phi
    
  } ## end synth_trace loop
  
  ## Scale importance density for each sample by its posterior probability (on log scale)
  log_imp_dens_over_posterior <- log_imp_dens_sample - log_posterior
  ## now scale by max term to facilitate summation
  max_log_term <- max(log_imp_dens_over_posterior)
  ## Y = X - max_X
  offset_log_imp_dens_over_posterior <- log_imp_dens_over_posterior - max_log_term
  ## Z = sum(exp(vec(Y)))
  offset_sum_imp_dens_over_posterior <- sum(exp(offset_log_imp_dens_over_posterior))
  log_sum_imp_dens_over_posterior <- log(offset_sum_imp_dens_over_posterior) + max_log_term
  ## ln(ML) = ln(n) -(Z + max_X)
  log_marg_lik <- log(n.samples) - log_sum_imp_dens_over_posterior
  ##marg_lik = 1.0/(log_inv_marg_lik/n.samples) # equation 9
  return(log_marg_lik)
}

#' Find and return list of optimal codons
#' 
#' \code{findOptimalCodon} extracrs the optimal codon for each amino acid.
#' 
#' @param csp a \code{data.frame} as returned by \code{getCSPEstimates}.
#'
#' @return A named list with with optimal codons for each amino acid.
#'
#' @examples 
#' genome_file <- system.file("extdata", "genome.fasta", package = "AnaCoDa")
#'
#' genome <- initializeGenomeObject(file = genome_file)
#' sphi_init <- 1
#' numMixtures <- 1
#' geneAssignment <- rep(1, length(genome))
#' parameter <- initializeParameterObject(genome = genome, sphi = sphi_init, 
#'                                        num.mixtures = numMixtures, 
#'                                        gene.assignment = geneAssignment, 
#'                                        mixture.definition = "allUnique")
#' model <- initializeModelObject(parameter = parameter, model = "ROC")
#' samples <- 2500
#' thinning <- 50
#' adaptiveWidth <- 25
#' mcmc <- initializeMCMCObject(samples = samples, thinning = thinning, 
#'                              adaptive.width=adaptiveWidth, est.expression=TRUE, 
#'                              est.csp=TRUE, est.hyper=TRUE, est.mix = TRUE) 
#' divergence.iteration <- 10
#' \dontrun{
#' runMCMC(mcmc = mcmc, genome = genome, model = model, 
#'         ncores = 4, divergence.iteration = divergence.iteration)
#' 
#' csp_mat <- getCSPEstimates(parameter, CSP="Selection")
#' opt_codons <- findOptimalCodon(csp_mat)
#' }

calculateMarginalLogLikelihood(parameter, mcmc, mixture, n.samples, divisor)

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

