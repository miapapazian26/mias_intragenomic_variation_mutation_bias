#this code analyzes the objects created in anacoda script
#converts model output to useful tibbles 

library(R.utils) ## need for gunzipping of fasta files.
## mike's path to the development version
## library(AnaCoDa, lib.loc="~/R/lib/4.4-AnaCoDa/upstream-develop/")
## path to use if using r-cran version
library(AnaCoDa)
library(tidyr)
library(dplyr)

#load previous workspace
species <- "candida_tenuis"
load(file.path("Results", paste0(species, "_run.info.Rda")))

#regenerate genome object 
gunzip(gz_file, remove=FALSE, temp=TRUE, skip = TRUE)  #unzips fasta to tempdir() output
genome <- initializeGenomeObject(file = file.path(tempdir(), fasta_file))

#load model output objects 
parameter <- loadParameterObject(param_file)
mcmc <- loadMCMCObject(mcmc_file)

##
# trace analysis
##
# get trace
trace <- parameter$getTraceObject()

# plot trace of gene to test previous code works ###
plot(x = trace, what = "Expression", mixture = 1, geneIndex = 1)

# phi trace (synthesis rate)
synth_trace_list <- trace$getSynthesisRateTrace()[[1]]
trace_to_tibble <- function(list_trace, id_names) {
  ## convert list of traces into a matrix
  matrix <- do.call(rbind, list_trace)
  sample_n <- 1:ncol(matrix)-1
  colnames(matrix) <- sample_n
  ## create trace tibble
  tibble <- as_tibble(matrix) %>% 
    bind_cols(tibble(id_names), .)
  return(tibble)
}

synth_trace <- trace_to_tibble(synth_trace_list,
                               id_names = getNames(genome))
head(synth_trace) # check first few rows
dim(synth_trace) # check number of columns

#USE SYNTH TRACE MIX IF INDEXING BY ROW/SAMPLE NUMBER
synth_trace_mix <- do.call(rbind, synth_trace_list[[1]])
head(synth_trace_mix) # check first few rows
length(synth_trace_mix) # check number of columns

##
# codon-specific parameter estimates
##

# get csp estimates
csp_est <- getCSPEstimates(parameter = parameter, mixture = 1, samples = 100)
head(csp_mat)

#filter out reference codons which have a mean of 0
csp_codons <- csp_est$Mutation %>%
  filter(Mean !=0) %>%
  select("Codon")

# get csp traces
csp_cat <- c("Mutation", "Selection")
csp_wide_list <- list()

for(i in 1:2) {
  param = csp_cat[i] ## don't need double [] because already atomic object
  trace_list <- trace$getCodonSpecificParameterTrace((i-1))[[1]]
  csp_wide_list[[i]] <- trace_to_tibble(trace_list, csp_codons)  %>%
    rowwise() %>%  #use rowwise because codonToAA isn't vectorized.
    mutate(AA = codonToAA(Codon), .before = "Codon")
  
}
names(csp_wide_list) <- csp_cat

csp_long_list <- list()
for(i in 1:2) {
  param = csp_cat[i] ## don't need double [] because already atomic object
  csp_long_list[[i]] <- pivot_longer(csp_wide_list[[i]],
                                     cols = !c("AA", "Codon"),
                                     values_to = param,
                                     names_to = "Sample") %>%
    mutate(Sample = as.numeric(Sample))
}
names(csp_long_list) <- csp_cat

csp_long <- left_join(csp_long_list$Mutation,
                      csp_long_list$Selection,
                      by = join_by(AA, Codon, Sample)
)
head(csp_long)

#code not being used currently 
##
# mutation/selection trace
##
mutationTrace <- trace$getCodonSpecificParameterTrace(0)
selectionTrace <- trace$getCodonSpecificParameterTrace(1)

##
# likelihood analysis(what i am comparing MY results to)
##
# convert log-likelihood trace to dataframe
logLikeTrace_df <- data.frame(iteration = 1:length(logLikeTrace), logLik = logLikeTrace)
head(logLikeTrace_df)

# calculate marginal likelihood
parameter <- loadParameterObject("parameter.Rda")
mcmc <- loadMCMCObject("mcmc.Rda")
calculate_marginal_likelihood(parameter, mcmc, mixture = 1, samples = 100, scaling = 1.5)

# calculate bayes factor for two models
parameter1 <- loadParameterObject("parameter1.Rda")
parameter2 <- loadParameterObject("parameter2.Rda")
mcmc1 <- loadMCMCObject("mcmc1.Rda")
mcmc2 <- loadMCMCObject("mcmc2.Rda")
mll1 <- calculate_marginal_likelihood(parameter1, mcmc1, mixture = 1, samples = 500, scaling = 1.5)
mll2 <- calculate_marginal_likelihood(parameter2, mcmc2, mixture = 1, samples = 500, scaling = 1.5)
cat("Bayes factor: ", mll1 - mll2, "\n")


#CODE: from the documentation 
#' Get Codon Counts For all Amino Acids
#' 
#' 
#' @param genome A genome object from which the counts of each
#' codon can be obtained.
#'  
#' @return Returns a data.frame storing the codon counts for each amino acid. 
#' 
#' @description provides the codon counts for a fiven amino acid across all genes
#' 
#' @details The returned matrix containes a row for each gene and a column 
#' for each synonymous codon of \code{aa}.
#' 
#' @examples 
#' 
#' genome_file <- system.file("extdata", "genome.fasta", package = "AnaCoDa")
#'  
#' ## reading genome
#' genome <- initializeGenomeObject(file = genome_file)
#' counts <- getCodonCounts(genome)
#' 
getCodonCounts <- function(genome){
  codons <- codons()
  ORF <- getNames(genome)
  codonCounts <- lapply(codons, function(codon) {
    codonCounts <- genome$getCodonCountsPerGene(codon)
  })
  codonCounts <- do.call("cbind", codonCounts)
  colnames(codonCounts) <- codons
  rownames(codonCounts) <- ORF
  return(as.data.frame(codonCounts,stringsAsFactors = F))
}

#' Get Codon Counts For a specific Amino Acid
#' 
#' @param aa One letter code of the amino acid for which the codon counts should be returned
#' 
#' @param genome A genome object from which the counts of each
#' codon can be obtained.
#'  
#' @return Returns a data.frame storing the codon counts for the specified amino acid. 
#' 
#' @description provides the codon counts for a fiven amino acid across all genes
#' 
#' @details The returned matrix containes a row for each gene and a coloumn 
#' for each synonymous codon of \code{aa}.
#' 
#' @examples 
#' 
#' genome_file <- system.file("extdata", "genome.fasta", package = "AnaCoDa")
#'  
#' ## reading genome
#' genome <- initializeGenomeObject(file = genome_file)
#' counts <- getCodonCountsForAA("A", genome)
#' 
getCodonCountsForAA <- function(aa, genome){
  # get codon count for aa
  codons <- AAToCodon(aa, F)
  codonCounts <- lapply(codons, function(codon){
    codonCounts <- genome$getCodonCountsPerGene(codon)
  })
  codonCounts <- do.call("cbind", codonCounts)
  return(codonCounts)
}

#' calculates the synonymous codon usage order (SCUO) 
#' 
#' \code{calculateSCUO} calulates the SCUO value for each gene in genome. Note that if a codon is absent, this will be treated as NA and will be skipped in final calculation
#' 
#' @param genome A genome object initialized with \code{\link{initializeGenomeObject}}.
#' 
#' @return returns the SCUO value for each gene in genome
#' 
#' @examples 
#' 
#' genome_file <- system.file("extdata", "genome.fasta", package = "AnaCoDa")
#'  
#' ## reading genome
#' genome <- initializeGenomeObject(file = genome_file)
#' scuo <- calculateSCUO(genome)
#' 
calculateSCUO <- function(genome)
{
  aas <- aminoAcids()
  aas <- aas[which(!aas %in% c("M","X","W","J"))]
  genes <- genome$getGenes(F)
  scuo.values <- data.frame(ORF=getNames(genome), SCUO=rep(NA, length(genome)))
  for(i in 1:length(genes))
  {
    g <- genes[[i]]
    total.aa.count <- g$length()/3
    
    scuo.per.aa <- unlist(lapply(X = aas, FUN = function(aa)
    {
      codon <- AAToCodon(aa = aa, focal = F)
      num.codons <- length(codon)
      aa.count <- g$getAACount(aa)
      if(aa.count == 0) return(0)
      
      codon.count <- unlist(lapply(codon, FUN = function(c){return(g$getCodonCount(c))}))
      codon.propotions <- codon.count / aa.count
      aa.entropy <- -1*sum(codon.propotions * log(codon.propotions),na.rm = T)
      max.entropy <- -log(1/num.codons)
      norm.entropy.diff <- (max.entropy - aa.entropy) / max.entropy
      
      comp.ratio <- aa.count / total.aa.count
      
      scuo.aa <- comp.ratio * norm.entropy.diff
      scuo.aa
    }))
    scuo.values[i,"SCUO"] <- sum(scuo.per.aa,na.rm = T)
  }
  return(scuo.values)
}

#' Length of Genome
#' 
#' \code{length} gives the length of a genome
#' 
#' @param x A genome object initialized with \code{\link{initializeGenomeObject}}.
#' 
#' @return returns the number of genes in a genome
#' 
#' @examples 
#' 
#' genome_file <- system.file("extdata", "genome.fasta", package = "AnaCoDa")
#'  
#' ## reading genome
#' genome <- initializeGenomeObject(file = genome_file)
#' length(genome) # 10
#' 
length.Rcpp_Genome <- function(x) {
  return(x$getGenomeSize(F))
}

#' Summary of Genome
#' 
#' \code{summary} summarizes the description of a genome, such as number of genes and average gene length.
#' 
#' @param object A genome object initialized with \code{\link{initializeGenomeObject}}.
#' 
#' @param ... Optional, additional arguments to be passed to the main summary function 
#' that affect the summary produced.
#'
#' @return This function returns by default an object of class c("summaryDefault", table").
summary.Rcpp_Genome <- function(object, ...) {
  # TODO output stuff like:
  # - no. of genes
  # - avg. gene length
  # - avg. A,C,G,T content
  # - avg. AA composition
  # - ...
  summary(object, ...)
}


#' Gene Names of Genome
#' 
#' 
#' @param genome A genome object initialized with \code{\link{initializeGenomeObject}}.
#' 
#' @param simulated A logical value denoting if the gene names to be listed are simulated or not.
#' The default value is FALSE.
#' 
#' @description returns the identifiers of the genes within the genome specified.
#' 
#' @return gene.names Returns the names of the genes as a vector of strings.
#' 
#' @examples 
#' 
#' genome_file <- system.file("extdata", "genome.fasta", package = "AnaCoDa")
#'  
#' ## reading genome
#' genome <- initializeGenomeObject(file = genome_file)
#'
#' ## return all gene ids for the genome
#' geneIDs <- getNames(genome, FALSE)
#' 
getNames <- function(genome, simulated = FALSE)
{
  genes <- genome$getGenes(simulated)
  gene.names <- unlist(lapply(1:length(genes), function(i){return(genes[[i]]$id)}))
  return(gene.names)
}

#' Get gene observed synthesis rates
#' 
#' \code{getObservedSynthesisRateSet} returns the observed 
#' synthesis rates of the genes within the genome specified.
#' 
#' @param genome A genome object initialized with \code{\link{initializeGenomeObject}}.
#' 
#' @param simulated A logical value denoting if the synthesis 
#' rates to be listed are simulated or not. The default value is FALSE.
#' 
#' @return Returns a data.frame with the observed expression values in genome
#' 
#' @examples 
#' 
#' genome_file <- system.file("extdata", "genome.fasta", package = "AnaCoDa")
#' expression_file <- system.file("extdata", "expression.csv", package = "AnaCoDa") 
#' ## reading genome
#' genome <- initializeGenomeObject(file = genome_file)
#' 
#'
#' ## return expression values as a data.frame with gene ids in the first column.
#' expressionValues <- getObservedSynthesisRateSet(genome = genome)
#' 
getObservedSynthesisRateSet <- function(genome, simulated = FALSE)
{
  genes <- genome$getGenes(simulated)
  expression <- lapply(1:length(genes), function(i){return(genes[[i]]$getObservedSynthesisRateValues())})
  ids <- getNames(genome, simulated)
  mat <- do.call(rbind,expression)
  return(cbind.data.frame(ids, mat,stringsAsFactors=F))
}


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


#found in parameterObject.R (line 2031)
numMutationCategories <- tempEnv$paramBase$numMut
numSelectionCategories <- tempEnv$paramBase$numSel
max <- tempEnv$paramBase$lastIteration + 1

if (i == 1){
  
  codonSpecificParameterTraceMut <- vector("list", length=numMutationCategories)
  for (j in 1:numMutationCategories) {
    codonSpecificParameterTraceMut[[j]] <- vector("list", length=length(tempEnv$mutationTrace[[j]]))
    for (k in 1:length(tempEnv$mutationTrace[[j]])){
      codonSpecificParameterTraceMut[[j]][[k]] <- tempEnv$mutationTrace[[j]][[k]][1:max]
      codonSpecificParameterTraceSel[[j]][[k]] <- tempEnv$selectionTrace[[j]][[k]][1:max]
    }
  }
  
  codonSpecificParameterTraceSel <- vector("list", length=numSelectionCategories)
  for (j in 1:numSelectionCategories) {
    codonSpecificParameterTraceSel[[j]] <- vector("list", length=length(tempEnv$selectionTrace[[j]]))
    for (k in 1:length(tempEnv$selectionTrace[[j]])){
      codonSpecificParameterTraceMut[[j]][[k]] <- tempEnv$mutationTrace[[j]][[k]][1:max]
      codonSpecificParameterTraceSel[[j]][[k]] <- tempEnv$selectionTrace[[j]][[k]][1:max]
    }
}}


#found in RibModelFramework/R/plotTraceObject.R
pdf("LKs5000trace_plot_mutation.pdf", width = 12, height = 16)
par(mar = c(1, 1, 1, 1))  # Set minimal margins
plot(trace, what = "Mutation", mixture = 1)  # replace with specific plot call
dev.off()

pdf("LKs5000trace_plot_selection.pdf", width = 12, height = 16)
par(mar = c(1, 1, 1, 1))  # Set minimal margins
plot(trace, what = "Selection", mixture = 1)  # replace with specific plot call
dev.off()

#documented code   
#Intended to combine 2D traces (vector of vectors) read in from C++. The firs
#element of the second trace is omited since it should be the same as the 
    #last value of the first trace.
    combineTwoDimensionalTrace <- function(trace1, trace2,start=2,end=NULL){
      if(start < 2)
      {
        print("Start must be at least 2 because the last element of first trace is first of second trace. Setting start = 2.")
      }
      if(is.null(end) || end <= start)
      {
        print(paste0("End must be greater than start. Setting end to length(trace2) = ", length(trace2)))
        end <- length(trace2)
      }
      for (size in 1:length(trace1))
      {
        trace1[[size]]<- c(trace1[[size]], trace2[[size]][start:end])
      }
      return(trace1)
    })

#Intended to combine 3D traces (vector of vectors of vectors) read in from C++. The first
    #element of the second trace is omited since it should be the same as the 
    #last value of the first trace.
    combineThreeDimensionalTrace <- function(trace1, trace2, max){
      
      for (size in 1:length(trace1)){
        for (sizeTwo in 1:length(trace1[[size]])){
          trace1[[size]][[sizeTwo]] <- c(trace1[[size]][[sizeTwo]], 
                                         trace2[[size]][[sizeTwo]][2:max])
        }
      }
    }


    
    
    #' Calculate the Effective Number of Codons
    #' 
    #' 
    #' \code{getNc} returns the Effective Number of Codons for a genome.
    #' 
    #' @param genome A genome object initialized with \code{\link{initializeGenomeObject}}.
    #' 
    #' @return Returns a named vector with the Effective Number of Codons for each gene
    #' 
    #' @examples 
    #' 
    #' genome_file <- system.file("extdata", "more_genes.fasta", package = "AnaCoDa")
    #' ## reading genome
    #' genome <- initializeGenomeObject(file = genome_file)
    #'
    #' nc <- getNc(genome)
    #' 
    getNc <- function(genome)
    {
      aa.vec <- aminoAcids()
      aa.vec <- aa.vec[-length(aa.vec)]
      
      f.mat <- matrix(0, ncol = 6, nrow = length(genome))
      division.counter <- matrix(0, ncol = 6, nrow = length(genome))
      
      for(aa in aa.vec)
      {
        if(aa == "M" || aa == "W") next # contribution of M and W is 2 in total
        
        codonCountForAA <- getCodonCountsForAA(aa, genome = genome)
        n <- rowSums(codonCountForAA)
        pi <- codonCountForAA / n
        
        f.vec <- ( ((n*rowSums(pi*pi))-1) / (n-1) )
        f.vec[!is.finite(f.vec)] <- 0
        
        ncodons <- length(AAToCodon(aa))
        f.mat[, ncodons] <- f.mat[, ncodons] + f.vec
        division.counter[n > 1, ncodons] <- division.counter[n > 1, ncodons] + 1
      }
      
      # adjusted number of AA with codons 2, 4, and 6 since we split Serine
      meanF <- data.frame(SF2=f.mat[,2]/division.counter[, 2], SF3=f.mat[,3], SF4=f.mat[,4]/division.counter[, 4], SF6=f.mat[,6]/division.counter[, 6])
      rare.Ile <- meanF$SF3 < 1
      meanF$SF3[rare.Ile] <- (meanF$SF2[rare.Ile] + meanF$SF4[rare.Ile])/2 # correcting for rare or mising Ile as suggested by Wright (1990, p25)
      
      #Wright (1990) Eqn. 3 adjusted for split serine
      Nc <- 2 + 10/meanF$SF2 + 1/meanF$SF3 + 6/meanF$SF4 + 2/meanF$SF6
      Nc[Nc > 61] <- 61 # revising Nc as suggested by Wright (1990, p25)
      names(Nc) <- getNames(genome, FALSE)
      return(Nc)  
    }
    
    #' Calculate the Effective Number of Codons for each Amino Acid
    #' 
    #' 
    #' \code{getNcAA} returns the Effective Number of Codons for each Amino Acid.
    #' 
    #' @param genome A genome object initialized with \code{\link{initializeGenomeObject}}.
    #' 
    #' @return Returns an object of type \code{data.frame} with the Effective Number of Codons
    #' for each amino acid in each gene.
    #' 
    #' @examples 
    #' 
    #' genome_file <- system.file("extdata", "more_genes.fasta", package = "AnaCoDa")
    #' ## reading genome
    #' genome <- initializeGenomeObject(file = genome_file)
    #'
    #' nc <- getNcAA(genome)
    #' 
    getNcAA <- function(genome)
    {
      aa.vec <- aminoAcids()
      aa.vec <- aa.vec[-length(aa.vec)]
      
      f.mat <- data.frame(matrix(NA, ncol = length(aa.vec), nrow = length(genome)))
      colnames(f.mat) <- aa.vec
      for(aa in aa.vec)
      {
        if(aa == "M" || aa == "W") next # contribution of M and W is 2 in total
        
        codonCountForAA <- getCodonCountsForAA(aa, genome = genome)
        n <- rowSums(codonCountForAA)
        pi <- codonCountForAA / n
        
        f.vec <- ( ((n*rowSums(pi*pi))-1) / (n-1) )
        f.vec <- 1/f.vec
        f.vec[!is.finite(f.vec)] <- NA
        
        f.mat[[aa]] <- f.vec
      }
      rownames(f.mat) <- getNames(genome, FALSE)
      return(f.mat)
    }    
    }

    

    
    