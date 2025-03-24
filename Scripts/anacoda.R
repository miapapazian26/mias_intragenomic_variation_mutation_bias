#ex: using codon data in the form of CDS in fasta format with one mixture (ROC)
#the following example illustrates how you would estimate parameters under the ROC model of a given set of protein coding genes, assuming the same mutation and selection regime for all genes.

# load necessary libraries
library(AnaCoDa)
library(tidyr)
library(dplyr)

##
# initialization
##
# initialize genome object
genome <- initializeGenomeObject(file = "fasta/revisit_cds_data/candida_tenuis.max.cds")

# initialize parameter object
parameter <- initializeParameterObject(genome = genome, sphi = 1, num.mixtures = 1, gene.assignment = rep(1, length(genome)))

# save and load parameter object
param_file <- file.path(tempdir(), "parameter.Rda")
writeParameterObject(parameter = parameter, file = param_file)
load(param_file)

# initialize mcmc object
mcmc <- initializeMCMCObject(samples = 100, thinning = 10, adaptive.width = 50)

# initialize model
model <- initializeModelObject(parameter = parameter, model = "ROC")

##
# run mcmc
##
# run mcmc (modifies 'mcmc' in place)
runMCMC(mcmc = mcmc, genome = genome, model = model)

# save updated mcmc object 
mcmc_file <- file.path(tempdir(), "mcmc_results.Rda")
writeMCMCObject(mcmc = mcmc, file = mcmc_file)

##
# trace analysis
##
# get trace
trace <- parameter$getTraceObject()

# plot trace of gene ###
plot(x = trace, what = "Expression", mixture = 1, geneIndex = 1)

# phi trace (synthesis rate)
synth_trace_list <- trace$getSynthesisRateTrace()
synth_trace_mix <- do.call(rbind, synth_trace_list[[1]])
head(synth_trace_mix) # check first few rows
length(synth_trace_mix) # check number of columns

##
# codon-specific parameter estimates
##
# get csp estimates
csp_est <- getCSPEstimates(parameter = parameter, mixture = 1, samples = 100)
head(csp_mat)

# get csp traces
csp_trace <- parameter$getTraceObject()
csp_trace_data <- csp_trace$getCodonSpecificParameterTrace(0)
csp_trace_1 <- csp_trace_data[[1]]
csp_trace_df1 <- as.data.frame(csp_trace_1)

# reshape the data into long format
csp_trace_long1 <- csp_trace_df1 %>%
  pivot_longer(cols = everything(), names_to = "Parameter", values_to = "Value")

head(csp_trace_long1)
summary(csp_trace_long1) # summary stats 

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

##
# codon usage analysis
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

# check codon counts for first gene
codon_usage_per_gene[[gene_names[1]]]

##
# log-likelihood loop analysis
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

