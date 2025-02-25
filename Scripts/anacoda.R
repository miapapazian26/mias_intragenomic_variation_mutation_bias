#ex: using codon data in the form of CDS in fasta format with one mixture (ROC)
#the following example illustrates how you would estimates parameters under the ROC model of a given set of protein coding genes, assuming the same mutation and selection regime for all genes.

library(AnaCoDa)
library(tidyr)
library(dplyr)

#initialize genome object
genome <- initializeGenomeObject(file = "fasta/revisit_cds_data/candida_tenuis.max.cds")

#initialize parameter object
parameter <- initializeParameterObject(genome = genome, sphi = 1, num.mixtures = 1, gene.assignment = rep(1, length(genome)))

#save parameter object
param_file <- file.path(tempdir(), "parameter.Rda")
writeParameterObject(parameter = parameter, file = param_file)

#initialize MCMC object
mcmc <- initializeMCMCObject(samples = 100, thinning = 10, adaptive.width = 50)

#initialize model
model <- initializeModelObject(parameter = parameter, model = "ROC")

#run MCMC (modifies 'mcmc' in place)
runMCMC(mcmc = mcmc, genome = genome, model = model)

#save updated MCMC object 
mcmc_file <- file.path(tempdir(), "mcmc_results.Rda")
writeMCMCObject(mcmc = mcmc, file = mcmc_file)

trace <- parameter$getTraceObject()

#plot trace of gene ###
plot(x = trace, what = "Expression", mixture = 1, geneIndex = 50)
#x axis is number of steps 
#y is estimated expression level for gene index ###

#synthesis trace 
synthesis_trace <- trace$getSynthesisRateTraceForGene(1)
long_df <- as.data.frame(synthesis_trace)
head(long_df)

synth_trace_list <- trace$getSynthesisRateTrace()
synth_trace_mix <- synth_trace_list[[1]] %>% tibble()

length(synth_trace_list)

#get csp
csp_mat <- getCSPEstimates(parameter = parameter, mixture = 1, samples = 100)
head(csp_mat)


