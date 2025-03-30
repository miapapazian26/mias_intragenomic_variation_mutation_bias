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

