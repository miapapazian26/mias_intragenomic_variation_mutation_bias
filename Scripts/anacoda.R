#ex: using codon data in the form of CDS in fasta format with one mixture (ROC)
#the following example illustrates how you would estimate parameters under the ROC model of a given set of protein coding genes, assuming the same mutation and selection regime for all genes.
#this code creates the relevant objects 

# load necessary libraries
library(R.utils) #need this for gunzipping of fasta files 
# mikes path to the development version
#library(AnaCoDa, lib.loc="~/R/lib/4.4-AnaCoDa/upstream-develop/")
# path to use if using r-cran version 
library(AnaCoDa)
library(tidyr)
library(dplyr)

##
# initialization
##

##
# initialize genome object (old)
genome <- initializeGenomeObject(file = "fasta/revisit_cds_data/candida_tenuis.max.cds")
##COPY FILES INTO THE REPO. THIS LINK ONLY WORKS FOR ME 

species <- "candida_tenuis"
fasta_file <- paste0(species, ".max.cds")
#best to use file.path to avoid issues with windows file systems differing from linux and macs
gz_file <- file.path("fasta", "revisit_cds_data",
                     paste0(fasta_file, ".gz"))

#unzip file
gunzip(gz_file, remove=FALSE, temp=TRUE, skip = TRUE)  #unzips fasta to tempdir() output

#initialize genome object
genome <- initializeGenomeObject(file = file.path(tempdir(), fasta_file))

# initialize parameter object
parameter <- initializeParameterObject(genome = genome, sphi = 1, num.mixtures = 1, gene.assignment = rep(1, length(genome)))

# initialize mcmc object
mcmc <- initializeMCMCObject(samples = 100, thinning = 10, adaptive.width = 50)

# initialize model
model <- initializeModelObject(parameter = parameter, model = "ROC")

##
# run mcmc
##

# run mcmc (modifies 'mcmc' in place)
runMCMC(mcmc = mcmc, genome = genome, model = model)

#save mcmc_file into results 
mcmc_file <- file.path("Results", paste0(species, "_mcmc.Rda"))
writeMCMCObject(mcmc = mcmc, file = mcmc_file)

# save and load parameter object - move to after run is complete
param_file <- file.path("Results", paste0(species, "_parameter.Rda"))
writeParameterObject(parameter = parameter, file = param_file)

#saving all non-model related objects 
save.image(file.path("Results", paste0(species, "_run.info.Rda")))

#end of model creation and saving of results 