#ex: using codon data in the form of CDS in fasta format with one mixture (ROC)
#the following example illustrates how you would estimates parameters under the ROC model of a given set of protein coding genes, assuming the same mutation and selection regime for all genes.
library(AnaCoDa)
genome <- initializeGenomeObject(file = "fasta/revisit_cds_data/candida_tenuis.max.cds")
parameter <- initializeParameterObject(genome = genome, sphi = 1, num.mixtures = 1, gene.assignment = rep(1, length(genome)))
mcmc <- initializeMCMCObject(samples = 2000, thinning = 10, adaptive.width = 50)
model <- initializeModelObject(parameter = parameter, model = "ROC")

mcmc_results <- runMCMC(mcmc = mcmc, genome = genome, model = model)

if (!is.null(mcmc_results)) {
  save(mcmc_results, file = "mcmc_results.RData")
  print("MCMC results saved successfully!")
  print(mcmc_results)
} else {
  print("MCMC did not return results.")
}

save(mcmc_results, file = "mcmc_results.RData")

print(mcmc_results)


#checking the data 
genome_names <- getNames(genome)
print(genome_names)
codon_counts <- getCodonCounts(genome)
print(codon_counts)
