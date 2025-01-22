#ex: using codon data in the form of CDS in fasta format with one mixture (ROC)
#the following example illustrates how you would estimates parameters under the ROC model of a given set of protein coding genes, assuming the same mutation and selection regime for all genes.
genome <- initializeGenomeObject(file = "fasta/revisit_cds_data/candida_tenuis.max.cds")
parameter <- initializeParameterObject(genome = genome, sphi = 1, num.mixtures = 1, gene.assignment = rep(1, length(genome)))
mcmc <- initializeMCMCObject(samples = 5000, thinning = 10, adaptive.width=50)
model <- initializeModelObject(parameter = parameter, model = "ROC")
runMCMC(mcmc = mcmc, genome = genome, model = model)