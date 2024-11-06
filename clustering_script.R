clustering_script <- function(input, output = "./", k = 2) {
  # Load necessary libraries
  library(seqinr)
  library(AnaCoDa)
  library(ca)
  library(cluster)
  library(argparse)
  
  # Set up argument parser (for standalone usage, can be ignored within function)
  # parser <- ArgumentParser()
  # parser$add_argument("-i","--input",help="FASTA file with protein-coding sequences",type="character")
  # parser$add_argument("-o","--output",help="Output directory",type="character",default="./")
  # parser$add_argument("-k","--num_clusters",help="Number of clusters to use in CLARA",type="integer",default=2)
  
  # Define paths
  genome.file <- file.path("Genomes", "Genomes", "cds_cleaned", input)
  
  # Initialize genome object
  genome <- initializeGenomeObject(genome.file)
  size <- length(genome)
  
  # Extract genes and IDs
  genes <- lapply(1:size, function(x) {
    tmp <- genome$getGeneByIndex(x, F)
    tmp$seq
  })
  genes <- lapply(genes, function(x) {
    unlist(strsplit(x, split = "", fixed = T))
  })
  
  ids <- lapply(1:size, function(x) {
    tmp <- genome$getGeneByIndex(x, F)
    tmp$id
  })
  df <- data.frame(genes = unlist(ids))
  df[words()] <- numeric(length = size)
  
  # Calculate codon counts and frequencies
  codon.counts <- lapply(genes, function(x) {
    tmp <- uco(x, frame = 0, index = "eff")
    tmp <- tmp / (length(x) / 3)
  })
  
  for (i in 1:size) {
    df[i, words()] <- codon.counts[[i]][words()]
  }
  df[is.na(df)] <- 0
  
  # Clean up data frame and remove stop codons
  rownames(df) <- df$genes
  df <- df[, 2:ncol(df)]
  df <- df[, which(!colnames(df) %in% c("taa", "tga", "tag"))]
  
  # Perform Correspondence Analysis
  fit <- ca(df)
  print("Beginning clustering...")
  
  # Perform CLARA clustering on the CA results
  clarax <- clara(x = fit$rowcoord[, 1:4], k = k, samples = 200, sampsize = length(genome) / 2, rngR = T, pamLike = T)
  print("Done...")
  
  # Output the clustering results
  clusters <- data.frame(Gene = names(clarax$clustering), Cluster = clarax$clustering, stringsAsFactors = F)
  write.table(clusters, file.path(output, input), sep = "\t", col.names = T, row.names = F, quote = F)
  
  # Generate PDF with plots
  pdf(file.path(output, paste0(input, "_", k, ".pdf")))
  plot(fit, labels = c(0, 2))
  plot(factoextra::fviz_cluster(clarax, labelsize = 0))
  dev.off()
  
  # Return the clustering results for testing purposes
  return(clusters)
}
