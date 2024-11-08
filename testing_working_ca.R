library(testthat)

#source clustering script
source("~/repositories/mias_intragenomic_variation_mutation_bias/clustering_script.R")  

#get a list of all CDS files in the Data folder
cds_files <- list.files("Data", pattern = "\\.cds$", full.names = TRUE)

#loop through each CDS file and run the tests
for (file in cds_files) {
  test_that(paste("Testing clustering for", basename(file)), {
    set.seed(42)  # Fix seed for reproducibility
    
    #run the script and capture the clustering output
    output1 <- clustering_script(input = file, output = "./", k = 2)
    
    #run the script a second time and capture output
    output2 <- clustering_script(input = file, output = "./", k = 2)
    
    #check that cluster assignments are identical across runs
    expect_equal(output1$Cluster, output2$Cluster)
    
    #test for expected number of clusters
    expect_equal(length(unique(output1$Cluster)), 2)
    
    # run the script with similar sequences and check for stable clustering
    # (modify this part based on your actual test data and IDs)
    output <- clustering_script(input = file, output = "./", k = 2)
    similar_cluster <- output$Cluster[output$Gene == "similar_gene_id"]  # replace with real gene ID
    expect_true(all(output$Cluster[output$Gene %in% similar_gene_ids] == similar_cluster))  # replace with real IDs
  })
}
