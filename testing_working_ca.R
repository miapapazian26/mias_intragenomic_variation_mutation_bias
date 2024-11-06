library(testthat)

# Load the clustering script
source("path/to/clustering_script.R")  # Update with the correct path

test_that("CLARA clustering is reproducible with controlled inputs", {
  set.seed(42)  # Fix seed for reproducibility
  
  # Run the script and capture the clustering output
  output1 <- clustering_script(input = "Data/test_sequences_1.cds", output = "./", k = 2)
  
  # Run the script a second time and capture output
  output2 <- clustering_script(input = "Data/test_sequences_1.cds", output = "./", k = 2)
  
  # Check that cluster assignments are identical across runs
  expect_equal(output1$Cluster, output2$Cluster)
})

test_that("Expected number of clusters", {
  # Run the script with controlled test data
  set.seed(42)
  output <- clustering_script(input = "Data/test_sequences_1.cds", output = "./", k = 2)
  
  # Verify that exactly two clusters were created
  expect_equal(length(unique(output$Cluster)), 2)
})

test_that("Cluster assignment stability with sequence changes", {
  # Run the script with similar sequences and check for stable clustering
  set.seed(42)
  output <- clustering_script(input = "Data/similar_sequences.cds", output = "./", k = 2)
  
  # Check that similar sequences cluster together (assuming IDs reflect similarity)
  similar_cluster <- output$Cluster[output$Gene == "similar_gene_id"]
  expect_true(all(output$Cluster[output$Gene %in% similar_gene_ids] == similar_cluster))
})


test_that("Cluster assignment stability with sequence changes", {
  #run the script with similar sequences and check for stable clustering
  set.seed(42)
  output <- clustering_script(input = "similar_sequences.fasta", output = "./", k = 2)
  
  #check that similar sequences cluster together (assuming IDs reflect similarity)
  similar_cluster <- output$Cluster[output$Gene == "similar_gene_id"]
  expect_true(all(output$Cluster[output$Gene %in% similar_gene_ids] == similar_cluster))
})
