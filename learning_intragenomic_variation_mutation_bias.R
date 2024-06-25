#LEARNING HOW TO MANUEVER THROUGH THE DATA
# Load required libraries (working)
library(seqinr)
library(data.table)

# Check the content of the FASTA file
fasta_file <- "Data/Expression/candida_albicans.max.cds"
fasta_content <- readLines(fasta_file)
head(fasta_content)

# After looking further, the file is not Fasta but CSV based on the output.
# Load required packages
library(data.table)

# Path to the .cds file
cds_file <- "Data/Expression/candida_albicans.max.cds"

# Read the .cds file using fread from data.table
cds_data <- fread(cds_file)

# Load data from .cds files
predicted_values <- fread("path/to/predicted_values.cds")
empirical_estimates <- fread("path/to/empirical_estimates.cds")

#Inspect the Data
head(cds_data)

# Print column names and sample rows
head(predicted_values)
head(empirical_estimates)

# View structure of dataframes
str(predicted_values)
str(empirical_estimates)

# Summary statistics
summary(predicted_values)
summary(empirical_estimates)

# Find the mean values of each column
mean_values <- sapply(predicted_values, mean, na.rm = TRUE)

# Print the mean values
print(mean_values)

# If you only want to find the mean values for specific columns, you can specify those columns in the sapply() function
# For example, if we are only looking for the mean of the "predicted_prob" column, you would use:
mean_predicted_prob <- mean(predicted_values$predicted_prob, na.rm = TRUE)
print(mean_predicted_prob)

# Creating .fasta files out of gtf and gff files//WORKING
# Example steps in R (conceptual)

# Load necessary libraries
library(tidyverse)  # For data manipulation
library(pROC)       # For ROC curve calculation

# Step 1: Read and process GTF/GFF file (example assumes GTF)
gtf_file <- "path/to/your/file.gtf"
gtf_data <- read.table(gtf_file, header = FALSE, sep = "\t")

# Step 2: Extract relevant features (e.g., gene names, biotypes, coordinates)
# Example: Extract gene names and biotypes
genes <- gtf_data$V9  # Assuming column with gene names
biotypes <- gtf_data$V10  # Assuming column with biotypes

# Extract relevant features cont. 
# Load necessary libraries
library(data.table)  # For efficient file reading
library(dplyr)       # For data manipulation

# Example: Read GTF file
gtf_file <- "path/to/your/file.gtf"
gtf_data <- fread(gtf_file, header = FALSE)  # Using fread from data.table for efficient reading

# Example: Extract relevant features (assuming basic columns in GTF format)
# Adjust column indices according to your specific GTF file structure
# For GTF format, columns typically include seqname, source, feature, start, end, score, strand, frame, attributes
# Example: Extract gene names and biotypes
genes <- gtf_data$V9  # Assuming column with gene names
biotypes <- gtf_data$V10  # Assuming column with biotypes

# Example: Print first few rows to verify data
head(gtf_data)

# Example: Print extracted features
cat("Gene names:", unique(genes), "\n")
cat("Biotypes:", unique(biotypes), "\n")

# Step 3: Define your classification task
# Example: Classify based on biotype, e.g., protein_coding vs non-coding
gtf_data$class <- ifelse(biotypes == "protein_coding", "functional", "non-functional")

# Step 4: Obtain true labels and predicted scores (hypothetical example)
# Hypothetical example: You have a method that predicts functionality based on features extracted
# For illustration, generate random scores and labels
set.seed(123)
gtf_data$predicted_score <- runif(nrow(gtf_data))
gtf_data$true_label <- ifelse(gtf_data$predicted_score > 0.5, "functional", "non-functional")

# Step 5: Calculate ROC curve
roc_curve <- roc(gtf_data$true_label, gtf_data$predicted_score)

# Step 6: Plot ROC curve
plot(roc_curve, main = "ROC Curve", col = "blue")

# Print AUC
# AUC close to 1 indicates a good classifier, while AUC close to 0.5 indicates a classifier performing no better than random
print(paste("AUC:", round(auc(roc_curve), 2)))



# SIMPLY RUNNING CODE

# Reading in TSV files 
# Define the file path to your TSV file
file_path <- "Data/Expression/Individual_datasets/candida_albicans.max.cds/SRR7685408_tpm_kallisto/abundance.tsv"

# Read the TSV file into a data frame
data <- read.delim(file_path, header = TRUE, sep = "\t")

# Print the first few rows of the data frame to check if it is loaded correctly
head(data)

# DATA MANIPULATION
head(data)  # View the first few rows
tail(data)  # View the last few rows
str(data)   # Display the structure of the data frame
summary(data)  # Provide summary statistics for each column
mean(data$Score) #calculate mean
median(data$Score) #calculate median
sd(data$Score) #calculate standard deviation 

# DATA CLEANING 
# handle missing values 
# Check for missing values
sum(is.na(data))

# Remove rows with missing values
data <- na.omit(data)

# Replace missing values with a specific value (e.g., mean of the column)
data$Age[is.na(data$Age)] <- mean(data$Age, na.rm = TRUE)

# Group by and summarize
# load libararies
library(dplyr)

# Group by a column and summarize, age is example
summary_data <- data %>%
  group_by(Age) %>%
  summarize(Average_Score = mean(Score, na.rm = TRUE))

print(summary_data)

# DATA VISUALIZATION/GRAPHING
# organize data 
# PERFORMING ROC ANALYSIS
# install packages & load 
library("pROC")
data <- read.table(file_path, header = TRUE, sep = "\t") #named data 

#check structure of data 
str(data) #structure
colnames(data) #column names
names(data) #list column names with index positions

# load libraries
library(ggplot2)
library(tidyr)

# reshape data
data_long <- gather(data, key = "Species", value = "MaxValue")
#assumes your data frame data has columns with species names (as listed in your names(data) output) and corresponding maximum values.

# Basic scatter plot
# Scatter plot of Age vs. Score (example variables)
ggplot(data, aes(x = Age, y = Score)) +
  geom_point() +
  labs(title = "Age vs. Score", x = "Age", y = "Score")

# Histogram
# (example variable)
ggplot(data, aes(x = eff_length)) +
  geom_histogram(binwidth = 5) +
  labs(title = "example", x = "eff_length", y = "length")

# Bar Plot
ggplot(data_long, aes(x = Species, y = MaxValue)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(title = "Maximum Expression Levels by Species",
       x = "Species", y = "Maximum Value")

# INDIVIDUAL DATA SETS
# install packages & load 
library("pROC")

#read the file into R, this one is tsv 
file_path <- "Data/Expression/Individual_datasets/candida_albicans.max.cds/SRR7685408_tpm_kallisto/abundance.tsv"
#read into data frame 
data <- read.delim(file_path, header = TRUE, sep = "\t")

# CREATING A SYMLINK 
# in the bash terminal, navigate to the directory where you would like the symlink to be located. 
cd repositories/mias_intragenomic_variation_mutation_bias/fasta
# create the link to the other data folder in the other repo
ln -s /~/revisit_Cope-and-Shah_2022_yeast-analysis/data/cds
# check to confirm that the link was created successfully 
ls -l 

#ENTERING AND EXITING PYTHON
exit() #to leave 