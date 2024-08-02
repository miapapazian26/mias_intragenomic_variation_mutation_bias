# LEARNING HOW TO MANUEVER THROUGH THE DATA
# Load required libraries (working)
library(seqinr)
library(data.table)

# Making a directory in R, use command 
dir.create()

# Checking content
fasta_file <- "Data/Expression/candida_albicans.max.cds"
fasta_content <- readLines(fasta_file)
head(fasta_content)

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
ls -la #list all files w details 
ls -a #list all files including hidden ones 
#check structure of data 
str(data) #structure
colnames(data) #column names
names(data) #list column names with index positions


# DATA CLEANING 
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


# PERFORMING ROC ANALYSIS
#correct input from runRoc_parallel.sh script to address the correct genome
#let the analysis finish running before continuing

#must unzip the .gz files to ensure the input is recognized correctly.
#in bash: 
gunzip mias_intragenomic_variation_mutation_bias/fasta/revisit_cds_data/*.gz
#use script runRoc_parallel and edit for which species needed to analyze 

#LOOPING 
#in bash
#loop was easily written in by coding: 
for SPECIES in "${SPECIES_LIST[@]}"; do 

INPUT="${INFOLDER}/${SPECIES}${SUFFIX}"


#geweke score: a stats test used to asses the convergence of MCMC simulations
#ensures that the posterior samples are representative of the true posterior distribution

#PLAYING WITH AnaCoDa's PLOTTING FUNCTIONS 

#analyzing a second genome 
#saved all results to its own personal folder 

#comparing parameter estimates 

#handling orthologous genes? 

#using blast 
#types 
blastn #compares a nucleotide query sequence against a nucleotide sequence database.
blastp #compares an amino acid query sequence against a protein sequence database.
blastx #compares a nucleotide query sequence translated in all reading frames against a protein sequence database.
tblastn #compares a protein query sequence against a nucleotide sequence database translated in all reading frames.
tblastx #compares the six-frame translations of a nucleotide query sequence against the six-frame translations of a nucleotide sequence database.

#RE-CREATING FIGURES 
#using the 'predicting_phi.Rmd file', we can recreate the figures from the paper. 
#saves results to an html file link, had to go through and debug. 

# PLOTTING
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
cd /repositories/mias_intragenomic_variation_mutation_bias/fasta
# create the link to the other data folder in the other repo
ln -s /~/revisit_Cope-and-Shah_2022_yeast-analysis/data/cds
# check to confirm that the link was created successfully 
ls -l 


# ENTERING AND EXITING PYTHON
exit() #to leave 

#MAKING COMMITS IN BASH
#to make commits: 
git status 
git add . #adds all changes  
git commit -m #makes commits and adds a message 


