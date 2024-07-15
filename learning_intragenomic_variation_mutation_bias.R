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
#help
#needed libraries? 
install.packages("ggplot2")
install.packages("plotROC")
library(plotROC)
shiny_plotROC()

#unzipping the .gz files maybe??
#in bash: 
gunzip mias_intragenomic_variation_mutation_bias/fasta/revisit_cds_data/*.gz


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




