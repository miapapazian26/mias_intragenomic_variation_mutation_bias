#!/usr/bin/env Rscript

#load necessary libraries
library(ggplot2)

#define the path to the results file
result_file <- "Results/ConstMut/candida_parapsilosis/results.csv"

#check if the results file exists
if (!file.exists(result_file)) {
  stop("Results file not found.")
}

#read the results file
data <- read.csv(result_file)

#plot ROC curve
ggplot(data, aes(x = FalsePositiveRate, y = TruePositiveRate)) +
  geom_line(color = "blue") +
  labs(title = "ROC Curve for candida_parapsilosis", 
       x = "False Positive Rate", 
       y = "True Positive Rate") +
  theme_minimal() +
  ggsave("ROC_Curve_candida_parapsilosis.png")

