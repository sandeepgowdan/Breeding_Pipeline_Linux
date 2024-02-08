#!/usr/bin/env Rscript

# Load required library

# Retrieve command-line arguments
args <- commandArgs(trailingOnly = TRUE)

## Check and create 'results' directory if not exists
if (!file.exists("results")) {
  dir.create("results")
}

# Check if there are at least two command-line arguments (input and output files)
if (length(args) < 2) {
  stop("Usage: Rscript correlation.R input_file output_file")
}

# Access input and output file names
input_file <- args[1]
output_file <- args[2]

# Load data from the input file (assuming it's a CSV file)
data <- read.csv(input_file)

data <- data[4:12]
# Compute correlation matrix
cor_matrix <- cor(data)

# Save the correlation matrix to the specified output file
file <- paste("results", paste("pearson correlation.csv", sep = "_"), sep = "/")
write.table(cor_matrix, file = file, sep = ",", row.names = TRUE, col.names = TRUE)