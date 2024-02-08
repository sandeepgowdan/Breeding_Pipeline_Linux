#!/usr/bin/env Rscript
## Usage: Rscript_vcf_to_Numeric.R <vcf_file> <pheno_file> <trait>
## Example: Rscript vcf-file-processing.R soynam_geno_sub.vcf  soyname_pheno_all.csv protein

### ANOVA, MIXED MODELS for BALANCED AND UNBALANCED DATA  RCBD DATA###
## Clear the environment and set working directory
rm(list = ls())
library(vcfR)
library(imputeTS)
# Read genotypic data
# Use of package vcfR to read vcf files, extract relevant information, 
# and convert it to formats usable by GAPIT
## Check if arguments are provided properly
# Check if command-line arguments are provided
if (length(commandArgs(trailingOnly = TRUE)) < 2) {
  stop("Usage: Rscript script.R <genofile.csv> <phenofile.csv>")
}

unlink("results",recursive=TRUE)
dir.create("results", showWarnings = TRUE)

geno_file <- commandArgs(trailingOnly = TRUE)[1]

#file reading
geno <- read.vcfR(geno_file)

#extract marker information from vcf file
geno1<- extract.gt(geno)


#show 10 rows and 10 columns, just to check
geno1[1:10,1:10]


# Converting marker formats to single digits (0 and 2 for homozygotes, 1 for heterozygotes)
geno2 <- ifelse(geno1 == "0/0", 0,
                ifelse(geno1 == "1/1", 2, 1))


#transpose matrix because GAPIT needs accessions in rows, markers in columns
geno3 <- t(geno2)

#conversion of file from a matrix to a data frame, to allow easy manipulation; 
#the result is creating a 
#new column with the row names
geno4 <- as.data.frame(cbind(rownames(geno3),
                             geno3))

#new column, number 1, named as "taxa"
colnames(geno4)[1] <- "taxa"

# Assuming your data frame is named 'geno_numeric2'
threshold <- 0.2  # Set the threshold for missing values (20%)

# Calculate the proportion of missing values in each column
missing_values_proportion <- colMeans(is.na(geno4), na.rm = TRUE)

# Subset the data frame to include only columns with less than or equal to 20% missing values
geno5 <- geno4[, missing_values_proportion <= threshold]


#in this case, marker data were read as text, not numbers; next line converts 
#marker columns to numbers
geno5[,2:ncol(geno5)] <- sapply(geno5[,2:ncol(geno5)],
                                as.numeric)

# Replace NA values with imputed values using linear interpolation
geno6 <- na_interpolation(geno5, option = "linear")

# Alternatively, you can use other imputation methods like 'spline', 'stine', etc.
# geno_numeric2_imputed <- na.interpolation(geno_numeric2, option = "spline")

# View the imputed data
write.csv(geno6, file = paste("results/genotype.csv"), row.names = TRUE)

#creation of new file with command from vcfR, which extracts relevant 
#information for each marker:
#identifier, chromosome, position
genome_map <- getFIX(geno)


#reorder columns as required by GAPIT: ID, chromosome, position, and rename 
#columns as required by GAPIT
genome_map1 <- as.data.frame(genome_map[,c(3,1,2)])


colnames(genome_map1) <- c("Name", "Chromosome", "Position")


# Read phenotypic data
#command read.csv2 reads csv filed with semi colon as column separator; path 
#may have to be adjusted

# Read phenotypic data using read.csv with comma as separator
pheno_file <- commandArgs(trailingOnly = TRUE)[2]
pheno_data <- read.csv(pheno_file, sep = ",")

#removes rows (genotypes) with missing data in any variable
pheno_data <- na.omit(pheno_data)

#create a new object with intersection of files geno6 and pheno_data, 
#with rows in which $id includes $taxa
geno7 <- geno6[geno6$taxa %in% pheno_data$RIL,]

#same as previous command, but opposite
pheno_data_2 <- pheno_data[pheno_data$RIL %in% geno6$taxa,]

#this is not needed, but explains the next command
summary(geno7$taxa == pheno_data_2$RIL)

#match order of files; GAPIT requires same order of genotypic and phenotypic data  
pheno_data_2 <- pheno_data_2[match(geno7$taxa, pheno_data_2$RIL),]

#change name of first column, to give the same name as recommended in the GAPIT manual
colnames(pheno_data_2)[1] <- c("Taxa")

# Identify common column names between genotpe file geno7 qnd map file 
common_columns <- intersect(colnames(geno7), genome_map1$Name)

# Subset genome_map1_map based on common column names
genome_map2 <- genome_map1[genome_map1$Name %in% common_columns, ]

write.csv(genome_map2, file = paste("results/map.csv"), row.names = TRUE)
write.csv(pheno_data_2 , file = paste("results/Phenotype.csv"), row.names = TRUE)


# Create a data frame with dimensions
dim_df <- data.frame(
  "Dataframe" = c("genotype", "map", "phenotype"),
  "Rows" = c(nrow(geno6), nrow(genome_map2), nrow(pheno_data_2)),
  "Columns" = c(ncol(geno6), ncol(genome_map2), ncol(pheno_data_2))
)

# Write the data frame to a CSV file
write.csv(dim_df, file = "results/dimensions.csv", row.names = FALSE)
