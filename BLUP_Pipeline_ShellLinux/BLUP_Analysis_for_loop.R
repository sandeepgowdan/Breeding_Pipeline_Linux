#!/usr/bin/env Rscript
## Usage: Rscript BLUP_analysis.R <working dir.> <input_file> <phenotypeName> <output_dir>
## Example: Rscript BLUP_analysis.R /home/sandeep/Desktop/Git/BLUP example_data/maize.csv  height results/GY

### ANOVA, MIXED MODELS for BALANCED AND UNBALANCED DATA  RCBD DATA###
## Clear the environment and set working directory
rm(list = ls())

library(lme4)
library(lmerTest)
library(lsmeans)

## Check if arguments are provided properly
# Check if command-line arguments are provided
if (length(commandArgs(trailingOnly = TRUE)) < 3) {
  stop("Usage: Rscript script.R <input_file> <trait> <output_dir>")
}

unlink(commandArgs(trailingOnly = TRUE)[3], recursive=TRUE)
dir.create(commandArgs(trailingOnly = TRUE)[3], showWarnings = TRUE)

# Get command-line arguments
input_file <- commandArgs(trailingOnly = TRUE)[1]
trait <- commandArgs(trailingOnly = TRUE)[2]
output_dir <- commandArgs(trailingOnly = TRUE)[3]

# Load the data from the input file
data <- read.csv(input_file)

# Convert REP, GEN, and ENV to factors
data$REP <- as.factor(data$REP)
data$GEN <- as.factor(data$GEN)
data$ENV <- as.factor(data$ENV)

datacheck <- table(data$GEN, data$ENV, data$REP) # Is the data balanced
write.table(datacheck, file = paste(output_dir, "/", trait, "_Check_data_balanced_or_unbalanced.csv", sep = ""), sep = ",", row.names = TRUE, col.names = TRUE)

# Plot the data location wise
# Open a PDF device to save the plots
pdf(paste(output_dir, "/", trait, "_location_wise_plots.pdf", sep = ""), width = 10)  # Adjust width as needed

# Plot the data location wise
par(mfrow = c(1,2))
plot1 <- boxplot(data[[trait]] ~ data$ENV) # effect of location
plot2 <- boxplot(data[[trait]] ~ data$REP + data$ENV) # effect of block within location

# Close the PDF device
dev.off()

### FIT THE MODEL ##

Fitted <- lmer(data[[trait]] ~ GEN + (1|ENV), data = data)

# Save the fitted model
# Capture the printed summary
summary_text <- capture.output(summary(Fitted))

# Write to a text file
writeLines(summary_text, paste(output_dir, "/", trait, "_fixed_effect_model.txt", sep = ""))

# Save the model summary coefficients to a data frame
model_summary <- data.frame(summary(Fitted)$coefficients)

# Write both the fitted model and its summary coefficients to a CSV file
write.csv(model_summary, file = paste(output_dir, "/", trait, "_fixed_model_summary_coefficients.csv", sep = ""), row.names = TRUE)

# Write model summary to a CSV file
write.csv(summary(Fitted)$coefficients, file = paste(output_dir, "/", trait, "_model_summary_fixed_effect.csv", sep = ""), row.names = TRUE)

# Extracting design matrix of fixed effects (Intercept and Genotype)
X <- model.matrix(~ GEN, data = data)
X <- X[!duplicated(X), ]

# Extracting the beta coefficients
Beta <- fixef(Fitted)

# Write the Beta coefficients to a CSV file
write.csv(Beta, file = paste(output_dir, "/", trait, "_Beta_coefficients.csv", sep = ""), row.names = TRUE)

# Obtaining the BLUEs of genotypes
BLUEs_Gen <- X %*% Beta

# Write BLUEs_Gen to a CSV file
write.csv(BLUEs_Gen, file = paste(output_dir, "/", trait, "_BLUEs_Gen.csv", sep = ""), row.names = FALSE)

# Estimating least square means by GEN
Lsmeans_Gen <- lsmeans(Fitted, ~ GEN)

# Creating the data frame of GEN and BLUEs
BLUEs_Gen1 <- data.frame(GID = levels(data$GEN),
                         BLUEs_lSmean = summary(Lsmeans_Gen)$lsmean,
                         BLUE_lmm = BLUEs_Gen
)

# BLUP of genotypes
Fitted2 <- lmer(data[[trait]] ~ (1|GEN) + (1|ENV) + (1|(GEN:ENV)), data = data)
# Save the fitted model
# Capture the printed summary
summary_text <- capture.output(summary(Fitted2))

# Write to a text file
writeLines(summary_text, paste(output_dir, "/", trait, "_random_model_summary.txt", sep = ""))

# Save the model summary coefficients to a data frame
model_summary2 <- data.frame(summary(Fitted2)$coefficients)

# Write both the fitted model and its summary coefficients to a CSV file
write.csv(model_summary2, file = paste(output_dir, "/", trait, "_random_model_summary_coefficients.csv", sep = ""), row.names = FALSE)

# Fixed effect=Intercept
Intercept <- fixef(Fitted2)
U_ref <- c(ranef(Fitted2)$GEN)

# BLUP of Genotypes
BLUP_Gen2 <- Intercept + U_ref$'(Intercept)'

# Creating a data frame of GEN and BLUPs
BLUPs_Gen1 <- data.frame(GID = levels(data$GEN),
                         BLUPs = BLUP_Gen2)

# Creating a data frame of GEN, BLUEs, and BLUPs
Combined_Gen <- data.frame(
  GID = levels(data$GEN),
  BLUEs = BLUEs_Gen1$BLUEs,
  BLUE_lmm = BLUEs_Gen,
  BLUPs = BLUP_Gen2
)

# Writing combined data to CSV file
write.csv(Combined_Gen, file = paste(output_dir, "/", trait, "_Combined_data.csv", sep = ""), row.names = FALSE)
