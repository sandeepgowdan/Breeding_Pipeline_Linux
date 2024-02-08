#!/usr/bin/env Rscript
## Usage: Rscript BLUP_analysis.R <working dir.> <input_file> <phenotypeName>
## Example: Rscript BLUP_analysis.R /home/sandeep/Desktop/Git/BLUP example_data/maize.csv  height

### ANOVA, MIXED MODELS for BALANCED AND UNBALANCED DATA  RCBD DATA###
## Clear the environment and set working directory
rm(list = ls())

library(lme4)
library(lmerTest)
library(lsmeans)

## Check if arguments are provided properly
# Check if command-line arguments are provided
if (length(commandArgs(trailingOnly = TRUE)) < 2) {
  stop("Usage: Rscript script.R <input_file> <trait>")
}

unlink("results", recursive=TRUE)
dir.create("results", showWarnings = TRUE)

# Get command-line arguments
input_file <- commandArgs(trailingOnly = TRUE)[1]
trait <- commandArgs(trailingOnly = TRUE)[2]

# Load the data from the input file
data <- read.csv(input_file)

# Convert REP, GEN, and ENV to factors
data$REP <- as.factor(data$REP)
data$GEN <- as.factor(data$GEN)
data$ENV <- as.factor(data$ENV)

datacheck <- table(data$GEN, data$ENV, data$REP) # Is the data balanced
write.table(datacheck, file = paste("results/", trait, "Check_data_balanced or unbalanced.csv"), sep = ",", row.names = TRUE, col.names = TRUE)

# Plot the data location wise
# Open a PDF device to save the plots
pdf(paste("results/", trait, "location_wise_plots.pdf"), width = 10)  # Adjust width as needed

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
writeLines(summary_text, paste("results/", trait, "fixed_effect_model.txt"))

# Save the model summary coefficients to a data frame
model_summary <- data.frame(summary(Fitted)$coefficients)

# Write both the fitted model and its summary coefficients to a CSV file
write.csv(model_summary, file = paste("results/", trait, "fixed_model_summary_oefficinets.csv"), row.names = TRUE)


# Write model summary to a CSV file
write.csv(summary(Fitted)$coefficients, file = paste("results/", trait, "model_summary_fixed_effect.csv"), row.names = TRUE)

# Extracting design matrix of fixed effects (Intercept and Genotype)
X <- model.matrix(~ GEN, data = data)
X <- X[!duplicated(X), ]

# Extracting the beta coefficients
Beta <- fixef(Fitted)

# Write the Beta coefficients to a CSV file
write.csv(Beta, file = paste("results/", trait, "Beta_coefficients.csv"), row.names = TRUE)


# Obtaining the BLUEs of genotypes
BLUEs_Gen <- X %*% Beta

# Write BLUEs_Gen to a CSV file
write.csv(BLUEs_Gen, file = paste("results/", trait, "BLUEs_Gen.csv"), row.names = TRUE)


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
writeLines(summary_text, paste("results/",trait,"random_model_summary.txt"))

# Save the model summary coefficients to a data frame
model_summary2 <- data.frame(summary(Fitted2)$coefficients)

# Write both the fitted model and its summary coefficients to a CSV file
write.csv(model_summary2, file = paste("results/", trait, "random_model_coefficients.csv"), row.names = FALSE)

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
write.csv(Combined_Gen, file = paste("results/", trait, "Comb_BLUEs_BLUPs.csv"), row.names = FALSE)
