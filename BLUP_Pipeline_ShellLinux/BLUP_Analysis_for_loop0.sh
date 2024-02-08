#!/bin/bash

# Define an array of trait names
traits=("GY" "HM")

# Loop over each trait
for trait in "${traits[@]}"; do
    echo "Processing $trait"
    
    # Run Rscript with trait as an argument
    Rscript --vanilla BLUP_Analysis.R data.csv "$trait"
done
