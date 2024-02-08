#!/bin/bash

# Define an array of trait names
traits=("GY" "HM")

# Loop over each trait
for trait in "${traits[@]}"; do
    echo "Processing $trait"
    
    # Create a directory for the current trait if it doesn't exist
    trait_dir="results/$trait"
    mkdir -p "$trait_dir"

    # Run Rscript with trait as an argument
    Rscript --vanilla etra.R data.csv "$trait" "$trait_dir"
done
