#!/bin/bash

# Create logs directory if it doesn't exist
mkdir -p logs

celltype=$1
assay=$2

# Define an array of groups
groups=("Infancy" "Early_Childhood" "Late_Childhood" "Adolescence" "Adulthood" "Late_Adulthood")

# Initialize a counter variable
counter=0

# Run pairwise comparisons
for i in "${!groups[@]}"; do
    for j in $(seq $((i + 1)) $((${#groups[@]} - 1))); do
        group1="${groups[$i]}"
        group2="${groups[$j]}"
        echo "Running Rscript for comparison: $group1 vs $group2"
        nohup Rscript DSS.R "$group1" "$group2" ${celltype} ${assay} &
        
        # Increment the counter
        counter=$((counter + 1))

        # Check if 3 jobs are running
        if [ $counter -eq 3 ]; then
            # Wait for all background jobs to complete
            wait
            # Reset the counter
            counter=0
        fi
    done
done

# Wait for any remaining background processes to complete
wait

echo "All Rscript jobs are done."
