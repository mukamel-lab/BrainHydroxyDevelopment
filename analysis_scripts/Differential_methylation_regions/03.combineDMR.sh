#!/bin/bash

# Combine DMR results from all chromosomes into a single file for each pairwise comparison
groups=("Infancy" "Early_Childhood" "Late_Childhood" "Adolescence" "Adulthood" "Late_Adulthood")
celltypes=("GABA" "GLU")
assays=("BS" "OXBS")


for celltype in "${celltypes[@]}"; do
    for Assay in "${assays[@]}"; do
        for i in "${!groups[@]}"; do
            for j in $(seq $((i + 1)) $((${#groups[@]} - 1))); do
                group1="${groups[$i]}"
                group2="${groups[$j]}"
                combined_file="${Assay}_DMR_results_ExSCZ_NoNorm/combined/${group1}_${group2}_${celltype}_${Assay}_DMR.csv"
                > "$combined_file"

                # Variable to keep track of whether the header has been copied
                header_copied=0
                echo ${combined_file}
                for chr in {1..22} X; do
                    file_name="${Assay}_DMR_results_ExSCZ_NoNorm/${group1}.${group2}.${celltype}_chr${chr}_DMR_unadjusted.csv"

                    # Check if the file exists
                    if [[ -f "$file_name" ]]; then
                        # If the header has not been copied yet, copy it
                        if [[ "$header_copied" -eq 0 ]]; then
                            head -n 1 "$file_name" > "$combined_file"
                            header_copied=1
                        fi
                        # Append the content without the header
                        tail -n +2 "$file_name" >> "$combined_file"
                    fi
                done
            done
        done
    done
done
