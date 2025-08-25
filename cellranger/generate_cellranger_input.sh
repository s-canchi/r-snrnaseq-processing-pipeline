#!/bin/bash

"""
Author: Saranya Canchi
Date Created: 2024-08
Description:
- Scans sample subdirectories to find and record fastq file paths.
- Produces a tab-delimited list linking each sample ID with its fastq directory for workflow input.
"""

# Set base dir containing the {id} folders
base_dir="/path/to/folder/"

# Output file to store {id} {path_to_fastq} pairs
output_file="/path/to/folder/fastq_input.txt"

# Create or clear the output file
> "$output_file"

# Loop through each {id}/fastq dir and append to the output file
for fastq_dir in "$base_dir"/*/fastq/; do
    # Check if it's a dir
    if [ -d "$fastq_dir" ]; then
        # Extract the ID from the dir name
        id=$(basename $(dirname "$fastq_dir"))
        
        # Append the id and path to the output file
        echo "$id $fastq_dir" >> "$output_file"
    fi
done

echo "Samples file created: $output_file"
