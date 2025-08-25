#!/bin/bash

"""
Author: Saranya Canchi
Date Created: 2024-08
Description:
- Submits parallel SLURM jobs to process and merge raw RNA-seq data for multiple dataset prefixes.
- Automates environment setup, logging, and execution of an R script for batch data processing.
"""

# Base dir
output_dir="/path/to/folder/"

# Create logs dir
mkdir -p ${output_dir}/slurm_logs

# Define dataset prefixes 
dataset_prefixes=("prefix1" "prefix2")

for prefix in "${dataset_prefixes[@]}"
do

  sbatch <<EOF
#!/bin/bash
#SBATCH --account=[ACCOUNT_NAME]
#SBATCH --job-name=${prefix}_merge_raw_rna 
#SBATCH --output=${output_dir}/slurm_logs/${prefix}_merge_raw_rna_%j.out 
#SBATCH --error=${output_dir}/slurm_logs/${prefix}_merge_raw_rna_%j.err  
#SBATCH --time=12:00:00                 
#SBATCH --partition=standard            
#SBATCH --ntasks=1                      
#SBATCH --cpus-per-task=6               
#SBATCH --mem=100G                       
#SBATCH --mail-type=END,FAIL            

# Set wd
cd ${output_dir}

# Activate conda env
source /path/to/conda/
conda activate [ENV_NAME]

# Run the R script
Rscript /path/to/folder/ ${prefix}

# Deactivate env
conda deactivate
EOF
done
