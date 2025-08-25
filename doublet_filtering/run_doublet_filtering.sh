#!/bin/bash

"""
Author: Saranya Canchi
Date Created: 2024-08
Description:
- Submits jobs to run doublet detection for multiple snRNA-seq datasets.
- Handles environment setup and logging for each dataset.  
"""


# Base dir
output_dir="/path/to/folder/"

# Define dataset prefixes 
dataset_prefixes=("prefix1" "prefix2")

for prefix in "${dataset_prefixes[@]}"
do

  sbatch <<EOF
#!/bin/bash
#SBATCH --account=[ACCOUNT_NAME]
#SBATCH --job-name=${prefix}_doublet
#SBATCH --output=${output_dir}/slurm_logs/${prefix}_doublet_%j.out 
#SBATCH --error=${output_dir}/slurm_logs/${prefix}_doublet_%j.err  
#SBATCH --time=24:00:00                 
#SBATCH --partition=standard            
#SBATCH --ntasks=1                      
#SBATCH --cpus-per-task=2               
#SBATCH --mem=120G                       
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
