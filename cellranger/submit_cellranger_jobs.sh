#!/bin/bash

"""
Author: Saranya Canchi
Date Created: 2024-08
Description:
- Automates the creation and submission of per-sample Cell Ranger job scripts from a list of sample IDs and FASTQ paths.
- Generates individualized job and log files to streamline batch processing of snRNA-seq data.
"""

#SBATCH --account=[ACCOUNT_NAME]
#SBATCH --job-name=submit_cellranger_jobs
#SBATCH --output=slurm_logs/submit_cellranger_jobs_%j.out
#SBATCH --error=slurm_logs/submit_cellranger_jobs_%j.err
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=2G
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=standard

# Inputs
input_file="/path/to/folder/"
job_script_dir="/path/to/folder/"
cellranger_output_dir="/path/to/folder/"   
slurm_log_dir="/path/to/folder/"                  

# Ensure dir exist
mkdir -p ${job_script_dir}
mkdir -p ${cellranger_output_dir}
mkdir -p ${slurm_log_dir}

# Path to the reference transcriptome
REFERENCE="/path/to/folder/"

# Read each line from the input file and submit a job
while IFS= read -r line; do
    id=$(echo $line | awk '{print $1}')
    fastq_path=$(echo $line | awk '{print $2}')
    
    # Remove trailing slash from fastq_path 
    fastq_path=$(echo ${fastq_path%/})

    # Create subdir for SLURM logs for each id
    sample_log_dir="${slurm_log_dir}/${id}"
    mkdir -p ${sample_log_dir}
    
    # Create a job script for each id
    job_script="${job_script_dir}/${id}_cellranger_job.slurm"
    
    cat <<EOT > ${job_script}
#!/bin/bash
#SBATCH --account=[ACCOUNT_NAME]
#SBATCH --job-name=${id}_cellranger
#SBATCH --output=${sample_log_dir}/cellranger_${id}_%j.out
#SBATCH --error=${sample_log_dir}/cellranger_${id}_%j.err
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=80g
#SBATCH --mail-type=END,FAIL
#SBATCH --partition=standard

# Load modules
module load Bioinformatics
module load cellranger/7.0.0

# Change directory
cd ${cellranger_output_dir}

# Run cellranger
cellranger count --id=${id} \
--transcriptome=${REFERENCE} \
--fastqs=${fastq_path} \
--include-introns=true \
--localcores=10 \
--localmem=80 \
--disable-ui 
EOT

    # Submit the job
    sbatch ${job_script}
done < "$input_file"
