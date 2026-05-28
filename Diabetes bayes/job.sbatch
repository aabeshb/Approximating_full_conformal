#!/bin/bash
#SBATCH --job-name=Tournament_bayes_diabetes  # A single job name for the array
#SBATCH --partition=caslake  #use caslake or br Partition
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --account=pi-rina
#SBATCH --time=04:00:00  # Running time in the format - D-HH:MM:SS
#SBATCH -o log_files/output_file_%A_%a.out # Standard output
#SBATCH -e log_files/error_file_%A_%a.err # Standard error
module load R/4.3.1
Rscript run.R