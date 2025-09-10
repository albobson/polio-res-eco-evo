#!/bin/bash

## Define the output directory for .o and .e files 
#$ -o ./log_snake/
#$ -e ./log_snake/
## Define resources
#$ -l mfree=2G
#$ -l h_rt=4:00:00:00
## Define working directory
#$ -cwd
## Send an email when complete
#$ -m bea
#$ -M <your_email@domain.edu>

# Activate polv_env environment
source <your_conda_env_path/envs/polv_env/bin/activate>

# Set a variable for your temporary directory using $TMPDIR
export TMPDIR="/jobs_folder"

## Run your Snakemake pipeline, specifying the temporary directory
snakemake --cluster "qsub -cwd -o {resources.log_loc} -e {resources.log_loc} -l mfree={resources.mfree} -l disk_free={resources.disk_free} -l h_rt={resources.cluster_time}" --jobs 20