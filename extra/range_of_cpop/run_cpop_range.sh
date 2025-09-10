#!/bin/bash

## Define the output directory for .o and .e files 
#$ -o ./log_cpop/
#$ -e ./log_cpop/
## Define resources
#$ -l mfree=3G
#$ -pe serial 20
#$ -l h_rt=7:00:00:00
## Define working directory
#$ -cwd
## Send an email when complete
#$ -m bea
#$ -M <your_email@domain.edu>

# Activate polv_env environment
source <your_conda_env_path/envs/polv_env/bin/activate>

# Set a variable for your temporary directory using $TMPDIR
export TMPDIR="/jobs_folder"

# Run R script

Rscript 250810_range_of_cpop/250810_cpop.R --verbose