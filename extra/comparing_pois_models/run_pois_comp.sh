#!/bin/bash

## Define the output directory for .o and .e files 
#$ -o ./log_pois/
#$ -e ./log_pois/
## Define resources
#$ -pe serial 20
#$ -l mfree=3G
#$ -l h_rt=2:00:00:00
## Define working directory
#$ -cwd
## Send an email when complete
#$ -m bea
#$ -M <your_email@domain.edu>

# Activate polv_env environment
source <path_to_your_conda_env_here/envs/polv_env/bin/activate>

# Set a variable for your temporary directory using $TMPDIR
export TMPDIR="/jobs_folder"

## Run the R script
Rscript 250818_comparing_pois_models.R