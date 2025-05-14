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
#$ -M alexrob@uw.edu

# Activate snakemake environment
source /net/feder/vol1/home/alexrob/mambaforge/envs/snakemake/bin/activate

# Set a variable for your temporary directory using $TMPDIR
export TMPDIR="/jobs_folder"

## Run your Snakemake pipeline, specifying the temporary directory
snakemake --cluster "qsub -cwd -o ./log_scr/ -e ./log_scr/ -l mfree=3G -l disk_free=1G -l h_rt=0:15:0" --jobs 20