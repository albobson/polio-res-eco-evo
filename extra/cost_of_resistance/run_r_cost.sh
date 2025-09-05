#!/bin/bash

## Define the output directory for .o and .e files 
#$ -o ./log_res/
#$ -e ./log_res/
## Define resources
#$ -pe serial 20
#$ -l mfree=3G
#$ -l h_rt=2:00:00:00
## Define working directory
#$ -cwd
## Send an email when complete
#$ -m bea
#$ -M alexrob@uw.edu

# Activate polv_env environment
source /net/feder/vol1/home/alexrob/mambaforge/activate polv_env

# Set a variable for your temporary directory using $TMPDIR
export TMPDIR="/jobs_folder"

## Run the R script
Rscript 250819_resistance_analysis.R