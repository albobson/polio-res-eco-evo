Last updated 2025-09-10.

This directory contains all of the scripts and data to run the analysis for *Intracellular interactions shape antiviral resistance outcomes in poliovirus via eco-evolutionary feedback*, by Robertson, Kerr and Feder (2025).

Directory structure:

-   cart: Illustrator/Affinity files for each figure that requires letters and/or cartoons.

-   dat: Reference data published by Tanner et al 2014 and Collett et al 2017

-   env: Contains the conda environment .yaml files used to run this analysis with Snakemake. This file contains all dependencies for the various R scripts to be run.

-   extra: Additional supplemental analyses are present in this directory. Each of these analyses has its own README.md to explain the order in which they should be run.

-   runs: Where the results from the pipeline are saved (This folder will contain specific sub folders based on the initial conditions of the run.) - dat_gen: All of the data generated during a pipeline run.

-   params: Data files corresponding to parameter estimates generated during the model run.

-   sims: .csv's of large simulations. These are saved here so that the simulations do not need to be run every time that the pipeline is started (given the parameters have not changed).

-   res: All of the results presented in the manuscript.

-   This directory contains a subdirectory for each main text figure and supplemental figure, as well as a .pdf of all main-text and supplemental figures together for easy viewing. Note, these will not contain the cartoons and letters. These were manually added using Affinity Designer, and those files are found in /cart/.

-   scr: A directory containing all relevant scripts to run the analysis. These will automatically be referenced by the snakemake pipeline.

In this top level directory there is a snakefile, which has all of the instructions for snakemake to execute the pipeline. The resources allocated for each job can be modified in this snakefile, depending on system requirements. The current pipeline was run with a maximum of 40 cores per job on a high performance computing cluster running Ubuntu 22.04.4. With the current parameters, a run of the pipeline takes about 8 hours to complete.

To run the pipeline, run the run_snakemake_pipeline.sh script in a terminal.

Dependencies:

-   The dependencies for R are listed in env/polv_env.ymal. This environment should be added to your conda environment list.

-   This code was executed using snakemake 6.15.1.

-   Before starting the pipeline, change the source of this conda env in run_snakemake_pipeline.sh. Likewise, change the email to alert on pipeline termination.
