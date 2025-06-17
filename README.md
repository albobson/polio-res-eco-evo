This directory contains all of the scripts and data to run the analysis for _Intracellular interactions shape antiviral resistance outcomes in poliovirus via eco-evolutionary feedback_, by Robertson, Kerr and Feder (2025).

Directory structure:
- cart: Illustrator/Affinity files for each figure that requires letters and/or cartoons. 
- dat: Reference data from Tanner et al 2014 and Collett et al 2017
- env: Contains the conda environment .yaml files used to run this analysis with Snakemake
- extra: Files for extra analyses that will be added to the Snakemake pipeline (eventually)
- runs: Where the results from the pipeline are saved (This folder will contain specific sub folders based on the initial conditions of the run.)
    - dat_gen: All of the data generated during a pipeline run. 
        - params: Parameter estimates
        - sims: .csv's of large simulations that take a long time to be saved.
    - res: All of the results presented in the manuscript
        - This directory contains a subdirectory for each figure, as well as a .pdf of all main-text and supplemental figures together.
- scr: Scripts containing the functions referenced in the main analysis scripts

In this main directory, there is a snakefile, which has all of the instructions for snakemake to execute the pipeline.

To run the pipeline, run the run_snakemake_pipeline.sh script in a terminal. 
