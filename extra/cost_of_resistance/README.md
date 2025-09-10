This analysis assesses the impact that an added cost of resistance would have on the model and is reported in a supplemental figure.

To run this script, run "run_r_cost.sh".

Dependencies:

-   This script relies on a previous pipeline run having been completed (the main snakemake pipeline).

<!-- -->

-   Inside of run_r_cost.sh, make sure to replace the path to the conda environment with your own path, and update the email notification settings if desired.

<!-- -->

-   If you have updated any parameters (including the date run) from the snakemake pipeline, those will need to be adjusted in the 250819_resistance_analysis.R file. It is currently set to run with 40 cores, so adjust as needed.

With 40 cores, this analysis takes about two hours to run. This will create individual folders containing all of the results, as well as a combined figure in this directory. Letters can be added to this figure with the .afdesign fild in this directory.
