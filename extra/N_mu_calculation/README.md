This analysis assesses extinction probability as a function of gamma X mu (a supplemental figure).

To run this script, run "snake_ne_analysis.sh". This will start a snakemake pipeline, where individual jobs will be submitted, and their results will be collated.

Individual scripts for this analysis can be found in ./scr/. The results will be deposited in ./results/.

Dependencies:

-   This script relies on a previous pipeline run having been completed (the main snakemake pipeline).

-   Inside of snake_ne_analysis.sh, make sure to replace the path to the conda environment with your own path, and update the email notification settings if desired.

-   The number of jobs submitted at one time can be modified in snake_ne_analysis.sh.

Running 20 jobs at a time (with a single core each), this analysis takes about four hours to run. This will output the result of each individual simulation, as well as the figure reported in the supplement of our paper.
