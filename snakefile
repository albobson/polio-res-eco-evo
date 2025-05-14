## 25-03-24
## TO DO

    ## Split generate figure 3A script into two pieces? Right now the fitness
    ## function is generated and the plot is generated
    
    ## Same as above for trajectory in fig 5C and fig 6C
    

#### Intro ####

## Central pipeline for running the analyses and generating figures for this
## manuscript. There are 2 input parameters that must be set at the beginning of
## the simulation:

  ## 1) fit_func_use - The fitness function to be used
  
    ## This is a descriptive title at the moment.
      
  ## 2) mu - The mutation rate of poliovirus
  
    ## This was inferred by assessing the mutations which confer resistance, and
    ## the possible ways to get there, given the WT amino acid. This was
    ## calculated to be 2x10-5 (see Supplemental Materials and Methods).
    
## Based on the input parameters, the model will run optimizations to determine
## the average burst size per cell and the P:PFU ratio, matching to Tanner et
## al. (2014) mixed cell culture data. Separately, it will optimize for the
## parameters that govern immune clearance, comparing to Collete et al. (2017). 

## It will then run the subsequent analyses and plot their outputs in their
## respective folders, based on the ordering of the figures in the manuscript.

## All of the figures and parameter optimization tables will be saved in a
## respective sub-folder in the /res/ folder. The sub-folder will be named
## according to the fitness function used and the mutation rate. A single file
## titled "model_out.csv" will contain the date that the model was run, the
## input parameters, and the optimized values.

#### Environment                                                            ####
## Set Conda environment globally
conda: "env/polv_env.yml"


#### Define input parameters                                                ####
## Get the current date
from datetime import datetime
## Store run date
# init_date = datetime.today().strftime('%Y-%m-%d')

## At the moment, this will need to be manually set, as I was having trouble
## when scripts took more than one day to complete.
init_date = "2025-05-09" 

#### Simulation parameters                                                  ####
## Fitness function to use
  ## This is a descriptive parameter to keep track of what sort of fitness
  ## function was used to generate the data.
init_fit_func = ["logistic"]

## Mutation rate to use
init_mu = [2e-05]


#### Plotting parameters                                                    ####
## The colors that will be used for resistant and susceptible:
rs_colors = ["#c94d4d", "#688cc8"]

## Axis text size
axis_text_size = [6]

## Legend text size
lege_text_size = [6]

## fig 2 width and height
fig2_dim = [8, 6.4]

## fig 3 width and height
fig3_dim = [8, 8]

## fig 4 width and height
fig4_dim = [18, 9]

## fig 5 width and height
fig5_dim = [11, 9]

## sfig 1 width and height
sfig1_dim = [8, 6.4]

## sfig 2 width and height
sfig2_dim = [18, 12]

## sfig 3 width and height
sfig3_dim = [11, 11]

## sfig 4 width and height
sfig4_dim = [9, 7]

## fig 5 width and height
sfig5_dim = [11, 9]

## fig 7 width and height
sfig6_dim = [11, 4.5]

## The name for the folder where this run will be stored:
## For spatial simulations
out_file_name = "runs/ddt_" + init_fit_func[0] + "_mu_" + str(init_mu[0]) + "_" + init_date + "/"

## Create this new directory
from pathlib import Path
Path(out_file_name).mkdir(parents=True, exist_ok=True)

## For the drug stringency and dominance analyses
## This changes how many trials are run
init_n_trials = [4]
## Number of individuals per trial
init_indiv = [100]

## The percent of individuals in each trial to be administered pocap at 24 hpi.
  ## Here, all are receiving drug at 24 hpi
init_perc_24 = [1]


#### Define Rules for Snakemake                                             ####
## Rule All: Tells Snakemake what we want in the end of the run
rule all:
    input:
        ## This is the collated file which I want in the end 
        out_file_name + "res/all_figs.pdf",
        out_file_name + "res/all_sfigs.pdf"

################################################################################
######################## Rules to run the optimizations ########################
################################################################################

## Rule fit_func -- Create the fitness function which will be used by all simulations
rule fit_func:
    input:
       ## These are the cell culture dataframes
        "dat/tanner_pure_culture.csv",
        "dat/tanner_mix_culture.csv"
    output:
        out_file_name + "dat_gen/params/logistic_fitness_function.csv",
        out_file_name + "dat_gen/params/tan_fit_df.csv",
        out_file_name + "dat_gen/params/logistic_fit_summary.csv"
    params:
        date = init_date,
        fit_func = init_fit_func,
        mu = init_mu
    resources:
        log_loc=str('./log_scr/'),
        mfree=str('3G'),
        cluster_time=str('0:00:10:00'),
        mem_mb=str('3G'),
        disk_free=str('2G'),
        runtime=str('1h'),
        n_cores=str('1')
    script:
        'scr/optimizations/gen_fit_func.R'

## Rule opt_model -- Based on the fitness function, find optimum model params
rule opt_model:
    input:
        out_file_name + "dat_gen/params/logistic_fitness_function.csv"
    output:
        out_file_name + "dat_gen/params/optim_params.csv"
    params:
        date = init_date,
        fit_func = init_fit_func,
        mu = init_mu
    resources:
        log_loc=str('./log_scr/'),
        mfree=str('3G'),
        cluster_time=str('0:02:00:00'),
        mem_mb=str('3G'),
        disk_free=str('4G'),
        runtime=str('2h'),
        n_cores=str('7') ## Only really need 7 because we're running 6 sims each
    script:
        'scr/optimizations/core_parameter_optimization.R'


## Rule optim_clin_trial_params -- Find the optimum params for the clin. trial
rule optim_clin_trial_params:
    input:
        "dat/collett_trial.csv", ## real clinical trial
        out_file_name + "dat_gen/params/optim_params.csv",
        out_file_name + "dat_gen/params/logistic_fitness_function.csv"
    output:
        out_file_name + "dat_gen/params/optim_clin_trial_params.rds",
        out_file_name + "dat_gen/params/optim_clin_trial_params.csv"
    params:
        date = init_date,
        fit_func = init_fit_func,
        mu = init_mu
    resources:
        log_loc=str('./log_scr/'),
        mfree=str('3G'),
        cluster_time=str('4:00:00:00'),
        mem_mb=str('3G'),
        disk_free=str('4G'),
        disk_mb=str('4G'),
        runtime=str('4d'),
        n_cores=str('20-25')
    script:
        'scr/optimizations/sim_clin_trial_param.R'

################################################################################
######################### Rules to run big simulations #########################
################################################################################

## Rule sim_clin_trial
rule sim_clin_trial:
    input:
        "dat/collett_trial.csv", ## real clinical trial
        out_file_name + "dat_gen/params/optim_params.csv",
        out_file_name + "dat_gen/params/logistic_fitness_function.csv",
        out_file_name + "dat_gen/params/optim_clin_trial_params.csv" ## parameters for sims
    output:
        out_file_name + "dat_gen/sims/sim_clinical_trial.csv"
    params:
        date = init_date,
        fit_func = init_fit_func,
        mu = init_mu
    resources:
        log_loc=str('./log_scr/'),
        mfree=str('3G'),
        cluster_time=str('0:02:00:00'),
        mem_mb=str('3G'),
        disk_free=str('4G'),
        disk_mb=str('4G'),
        runtime=str('2h'),
        n_cores=str('10')
    script:
        'scr/sims/simulated_clinical_trial.R'
        

## Rule stringency_trials -- Run clinical trials scross n_trial stringencies
rule stringency_trials:
    input:
        out_file_name + "dat_gen/params/optim_params.csv",
        out_file_name + "dat_gen/params/logistic_fitness_function.csv",
        out_file_name + "dat_gen/params/optim_clin_trial_params.csv" ## parameters for sims
    output:
        out_file_name + "dat_gen/sims/full_stringency_trials.csv"
    params:
        date = init_date,
        fit_func = init_fit_func,
        mu = init_mu,
        n_trials = init_n_trials, ## This changes how many trials are run
        n_indiv = init_indiv,  ## Number of individuals per trial
        perc_24 = init_perc_24
    resources:
        log_loc=str('./log_scr/'),
        mfree=str('3G'),
        cluster_time=str('3:00:00:00'),
        mem_mb=str('3G'),
        disk_free=str('4G'),
        disk_mb=str('4G'),
        runtime=str('3d'),
        n_cores=str('20')
    script:
        'scr/sims/drug_stringency_clin_runs.R'
        
## Rule all_stringency_trials -- Run clinical trials along all stringencies
rule all_stringency_trials:
    input:
        out_file_name + "dat_gen/params/optim_params.csv",
        out_file_name + "dat_gen/params/logistic_fitness_function.csv",
        out_file_name + "dat_gen/params/optim_clin_trial_params.csv" ## parameters for sims
    output:
        out_file_name + "dat_gen/sims/all_full_stringency_trials.csv"
    params:
        date = init_date,
        fit_func = init_fit_func,
        mu = init_mu
    resources:
        log_loc=str('./log_scr/'),
        mfree=str('3G'),
        cluster_time=str('3:00:00:00'),
        mem_mb=str('3G'),
        disk_free=str('4G'),
        disk_mb=str('4G'),
        runtime=str('3d'),
        n_cores=str('20')
    script:
        'scr/sims/full_stringency_clin_runs.R'
        
## Rule dominance_trials -- Run clinical trials with different drug stringencies
rule dominance_trials:
    input:
        out_file_name + "dat_gen/params/optim_params.csv",
        out_file_name + "dat_gen/params/logistic_fitness_function.csv",
        out_file_name + "dat_gen/params/optim_clin_trial_params.csv" ## parameters for sims
    output:
        out_file_name + "dat_gen/sims/full_dominance_trials.csv"
    params:
        date = init_date,
        fit_func = init_fit_func,
        mu = init_mu,
        n_trials = init_n_trials, ## This changes how many trials are run
        n_indiv = init_indiv,  ## Number of individuals per trial
        perc_24 = init_perc_24
    resources:
        log_loc=str('./log_scr/'),
        mfree=str('3G'),
        cluster_time=str('3:00:00:00'),
        mem_mb=str('3G'),
        disk_free=str('4G'),
        disk_mb=str('4G'),
        runtime=str('3d'),
        n_cores=str('20')
    script:
        'scr/sims/drug_dominance_clin_runs.R'

################################################################################
##################### Rules to generate plots for the data #####################
################################################################################

#### Figure 2                                                               ####
## Rule fig2A -- Generate figure 2A - Partial cell culture
rule fig2A:
    input:
        "dat/tanner_pure_culture.csv", ## These are the cell culture dataframes
        "dat/tanner_mix_culture.csv"
    output:
        out_file_name + "res/fig2/A.rds"
    params:
        date = init_date,
        fit_func = init_fit_func,
        mu = init_mu,
        rs_colors = rs_colors,
        axis_text_size = axis_text_size,
        lege_text_size = lege_text_size,
        fig2_dim = fig2_dim
    resources:
        log_loc=str('./log_scr/'),
        mfree=str('3G'),
        cluster_time=str('0:00:10:00'),
        mem_mb=str('3G'),
        disk_free=str('4G'),
        disk_mb=str('4G'),
        runtime=str('10m'),
        n_cores=str('1')
    script:
        'scr/panels/figure2A_partial_cell_culture_data.R'

## Rule fig2B -- Generature figure 2B - Prev. Clinical Trial
rule fig2B:
    input:
        "dat/collett_trial.csv"
    output:
        out_file_name + "res/fig2/B.rds"
    params:
        date = init_date,
        fit_func = init_fit_func,
        mu = init_mu,
        rs_colors = rs_colors,
        axis_text_size = axis_text_size,
        lege_text_size = lege_text_size,
        fig2_dim = fig2_dim
    resources:
        log_loc=str('./log_scr/'),
        mfree=str('3G'),
        cluster_time=str('0:00:10:00'),
        mem_mb=str('3G'),
        disk_free=str('4G'),
        disk_mb=str('4G'),
        runtime=str('10m'),
        n_cores=str('1')
    script:
        'scr/panels/figure2B_clinical_trial_data.R'


#### Figure 3                                                               ####
## Rule fig3B -- Generate figure 3B - Fitness Function
rule fig3B:
    input:
      ## This is the fitness func. csv
        out_file_name + "dat_gen/params/logistic_fitness_function.csv",
      ## This is the Tanner data
        out_file_name + "dat_gen/params/tan_fit_df.csv"
    output:
        out_file_name + "res/fig3/B.rds"
    params:
        date = init_date,
        fit_func = init_fit_func,
        mu = init_mu,
        rs_colors = rs_colors,
        axis_text_size = axis_text_size,
        lege_text_size = lege_text_size,
        fig3_dim = fig3_dim
    resources:
        log_loc=str('./log_scr/'),
        mfree=str('3G'),
        cluster_time=str('0:00:10:00'),
        mem_mb=str('3G'),
        disk_free=str('4G'),
        disk_mb=str('4G'),
        runtime=str('10m'),
        n_cores=str('1')
    script:
        'scr/panels/figure3B_fitness_function.R'


## Rule fig3C -- Generate figure 3C - Full Cell Culture
rule fig3C:
    input:
        "dat/tanner_pure_culture.csv", ## These are the cell culture dataframes
        "dat/tanner_mix_culture.csv"
    output:
        out_file_name + "res/fig3/C.rds"
    params:
        date = init_date,
        fit_func = init_fit_func,
        mu = init_mu,
        rs_colors = rs_colors,
        axis_text_size = axis_text_size,
        lege_text_size = lege_text_size,
        fig3_dim = fig3_dim
    resources:
        log_loc=str('./log_scr/'),
        mfree=str('3G'),
        cluster_time=str('0:00:10:00'),
        mem_mb=str('3G'),
        disk_free=str('4G'),
        disk_mb=str('4G'),
        runtime=str('10m'),
        n_cores=str('1')
    script:
        'scr/panels/figure3C_full_cell_culture.R'

## Rule fig3D -- Generate figure 3D - Simulated cell culture
rule fig3D:
    ## This requires that both the fitness function and the optimized
    ## parameters have been generated
    input:
        out_file_name + "dat_gen/params/optim_params.csv",
        out_file_name + "dat_gen/params/logistic_fitness_function.csv"
    output:
        out_file_name + "res/fig3/D.rds"
    params:
        date = init_date,
        fit_func = init_fit_func,
        mu = init_mu,
        rs_colors = rs_colors,
        axis_text_size = axis_text_size,
        lege_text_size = lege_text_size,
        fig3_dim = fig3_dim
    resources:
        log_loc=str('./log_scr/'),
        mfree=str('3G'),
        cluster_time=str('0:00:10:00'),
        mem_mb=str('3G'),
        disk_free=str('4G'),
        disk_mb=str('4G'),
        runtime=str('10m'),
        n_cores=str('7')
    script:
        'scr/panels/figure3D_simulated_cell_culture.R'


#### Figure 4                                                               ####
## Rule fig4A -- Generate figure 4A - Res. suppressed at high MOI
rule fig4A:
    ## This requires that both the fitness function and the optimized
    ## parameters have been generated
    input:
        out_file_name + "dat_gen/params/optim_params.csv",
        out_file_name + "dat_gen/params/logistic_fitness_function.csv"
    output:
        out_file_name + "res/fig4/A.rds"
    params:
        date = init_date,
        fit_func = init_fit_func,
        mu = init_mu,
        rs_colors = rs_colors,
        axis_text_size = axis_text_size,
        lege_text_size = lege_text_size,
        fig4_dim = fig4_dim
    resources:
        log_loc=str('./log_scr/'),
        mfree=str('3G'),
        cluster_time=str('0:00:30:00'),
        mem_mb=str('3G'),
        disk_free=str('4G'),
        disk_mb=str('4G'),
        runtime=str('30m'),
        n_cores=str('10')
    script:
        'scr/panels/figure4A_res_suppressed_at_high_moi.R'

## Rule fig4D -- Generate figure 4D - Trajectory plot
rule fig4D:
    input:
        out_file_name + "dat_gen/params/optim_params.csv",
        out_file_name + "dat_gen/params/logistic_fitness_function.csv"
    output:
        out_file_name + "res/fig4/D.rds",
        out_file_name + "dat_gen/sims/traj_for_phase_plane.csv" ## Traj. path in phase plane
    params:
        date = init_date,
        fit_func = init_fit_func,
        mu = init_mu,
        rs_colors = rs_colors,
        axis_text_size = axis_text_size,
        lege_text_size = lege_text_size,
        fig4_dim = fig4_dim
    resources:
        log_loc=str('./log_scr/'),
        mfree=str('3G'),
        cluster_time=str('0:00:10:00'),
        mem_mb=str('3G'),
        disk_free=str('4G'),
        disk_mb=str('4G'),
        runtime=str('10m'),
        n_cores=str('1')
    script:
        'scr/panels/figure4D_trajectory_plot.R'

## Rule fig4E -- Generate figure 4E - Phase plane
rule fig4E:
    input:
        out_file_name + "dat_gen/params/optim_params.csv",
        out_file_name + "dat_gen/params/logistic_fitness_function.csv",
        out_file_name + "dat_gen/sims/traj_for_phase_plane.csv" 
            ## This is the trajectory path through the phase plane from above
    output:
        out_file_name + "res/fig4/E.rds"
    params:
        date = init_date,
        fit_func = init_fit_func,
        mu = init_mu,
        rs_colors = rs_colors,
        axis_text_size = axis_text_size,
        lege_text_size = lege_text_size,
        fig4_dim = fig4_dim
    resources:
        log_loc=str('./log_scr/'),
        mfree=str('3G'),
        cluster_time=str('0:00:45:00'),
        mem_mb=str('3G'),
        disk_free=str('4G'),
        disk_mb=str('4G'),
        runtime=str('45m'),
        n_cores=str('10')
    script:
        'scr/panels/figure4E_phase_plane.R'


## Rule fig4F -- Generate figure 4F - Simulated clinical trial vs real trial
rule fig4F:
    input:
        "dat/collett_trial.csv", ## real clinical trial
        out_file_name + "dat_gen/params/optim_params.csv",
        out_file_name + "dat_gen/params/logistic_fitness_function.csv",
        out_file_name + "dat_gen/sims/sim_clinical_trial.csv" ## Simulated clinical trial
    output:
        out_file_name + "res/fig4/F.rds"
    params:
        date = init_date,
        fit_func = init_fit_func,
        mu = init_mu,
        rs_colors = rs_colors,
        axis_text_size = axis_text_size,
        lege_text_size = lege_text_size,
        fig4_dim = fig4_dim
    resources:
        log_loc=str('./log_scr/'),
        mfree=str('3G'),
        cluster_time=str('0:00:10:00'),
        mem_mb=str('3G'),
        disk_free=str('4G'),
        disk_mb=str('4G'),
        runtime=str('10m'),
        n_cores=str('1')
    script:
        'scr/panels/figure4F_simulated_trial_vs_real.R'


## Rule fig4G -- Generate figure 4G - Late clearer in the clinical trial
rule fig4GH:
    input:
        "dat/collett_trial.csv", ## real clinical trial
        out_file_name + "dat_gen/params/optim_params.csv",
        out_file_name + "dat_gen/params/logistic_fitness_function.csv",
        out_file_name + "dat_gen/sims/sim_clinical_trial.csv" ## Simulated clinical trial
    output:
        out_file_name + "res/fig4/GH.rds"
    params:
        date = init_date,
        fit_func = init_fit_func,
        mu = init_mu,
        rs_colors = rs_colors,
        axis_text_size = axis_text_size,
        lege_text_size = lege_text_size,
        fig4_dim = fig4_dim
    resources:
        log_loc=str('./log_scr/'),
        mfree=str('3G'),
        cluster_time=str('0:00:10:00'),
        mem_mb=str('3G'),
        disk_free=str('4G'),
        disk_mb=str('4G'),
        runtime=str('10m'),
        n_cores=str('1')
    script:
        'scr/panels/figure4GH_early_late_traj.R'



#### Figure 5                                                               ####
## Rule fig5A -- Generate 5A - The range of less stringent drugs
rule fig5A:
    input:
        out_file_name + "dat_gen/params/optim_params.csv",
        out_file_name + "dat_gen/params/logistic_fitness_function.csv"
    output:
        out_file_name + "res/fig5/A.rds"
    params:
        date = init_date,
        fit_func = init_fit_func,
        mu = init_mu,
        rs_colors = rs_colors,
        axis_text_size = axis_text_size,
        lege_text_size = lege_text_size,
        fig5_dim = fig5_dim,
        n_trials = init_n_trials ## This changes how many trials are run
    resources:
        log_loc=str('./log_scr/'),
        mfree=str('3G'),
        cluster_time=str('0:00:10:00'),
        mem_mb=str('3G'),
        disk_free=str('4G'),
        disk_mb=str('4G'),
        runtime=str('10m'),
        n_cores=str('1')
    script:
        'scr/panels/figure5A_range_of_worse.R'


## Rule fig5B -- Generate figures 5B - Trajectory of less stringent drug
rule fig5B:
    input:
        out_file_name + "dat_gen/params/optim_params.csv",
        out_file_name + "dat_gen/params/logistic_fitness_function.csv"
    output:
        out_file_name + "res/fig5/B.rds"
    params:
        date = init_date,
        fit_func = init_fit_func,
        mu = init_mu,
        rs_colors = rs_colors,
        axis_text_size = axis_text_size,
        lege_text_size = lege_text_size,
        fig5_dim = fig5_dim
    resources:
        log_loc=str('./log_scr/'),
        mfree=str('3G'),
        cluster_time=str('0:00:10:00'),
        mem_mb=str('3G'),
        disk_free=str('4G'),
        disk_mb=str('4G'),
        runtime=str('10m'),
        n_cores=str('1')
    script:
        'scr/panels/figure5B_trajectory.R'


## Rule fig5CD -- Generate figures 5CD - All less stringent drugs trials
rule fig5CD:
    input:
        out_file_name + "dat_gen/params/optim_params.csv",
        out_file_name + "dat_gen/params/logistic_fitness_function.csv",
        out_file_name + "dat_gen/sims/full_stringency_trials.csv"
    output:
        out_file_name + "res/fig5/C.rds",
        out_file_name + "res/fig5/D.rds"
    params:
        date = init_date,
        fit_func = init_fit_func,
        mu = init_mu,
        rs_colors = rs_colors,
        axis_text_size = axis_text_size,
        lege_text_size = lege_text_size,
        fig5_dim = fig5_dim
    resources:
        log_loc=str('./log_scr/'),
        mfree=str('3G'),
        cluster_time=str('0:00:10:00'),
        mem_mb=str('3G'),
        disk_free=str('4G'),
        disk_mb=str('4G'),
        runtime=str('10m'),
        n_cores=str('1')
    script:
        'scr/panels/figure5CD_less_str_drug_vs_pocap.R'

################################################################################
############## Rules to generate full figures from previous plots ##############
################################################################################
rule gen_fig2:
    input:
        out_file_name + "res/fig2/A.rds",
        out_file_name + "res/fig2/B.rds"
    output:
        out_file_name + "res/fig2/fig2.pdf"
    params:
        date = init_date,
        fit_func = init_fit_func,
        mu = init_mu,
        rs_colors = rs_colors,
        axis_text_size = axis_text_size,
        lege_text_size = lege_text_size,
        fig2_dim = fig2_dim
    resources:
        log_loc=str('./log_scr/'),
        mfree=str('3G'),
        cluster_time=str('0:00:10:00'),
        mem_mb=str('3G'),
        disk_free=str('4G'),
        disk_mb=str('4G'),
        runtime=str('10m'),
        n_cores=str('1')
    script:
        'scr/full_figures/gen_fig_2.R'

rule gen_fig3:
    input:
        out_file_name + "res/fig3/B.rds",
        out_file_name + "res/fig3/C.rds",
        out_file_name + "res/fig3/D.rds"
    output:
        out_file_name + "res/fig3/fig3.pdf"
    params:
        date = init_date,
        fit_func = init_fit_func,
        mu = init_mu,
        rs_colors = rs_colors,
        axis_text_size = axis_text_size,
        lege_text_size = lege_text_size,
        fig3_dim = fig3_dim
    resources:
        log_loc=str('./log_scr/'),
        mfree=str('3G'),
        cluster_time=str('0:00:10:00'),
        mem_mb=str('3G'),
        disk_free=str('4G'),
        disk_mb=str('4G'),
        runtime=str('10m'),
        n_cores=str('1')
    script:
        'scr/full_figures/gen_fig_3.R'

rule gen_fig4:
    input:
        out_file_name + "res/fig4/A.rds",
        out_file_name + "res/fig4/D.rds",
        out_file_name + "res/fig4/E.rds",
        out_file_name + "res/fig4/F.rds",
        out_file_name + "res/fig4/GH.rds"
    output:
        out_file_name + "res/fig4/fig4.pdf"
    params:
        date = init_date,
        fit_func = init_fit_func,
        mu = init_mu,
        rs_colors = rs_colors,
        axis_text_size = axis_text_size,
        lege_text_size = lege_text_size,
        fig4_dim = fig4_dim
    resources:
        log_loc=str('./log_scr/'),
        mfree=str('3G'),
        cluster_time=str('0:00:10:00'),
        mem_mb=str('3G'),
        disk_free=str('4G'),
        disk_mb=str('4G'),
        runtime=str('10m'),
        n_cores=str('1')
    script:
        'scr/full_figures/gen_fig_4.R'

rule gen_fig5:
    input:
        out_file_name + "res/fig5/A.rds",
        out_file_name + "res/fig5/B.rds",
        out_file_name + "res/fig5/C.rds",
        out_file_name + "res/fig5/D.rds"
    output:
        out_file_name + "res/fig5/fig5.pdf"
    params:
        date = init_date,
        fit_func = init_fit_func,
        mu = init_mu,
        rs_colors = rs_colors,
        axis_text_size = axis_text_size,
        lege_text_size = lege_text_size,
        fig5_dim = fig5_dim
    resources:
        log_loc=str('./log_scr/'),
        mfree=str('3G'),
        cluster_time=str('0:00:10:00'),
        mem_mb=str('3G'),
        disk_free=str('4G'),
        disk_mb=str('4G'),
        runtime=str('10m'),
        n_cores=str('1')
    script:
        'scr/full_figures/gen_fig_5.R'


#### Pull all figures into one document                                     ####
rule collate_all:
    input:
        out_file_name + "res/fig2/fig2.pdf",
        out_file_name + "res/fig3/fig3.pdf",
        out_file_name + "res/fig4/fig4.pdf",
        out_file_name + "res/fig5/fig5.pdf"
    output:
        out_file_name + "res/all_figs.pdf"
    params:
        date = init_date,
        fit_func = init_fit_func,
        mu = init_mu
    resources:
        log_loc=str('./log_scr/'),
        mfree=str('3G'),
        cluster_time=str('0:00:10:00'),
        mem_mb=str('3G'),
        disk_free=str('4G'),
        disk_mb=str('4G'),
        runtime=str('10m'),
        n_cores=str('1')
    script:
        'scr/full_figures/collate_plots.R'


################################################################################
########################### Supplemental figures ###############################
################################################################################

## Rule sfig1 -- Generate sup. fig 1 - MOI controls VPC
rule sfig1:
    input:
    output:
        out_file_name + "res/sup/sfig1/sfig1.pdf"
    params:
        date = init_date,
        fit_func = init_fit_func,
        mu = init_mu,
        rs_colors = rs_colors,
        axis_text_size = axis_text_size,
        lege_text_size = lege_text_size,
        sfig1_dim = sfig1_dim
    resources:
        log_loc=str('./log_scr/'),
        mfree=str('3G'),
        cluster_time=str('0:00:10:00'),
        mem_mb=str('3G'),
        disk_free=str('4G'),
        disk_mb=str('4G'),
        runtime=str('10m'),
        n_cores=str('1')
    script:
        'scr/supplemental/sfigure1_moi_coinfection.R'

## Rule sfig2 -- Generate sup. fig 2- MOI controls VPC
rule sfig2:
    input:
    output:
        out_file_name + "res/sup/sfig2/sfig2.pdf"
    params:
        date = init_date,
        fit_func = init_fit_func,
        mu = init_mu,
        rs_colors = rs_colors,
        axis_text_size = axis_text_size,
        lege_text_size = lege_text_size,
        sfig2_dim = sfig2_dim
    resources:
        log_loc=str('./log_scr/'),
        mfree=str('3G'),
        cluster_time=str('0:00:10:00'),
        mem_mb=str('3G'),
        disk_free=str('4G'),
        disk_mb=str('4G'),
        runtime=str('10m'),
        n_cores=str('10')
    script:
        'scr/supplemental/sfigure2_capsid_sub_comp_single.R'

## Rule sfig3 -- Generate sup. fig 3 - MOI controls VPC
rule sfig3:
    input:
    output:
        out_file_name + "res/sup/sfig3/sfig3.pdf"
    params:
        date = init_date,
        fit_func = init_fit_func,
        mu = init_mu,
        rs_colors = rs_colors,
        axis_text_size = axis_text_size,
        lege_text_size = lege_text_size,
        sfig3_dim = sfig3_dim
    resources:
        log_loc=str('./log_scr/'),
        mfree=str('3G'),
        cluster_time=str('0:00:10:00'),
        mem_mb=str('3G'),
        disk_free=str('4G'),
        disk_mb=str('4G'),
        runtime=str('10m'),
        n_cores=str('1')
    script:
        'scr/supplemental/sfigure3_capsid_sub_comp_multi.R'
        
## Rule sfig4 -- Generate sup. fig 4 - Clin Trial Trajectories
rule sfig4:
    input:
        out_file_name + "dat_gen/params/logistic_fitness_function.csv",
        out_file_name + "dat_gen/sims/sim_clinical_trial.csv"
    output:
        out_file_name + "res/sup/sfig4/sfig4.pdf"
    params:
        date = init_date,
        fit_func = init_fit_func,
        mu = init_mu,
        rs_colors = rs_colors,
        axis_text_size = axis_text_size,
        lege_text_size = lege_text_size,
        sfig4_dim = sfig4_dim
    resources:
        log_loc=str('./log_scr/'),
        mfree=str('3G'),
        cluster_time=str('0:00:10:00'),
        mem_mb=str('3G'),
        disk_free=str('4G'),
        disk_mb=str('4G'),
        runtime=str('10m'),
        n_cores=str('1')
    script:
        'scr/supplemental/sfigure4_clin_trial_trajectories.R'

#### Supplemental figure 5                                                  ####
## Rule sfig5A -- Generate sfig 5A - Alternative full dominant fitness curve####
rule sfig5A:
    input:
        out_file_name + "dat_gen/params/optim_params.csv",
        out_file_name + "dat_gen/params/logistic_fitness_function.csv"
    output:
        out_file_name + "res/sup/sfig5/A.rds"
    params:
        date = init_date,
        fit_func = init_fit_func,
        mu = init_mu,
        rs_colors = rs_colors,
        axis_text_size = axis_text_size,
        lege_text_size = lege_text_size,
        sfig5_dim = sfig5_dim
    resources:
        log_loc=str('./log_scr/'),
        mfree=str('3G'),
        cluster_time=str('0:00:10:00'),
        mem_mb=str('3G'),
        disk_free=str('4G'),
        disk_mb=str('4G'),
        runtime=str('10m'),
        n_cores=str('1')
    script:
        'scr/supplemental/sfigure5A_alt_fit_curves.R'

## Rule sfig5B -- Generate sfigure 5B - Resistance change one step more dom ####
rule sfig5B:
    input:
        out_file_name + "dat_gen/params/optim_params.csv",
        out_file_name + "dat_gen/params/logistic_fitness_function.csv"
    output:
        out_file_name + "res/sup/sfig5/B.rds"
    params:
        date = init_date,
        fit_func = init_fit_func,
        mu = init_mu,
        rs_colors = rs_colors,
        axis_text_size = axis_text_size,
        lege_text_size = lege_text_size,
        sfig5_dim = sfig5_dim
    resources:
        log_loc=str('./log_scr/'),
        mfree=str('3G'),
        cluster_time=str('0:00:10:00'),
        mem_mb=str('3G'),
        disk_free=str('4G'),
        disk_mb=str('4G'),
        runtime=str('10m'),
        n_cores=str('10')
    script:
        'scr/supplemental/sfigure5B_resistance_change.R'

## Rule sfig5C -- Generate sfigures 5C - Trajectory of more dominant drug   ####
rule sfig5C:
    input:
        out_file_name + "dat_gen/params/optim_params.csv",
        out_file_name + "dat_gen/params/logistic_fitness_function.csv"
    output:
        out_file_name + "res/sup/sfig5/C.rds"
    params:
        date = init_date,
        fit_func = init_fit_func,
        mu = init_mu,
        rs_colors = rs_colors,
        axis_text_size = axis_text_size,
        lege_text_size = lege_text_size,
        sfig5_dim = sfig5_dim
    resources:
        log_loc=str('./log_scr/'),
        mfree=str('3G'),
        cluster_time=str('0:00:10:00'),
        mem_mb=str('3G'),
        disk_free=str('4G'),
        disk_mb=str('4G'),
        runtime=str('10m'),
        n_cores=str('1')
    script:
        'scr/supplemental/sfigure5C_trajectory.R'

## Rule sfig5D -- Generate 5D - The range of less stringent drugs           ####
rule sfig5D:
    input:
        out_file_name + "dat_gen/params/optim_params.csv",
        out_file_name + "dat_gen/params/logistic_fitness_function.csv",
    output:
        out_file_name + "res/sup/sfig5/D.rds",
    params:
        date = init_date,
        fit_func = init_fit_func,
        mu = init_mu,
        rs_colors = rs_colors,
        axis_text_size = axis_text_size,
        lege_text_size = lege_text_size,
        sfig5_dim = sfig5_dim
    resources:
        log_loc=str('./log_scr/'),
        mfree=str('3G'),
        cluster_time=str('0:00:10:00'),
        mem_mb=str('3G'),
        disk_free=str('4G'),
        disk_mb=str('4G'),
        runtime=str('10m'),
        n_cores=str('1')
    script:
        'scr/supplemental/sfigure5D_range_of_dom.R'

## Rule sfig5EF -- Generate sfigures 5EF - All less stringent drugs trials ####
rule sfig5EF:
    input:
        out_file_name + "dat_gen/params/optim_params.csv",
        out_file_name + "dat_gen/params/logistic_fitness_function.csv",
        out_file_name + "dat_gen/sims/full_dominance_trials.csv"
    output:
        out_file_name + "res/sup/sfig5/E.rds",
        out_file_name + "res/sup/sfig5/F.rds"
    params:
        date = init_date,
        fit_func = init_fit_func,
        mu = init_mu,
        rs_colors = rs_colors,
        axis_text_size = axis_text_size,
        lege_text_size = lege_text_size,
        sfig5_dim = sfig5_dim
    resources:
        log_loc=str('./log_scr/'),
        mfree=str('3G'),
        cluster_time=str('0:00:10:00'),
        mem_mb=str('3G'),
        disk_free=str('4G'),
        disk_mb=str('4G'),
        runtime=str('10m'),
        n_cores=str('1')
    script:
        'scr/supplemental/sfigure5EF_more_dom_drug_vs_pocap.R'

## Pull all of the sfigure 5 plots together                                 ####
rule gen_sfig5:
    input:
        out_file_name + "res/sup/sfig5/A.rds",
        out_file_name + "res/sup/sfig5/B.rds",
        out_file_name + "res/sup/sfig5/C.rds",
        out_file_name + "res/sup/sfig5/D.rds",
        out_file_name + "res/sup/sfig5/E.rds",
        out_file_name + "res/sup/sfig5/F.rds"
    output:
        out_file_name + "res/sup/sfig5/sfig5.pdf"
    params:
        date = init_date,
        fit_func = init_fit_func,
        mu = init_mu,
        rs_colors = rs_colors,
        axis_text_size = axis_text_size,
        lege_text_size = lege_text_size,
        sfig5_dim = sfig5_dim
    resources:
        log_loc=str('./log_scr/'),
        mfree=str('3G'),
        cluster_time=str('0:00:10:00'),
        mem_mb=str('3G'),
        disk_free=str('4G'),
        disk_mb=str('4G'),
        runtime=str('10m'),
        n_cores=str('1')
    script:
        'scr/supplemental/gen_sfig_5.R'



## Rule sfig6a -- Generate Sfig 6A - The range of less stringent drugs
rule sfig6A:
    input:
        out_file_name + "dat_gen/params/optim_params.csv",
        out_file_name + "dat_gen/params/logistic_fitness_function.csv",
    output:
        out_file_name + "res/sup/sfig6/A.rds"
    params:
        date = init_date,
        fit_func = init_fit_func,
        mu = init_mu,
        n_trials = init_n_trials, ## This changes how many trials are run
        rs_colors = rs_colors,
        axis_text_size = axis_text_size,
        lege_text_size = lege_text_size,
        sfig6_dim = sfig6_dim
    resources:
        log_loc=str('./log_scr/'),
        mfree=str('3G'),
        cluster_time=str('0:00:10:00'),
        mem_mb=str('3G'),
        disk_free=str('4G'),
        disk_mb=str('4G'),
        runtime=str('10m'),
        n_cores=str('4')
    script:
        'scr/supplemental/sfigure6A_change_in_res_freq_by_weakness.R'

## Rule sfig6BC -- Generate figures S6BC - All less stringent drugs trials
rule sfig6BC:
    input:
        out_file_name + "dat_gen/params/optim_params.csv",
        out_file_name + "dat_gen/params/logistic_fitness_function.csv",
        out_file_name + "dat_gen/sims/all_full_stringency_trials.csv"
    output:
        out_file_name + "res/sup/sfig6/B.rds",
        out_file_name + "res/sup/sfig6/C.rds"
    params:
        date = init_date,
        fit_func = init_fit_func,
        mu = init_mu,
        rs_colors = rs_colors,
        axis_text_size = axis_text_size,
        lege_text_size = lege_text_size,
        sfig6_dim = sfig6_dim
    resources:
        log_loc=str('./log_scr/'),
        mfree=str('3G'),
        cluster_time=str('0:00:10:00'),
        mem_mb=str('3G'),
        disk_free=str('4G'),
        disk_mb=str('4G'),
        runtime=str('10m'),
        n_cores=str('1')
    script:
        'scr/supplemental/sfigure6BC_less_str_drug_vs_pocap.R'

## Pull all of the sfigure 6 plots together                                 ####
rule gen_sfig6:
    input:
        out_file_name + "res/sup/sfig6/A.rds",
        out_file_name + "res/sup/sfig6/B.rds",
        out_file_name + "res/sup/sfig6/C.rds",
    output:
        out_file_name + "res/sup/sfig6/sfig6.pdf"
    params:
        date = init_date,
        fit_func = init_fit_func,
        mu = init_mu,
        rs_colors = rs_colors,
        axis_text_size = axis_text_size,
        lege_text_size = lege_text_size,
        sfig6_dim = sfig6_dim
    resources:
        log_loc=str('./log_scr/'),
        mfree=str('3G'),
        cluster_time=str('0:00:10:00'),
        mem_mb=str('3G'),
        disk_free=str('4G'),
        disk_mb=str('4G'),
        runtime=str('10m'),
        n_cores=str('1')
    script:
        'scr/supplemental/gen_sfig_6.R'


#### Pull all figures into one document                                     ####
rule collate_all_sup:
    input:
        out_file_name + "res/sup/sfig1/sfig1.pdf",
        out_file_name + "res/sup/sfig2/sfig2.pdf",
        out_file_name + "res/sup/sfig3/sfig3.pdf",
        out_file_name + "res/sup/sfig4/sfig4.pdf",
        out_file_name + "res/sup/sfig5/sfig5.pdf",
        out_file_name + "res/sup/sfig6/sfig6.pdf"
    output:
        out_file_name + "res/all_sfigs.pdf"
    params:
        date = init_date,
        fit_func = init_fit_func,
        mu = init_mu
    resources:
        log_loc=str('./log_scr/'),
        mfree=str('3G'),
        cluster_time=str('0:00:10:00'),
        mem_mb=str('3G'),
        disk_free=str('4G'),
        disk_mb=str('4G'),
        runtime=str('10m'),
        n_cores=str('1')
    script:
        'scr/supplemental/collate_sup_plots.R'
