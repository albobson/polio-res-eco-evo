
######## snakemake preamble start (automatically inserted, do not edit) ########
library(methods)
Snakemake <- setClass(
    "Snakemake",
    slots = c(
        input = "list",
        output = "list",
        params = "list",
        wildcards = "list",
        threads = "numeric",
        log = "list",
        resources = "list",
        config = "list",
        rule = "character",
        bench_iteration = "numeric",
        scriptdir = "character",
        source = "function"
    )
)
snakemake <- Snakemake(
    input = list('runs/ddt_logistic_mu_2e-05_2025-03-24/dat_gen/params/optim_params.csv', 'runs/ddt_logistic_mu_2e-05_2025-03-24/dat_gen/params/logistic_fitness_function.csv', 'runs/ddt_logistic_mu_2e-05_2025-03-24/dat_gen/params/optim_clin_trial_params.csv'),
    output = list('runs/ddt_logistic_mu_2e-05_2025-03-24/dat_gen/sims/full_stringency_trials.csv'),
    params = list('2025-03-24', c('logistic'), c(2e-05), c(5000), c(1), c(1), "date" = '2025-03-24', "fit_func" = c('logistic'), "mu" = c(2e-05), "n_trials" = c(5000), "n_indiv" = c(1), "perc_24" = c(1)),
    wildcards = list(),
    threads = 1,
    log = list(),
    resources = list('mem_mb', 'disk_mb', 'tmpdir', 'log_loc', 'mfree', 'cluster_time', 'disk_free', 'runtime', 'n_cores', "mem_mb" = '3G', "disk_mb" = '4G', "tmpdir" = '/tmp/12722435.1.feder-short.q', "log_loc" = './log_scr/', "mfree" = '3G', "cluster_time" = '3:00:00:00', "disk_free" = '4G', "runtime" = 4320, "n_cores" = '20'),
    config = list(),
    rule = 'stringency_trials',
    bench_iteration = as.numeric(NA),
    scriptdir = '/net/feder/vol1/home/alexrob/projects/dominant_drug_targets/pipeline/scr/sims',
    source = function(...){
        wd <- getwd()
        setwd(snakemake@scriptdir)
        source(...)
        setwd(wd)
    }
)


######## snakemake preamble end #########
########## Clinical trial simulations with variable drug stringencies ##########

## Reason: 

## This script will run multiple clinical trials across the range of potential
## fitness functions. The data will be used to generate figure 6C-H

#### Set up ####

## Set up where to save
## Date of the run
init_date <- as.character(snakemake@params[["date"]])

## Fitness function
init_fit_func <- as.character(snakemake@params[["fit_func"]])

## Mutation rate
init_mu <- as.numeric(snakemake@params[["mu"]])

## Create the filepath where things will be saved
filepath <- paste0("runs/ddt_", init_fit_func, "_mu_", init_mu, "_", init_date, "/")

## Read in fitness function and parameters
fit_func <- read.csv(paste0(filepath, "dat_gen/params/logistic_fitness_function.csv"))
optim_params <- read.csv(paste0(filepath, "dat_gen/params/optim_params.csv"))
clin_trial_params <- read.csv(paste0(filepath, "dat_gen/params/optim_clin_trial_params.csv"))

## Source functions to run simulations
source("scr/polv_DDT_functions.R")  ## Main functions 

## Tell R to not use scientific notation:
options(scipen = 999)

## Try to create the directory (just in case it's not already created)
suppressWarnings(dir.create(paste0(filepath, "dat_gen/sims")))



#### Run simulations                                                        ####
## Find how many variations we want to run 
n_trials <- as.numeric(snakemake@params[["n_trials"]])

## Find how many individuals per trial 
n_indiv <- as.numeric(snakemake@params[["n_indiv"]])

## How many participants to be run at 24 hours (as opposed to 72)?
## 70/93 were at 24 in the original trial.
per_24 <- as.numeric(snakemake@params[["perc_24"]])

## Find the minimum p(WT survival for each trial), evenly spaced across log scale
min_surv <- 10^seq(from = log10(fit_func$prob_surv[1]), to = 0, length.out = n_trials)

## Create the cluster
print(paste0("Cores detected: ", detectCores()))

## Register cluster
cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl = cl)

## Make sure that it has access to the doParallel package
# clusterEvalQ(cl, .libPaths("/net/feder/vol1/home/alexrob/R/x86_64-pc-linux-gnu-library/4.3"))


#### Run the sims                                                           ####
## First loop will adjust the fitness function. Second loop will run the sims
mixed_sims <- foreach(i = 1:n_trials, 
                      .combine = 'rbind', 
                      .packages = c('dplyr')) %do% {
                        ## Create null DFs
                        trial_24_h <- trial_72_h <- NULL
                        
                        ## Create scaling function
                        scaled_vector <- function(og_probs, new_offset) {
                          p <- new_offset + ((og_probs - og_probs[1]) * (1 - new_offset) / (max(og_probs) - og_probs[1]))
                          p <- pmin(p, 1)  # Ensure values do not exceed 1
                          return(p)
                        }
                        
                        fit_func$prob_surv_w <- scaled_vector(og_probs = fit_func$prob_surv, new_offset = min_surv[i])
                        
                        print(paste0("Trial # ", i))
                        
                        ## Generate all of the seeds to be used for these sims
                        set.seed(12345)
                        
                        ## Row: trial number, Column: participant ID
                        seed_vec <- round(runif(n_trials * n_indiv, 1, 1000000))
                        
                        seed_matx <- matrix(data = seed_vec, nrow = n_trials, ncol = n_indiv)
                        
                        
                        ## 24 hours
                        trial_24_h <- foreach(k = 1:round((n_indiv*(per_24))), ## Percentage that were administered at 24 hours
                                .combine = 'rbind', 
                                .packages = c('dplyr')) %dopar% {
                                  
                                  (stoch_polv(n = 300, 
                                              moi_wt_start = 1, 
                                              moi_mut_start = 0, 
                                              t_pocap = 3, 
                                              id = k,
                                              imm_delay = clin_trial_params$optim_t_imm,
                                              c_pop = clin_trial_params$optim_c_pop, 
                                              imm_m = clin_trial_params$optim_imm_m, 
                                              imm_sd = clin_trial_params$optim_imm_sd,
                                              fit_func_in = fit_func$prob_surv_w,
                                              v_prog = optim_params$optim_v_prog,
                                              p2pfu = optim_params$optim_p2pfu,
                                              seed_in = seed_matx[i, k]
                                  ) %>% 
                                    mutate(wt_rep_abil = fit_func$prob_surv_w[1]))
                                }
                        
                        
                        ## check to see that there are some 72 hour administrations to be run
                        if(per_24 != 1) {
                            ## 72 hours
                            trial_72_h <- foreach(k = round((n_indiv*(per_24)+1)):n_indiv,  ## Percentage that were administered at 72 hours
                                                  .combine = 'rbind', 
                                                  .packages = c('dplyr')) %dopar% {
                                                    
                                                    (stoch_polv(n = 300, 
                                                                moi_wt_start = 1, 
                                                                moi_mut_start = 0, 
                                                                t_pocap = 9, 
                                                                id = k,
                                                                imm_delay = clin_trial_params$optim_t_imm,
                                                                c_pop = clin_trial_params$optim_c_pop, 
                                                                imm_m = clin_trial_params$optim_imm_m, 
                                                                imm_sd = clin_trial_params$optim_imm_sd,
                                                                fit_func_in = fit_func$prob_surv_w,
                                                                v_prog = optim_params$optim_v_prog,
                                                                p2pfu = optim_params$optim_p2pfu,
                                                                seed_in = seed_matx[i, k]
                                                    ) %>% 
                                                      mutate(wt_rep_abil = fit_func$prob_surv_w[1]))
                                                  }
                        }
                        
                        
                        trial <- rbind(trial_24_h, trial_72_h)
                        
                        return(trial)
}

## Important to stop the cluster
stopImplicitCluster()
stopCluster(cl)

print("sims done")

## Save the file as a CSV
write.csv(mixed_sims, paste0(filepath, "dat_gen/sims/full_stringency_trials.csv"), row.names = FALSE)

print("csv saved")