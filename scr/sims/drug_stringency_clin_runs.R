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


## Make sure do parallel is loaded
library(doParallel)


#### Run simulations                                                        ####
## Find how many variations we want to run 
n_trials <- as.numeric(snakemake@params[["n_trials"]])

## Find how many individuals per trial 
n_indiv <- as.numeric(snakemake@params[["n_indiv"]])

## How many participants to be run at 24 hours (as opposed to 72)?
## 70/93 were at 24 in the original trial.
per_24 <- as.numeric(snakemake@params[["perc_24"]])

## Find the minimum p(WT survival for each trial), evenly spaced across log scale
# min_surv <- 10^seq(from = log10(fit_func$prob_surv[1]), to = 0, length.out = n_trials)

## Specifically finding 1, 10, 100 and 1000 fold weaker
scaled_vector <- function(og_probs, scaler) {
  ## New vector for storing
  new_vec <- NULL
  
  new_vec <- ((1- scaler * og_probs[1])/(1-og_probs[1])) * (og_probs - og_probs[1]) + scaler * og_probs[1]
  
  return(new_vec)
}

## Define the range for which we will scale the function
# min_scale <- 1
# max_scale <- 1000
# 
# scale_range <- 10^log10(seq(min_scale, max_scale, length.out = n_trials))

min_scale <- 0
max_scale <- n_trials-1

scale_range <- 10^c(min_scale:max_scale)

## Create the cluster
print(paste0("Cores detected: ", detectCores()))

## Register cluster
cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl = cl)

## Make sure that it has access to the doParallel package
# clusterEvalQ(cl, .libPaths("/net/feder/vol1/home/alexrob/R/x86_64-pc-linux-gnu-library/4.3"))


#### Run the sims                                                           ####
## First loop will adjust the fitness function. Second loop will run the sims
mixed_sims <- foreach(i = 1:length(scale_range), 
                      .combine = 'rbind', 
                      .packages = c('dplyr')) %:%
                          foreach(k = 1:n_indiv, 
                                  .combine = 'rbind', 
                                  .packages = c('dplyr')) %dopar% {
                        ## Create null DFs
                        trial_24_h <- trial_72_h <- NULL
                        
                        ## Create scaling function
                        # scaled_vector <- function(og_probs, new_offset) {
                        #   p <- new_offset + ((og_probs - og_probs[1]) * (1 - new_offset) / (max(og_probs) - og_probs[1]))
                        #   p <- pmin(p, 1)  # Ensure values do not exceed 1
                        #   return(p)
                        # }
                        # 
                        fit_func$prob_surv_w <- scaled_vector(og_probs = fit_func$prob_surv, scaler = scale_range[i])
                        
                        ## Generate all of the seeds to be used for these sims
                        set.seed(12345)
                        
                        ## Row: trial number, Column: participant ID
                        seed_vec <- round(runif(n_trials * n_indiv, 1, 1000000))
                        
                        ## 24 hours
                        trial_24_h <- stoch_polv(n = 300, 
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
                                                 seed_in = seed_vec[k]
                        ) %>% 
                          mutate(wt_rep_abil = fit_func$prob_surv_w[1],
                                 fold_strin = scale_range[i])
                        
                        trial <- trial_24_h
                        
                        return(trial)
}

## Important to stop the cluster
stopImplicitCluster()
stopCluster(cl)

print("sims done")

## Save the file as a CSV
write.csv(mixed_sims, paste0(filepath, "dat_gen/sims/full_stringency_trials.csv"), row.names = FALSE)

print("csv saved")
