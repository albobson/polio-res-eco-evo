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
n_trials <- 2000

## Find how many individuals per trial 
n_indiv <- as.numeric(snakemake@params[["n_indiv"]])

## How many participants to be run at 24 hours (as opposed to 72)?
## 70/93 were at 24 in the original trial.
per_24 <- as.numeric(snakemake@params[["perc_24"]])


#### Set up less stringent drugs                                            ####
## Set the range of dominance
## Define the range for which we will scale the function
min_scale <- 0
max_scale <- 60

scale_range <- seq(min_scale, max_scale, length.out = n_trials)

## Storage df
scale_store <- NULL

## generate functions across scale range
for(i in 1:length(scale_range)) {
  temp_scale <- NULL
  temp_df <- NULL
  
  temp_scale <- fit_func_logistic(k = 100, 
                                  min_val = min(fit_func$prob_surv),
                                  subunits = 60, mid = i)
  
  ## Set our minimum and maximum values to the values that we want
  temp_scale[1] <- min(fit_func$prob_surv)
  temp_scale[length(temp_scale)] <- max(fit_func$prob_surv)
  
  temp_df <- data.frame(subunits = 0:60, prob_surv = temp_scale, fit_type = scale_range[i])
  
  scale_store <- rbind(scale_store, temp_df)
}

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
                      .packages = c('dplyr')) %dopar% {
                        ## Create null DFs
                        trial_24_h <- NULL
                        
                        ## Find our fitness function for this run
                        new_fit_func <- NULL
                        new_fit_func <- scale_store %>%
                          filter(fit_type == scale_range[i])
                        
                        # print(paste0("Trial # ", i))
                        
                        ## Generate all of the seeds to be used for these sims
                        set.seed(12345)
                        
                        ## Row: trial number, Column: participant ID
                        seed_vec <- round(runif(n_trials * n_indiv, 1, 1000000))
                        
                        
                        ## 24 hours
                        trial_24_h <- stoch_polv(n = 300, 
                                                  moi_wt_start = 1, 
                                                  moi_mut_start = 0, 
                                                  t_pocap = 3, 
                                                  id = i,
                                                  imm_delay = clin_trial_params$optim_t_imm,
                                                  c_pop = clin_trial_params$optim_c_pop, 
                                                  imm_m = clin_trial_params$optim_imm_m, 
                                                  imm_sd = clin_trial_params$optim_imm_sd,
                                                  fit_func_in = new_fit_func$prob_surv,
                                                  v_prog = optim_params$optim_v_prog,
                                                  p2pfu = optim_params$optim_p2pfu,
                                                  seed_in = seed_vec[i]
                                                ) %>% 
                                                  mutate(dom_val = scale_range[i])
                        
                        trial <- trial_24_h
                        
                        return(trial)
                      }

## Important to stop the cluster
stopImplicitCluster()
stopCluster(cl)

print("sims done")

## Save the file as a CSV
write.csv(mixed_sims, paste0(filepath, "dat_gen/sims/full_dominance_trials.csv"), row.names = FALSE)

print("csv saved")