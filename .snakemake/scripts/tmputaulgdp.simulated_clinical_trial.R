
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
    input = list('dat/collett_trial.csv', 'runs/ddt_logistic_mu_2e-05_2025-03-24/dat_gen/params/optim_params.csv', 'runs/ddt_logistic_mu_2e-05_2025-03-24/dat_gen/params/logistic_fitness_function.csv', 'runs/ddt_logistic_mu_2e-05_2025-03-24/dat_gen/params/optim_clin_trial_params.csv'),
    output = list('runs/ddt_logistic_mu_2e-05_2025-03-24/dat_gen/sims/sim_clinical_trial.csv'),
    params = list('2025-03-24', c('logistic'), c(2e-05), "date" = '2025-03-24', "fit_func" = c('logistic'), "mu" = c(2e-05)),
    wildcards = list(),
    threads = 1,
    log = list(),
    resources = list('mem_mb', 'disk_mb', 'tmpdir', 'log_loc', 'mfree', 'cluster_time', 'disk_free', 'runtime', 'n_cores', "mem_mb" = '3G', "disk_mb" = '4G', "tmpdir" = '/tmp/12723306.1.feder-short.q', "log_loc" = './log_scr/', "mfree" = '3G', "cluster_time" = '0:02:00:00', "disk_free" = '4G', "runtime" = 120, "n_cores" = '10'),
    config = list(),
    rule = 'sim_clin_trial',
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
######################### Clinical trial simulation ############################

## Reason:

## This script will run a clinical trial in accordance with the clinical trial
## run in Collett et al. (2017)

#### Set Up                                                                 ####
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


#### Running simulations                                                    ####
cl <- makeCluster(detectCores() - 2)
registerDoParallel(cl)

set.seed(1249128)
trial_72_h <- foreach(k = 1:70, .combine = "rbind", 
                      .packages = c('dplyr')) %dopar% {
                        stoch_polv(n = 300, 
                                   moi_wt_start = 1, 
                                   moi_mut_start = 0, 
                                   t_pocap = 9, 
                                   id = k,
                                   imm_delay = clin_trial_params$optim_t_imm,
                                   c_pop = clin_trial_params$optim_c_pop, 
                                   imm_m = clin_trial_params$optim_imm_m, 
                                   imm_sd = clin_trial_params$optim_imm_sd,
                                   fit_func_in = fit_func$prob_surv,
                                   v_prog = optim_params$optim_v_prog,
                                   p2pfu = optim_params$optim_p2pfu,
                                   seed_in = k*2
                        ) %>%
                          mutate(treatment = paste0("Simulated Pocapavir")) %>%
                          mutate(time = time/3) ## Changing from replications to days
                      }


set.seed(12491284)
trial_24_h <- foreach(k = 71:93, .combine = "rbind") %dopar% {
  stoch_polv(n = 300, 
             moi_wt_start = 1, 
             moi_mut_start = 0, 
             t_pocap = 3, 
             id = k,
             imm_delay = clin_trial_params$optim_t_imm,
             c_pop = clin_trial_params$optim_c_pop, #usually 8*10^3
             imm_m = clin_trial_params$optim_imm_m, # Usually -1.6
             imm_sd = clin_trial_params$optim_imm_sd, # Usually 0.5
             fit_func_in = fit_func$prob_surv,
             v_prog = optim_params$optim_v_prog,
             p2pfu = optim_params$optim_p2pfu,
             seed_in = k*2
  ) %>% 
    mutate(treatment = paste0("Simulated Pocapavir")) %>%
    mutate(time = time/3) ## Changing from replications to days
}

stopCluster(cl)
stopImplicitCluster()

trial <- rbind(trial_72_h, trial_24_h)

#### Saving                                                                 ####

write.csv(trial, file = paste0(filepath, "dat_gen/sims/sim_clinical_trial.csv"))