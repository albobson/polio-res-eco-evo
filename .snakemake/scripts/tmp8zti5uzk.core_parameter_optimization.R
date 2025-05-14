
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
    input = list('runs/ddt_logistic_mu_2e-05_2025-03-24/dat_gen/params/logistic_fitness_function.csv'),
    output = list('runs/ddt_logistic_mu_2e-05_2025-03-24/dat_gen/params/optim_params.csv'),
    params = list('2025-03-24', c('logistic'), c(2e-05), "date" = '2025-03-24', "fit_func" = c('logistic'), "mu" = c(2e-05)),
    wildcards = list(),
    threads = 1,
    log = list(),
    resources = list('mem_mb', 'disk_mb', 'disk_mib', 'tmpdir', 'log_loc', 'mfree', 'cluster_time', 'disk_free', 'runtime', 'n_cores', "mem_mb" = '3G', "disk_mb" = 1000, "disk_mib" = 954, "tmpdir" = '/tmp/12666758.1.feder-short.q', "log_loc" = './log_scr/', "mfree" = '3G', "cluster_time" = '0:02:00:00', "disk_free" = '4G', "runtime" = 120, "n_cores" = '7'),
    config = list(),
    rule = 'opt_model',
    bench_iteration = as.numeric(NA),
    scriptdir = '/net/feder/vol1/home/alexrob/projects/dominant_drug_targets/pipeline/scr',
    source = function(...){
        wd <- getwd()
        setwd(snakemake@scriptdir)
        source(...)
        setwd(wd)
    }
)


######## snakemake preamble end #########
####################### Function Parameter Optimization ########################

## Reason:

## Based on the fitness curve, we need to optimize our unknown parameters. For
## this model, our unknown paramters are burst size and particle to plaque
## forming unit ratio. We need to have generated the fitness curve already in
## order to run this script.

#### Set Up                                                                 ####

## Source simulation functions
source("scr/polv_DDT_functions.R")  ## Main functions 
source("scr/model_optimization_funcs.R")  ## Optimization functions

## Read in Tanner cell culture data for comparison
tan_pure <- read.csv("dat/tanner_pure_culture.csv") ## Mono culture data
tan_mix <- read.csv("dat/tanner_mix_culture.csv")   ## Co-culture data

## Set up where to save
## Date of the run
init_date <- as.character(snakemake@params[["date"]])

## Fitness function
init_fit_func <- as.character(snakemake@params[["fit_func"]])

## Mutation rate
init_mu <- as.numeric(snakemake@params[["mu"]])

## Create the filepath where things will be saved
filepath <- paste0("runs/ddt_", init_fit_func, "_mu_", init_mu, "_", init_date, "/")

## Create this filepath
suppressWarnings(dir.create(paste0(filepath)))

## Read in fitness function dataframe
fit_func_df <- read.csv(paste0(filepath, "dat_gen/params/logistic_fitness_function.csv"))

## Isolating just the probabilities from this dataframe
fit_func <- as.vector(fit_func_df$prob_surv)

#### Clean up Tanner data for comparison                                    ####
## For this analysis, we will optimize on VP3-A24V
## Calculate the probability of survival for capsids with 0 resistant subunits
mah_a24v_sus_no_p <- tan_pure %>%
  filter(mutation == "VP3-A24V", 
         strain == "Mahoney", 
         sample == "sus", 
         pocap == "no")

mah_a24v_sus_w_p <- tan_pure %>%
  filter(mutation == "VP3-A24V", 
         strain == "Mahoney", 
         sample == "sus", 
         pocap == "yes")

## Calculate probability of survival
mah_a24v_yint <- mah_a24v_sus_w_p$pfu/mah_a24v_sus_no_p$pfu

## Filter the co-culture data
mah_a24v <- tan_mix %>%
  filter(mutation == "VP3-A24V", 
         strain == "Mahoney", 
         pocap == "yes")

#### Run optimization                                                       ####
## The initial values to start
params_in <- c(
  log(1000),    ## Burst size -- Logging this to explore log range
  log(100),     ## p2pfu -- Logging this to explore log range
  500           ## max_vpc -- Maximum viruses per cell 
) 

cc <- 0         ## This will count the cycles of optim
vals <- list()  ## And this will store the values of optim

## Number of cores to use
num_cores <-  detectCores() - 1 # Save one core

## Register parallel backend
cl <- makeCluster(num_cores)
registerDoParallel(cl)

## Running the function
system.time(
  optim_params <- optim(data = mah_a24v, 
                          par = params_in, 
                          fn = func_to_optim, ## This comes from the model_optimization.R script
                          fit_func_in = fit_func,
                          method = "Nelder-Mead"
  )
)

stopCluster(cl)
stopImplicitCluster()

#### . . Output ####
## Print the new parameters
# print(paste0(
#   "optim v_prog: ", exp(optim_params$par[1]), 
#   ", optim p2pfu: ", exp(optim_params$par[2]))
# )

## Saving these parameters in a new file for later use
optim_df <- data.frame(fit_type = fit_func_df$fit_type[1],      ## Type of fitness function used
                       optim_v_prog = exp(optim_params$par[1]), ## v_prog parameter value
                       optim_p2pfu = exp(optim_params$par[2]),  ## p2pfu parameter value
                       optim_max_vpc = optim_params$par[3],     ## max vpc value
                       minimum_fit_diff = optim_params$value)   ## Best fit difference

write.csv(optim_df, file = paste0(filepath, "dat_gen/params/optim_params.csv"), row.names = F)
