
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
    input = list('dat/collett_trial.csv', 'runs/ddt_logistic_mu_2e-05_2025-03-18/dat_gen/params/optim_params.csv', 'runs/ddt_logistic_mu_2e-05_2025-03-18/dat_gen/params/logistic_fitness_function.csv'),
    output = list('runs/ddt_logistic_mu_2e-05_2025-03-18/dat_gen/params/optim_clin_trial_params.rds', 'runs/ddt_logistic_mu_2e-05_2025-03-18/dat_gen/params/optim_clin_trial_params.csv'),
    params = list('2025-03-18', c('logistic'), c(2e-05), "date" = '2025-03-18', "fit_func" = c('logistic'), "mu" = c(2e-05)),
    wildcards = list(),
    threads = 1,
    log = list(),
    resources = list('mem_mb', 'disk_mb', 'tmpdir', 'log_loc', 'mfree', 'cluster_time', 'disk_free', 'runtime', 'n_cores', "mem_mb" = '3G', "disk_mb" = '4G', "tmpdir" = '/tmp/12560368.1.feder-short.q', "log_loc" = './log_scr/', "mfree" = '3G', "cluster_time" = '4:00:00:00', "disk_free" = '4G', "runtime" = 5760, "n_cores" = '20-25'),
    config = list(),
    rule = 'optim_clin_trial_params',
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
######################## Clinical trial Immune Params ##########################

## Reason:

## Here, we will be finding the optimum cell population and immune parameters
## for generating the placebo outcomes observed in the clinical trial.
## Originally, I was basing this off of the pocapavir (experimental) group,
## however, that distribution was influenced by the drug. We are interested in
## finding the optimum parameters to give us the outcomes observed without the
## drug.

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

## Read in Collett clinical trial data for comparison
coldat <- read.csv("dat/collett_trial.csv")

## Source functions to run simulations
source("scr/polv_DDT_functions.R")  ## Main functions 

## Number of cores to run the model
n_cores <- detectCores()-1

print(paste0("Number of cores detected: ", n_cores))


#### Running simulations                                                    ####

## For comparison, we will need to find the number of Collett clearers for each
## day from the placebo group
coldat_obs <- filter(coldat, grepl("Placebo", treatment, fixed = TRUE)) %>%
  mutate(time = clearange) %>%
  select(time) %>%
  group_by(time) %>%
  summarize(n_per_d = n())

## Number of simulations to run to generate probability distribution to compute
## log likelihood. Using 10 times the observed data
sims_to_run <- sum(coldat_obs$n_per_d) * 10

#### Function to optimize                                                   ####
func_optim_imm_param <- function(data, par, fit_func_in, ...) {
  ## Assigning parameters
  c_pop <- exp(par[1])      ## Params[1] is the log(c_pop)
  t_imm <- par[2]           ## Params[2] is the first day of the immune response
  imm_m <- par[3]           ## Params[3] is the mean of the immune function
  imm_sd <- par[4]          ## Params[4] is the sd of the immune function
  
  ## Reading in the other variables
  args <- list(...)
  ncores <- args$ncores
  sims_to_run <- args$sims_to_run
  optim_params <- args$optim_params
  
  ## changing the name of data to a different object
  coldat_obs <- data
  
  ## Generate a dataframe to store data and a df to bind with
  df_test <- NULL
  df_it <- NULL
  
  ## Creating a cluster to make this run faster
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  
  ## Placebo group
  trial_raw <- foreach(k = 1:sims_to_run, 
                         .combine = "rbind", 
                         .export = 'stoch_polv',
                     .packages = c('dplyr')) %dopar% {
                       stoch_polv(n = 300, 
                                  moi_wt_start = 1, 
                                  moi_mut_start = 0, 
                                  t_pocap = 0, 
                                  id = k,
                                  imm_delay = t_imm,
                                  c_pop = c_pop,
                                  imm_m = imm_m,
                                  imm_sd = imm_sd, 
                                  fit_func_in = fit_func_in,
                                  v_prog = optim_params$optim_v_prog,
                                  p2pfu = optim_params$optim_p2pfu,
                                  seed_in = k*2
                       ) %>% 
                         mutate(time = time/3) ## Changing from replications to days
                     }
  
  stopCluster(cl)
  stopImplicitCluster()
  
  
  ## Clean the data to find the dates of each clearance
  simsclean <- trial_raw %>% group_by(id) %>% 
    filter(type == "resistant") %>%
    select(id, time) %>% 
    filter(time == max(time)) %>%
    mutate(time = ceiling(time))
  
  ## Fix the fact that not every day was sampled
  simsclean$time <- ifelse(simsclean$time %in%  16:17, 18,
                           ifelse(simsclean$time %in% 19:21, 22,
                                  ifelse(simsclean$time %in% 23:28, 29,
                                         ifelse(simsclean$time %in% 30:42, 43,
                                                ifelse(simsclean$time > 43, 100, ## If the date observed was greater than 43, put 100. Must reject this
                                                       simsclean$time)))))
  
  ## If there are any times == 100, make the value very high
  if(any(simsclean$time == 100)) {
    
    tot_logL <- 1000
    
    ## Saving the parameters and optimization variable to visualize later
    ## First updating counter
    cc <<- cc+1
    
    ## Now recoding the values
    vals[[cc]] <<- c(cc, par, tot_logL)
  
  print(paste0("iteration: ", cc, ", Val: ", tot_logL))
  ## I need to flip the sign here, since optim tries to minimize values
  return(tot_logL)
  
  } else { ## If there are not any times == 100, record as usual
    ## Now calculate the PDF of this distribution
    trial_pdf <- simsclean %>%
      group_by(time) %>%
      summarize(freq = n()/sims_to_run)
    
    ## Now, we need to fill in the times that were not observed with 0s
    ## Note: the only problem with this is that it potentially ruins the likelihood function. An observation where there was 0 probability will result in a 0 likelihood.
    dates_sampled <- data.frame(time = c(seq(0, 15), 18, 22, 29, 43))
    
    ## We will combine this with the trial PDF
    trial_pdf <- merge(trial_pdf, dates_sampled, all = TRUE)
    
    ## Add a very small value to anything that was not observed. This way we do
    ## not deal with negative infinity due to log(0) error. I'm going to reweight
    ## with a value two orders of magnitude smaller than the number of simulations
    ## that were ran.
    trial_pdf[is.na(trial_pdf)] <- 1/(sims_to_run*100)
    
    ## Reweight the probabilities now that we have added a small number
    trial_pdf$freq <- trial_pdf$freq/sum(trial_pdf$freq)
    
    ## Now we need to calculate the log(likelihood). We will do this by summing
    ## the log probability of each observation. For ease, I will combine datasets.
    tot_df <- merge(coldat_obs, trial_pdf, all  = TRUE)
    
    ## Record where 0 observations were seen
    tot_df[is.na(tot_df)] <- 0
    
    tot_df$log_freq <- log(tot_df$freq)
    
    tot_df$logL <- tot_df$n_per_d*tot_df$log_freq
    
    tot_logL <- sum(tot_df$logL)
    
    ## Saving the parameters and optimization variable to visualize later
    ## First updating counter
    cc <<- cc+1
    
    ## Now recoding the values
    vals[[cc]] <<- c(cc, par, tot_logL)
    
    print(paste0("iteration: ", cc, ", Val: ", -tot_logL))
    ## I need to flip the sign here, since optim tries to minimize values
    return(-tot_logL)
  }
}

#### Run optimization                                                       ####
## The initial values to start
params_in <- c(
  log(8*10^3),      ## Params[1] is the log(c_pop)
  9,                ## Params[2] is the first time step of the immune response
  -1.6,             ## Params[3] is the mean of the immune function
  0.5               ## Params[4] is the sd of the immune function
) 

cc <- 0         ## This will count the cycles of optim
vals <- list()  ## And this will store the values of optim


## We need to generate a new fitness function, since this is comparing to the
## placebo group. I'm setting the fitness function == 1.
placebo_fit_func <- rep(1, 61)

## Running the function
system.time(
  optim_clin_trial_params <- optim(data = coldat_obs, 
                                  par = params_in, 
                                  fn = func_optim_imm_param, 
                                  fit_func_in = placebo_fit_func,
                                  ncores = n_cores,
                                  sims_to_run = sims_to_run,
                                  optim_params = optim_params,
                                  method = "Nelder-Mead"
  )
)

#### Output                                                                 ####
## Save the exact optim output object
suppressWarnings(dir.create(paste0(filepath, "dat_gen")))
suppressWarnings(dir.create(paste0(filepath, "dat_gen/params")))
saveRDS(optim_clin_trial_params, file = paste0(filepath, "dat_gen/params/optim_clin_trial_params.rds"))

## Save the values that the optim function iterated over
saveRDS(vals, file = paste0(filepath, "dat_gen/params/vals_optim_clin_trial_vals.rds"))

## Saving these parameters in a CSV forlater use
c_pop <- exp(optim_clin_trial_params$par[1])      ## Params[1] is the log(c_pop)
t_imm <- optim_clin_trial_params$par[2]           ## Params[2] is the first day of the immune response
imm_m <- optim_clin_trial_params$par[3]           ## Params[3] is the mean of the immune function
imm_sd <- optim_clin_trial_params$par[4]          ## Params[4] is the sd of the immune function

optim_df <- data.frame(fit_type = "Placebo",        ## Type of fitness function used
                       optim_c_pop = exp(optim_clin_trial_params$par[1]), ## c_pop parameter value
                       optim_t_imm = optim_clin_trial_params$par[2],  ## t_imm parameter value
                       optim_imm_m = optim_clin_trial_params$par[3], ## imm_m parameter value
                       optim_imm_sd = optim_clin_trial_params$par[4],  ## imm_sd parameter value
                       minimum_lLik = -optim_clin_trial_params$value)   ## Best fit difference


write.csv(optim_df, file = paste0(filepath, "dat_gen/params/optim_clin_trial_params.csv"), row.names = F)
