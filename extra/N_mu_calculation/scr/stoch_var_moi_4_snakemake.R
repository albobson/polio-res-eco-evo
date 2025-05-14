## 230830

#### Introduction ####
## Running version 2 of the stochastic model

#### Set up ####
## Tell R to display all scientific-nodated numbers as their full number
options(scipen=999)

## Get wildcards 
n_u <- as.numeric(snakemake@wildcards[["nu"]])
iter <- as.numeric(snakemake@wildcards[["run"]])
gens <- as.numeric(snakemake@params[["generations"]])


## Create the filepath where things will be saved
filepath <- paste0("../../pipeline/runs/ddt_logistic_mu_2e-5_2025-05-09/")

## geom_text() is in terms of mm, where element_text() is in terms of pt. I want
## my geom_text and axis text to be the same. Here, I am defining a conversion
## which will keep the text sizes the same between text types.
geom_text_conv = 0.3528


## Read in fitness function and parameters
fit_func <- read.csv(paste0("../../pipeline/runs/ddt_logistic_mu_2e-05_2025-05-09/dat_gen/params/logistic_fitness_function.csv"))
optim_params <- read.csv(paste0("../../pipeline/runs/ddt_logistic_mu_2e-05_2025-05-09/dat_gen/params/optim_params.csv"))


## Source functions to run simulations
source("../../pipeline/scr/polv_DDT_functions.R")  ## Main functions 
source("../../pipeline/scr/optimizations/model_optimization_funcs.R")  ## Optimization functions

## From the n_u, find the initial parameters to use
moi_res <- 0
moi_sus <- 10 ## MOI at 10 so that nearly every cell is infected
mut_rate <- 2e-5
cells <- n_u/mut_rate


#### Run sims ####
sim_stochv2 <- stoch_polv(n = gens,
                          imm_delay = gens + 1,
                          t_pocap = 0,
                           c_pop = cells,
                           moi_mut_start = moi_res, 
                           moi_wt_start = moi_sus,
                           bg_mutat = 1/mut_rate,
                           max_vpc = optim_params$optim_max_vpc,
                           v_prog = optim_params$optim_v_prog,
                           p2pfu = optim_params$optim_p2pfu,
                           fit_func_in = fit_func$prob_surv,
                          id = iter
)

df_length <- length(sim_stochv2$time)

n_u_data <- data.frame(nu = rep(n_u, df_length))

## Bind the n_u data to the dataframe
sim_stochv2 <- cbind(sim_stochv2, n_u_data)

#### Wrap-up ####
## Create the directory for exporting data
dir.create(file.path("sim_data/"), showWarnings = FALSE)

## Create the names of the files
file_name <- paste0("sim_data/n-u-", format(n_u, nsmall = 8), "_it-", iter, ".csv")

## Write the .csv's
write.csv(sim_stochv2, file_name, row.names = FALSE)