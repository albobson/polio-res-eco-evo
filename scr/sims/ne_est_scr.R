############################### Sup. Fig 7 #####################################

## Reason:

## Generate supplemental figure 7 (really 6, now), data for plotting extinction
## probability by Ne.

#### Set Up                                                                 ####

## Set up where to save
## Date of the run
init_date <- as.character(snakemake@params[["date"]])

## Fitness function
init_fit_func <- as.character(snakemake@params[["fit_func"]])

## Mutation rate
init_mu <- as.numeric(snakemake@params[["mu"]])

## Color scheme for resistant and susceptible
color_in <- snakemake@params[["rs_colors"]]

## Read in axis text size
axis_text_size <- as.numeric(snakemake@params[["axis_text_size"]])

## Legend text size
legend_text_size <- as.numeric(snakemake@params[["lege_text_size"]])

## This figure's size
width <-  as.numeric(snakemake@params[["sfig7_dim"]])[1]
height <- as.numeric(snakemake@params[["sfig7_dim"]])[2]

## Create the filepath where things will be saved
filepath <- paste0("runs/ddt_", init_fit_func, "_mu_", init_mu, "_", init_date, "/")

## Read in fitness function and parameters
fit_func <- read.csv(paste0(filepath, "dat_gen/params/logistic_fitness_function.csv"))
optim_params <- read.csv(paste0(filepath, "dat_gen/params/optim_params.csv"))

## Source functions to run simulations
source("scr/polv_DDT_functions.R")  ## Main functions 

## geom_text() is in terms of mm, where element_text() is in terms of pt. I want
## my geom_text and axis text to be the same. Here, I am defining a conversion
## which will keep the text sizes the same between text types.
geom_text_conv = 0.3528



########################## Run the simulations #################################
## Read in the parameters to use
start_in <-  as.numeric(snakemake@params[["point_start"]])
end_in <-  as.numeric(snakemake@params[["point_end"]])
n_points <-  as.numeric(snakemake@params[["n_points"]])

## Find the n_u values to run
n_u <- 10^seq()







## Note: saving to "filepath "dat_gen/sims/ne_sims.csv""