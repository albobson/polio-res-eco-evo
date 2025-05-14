
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
    input = list('dat/tanner_pure_culture.csv', 'dat/tanner_mix_culture.csv'),
    output = list('runs/ddt_logistic_mu_2e-05_2025-03-24/dat_gen/params/logistic_fitness_function.csv', 'runs/ddt_logistic_mu_2e-05_2025-03-24/dat_gen/params/tan_fit_df.csv', 'runs/ddt_logistic_mu_2e-05_2025-03-24/dat_gen/params/logistic_fit_summary.csv'),
    params = list('2025-03-24', c('logistic'), c(2e-05), "date" = '2025-03-24', "fit_func" = c('logistic'), "mu" = c(2e-05)),
    wildcards = list(),
    threads = 1,
    log = list(),
    resources = list('mem_mb', 'disk_mb', 'disk_mib', 'tmpdir', 'log_loc', 'mfree', 'cluster_time', 'disk_free', 'runtime', 'n_cores', "mem_mb" = '3G', "disk_mb" = 1000, "disk_mib" = 954, "tmpdir" = '/tmp/12666695.1.feder-short.q', "log_loc" = './log_scr/', "mfree" = '3G', "cluster_time" = '0:00:10:00', "disk_free" = '2G', "runtime" = 60, "n_cores" = '1'),
    config = list(),
    rule = 'fit_func',
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
########################### Generate Fit Func ##################################

## Reason:

## This figure pannel shows what the fitness function looks like for our
## analysis. For a comparison of other fitness functions, please look at
## supplemental figure 1.

#### Set Up                                                                 ####

## libraries
require(dplyr)
require(ggplot2)
require(minpack.lm)

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

## Read in Tanner cell culture data
tan_pure <- read.csv("dat/tanner_pure_culture.csv") ## Mono culture data
tan_mix <- read.csv("dat/tanner_mix_culture.csv")   ## Co-culture data

#### Survival probabilities ####
## For this analysis, we are using the VP3-A24V mutation in Mahoney strain 
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

## Using co-culture data, define the fitness curve
fit_df_wo0 <- data.frame(subunits = 60*mah_a24v$moi_mut/(mah_a24v$moi_mut+mah_a24v$moi_wt),
                         prob_surv = mah_a24v$pfu/mah_a24v$pfu[1])

## Incorporate the 0 resistant subunit survival
fit_df <- rbind(fit_df_wo0, data.frame(subunits = 0, prob_surv = mah_a24v_yint))

## Save this
write.csv(x = fit_df, file = paste0(filepath, "dat_gen/params/tan_fit_df.csv"))


###### Generate logistic fitness function                                   ####
x <- fit_df$subunits
y <- fit_df$prob_surv

## Define the endpoints using the desired x values.
## (Here we assume that there is a data point at x = 0 and one at x = 60.)
y0   <- fit_df$prob_surv[fit_df$subunits == 0]
y60  <- fit_df$prob_surv[fit_df$subunits == 60]

## Define the scaled logistic function
scaled_logistic <- function(x, k, x0) {
  L <- function(x) 1 / (1 + exp(-k * (x - x0)))
  ## Compute scaled value ensuring f(0) = y0 and f(60) = y60
  y0 + (y60 - y0) * ((L(x) - L(0)) / (L(60) - L(0)))
}

## Provide starting estimates. 
start_k <- (max(y) - min(y)) / (max(x) - min(x)) # Alternative slope estimate
start_x0 <- 30 # Using 40, as this was estimated in Tanner (2014)

## Fit the data
fit <- nlsLM(prob_surv ~ scaled_logistic(subunits, k, x0),
             data = fit_df,
             start = list(k = start_k, x0 = start_x0),
             lower = c(k = 0.0001, x0 = 0),
             upper = c(k = 100, x0 = 60)
)

## Summarize the model
summary(fit)

## Extract estimated parameters
summ_fit <- summary(fit)
new_k  <- summ_fit$coefficients[1, 1]
new_x0 <- summ_fit$coefficients[2, 1]

## Create a new data frame for plotting the fitted curve over a dense grid
fix_logistic <- data.frame(subunits = seq(0, 60, length.out = 61))
fix_logistic$prob_surv <- scaled_logistic(x = fix_logistic$subunits, k = new_k, x0 = new_x0)
fix_logistic$fit_type <- "logistic"

#### Save the fitness function for later use                                ####
suppressWarnings(dir.create(paste0(filepath, "dat_gen/")))
suppressWarnings(dir.create(paste0(filepath, "dat_gen/params")))
write.csv(x = fix_logistic, file = paste0(filepath, "dat_gen/params/logistic_fitness_function.csv"), row.names = F)

## Save the fit summary parameters
write.csv(x = summ_fit, file = paste0(filepath, "dat_gen/params/logistic_fit_summary.csv"), row.names = F)