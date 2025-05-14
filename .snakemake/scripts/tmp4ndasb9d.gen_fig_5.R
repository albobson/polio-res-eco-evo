
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
    input = list('runs/ddt_logistic_mu_2e-05_2025-03-24/res/fig5/A.rds', 'runs/ddt_logistic_mu_2e-05_2025-03-24/res/fig5/B.rds', 'runs/ddt_logistic_mu_2e-05_2025-03-24/res/fig5/C.rds', 'runs/ddt_logistic_mu_2e-05_2025-03-24/res/fig5/D.rds', 'runs/ddt_logistic_mu_2e-05_2025-03-24/res/fig5/E.rds', 'runs/ddt_logistic_mu_2e-05_2025-03-24/res/fig5/F.rds'),
    output = list('runs/ddt_logistic_mu_2e-05_2025-03-24/res/fig5/fig5.pdf'),
    params = list('2025-03-24', c('logistic'), c(2e-05), "date" = '2025-03-24', "fit_func" = c('logistic'), "mu" = c(2e-05)),
    wildcards = list(),
    threads = 1,
    log = list(),
    resources = list('mem_mb', 'disk_mb', 'tmpdir', 'log_loc', 'mfree', 'cluster_time', 'disk_free', 'runtime', 'n_cores', "mem_mb" = '3G', "disk_mb" = '4G', "tmpdir" = '/tmp/12724658.1.feder-short.q', "log_loc" = './log_scr/', "mfree" = '3G', "cluster_time" = '0:00:10:00', "disk_free" = '4G', "runtime" = 10, "n_cores" = '1'),
    config = list(),
    rule = 'gen_fig5',
    bench_iteration = as.numeric(NA),
    scriptdir = '/net/feder/vol1/home/alexrob/projects/dominant_drug_targets/pipeline/scr/full_figures',
    source = function(...){
        wd <- getwd()
        setwd(snakemake@scriptdir)
        source(...)
        setwd(wd)
    }
)


######## snakemake preamble end #########
############################ Generate figure 5 #################################

## This script will generate figure 5 from the disparate panels created in the
## specific panel scripts

#### Set up                                                                 ####
## Libraries 
require(dplyr)
require(ggplot2)
require(cowplot)
require(tidyr)

## Set up where to save
## Date of the run
init_date <- as.character(snakemake@params[["date"]])

## Fitness function
init_fit_func <- as.character(snakemake@params[["fit_func"]])

## Mutation rate
init_mu <- as.numeric(snakemake@params[["mu"]])

## Create the filepath where things will be saved
filepath <- paste0("runs/ddt_", init_fit_func, "_mu_", init_mu, "_", init_date, "/")


## Grab the .rds files
## Create the filepath where things will be saved
rds_files <- list.files(path = paste0(filepath, "res/fig5"), pattern = "\\.rds$", 
                        recursive = TRUE, full.names = TRUE)
plot_list <- lapply(rds_files, readRDS)


#### Plotting together                                                      ####
## Rebuilding figure 6A
fig5a <- plot_list[[1]]

## Rebuilding figure 6B
fig5b <- plot_list[[2]]

## Rebuilding figure 6C
fig5c <- plot_list[[3]]

## Creating the first row
top_row <- plot_grid(fig5a, fig5b, fig5c,
                     ncol = 3, axis = "l", align = "v")
top_row

## Rebuilding fig 6D
fig5d <- plot_list[[4]]

fig5d

## Rebuilding fig 6E
fig5e <- plot_list[[5]]

fig5e

## Rebuilding fig 6F
fig5f <- plot_list[[6]]

fig5f

## Building the middle row
mid_row <- plot_grid(fig5d, fig5e, fig5f,
                     ncol = 3, axis = "l", align = "v")

mid_row


## All together
fig5 <- plot_grid(top_row, mid_row,
                  nrow = 2, 
                  axis = "l", 
                  align = "v"
)


ggsave(paste0(filepath, "res/fig5/fig5.png"), fig5, h = 22*(2/3), w = 18, unit = "cm", bg = "transparent", dpi = 300)
ggsave(paste0(filepath, "res/fig5/fig5.pdf"), fig5, h = 22*(2/3), w = 18, unit = "cm", bg = "transparent", dpi = 300)
ggsave(paste0(filepath, "res/fig5/fig5.svg"), fig5, h = 22*(2/3), w = 18, unit = "cm", bg = "transparent", dpi = 300, device = "svg")
# saveRDS(fig5, paste0(filepath, "res/fig5/fig5.rds"))
