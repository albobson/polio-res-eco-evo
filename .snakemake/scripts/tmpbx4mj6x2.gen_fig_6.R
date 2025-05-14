
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
    input = list('runs/ddt_logistic_mu_2e-05_2025-03-24/res/fig6/A.rds', 'runs/ddt_logistic_mu_2e-05_2025-03-24/res/fig6/B.rds', 'runs/ddt_logistic_mu_2e-05_2025-03-24/res/fig6/C.rds', 'runs/ddt_logistic_mu_2e-05_2025-03-24/res/fig6/D.rds', 'runs/ddt_logistic_mu_2e-05_2025-03-24/res/fig6/E.rds', 'runs/ddt_logistic_mu_2e-05_2025-03-24/res/fig6/F.rds'),
    output = list('runs/ddt_logistic_mu_2e-05_2025-03-24/res/fig6/fig6.pdf'),
    params = list('2025-03-24', c('logistic'), c(2e-05), "date" = '2025-03-24', "fit_func" = c('logistic'), "mu" = c(2e-05)),
    wildcards = list(),
    threads = 1,
    log = list(),
    resources = list('mem_mb', 'disk_mb', 'tmpdir', 'log_loc', 'mfree', 'cluster_time', 'disk_free', 'runtime', 'n_cores', "mem_mb" = '3G', "disk_mb" = '4G', "tmpdir" = '/tmp/12724610.1.feder-short.q', "log_loc" = './log_scr/', "mfree" = '3G', "cluster_time" = '0:00:10:00', "disk_free" = '4G', "runtime" = 10, "n_cores" = '1'),
    config = list(),
    rule = 'gen_fig6',
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
############################ Generate figure 6 #################################

## This script will generate figure 6 from the disparate panels created in the
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
rds_files <- list.files(path = paste0(filepath, "res/fig6"), pattern = "\\.rds$", 
                        recursive = TRUE, full.names = TRUE)
plot_list <- lapply(rds_files, readRDS)


#### Plotting together                                                      ####
## Rebuilding figure 6A
fig6a <- plot_list[[1]]

## Rebuilding figure 6B
fig6b <- plot_list[[2]]

## Rebuilding figure 6C
fig6c <- plot_list[[3]]

## Creating the first row
top_row <- plot_grid(fig6a, fig6b, fig6c,
                     ncol = 3, axis = "l", align = "v")
top_row

## Rebuilding fig 6D
fig6d <- plot_list[[4]]

fig6d

## Rebuilding fig 6E
fig6e <- plot_list[[5]]

fig6e

## Rebuilding fig 6F
fig6f <- plot_list[[6]]

fig6f

## Building the middle row
mid_row <- plot_grid(fig6d, fig6e, fig6f,
                     ncol = 3, axis = "l", align = "v")

mid_row


## All together
fig6 <- plot_grid(top_row, mid_row,
                  nrow = 2, 
                  axis = "l", 
                  align = "v"
)


ggsave(paste0(filepath, "res/fig6/fig6.png"), fig4, h = 22*(2/3), w = 18, unit = "cm", bg = "transparent", dpi = 300)
ggsave(paste0(filepath, "res/fig6/fig6.pdf"), fig4, h = 22*(2/3), w = 18, unit = "cm", bg = "transparent", dpi = 300)
ggsave(paste0(filepath, "res/fig6/fig6.svg"), fig4, h = 22*(2/3), w = 18, unit = "cm", bg = "transparent", dpi = 300, device = "svg")
# saveRDS(fig6, paste0(filepath, "res/fig6/fig6.rds"))
