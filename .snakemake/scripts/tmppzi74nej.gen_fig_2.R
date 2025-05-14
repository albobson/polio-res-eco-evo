
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
    input = list('runs/ddt_logistic_mu_2e-05_2025-03-24/res/fig2/A.rds', 'runs/ddt_logistic_mu_2e-05_2025-03-24/res/fig2/B.rds'),
    output = list('runs/ddt_logistic_mu_2e-05_2025-03-24/res/fig2/fig2.pdf'),
    params = list('2025-03-24', c('logistic'), c(2e-05), "date" = '2025-03-24', "fit_func" = c('logistic'), "mu" = c(2e-05)),
    wildcards = list(),
    threads = 1,
    log = list(),
    resources = list('mem_mb', 'disk_mb', 'tmpdir', 'log_loc', 'mfree', 'cluster_time', 'disk_free', 'runtime', 'n_cores', "mem_mb" = '3G', "disk_mb" = '4G', "tmpdir" = '/tmp/12724653.1.feder-short.q', "log_loc" = './log_scr/', "mfree" = '3G', "cluster_time" = '0:00:10:00', "disk_free" = '4G', "runtime" = 10, "n_cores" = '1'),
    config = list(),
    rule = 'gen_fig2',
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
############################ Generate figure 2 #################################

## This script will generate figure 2 from the disparate panels created in the
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

## Making the color map for plotting
color_in <- c("#c94d4d", "#688cc8")

## Grab the .rds files
## Create the filepath where things will be saved
rds_files <- list.files(path = paste0(filepath, "res/fig2"), pattern = "\\.rds$", 
                        recursive = TRUE, full.names = TRUE)
plot_list <- lapply(rds_files, readRDS)


#### Plotting together                                                      ####
## Rebuilding figure 2A
tan_pocap <- plot_list[1][[1]][[1]]

tan_trunc_p <- plot_list[1][[1]][[2]]

tan_mois <- plot_list[1][[1]][[3]]


full_tan_trunc <- plot_grid(tan_pocap, NULL, tan_trunc_p, NULL, tan_mois, ## Adding NULL to move plots closer together
                            rel_heights = c(0.3, -0.1, 1.2, -0.05, 0.3),
                            ncol = 1, axis = "l", align = "v")

full_tan_trunc


## Rebuilding figure 2B
tend_plot <- plot_list[[2]][[1]]

tend_plot_cols <- plot_list[[2]][[2]]


tend_p <- plot_grid(NULL, NULL, tend_plot, NULL, tend_plot_cols, ## Adding NULL to move plots closer together
                            rel_heights = c(0.3, -0.1, 1.2, -0.05, 0.3),
                            ncol = 1, axis = "l", align = "v")

tend_p


## All together
full_fig2 <- plot_grid(
  full_tan_trunc, NULL, tend_p, 
  rel_widths = c(
    1, 0.05, 0.9
  ),
  ncol = 3,
  axis = "l", 
  align = "v"
)


save_plot(plot = full_fig2, filename = paste0(filepath, "res/fig2/fig2.png"), base_height = 3.5, base_width = 7)#+3.5*0.05)
save_plot(plot = full_fig2, filename = paste0(filepath, "res/fig2/fig2.pdf"), base_height = 3.5, base_width = 7) #+3.5*0.05)
