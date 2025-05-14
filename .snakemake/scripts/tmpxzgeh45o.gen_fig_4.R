
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
    input = list('runs/ddt_logistic_mu_2e-05_2025-03-24/res/fig4/A.rds', 'runs/ddt_logistic_mu_2e-05_2025-03-24/res/fig4/D.rds', 'runs/ddt_logistic_mu_2e-05_2025-03-24/res/fig4/E.rds', 'runs/ddt_logistic_mu_2e-05_2025-03-24/res/fig4/F.rds', 'runs/ddt_logistic_mu_2e-05_2025-03-24/res/fig4/G.rds', 'runs/ddt_logistic_mu_2e-05_2025-03-24/res/fig4/H.rds'),
    output = list('runs/ddt_logistic_mu_2e-05_2025-03-24/res/fig4/fig4.pdf'),
    params = list('2025-03-24', c('logistic'), c(2e-05), "date" = '2025-03-24', "fit_func" = c('logistic'), "mu" = c(2e-05)),
    wildcards = list(),
    threads = 1,
    log = list(),
    resources = list('mem_mb', 'disk_mb', 'tmpdir', 'log_loc', 'mfree', 'cluster_time', 'disk_free', 'runtime', 'n_cores', "mem_mb" = '3G', "disk_mb" = '4G', "tmpdir" = '/tmp/12724659.1.feder-short.q', "log_loc" = './log_scr/', "mfree" = '3G', "cluster_time" = '0:00:10:00', "disk_free" = '4G', "runtime" = 10, "n_cores" = '1'),
    config = list(),
    rule = 'gen_fig4',
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
############################ Generate figure 4 #################################

## This script will generate figure 4 from the disparate panels created in the
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
rds_files <- list.files(path = paste0(filepath, "res/fig4"), pattern = "\\.rds$", 
                        recursive = TRUE, full.names = TRUE)
plot_list <- lapply(rds_files, readRDS)


#### Plotting together                                                      ####
## Rebuilding figure 4A
fig4a <- plot_list[[1]]

## Creating the first row
top_row <- plot_grid(fig4a, NULL, NULL,
                    ncol = 3, axis = "l", align = "v")
top_row

## Rebuilding fig 4D
fig4d <- plot_list[[2]]

fig4d

## rebuilding fig 4E
fig4e <- plot_list[[3]]

fig4e

## Building the middle row
mid_row <- plot_grid(fig4d, fig4e,
                     ncol = 2, axis = "l", align = "v")

mid_row


## Building fig 4F
fig4f <- plot_list[[4]]

fig4f

## Building fig 4G
fig4g_top <- plot_list[[5]][[1]]

fig4g_bottom <- plot_list[[5]][[2]]

fig4g <- plot_grid(fig4g_top, NULL, fig4g_bottom,
                   rel_heights = c(0.3, -0.1, 1),
                   ncol = 1, axis = "l", align = "v")

fig4g

## Building fig 4H
fig4h_top <- plot_list[[6]][[1]]

fig4h_bottom <- plot_list[[6]][[2]]

fig4gh_xaxis <- plot_list[[6]][[3]]

fig4h <- plot_grid(fig4h_top, NULL, fig4h_bottom, NULL, fig4gh_xaxis,
                   rel_heights = c(0.3, -0.1, 1, -0.15, 0.3),
                   ncol = 1, axis = "l", align = "v")

fig4h

## Stacking G and H and adding the axis label

fig4gh <- plot_grid(fig4g, NULL, fig4h,
                   rel_heights = c(0.9, -0.1, 1),
                   ncol = 1, 
                   axis = "l", 
                   align = "v")

fig4gh


## Creating the bottom row
bot_row <- plot_grid(NULL, fig4f, NULL, fig4gh,
                     rel_widths = c(0.065, 1, -0.1, 1.155),
                     nrow = 1, 
                     axis = "l", 
                     align = "v")

bot_row


## All together
fig4 <- plot_grid(top_row, mid_row, bot_row,
                  ncol = 1, 
                  axis = "l", 
                  align = "v"
                  )


ggsave(paste0(filepath, "res/fig4/fig4.png"), fig4, h = 22, w = 18, unit = "cm", bg = "transparent", dpi = 300)
ggsave(paste0(filepath, "res/fig4/fig4.pdf"), fig4, h = 22, w = 18, unit = "cm", bg = "transparent", dpi = 300)
ggsave(paste0(filepath, "res/fig4/fig4.svg"), fig4, h = 22, w = 18, unit = "cm", bg = "transparent", dpi = 300, device = "svg")
# saveRDS(fig4, paste0(filepath, "res/fig4/fig4.rds"))
