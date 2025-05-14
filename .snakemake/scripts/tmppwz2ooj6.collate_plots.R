
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
    input = list('runs/ddt_logistic_mu_2e-05_2025-03-24/res/fig2/A.rds', 'runs/ddt_logistic_mu_2e-05_2025-03-24/res/fig2/B.rds', 'runs/ddt_logistic_mu_2e-05_2025-03-24/res/fig3/B.rds', 'runs/ddt_logistic_mu_2e-05_2025-03-24/res/fig3/C.rds', 'runs/ddt_logistic_mu_2e-05_2025-03-24/res/fig3/D.rds', 'runs/ddt_logistic_mu_2e-05_2025-03-24/res/fig4/A.rds', 'runs/ddt_logistic_mu_2e-05_2025-03-24/res/fig4/B.rds', 'runs/ddt_logistic_mu_2e-05_2025-03-24/res/fig4/E.rds', 'runs/ddt_logistic_mu_2e-05_2025-03-24/res/fig5/B.rds', 'runs/ddt_logistic_mu_2e-05_2025-03-24/res/fig5/C.rds', 'runs/ddt_logistic_mu_2e-05_2025-03-24/res/fig6/A.rds', 'runs/ddt_logistic_mu_2e-05_2025-03-24/res/fig6/B.rds', 'runs/ddt_logistic_mu_2e-05_2025-03-24/res/fig6/C.rds', 'runs/ddt_logistic_mu_2e-05_2025-03-24/res/fig6/D.rds', 'runs/ddt_logistic_mu_2e-05_2025-03-24/res/fig6/E.rds'),
    output = list('runs/ddt_logistic_mu_2e-05_2025-03-24/res/all_res.pdf'),
    params = list('2025-03-24', c('logistic'), c(2e-05), "date" = '2025-03-24', "fit_func" = c('logistic'), "mu" = c(2e-05)),
    wildcards = list(),
    threads = 1,
    log = list(),
    resources = list('mem_mb', 'disk_mb', 'tmpdir', 'log_loc', 'mfree', 'cluster_time', 'disk_free', 'runtime', 'n_cores', "mem_mb" = '3G', "disk_mb" = '4G', "tmpdir" = '/tmp/12680209.1.feder-short.q', "log_loc" = './log_scr/', "mfree" = '3G', "cluster_time" = '0:02:00:00', "disk_free" = '4G', "runtime" = 120, "n_cores" = '1'),
    config = list(),
    rule = 'collate_plots',
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
## A script to bring all of the different plots together into their separate
## figures and into one core document

#### Set Up                                                                 ####

## Libraries to load
require(dplyr)
require(ggplot2)

## Set up where to save
## Date of the run
init_date <- as.character(snakemake@params[["date"]])

## Fitness function
init_fit_func <- as.character(snakemake@params[["fit_func"]])

## Mutation rate
init_mu <- as.numeric(snakemake@params[["mu"]])

## Create the filepath where things will be saved
filepath <- paste0("runs/ddt_", init_fit_func, "_mu_", init_mu, "_", init_date, "/")

## Recursively list .rds files under /res/
rds_files <- list.files(path = paste0(filepath, "res/"), pattern = "\\.rds$", 
                        recursive = TRUE, full.names = TRUE)
plot_list <- lapply(rds_files, readRDS)

print(plot_list)

## Sort the panels by their respective figures
fig1_panels <- plot_list[grepl("Fig1", rds_files)]
fig2_panels <- plot_list[grepl("Fig2", rds_files)]
fig3_panels <- plot_list[grepl("Fig3", rds_files)]
fig4_panels <- plot_list[grepl("Fig4", rds_files)]
fig5_panels <- plot_list[grepl("Fig5", rds_files)]
fig6_panels <- plot_list[grepl("Fig6", rds_files)]

## Temp way of sorting all of the pannels together
require(patchwork)
figure1 <- wrap_plots(figure1_panels, ncol = 2)
figure2 <- wrap_plots(figure2_panels, ncol = 2)