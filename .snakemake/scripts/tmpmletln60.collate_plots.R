
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
    input = list('runs/ddt_logistic_mu_2e-05_2025-03-24/res/fig2/fig2.pdf', 'runs/ddt_logistic_mu_2e-05_2025-03-24/res/fig3/fig3.pdf', 'runs/ddt_logistic_mu_2e-05_2025-03-24/res/fig4/fig4.pdf', 'runs/ddt_logistic_mu_2e-05_2025-03-24/res/fig5/fig5.pdf', 'runs/ddt_logistic_mu_2e-05_2025-03-24/res/fig6/fig6.pdf'),
    output = list('runs/ddt_logistic_mu_2e-05_2025-03-24/res/all_figs.pdf'),
    params = list('2025-03-24', c('logistic'), c(2e-05), "date" = '2025-03-24', "fit_func" = c('logistic'), "mu" = c(2e-05)),
    wildcards = list(),
    threads = 1,
    log = list(),
    resources = list('mem_mb', 'disk_mb', 'tmpdir', 'log_loc', 'mfree', 'cluster_time', 'disk_free', 'runtime', 'n_cores', "mem_mb" = '3G', "disk_mb" = '4G', "tmpdir" = '/tmp/12724805.1.feder-short.q', "log_loc" = './log_scr/', "mfree" = '3G', "cluster_time" = '0:00:10:00', "disk_free" = '4G', "runtime" = 10, "n_cores" = '1'),
    config = list(),
    rule = 'collate_all',
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
## A script to bring all of the different plots together into their separate
## figures and into one core document

#### Set Up                                                                 ####

## Libraries to load
require(dplyr)
require(ggplot2)
require(qpdf)

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
pdf_files <- list.files(
  path = paste0(filepath, "res/"),
  pattern = "^fig[0-9]\\.pdf$",
  full.names = TRUE,
  recursive = TRUE
)

pdf_combine(input = pdf_files, output = paste0(filepath, "res/all_figs.pdf"))
