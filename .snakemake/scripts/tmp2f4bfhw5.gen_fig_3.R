
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
    input = list('runs/ddt_logistic_mu_2e-05_2025-03-24/res/fig3/B.rds', 'runs/ddt_logistic_mu_2e-05_2025-03-24/res/fig3/C.rds', 'runs/ddt_logistic_mu_2e-05_2025-03-24/res/fig3/D.rds'),
    output = list('runs/ddt_logistic_mu_2e-05_2025-03-24/res/fig3/fig3.pdf'),
    params = list('2025-03-24', c('logistic'), c(2e-05), "date" = '2025-03-24', "fit_func" = c('logistic'), "mu" = c(2e-05)),
    wildcards = list(),
    threads = 1,
    log = list(),
    resources = list('mem_mb', 'disk_mb', 'tmpdir', 'log_loc', 'mfree', 'cluster_time', 'disk_free', 'runtime', 'n_cores', "mem_mb" = '3G', "disk_mb" = '4G', "tmpdir" = '/tmp/12724654.1.feder-short.q', "log_loc" = './log_scr/', "mfree" = '3G', "cluster_time" = '0:00:10:00', "disk_free" = '4G', "runtime" = 10, "n_cores" = '1'),
    config = list(),
    rule = 'gen_fig3',
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
############################ Generate figure 3 #################################

## This script will generate figure 3 from the disparate panels created in the
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
rds_files <- list.files(path = paste0(filepath, "res/fig3"), pattern = "\\.rds$", 
                        recursive = TRUE, full.names = TRUE)
plot_list <- lapply(rds_files, readRDS)



#### Plotting together                                                      ####
## Rebuilding figure 3B
fit_func_p <- plot_list[1][[1]]

## Rebuilding figure 3C
tan_pocap <- plot_list[2][[1]][[1]]

tan_p <- plot_list[2][[1]][[2]]

tan_mois <- plot_list[2][[1]][[3]]

full_tan_p <- plot_grid(tan_pocap, NULL, tan_p, NULL, tan_mois, ## Adding NULL to move plots closer together
                            rel_heights = c(0.3, -0.1, 1.2, -0.05, 0.3),
                            ncol = 1, axis = "l", align = "v")

full_tan_p

## Redbuilding figure 3D
my_pocap <- plot_list[3][[1]][[1]]

my_p <- plot_list[3][[1]][[2]]

my_mois <- plot_list[3][[1]][[3]]

full_my_p <- plot_grid(my_pocap, NULL, my_p, NULL, my_mois, ## Adding NULL to move plots closer together
                        rel_heights = c(0.3, -0.1, 1.2, -0.05, 0.3),
                        ncol = 1, axis = "l", align = "v")

full_my_p

#### All together

# ## All together
full_fig3 <- plot_grid(
  NULL, fit_func_p,
  full_tan_p, full_my_p,
  rel_heights = c(0.72, 1), ## Need to change the relative height of the top row to match C and D's bounds
  rel_widths = c(1),
  # rel_widths = c(),
  ncol = 2,
  axis = "l",
  align = "v"
)


full_fig3

## All together
# full_fig3 <- plot_grid(
#   NULL, NULL, NULL, NULL, NULL, ## First panel
#   NULL, NULL, fit_func_p, NULL, NULL, ## Second panel
#   tan_pocap, NULL, tan_p, NULL, tan_mois, ## Third pannel
#   my_pocap, NULL, my_p, NULL, my_mois, ## Fourth pannel
#   rel_heights = c(
#     0.3, -0.1, 1.2, -0.05, 0.3,
#     0.3, -0.1, 1.2, -0.05, 0.3,
#     0.3, -0.1, 1.2, -0.05, 0.3,
#     0.3, -0.1, 1.2, -0.05, 0.3
#     ),
#   ncol = 4,
#   axis = "l", 
#   align = "v"
# )
# 
# 
# full_fig3

save_plot(plot = full_fig3, filename = paste0(filepath, "res/fig3/fig3.png"), base_height = 7, base_width = 7)#+3.5*0.05)
save_plot(plot = full_fig3, filename = paste0(filepath, "res/fig3/fig3.pdf"), base_height = 7, base_width = 7) #+3.5*0.05)
ggsave(paste0(filepath, "res/fig3/fig3.svg"), full_fig3, h = 22, w = 18, unit = "cm", bg = "transparent", dpi = 300, device = "svg")