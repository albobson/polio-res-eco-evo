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

## Color scheme for resistant and susceptible
color_in <- snakemake@params[["rs_colors"]]

## Read in axis text size
axis_text_size <- as.numeric(snakemake@params[["axis_text_size"]])

## Legend text size
legend_text_size <- as.numeric(snakemake@params[["lege_text_size"]])

## This figure's size
width <-  as.numeric(snakemake@params[["fig2_dim"]])[1]
height <- as.numeric(snakemake@params[["fig2_dim"]])[2]

## Create the filepath where things will be saved
filepath <- paste0("runs/ddt_", init_fit_func, "_mu_", init_mu, "_", init_date, "/")

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
                            rel_heights = c(0.3, -0.15, 1.2, -0.1, 0.3),
                            ncol = 1, axis = "l", align = "v")


## Rebuilding figure 2B
tend_plot <- plot_list[[2]][[1]]

tend_plot_cols <- plot_list[[2]][[2]]


tend_p <- plot_grid(NULL, NULL, tend_plot, NULL, tend_plot_cols, ## Adding NULL to move plots closer together
                            rel_heights = c(0.3, -0.15, 1.2, -0.1, 0.3),
                            ncol = 1, axis = "l", align = "v")

## All together
full_fig2 <- plot_grid(
  full_tan_trunc, NULL, tend_p, 
  rel_widths = c(
    1, 0.05, 1
  ),
  ncol = 3,
  axis = "l", 
  align = "v"
)


save_plot(plot = full_fig2, filename = paste0(filepath, "res/fig2/fig2.png"), base_height = height, base_width = width, unit = "cm", dpi = 1200)
save_plot(plot = full_fig2, filename = paste0(filepath, "res/fig2/fig2.pdf"), base_height = height, base_width = width, unit = "cm", dpi = 1200)

