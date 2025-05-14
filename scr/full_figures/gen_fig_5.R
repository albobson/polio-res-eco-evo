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

## Color scheme for resistant and susceptible
color_in <- snakemake@params[["rs_colors"]]

## Read in axis text size
axis_text_size <- as.numeric(snakemake@params[["axis_text_size"]])

## Legend text size
legend_text_size <- as.numeric(snakemake@params[["lege_text_size"]])

## This figure's size
width <-  as.numeric(snakemake@params[["fig5_dim"]])[1]
height <- as.numeric(snakemake@params[["fig5_dim"]])[2]

## Create the filepath where things will be saved
filepath <- paste0("runs/ddt_", init_fit_func, "_mu_", init_mu, "_", init_date, "/")


## Grab the .rds files
## Create the filepath where things will be saved
rds_files <- list.files(path = paste0(filepath, "res/fig5"), pattern = "\\.rds$", 
                        recursive = TRUE, full.names = TRUE)
plot_list <- lapply(rds_files, readRDS)


#### Plotting together                                                      ####
## Rebuilding figure 5A
fig5a <- plot_list[[1]]

## Rebuilding figure 5B
fig5b <- plot_list[[2]]


## Creating the first row
top_row <- plot_grid(fig5a, fig5b,
                     ncol = 2, axis = "l", align = "v")


## Rebuilding figure 5C
fig5c <- plot_list[[3]]

## Rebuilding fig 5D
fig5d <- plot_list[[4]]


## Building the middle row
mid_row <- plot_grid(fig5c, fig5d,
                     ncol = 2, axis = "l", align = "v")


## All together
fig5 <- plot_grid(top_row, mid_row,
                  nrow = 2, 
                  axis = "l", 
                  align = "v"
)


save_plot(paste0(filepath, "res/fig5/fig5.png"), fig5,  base_height = height, base_width = width-width*0.3, unit = "cm", bg = "transparent", dpi = 1200)
save_plot(paste0(filepath, "res/fig5/fig5.pdf"), fig5,  base_height = height, base_width = width-width*0.3, unit = "cm", bg = "transparent", dpi = 1200)
# ggsave(paste0(filepath, "res/fig5/fig5.svg"), fig5,  base_height = height, base_width = width, unit = "cm", bg = "transparent", dpi = 1200, device = "svg")
# saveRDS(fig5, paste0(filepath, "res/fig5/fig5.rds"))
