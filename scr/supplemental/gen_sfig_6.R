############################ Generate sfigure 5 ################################

## This script will generate sfigure 5 from the disparate panels created in the
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
width <-  as.numeric(snakemake@params[["sfig6_dim"]])[1]
height <- as.numeric(snakemake@params[["sfig6_dim"]])[2]

## Create the filepath where things will be saved
filepath <- paste0("runs/ddt_", init_fit_func, "_mu_", init_mu, "_", init_date, "/")


## Grab the .rds files
## Create the filepath where things will be saved
rds_files <- list.files(path = paste0(filepath, "res/sup/sfig6"), pattern = "\\.rds$", 
                        recursive = TRUE, full.names = TRUE)
plot_list <- lapply(rds_files, readRDS)


#### Plotting together                                                      ####
## Rebuilding figure 6A
sfig6a <- plot_list[[1]]

## Rebuilding figure 6B
sfig6b <- plot_list[[2]]

## Rebuilding figure 6C
sfig6c <- plot_list[[3]]

## Creating the first row
sfig6 <- plot_grid(sfig6a, sfig6b, sfig6c,
                     ncol = 3, axis = "l", align = "v")
# 
# ## Rebuilding fig 6D
# sfig5d <- plot_list[[4]]
# 
# ## Rebuilding fig 6E
# sfig5e <- plot_list[[5]]
# 
# ## Rebuilding fig 6F
# sfig5f <- plot_list[[6]]
# 
# ## Building the middle row
# mid_row <- plot_grid(sfig5d, sfig5e, sfig5f,
#                      ncol = 3, axis = "l", align = "v")
# 
# ## All together
# sfig6 <- plot_grid(top_row, mid_row,
#                   nrow = 2, 
#                   axis = "l", 
#                   align = "v"
# )


save_plot(paste0(filepath, "res/sup/sfig6/sfig6.png"), sfig6, base_height = height, base_width = width, unit = "cm", bg = "transparent", dpi = 1200)
save_plot(paste0(filepath, "res/sup/sfig6/sfig6.pdf"), sfig6, base_height = height, base_width = width, unit = "cm", bg = "transparent", dpi = 1200)
# save_plot(paste0(filepath, "res/sfig5/sfig5.svg"), sfig5, base_height = height, base_width = width, unit = "cm", bg = "transparent", dpi = 1200, device = "svg")
# saveRDS(sfig5, paste0(filepath, "res/sfig5/sfig5.rds"))