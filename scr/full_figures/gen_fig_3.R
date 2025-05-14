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

## Color scheme for resistant and susceptible
color_in <- snakemake@params[["rs_colors"]]

## Read in axis text size
axis_text_size <- as.numeric(snakemake@params[["axis_text_size"]])

## Legend text size
legend_text_size <- as.numeric(snakemake@params[["lege_text_size"]])

## This figure's size
width <-  as.numeric(snakemake@params[["fig3_dim"]])[1]
height <- as.numeric(snakemake@params[["fig3_dim"]])[2]

## Create the filepath where things will be saved
filepath <- paste0("runs/ddt_", init_fit_func, "_mu_", init_mu, "_", init_date, "/")


## Grab the .rds files
## Create the filepath where things will be saved
rds_files <- list.files(path = paste0(filepath, "res/fig3"), pattern = "\\.rds$", 
                        recursive = TRUE, full.names = TRUE)
plot_list <- lapply(rds_files, readRDS)



#### Individual plots                                                       ####

## Rebuilding figure 3B
fit_func_p <- plot_list[1][[1]]

## Rebuilding figure 3C
tan_pocap <- plot_list[2][[1]][[1]]

tan_p <- plot_list[2][[1]][[2]]

tan_mois <- plot_list[2][[1]][[3]]

rel_cd_heights <- c(0.3, -0.15, 1.2, -0.1, 0.4)

full_tan_p <- plot_grid(tan_pocap, NULL, tan_p, NULL, tan_mois, ## Adding NULL to move plots closer together
                            rel_heights = rel_cd_heights,
                            ncol = 1, axis = "l", align = "v")

## Redbuilding figure 3D
my_pocap <- plot_list[3][[1]][[1]]

my_p <- plot_list[3][[1]][[2]]

my_mois <- plot_list[3][[1]][[3]]

full_my_p <- plot_grid(my_pocap, NULL, my_p, NULL, my_mois, ## Adding NULL to move plots closer together
                        rel_heights = rel_cd_heights,
                        ncol = 1, axis = "l", align = "v")



#### All together                                                           ####
top_row <- plot_grid(
  NULL, fit_func_p,
  ncol = 2
)

bottom_row <- plot_grid(
  full_tan_p, NULL, full_my_p,
  rel_widths = c(1, -0.035, 1),
  ncol = 3,
  axis = "l",
  align = "v"
)


# ## All together
full_fig3 <- plot_grid(
  NULL,
  top_row,
  NULL,
  bottom_row,
  rel_heights = c(0.05, 0.78, -0.04, 1), ## Need to change the relative height of the top row to match C and D's bounds
  nrow = 4,
  axis = "l",
  align = "v"
)


full_fig3


save_plot(plot = full_fig3, filename = paste0(filepath, "res/fig3/fig3.png"), base_height = height, base_width = width, units = "cm", bg = "transparent", dpi = 1200)
save_plot(plot = full_fig3, filename = paste0(filepath, "res/fig3/fig3.pdf"), base_height = height, base_width = width, units = "cm", bg = "transparent", dpi = 1200)
ggsave(paste0(filepath, "res/fig3/fig3.svg"), full_fig3, h = height, w = width, units = "cm", bg = "transparent", dpi = 1200, device = "svg")

