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

## Color scheme for resistant and susceptible
color_in <- snakemake@params[["rs_colors"]]

## Read in axis text size
axis_text_size <- as.numeric(snakemake@params[["axis_text_size"]])

## Legend text size
legend_text_size <- as.numeric(snakemake@params[["lege_text_size"]])

## This figure's size
width <-  as.numeric(snakemake@params[["fig4_dim"]])[1]
height <- as.numeric(snakemake@params[["fig4_dim"]])[2]

## Create the filepath where things will be saved
filepath <- paste0("runs/ddt_", init_fit_func, "_mu_", init_mu, "_", init_date, "/")


## Grab the .rds files
## Create the filepath where things will be saved
rds_files <- list.files(path = paste0(filepath, "res/fig4"), pattern = "\\.rds$", 
                        recursive = TRUE, full.names = TRUE)
plot_list <- lapply(rds_files, readRDS)


#### Plotting together                                                      ####
## Rebuilding figure 4A
# fig4a <- plot_list[[1]]
# 
# ## Rebuilding fig 4D
# fig4d <- plot_list[[2]]
# 
# ## Creating the first row
# top_row <- plot_grid(fig4a, NULL, fig4d,
#                      ncol = 3, axis = "l", align = "v")
# 
# ## rebuilding fig 4E
# fig4e <- plot_list[[3]]
# 
# ## Building fig 4F
# fig4f <- plot_list[[4]]
# 
# ## Building fig 4G
# fig4g_top <- plot_list[[5]][[1]]
# 
# fig4g_bottom <- plot_list[[5]][[2]]
# 
# fig4g <- plot_grid(fig4g_top, NULL, fig4g_bottom,
#                    rel_heights = c(0.3, -0.1, 1),
#                    ncol = 1, axis = "l", align = "v")
# 
# ## Building fig 4H
# fig4h_top <- plot_list[[6]][[1]]
# 
# fig4h_bottom <- plot_list[[6]][[2]]
# 
# fig4gh_xaxis <- plot_list[[6]][[3]]
# 
# fig4h <- plot_grid(fig4h_top, NULL, fig4h_bottom, NULL, fig4gh_xaxis,
#                    rel_heights = c(0.3, -0.1, 1, -0.15, 0.3),
#                    ncol = 1, axis = "l", align = "v")
# 
# ## Stacking G and H and adding the axis label
# fig4gh <- plot_grid(fig4g, NULL, fig4h,
#                     rel_heights = c(0.9, -0.1, 1),
#                     ncol = 1,
#                     axis = "l",
#                     align = "v")
# 
# 
# 
# ## Building the bottom row
# bot_row <- plot_grid(fig4e, fig4f, NULL, fig4gh,
#                         ncol = 4,
#                      rel_widths = c(1, 1, -0.1, 1),
#                         axis = "l",
#                         align = "v")
# 
# ## All together
# fig4 <- plot_grid(top_row, bot_row,
#                   nrow = 2,
#                   axis = "l",
#                   align = "v"
#                   )
# 
# 
# save_plot(plot = fig4, filename = paste0(filepath, "res/fig4/fig4.png"), base_height = height, base_width = width, units = "cm", bg = "transparent", dpi = 1200)
# save_plot(plot = fig4, filename = paste0(filepath, "res/fig4/fig4.pdf"), base_height = height, base_width = width, units = "cm", bg = "transparent", dpi = 1200)
# ggsave(paste0(filepath, "res/fig4/fig4.svg"), fig4, h = height, w = width, units = "cm", bg = "transparent", dpi = 1200, device = "svg")



## Rebuilding figure 4A
fig4a <- plot_list[[1]]



## rebuilding fig 4E
fig4e <- plot_list[[3]]

## Creating first column
col1 <- plot_grid(fig4a, fig4e,
                  ncol = 1,
                  nrow = 2,
                  # rel_heights = c(0.25, 0.25, 0.5),
                  axis = "tblr", 
                  align = "h")

## Building fig 4F
fig4f <- plot_list[[4]]

## Creating second column
col2 <- plot_grid(NULL, fig4f,
                  ncol = 1,
                  nrow = 2,
                  # rel_heights = c(1.0, 0.93, 0.07),
                  axis = "l", 
                  align = "v")

## Building fig 4G
# fig4g_top <- plot_list[[5]][[1]]
# 
# fig4g_bottom <- plot_list[[5]][[2]]

# fig4g <- plot_grid(fig4g_top, NULL, fig4g_bottom,
#                    rel_heights = c(0.5, -0.25, 1),
#                    ncol = 1, axis = "l", align = "v")

## Rebuilding fig 4D
fig4d <- plot_list[[2]]

fig4d_aligned <- plot_grid(NULL, fig4d, NULL,
                           ncol = 3,
                           # nrow = 2,
                           # rel_heights = c(1, 6),
                           rel_widths = c(-0.00, 0.95, 0.05),
                           axis = "lr", 
                           align = "h")

## Building fig 4G
fig4gh <- plot_list[[5]]

## Building the third column
col3 <- plot_grid(fig4d_aligned, fig4gh,
                  # ncol = 2,
                  nrow = 2,
                  # rel_heights = c(1, 6),
                  # rel_widths = c(0.8, 0.2),
                  axis = "lr", 
                  align = "h")

col3

## All together
fig4 <- plot_grid(col1, col2, col3,
                  ncol = 3, 
                  axis = "tblr", 
                  align = "h"
)


# save_plot(plot = fig4, filename = paste0(filepath, "res/fig4/fig4.png"), base_height = height, base_width = width, units = "cm", bg = "transparent", dpi = 1200)
save_plot(plot = fig4, filename = paste0(filepath, "res/fig4/fig4.pdf"), base_height = height, base_width = width, units = "cm", bg = "transparent", dpi = 1200)
# ggsave(paste0(filepath, "res/fig4/fig4.svg"), fig4, h = height, w = width, units = "cm", bg = "transparent", dpi = 1200, device = "svg")

