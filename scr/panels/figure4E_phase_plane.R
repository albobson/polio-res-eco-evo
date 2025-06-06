############################### Figure 4 E #####################################

## Reason:

## This plot shows the outcome of every possible initial MOI and resistance
## frequency

#### Set Up                                                                 ####

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

## Read in fitness function and parameters
fit_func <- read.csv(paste0(filepath, "dat_gen/params/logistic_fitness_function.csv"))
optim_params <- read.csv(paste0(filepath, "dat_gen/params/optim_params.csv"))

## Read in the trajectory plot data
    ## Data generated in figure4B_trajectory_plot.R
phase_plane_var <- read.csv(paste0(filepath, "dat_gen/sims/traj_for_phase_plane.csv"))

## Source functions to run simulations
source("scr/polv_DDT_functions.R")  ## Main functions 

## geom_text() is in terms of mm, where element_text() is in terms of pt. I want
## my geom_text and axis text to be the same. Here, I am defining a conversion
## which will keep the text sizes the same between text types.
geom_text_conv = 0.3528

#### Running simulations                                                    ####
mois <- 10^seq(from = -2, to = 2.35, length.out = 11)

# r_freq <- 10^seq(from = -4, to = 0, length.out = 11)
r_freq <- seq(from = 0, to = 1, length.out = 11)

## Set up parallel environment:
cl <- makeCluster((detectCores() - 2))
registerDoParallel(cl)

## Run simulation
sims <- foreach(i = mois, .combine = 'rbind', .packages = c('dplyr')) %:%
  foreach(j = r_freq, .combine = 'rbind', .packages = c('dplyr')) %dopar% {
    (determ_polv(n = 1,
                 moi_wt_start = (1-j)*i,
                 moi_mut_start = (j)*i,
                 fit_func_in = fit_func$prob_surv,
                 v_prog = optim_params$optim_v_prog,
                 p2pfu = optim_params$optim_p2pfu,
                 ))
  }

stopImplicitCluster()
stopCluster(cl)

## Save the simulations
suppressWarnings(dir.create(paste0(filepath, "dat_gen/sims")))
write.csv(sims, file = paste0(filepath, "dat_gen/sims/phase_plane_sims.csv"))

## Clean simulation
sims_clean <- sims %>%
  select(type, time, moi_res, moi_wt, init_wt, init_mut, pop_prop) %>%
  filter(type == "resistant") %>%
  group_by(init_mut, init_wt) %>%
  arrange(time) %>%
  mutate(pop_prop = as.numeric(pop_prop),
         init_mut = as.numeric(init_mut),
         init_wt = as.numeric(init_wt),
         moi_res = as.numeric(moi_res),
         moi_wt = as.numeric(moi_wt)) %>%
  mutate(delta_p = (pop_prop - lag(pop_prop, default = first(pop_prop))),
         delta_moi_tot = ((moi_res+moi_wt) - lag((moi_res+moi_wt), default = first((moi_res+moi_wt)))),
         init_moi_tot = (init_wt+init_mut),
         init_freq_mut = init_mut/(init_wt+init_mut),
         dist_col = abs(delta_moi_tot) + abs(delta_p)) %>%
  mutate(x_end = (init_moi_tot+cos(atan(delta_p/delta_moi_tot))),
         y_end = (init_freq_mut+sin(atan(delta_p/delta_moi_tot)))) %>%
  filter(time == 1) %>%
  ungroup()




#### Plotting an unscaled version                                           ####
plot4e <- ggplot(sims_clean, aes(x = init_moi_tot, 
                                 y = init_freq_mut
                                 )) +
  geom_segment(aes(
    xend = (init_moi_tot+(delta_moi_tot/75)), 
    yend = (init_freq_mut+(delta_p/10))
    ),
    arrow = arrow(length = unit(0.1, "cm")), 
    linewidth = 0.4) +
  # theme_light() +
  ylim(c(0,1)) +
  scale_x_continuous(trans="log10", breaks=10^(-2:3),
                     labels = sapply(c(-2:3),function(i){parse(text = sprintf("10^%d",i))}),
                     limits = c(10^(-2.1), 10^(2.6))) +
  # scale_y_continuous(trans="log10", 
  #                    breaks=10^(-4:0),
  #                    labels = sapply(c(-4:0),function(i){parse(text = sprintf("10^%d",i))}),
  #                    limits = c(10^(-4), 10^(0))
  # ) +
  geom_segment(phase_plane_var, 
               mapping = aes(
                 x = moi_tot, 
                 y = pop_prop, 
                 xend = (moi_tot+(delta_moi)), 
                 yend = (pop_prop+(delta_p)),
                 # color = factor(moi_tot, levels = moi_tot)
                 ),
               arrow = arrow(length = unit(0.15, "cm")), 
               linewidth = 0.5, 
               color = "#717171") +
  # labs(color = "Total MOI") +
  # scale_colour_manual(values = c("#FFA61C", "#64A61C", rep("#000", 5))) +
  xlab("Total MOI") + 
  ylab(expression(f[Res])) + # paste0("\U0394 Freq. Resistance")) +
  theme(    
    ## Text size
    text = element_text(size= axis_text_size), 
    axis.text = element_text(size= axis_text_size),
    ## Text color
    axis.text.x = element_text(colour="black"),
    axis.text.y = element_text(colour="black"),
    axis.title.y = element_text(size=axis_text_size + axis_text_size * .4),
    axis.title.x = element_text(),
    ## Changing the lines in the plot 
    panel.grid.major.x = element_line(color = "lightgrey"), ## Major x lines
    panel.grid.minor.x = element_blank(),                   ## Minor x lines
    panel.grid.major.y = element_line(color = "lightgrey"), ## Major y lines
    panel.grid.minor.y = element_blank(), ## Minor y lines
    ## For faceted plots: The strip is the top of the facet
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    ## Remove the background of the plot and panel
    panel.background = element_rect(fill='transparent'), 
    plot.background = element_rect(fill='transparent', color=NA),
    ## Legend stuff
    # legend.position = "none",
    legend.position = "inside",
    legend.position.inside = c(0.85,0.15),
    legend.text = element_text(size = legend_text_size),
    legend.title = element_blank(),
    legend.background = element_rect(fill = "transparent"),
    legend.key = element_rect(fill = "transparent"),
    legend.key.width = unit(0, 'cm'),
    legend.key.size = unit(0.3, "cm")
        )

ggsave(paste0(filepath, "res/fig4/E.png"), plot4e, h = height/2, w = width/3, units = "cm", bg = "transparent", dpi = 1200)
ggsave(paste0(filepath, "res/fig4/E.pdf"), plot4e, h = height/2, w = width/3, units = "cm", bg = "transparent", dpi = 1200)
ggsave(paste0(filepath, "res/fig4/E.svg"), plot4e, h = height/2, w = width/3, units = "cm", bg = "transparent", dpi = 1200, device = 'svg')
saveRDS(plot4e, paste0(filepath, "res/fig4/E.rds"))


# #### Plotting an unscaled version                                           ####
# plot4e_unscaled <- ggplot(sims_clean, aes(x = init_moi_tot, 
#                                           y = init_freq_mut
# )) +
#   geom_segment(aes(
#     xend = (init_moi_tot+(delta_moi_tot)), 
#     yend = (init_freq_mut+(delta_p))
#   ),
#   arrow = arrow(length = unit(0.2, "cm")), 
#   linewidth = 0.6) +
#   # theme_light() +
#   ylim(c(0,1)) +
#   scale_x_continuous(trans="log10", breaks=10^(-2:3),
#                      labels = sapply(c(-2:3),function(i){parse(text = sprintf("10^%d",i))}),
#                      limits = c(10^(-2.1), 10^(2.6))) +
#   geom_segment(phase_plane_var, 
#                mapping = aes(
#                  x = moi_tot, 
#                  y = pop_prop, 
#                  xend = (moi_tot+(delta_moi)), 
#                  yend = (pop_prop+(delta_p)),
#                  # color = factor(moi_tot, levels = moi_tot)
#                ),
#                arrow = arrow(length = unit(0.3, "cm")), 
#                linewidth = 1.4) +
#   # labs(color = "Total MOI") +
#   # scale_colour_manual(values = c("#FFA61C", "#64A61C", rep("#000", 5))) +
#   xlab("Total MOI") + 
#   ylab("Resistance Frequency") +
#   theme(    
#     ## Text size
#     text = element_text(size= axis_text_size), 
#     axis.text = element_text(size= axis_text_size),
#     ## Text color
#     axis.text.x = element_text(colour="black"),
#     axis.text.y = element_text(colour="black"),
#     axis.title.y = element_text(),
#     axis.title.x = element_text(),
#     ## Changing the lines in the plot 
#     panel.grid.major.x = element_line(color = "lightgrey"), ## Major x lines
#     panel.grid.minor.x = element_blank(),                   ## Minor x lines
#     panel.grid.major.y = element_line(color = "lightgrey"), ## Major y lines
#     panel.grid.minor.y = element_blank(), ## Minor y lines
#     ## For faceted plots: The strip is the top of the facet
#     strip.background = element_blank(),
#     strip.text.x = element_blank(),
#     ## Remove the background of the plot and panel
#     panel.background = element_rect(fill='transparent'), 
#     plot.background = element_rect(fill='transparent', color=NA),
#     ## Legend stuff
#     # legend.position = "none",
#     legend.position = "inside",
#     legend.position.inside = c(0.85,0.15),
#     legend.text = element_text(size = legend_text_size),
#     legend.title = element_blank(),
#     legend.background = element_rect(fill = "transparent"),
#     legend.key = element_rect(fill = "transparent"),
#     legend.key.width = unit(0, 'cm'),
#     legend.key.size = unit(0.3, "cm")
#   )
# 
# plot4e_unscaled
# 
# ggsave(paste0(filepath, "res/fig4/E_unscaled.png"), plot4e_unscaled, h = height/2/2, w = width/3, units = "cm", bg = "transparent", dpi = 1200)
# ggsave(paste0(filepath, "res/fig4/E_unscaled.pdf"), plot4e_unscaled, h = height/2/2, w = width/3, units = "cm", bg = "transparent", dpi = 1200)
# ggsave(paste0(filepath, "res/fig4/E_unscaled.svg"), plot4e_unscaled, h = height/2/2, w = width/3, units = "cm", bg = "transparent", dpi = 1200, device = 'svg')
# # saveRDS(plot4e_unscaled, paste0(filepath, "res/fig4/E.rds"))
