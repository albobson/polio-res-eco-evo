#### Script to test for the effect of resistance on the model outcomes ####

## Set up
## Read in original functions
source("../scr/polv_DDT_functions.R")

## Set number of cores
n_cores <- 20*2

## Read in values for simulations
#### Read in the previous parameters
init_date <- "2025-09-02" 

## Fitness function
init_fit_func <- "logistic"

## Mutation rate
init_mu <- "2e-05"

color_in <- c("#c94d4d", "#688cc8")

## Axis text size
axis_text_size = 6

## Legend text size
legend_text_size = 6

height = 4*4

width = 3*4


## Create the filepath where things will be saved
filepath <- paste0("../runs/ddt_", init_fit_func, "_mu_", init_mu, "_", init_date, "/")

## Read in fitness function and parameters
fit_func <- read.csv(paste0(filepath, "dat_gen/params/logistic_fitness_function.csv"))
optim_params <- read.csv(paste0(filepath, "dat_gen/params/optim_params.csv"))
clin_trial_params <- read.csv(paste0(filepath, "dat_gen/params/optim_clin_trial_params.csv"))
# clin_trial_params <- NULL
# clin_trial_params$optim_t_imm <- 9.39
# clin_trial_params$optim_c_pop <- 51087
# clin_trial_params$optim_imm_m <- -1.32
# clin_trial_params$optim_imm_sd <- 0.51

#### Define the first three resistance functions                            ####
sub_vec <- 0:60

suppressWarnings(dir.create(paste0("fitness_functions/")))

## full resistance cost 0.1 -- fully resistant virion will have a fitness of 1-0.1
## Slope
res_fit9 <- 0.9
m_rfun9 <- (res_fit9-1)/(60-0)
rfun9 <- (m_rfun9*sub_vec) + 1
# plot(sub_vec, rfun9)

fit_func_9 <- data.frame(subunits = sub_vec, 
                         prob_surv = rfun9*fit_func$prob_surv)

fit_func_9_p <- ggplot(fit_func_9, aes(x = sub_vec, y = prob_surv)) +
  geom_line(
    # linewidth = 0.7
  ) +
  xlab("Resistant subunits") +
  ylab("Probability of survival") +
  scale_y_continuous(trans="log10",
                     limits = c(0.0001,1),
                     breaks=10^(-4:0),
                     labels = c(sapply(-4:-1,function(i){parse(text = sprintf("10^%d",i))}), "1 ")) +
  theme_light() +
  ggtitle("Res. capsid fitness = 0.9") +
  theme(
    ## Text size
    text = element_text(size= axis_text_size),
    axis.text = element_text(size= axis_text_size),
    plot.title = element_text(size = axis_text_size, hjust = 0.5),
    ## Text color
    axis.text.x = element_text(colour="black"),
    axis.text.y = element_text(colour="black"),
    axis.title.y = element_text(),
    axis.title.x = element_text(),
    ## Changing the lines in the plot 
    panel.grid.major.x = element_line(color = "lightgrey"), ## Major x lines
    # panel.grid.minor.x = element_blank(),                   ## Minor x lines
    panel.grid.major.y = element_line(color = "lightgrey"), ## Major y lines
    panel.grid.minor.y = element_line(color = "lightgrey"), ## Minor y lines
    ## For faceted plots: The strip is the top of the facet
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    ## Remove the background of the plot and panel
    panel.background = element_rect(fill='transparent'), 
    plot.background = element_rect(fill='transparent', color=NA),
    ## Legend stuff
    legend.position = "none",
  ) 

fit_func_9_p

ggsave("fitness_functions/fit_func_9_p.pdf", fit_func_9_p, w = width/3, h = height/4, units = "cm")
ggsave("fitness_functions/fit_func_9_p.jpg", fit_func_9_p, w = width/3, h = height/4, units = "cm")


## full resistance cost 0.9 -- fully resistant virion will have a fitness of 1-0.1
## Slope
res_fit1 <- 0.1
m_rfun1 <- (res_fit1-1)/(60-0)
rfun1 <- (m_rfun1*sub_vec) + 1
plot(sub_vec, rfun1)

fit_func_1 <- data.frame(subunits = sub_vec, 
                         prob_surv = rfun1*fit_func$prob_surv)

fit_func_1_p <- ggplot(fit_func_1, aes(x = sub_vec, y = prob_surv)) +
  geom_line(
    # linewidth = 0.7
  ) +
  xlab("Resistant subunits") +
  ylab("Probability of survival") +
  scale_y_continuous(trans="log10",
                     limits = c(0.0001,1),
                     breaks=10^(-4:0),
                     labels = c(sapply(-4:-1,function(i){parse(text = sprintf("10^%d",i))}), "1 ")) +
  theme_light() +
  ggtitle("Res. capsid fitness = 0.1") +
  theme(
    ## Text size
    text = element_text(size= axis_text_size),
    axis.text = element_text(size= axis_text_size),
    plot.title = element_text(size = axis_text_size, hjust = 0.5),
    ## Text color
    axis.text.x = element_text(colour="black"),
    axis.text.y = element_text(colour="black"),
    axis.title.y = element_text(),
    axis.title.x = element_text(),
    ## Changing the lines in the plot 
    panel.grid.major.x = element_line(color = "lightgrey"), ## Major x lines
    # panel.grid.minor.x = element_blank(),                   ## Minor x lines
    panel.grid.major.y = element_line(color = "lightgrey"), ## Major y lines
    panel.grid.minor.y = element_line(color = "lightgrey"), ## Minor y lines
    ## For faceted plots: The strip is the top of the facet
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    ## Remove the background of the plot and panel
    panel.background = element_rect(fill='transparent'), 
    plot.background = element_rect(fill='transparent', color=NA),
    ## Legend stuff
    legend.position = "none",
  ) 

fit_func_1_p

ggsave("fitness_functions/fit_func_1_p.pdf", fit_func_1_p, w = width/3, h = height/4, units = "cm")
ggsave("fitness_functions/fit_func_1_p.jpg", fit_func_1_p, w = width/3, h = height/4, units = "cm")


## full resistance cost 0.99 -- fully resistant virion will have a fitness of 1-0.1
## Slope
res_fit01 <- 0.01
m_rfun01 <- (res_fit01-1)/(60-0)
rfun01 <- (m_rfun01*sub_vec) + 1
plot(sub_vec, rfun01)

fit_func_01 <- data.frame(subunits = sub_vec, 
                         prob_surv = rfun01*fit_func$prob_surv)

fit_func_01_p <- ggplot(fit_func_01, aes(x = sub_vec, y = prob_surv)) +
  geom_line(
    # linewidth = 0.7
  ) +
  xlab("Resistant subunits") +
  ylab("Probability of survival") +
  scale_y_continuous(trans="log10",
                     limits = c(0.0001,1),
                     breaks=10^(-4:0),
                     labels = c(sapply(-4:-1,function(i){parse(text = sprintf("10^%d",i))}), "1 ")) +
  theme_light() +
  ggtitle("Res. capsid fitness = 0.01") +
  theme(
    ## Text size
    text = element_text(size= axis_text_size),
    axis.text = element_text(size= axis_text_size),
    plot.title = element_text(size = axis_text_size, hjust = 0.5),
    ## Text color
    axis.text.x = element_text(colour="black"),
    axis.text.y = element_text(colour="black"),
    axis.title.y = element_text(),
    axis.title.x = element_text(),
    ## Changing the lines in the plot 
    panel.grid.major.x = element_line(color = "lightgrey"), ## Major x lines
    # panel.grid.minor.x = element_blank(),                   ## Minor x lines
    panel.grid.major.y = element_line(color = "lightgrey"), ## Major y lines
    panel.grid.minor.y = element_line(color = "lightgrey"), ## Minor y lines
    ## For faceted plots: The strip is the top of the facet
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    ## Remove the background of the plot and panel
    panel.background = element_rect(fill='transparent'), 
    plot.background = element_rect(fill='transparent', color=NA),
    ## Legend stuff
    legend.position = "none",
  ) 

fit_func_01_p

ggsave("fitness_functions/fit_func_01_p.pdf", fit_func_01_p, w = width/3, h = height/4, units = "cm")
ggsave("fitness_functions/fit_func_01_p.jpg", fit_func_01_p, w = width/3, h = height/4, units = "cm")

################################################################################
########################### Run the simulations ################################
################################################################################

suppressWarnings(dir.create(paste0("trajectories/")))

## Deterministic 9 trajectory

determ_run_9 <- determ_polv(n = 6,
                         moi_mut_start = 100*0.0001, 
                         moi_wt_start = 100*(1-0.0001),
                         fit_func_in = fit_func_9$prob_surv,
                         v_prog = optim_params$optim_v_prog,
                         p2pfu = optim_params$optim_p2pfu
) %>%
  mutate(fit_type = "Full resistant fitness = 0.9")

#### Plotting 9 trajectory                                                  ####
determ_9_p <- ggplot(determ_run_9, aes(x = time, y = moi_type, color = type, linetype = fit_type)) +
  # geom_point(size = 3) + 
  geom_line(linewidth = 0.7) +
  scale_color_manual(labels = c("Res.", "Sus."),
                     values = color_in,
                     name = "Genotype") +
  theme_light() + 
  scale_linetype(guide = "none") +
  xlab("Passages") + 
  ylab("Viral density (genomes/cell)") + 
  scale_y_log10(breaks=10^(-6:10),
                labels = sapply(-6:10,function(i){parse(text = sprintf("10^%d",i))}),
                limits = c(10^-6.2, 10^2.4)) +
  scale_x_continuous(breaks = c(0:max(determ_run_9$time))) +
  # guides(
  #   color = "none"
  # ) +
  # ggtitle(paste0(unique(determ_run_9$fit_type))) +
  theme(
    ## Text size
    text = element_text(size= axis_text_size),
    axis.text = element_text(size= axis_text_size),
    ## Text color
    axis.text.x = element_text(colour="black"),
    axis.text.y = element_text(colour="black"),
    axis.title.y = element_text(),
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
    legend.position.inside = c(0.82, 0.11),
    legend.text = element_text(size = legend_text_size),
    legend.title = element_blank(),
    # legend.key.width = unit(0, 'cm'),
    legend.key.size = unit(0.3, "cm"),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.2),  # Fills legend background & adds a border
    legend.box.background = element_rect(color = "black", linewidth = 0),  # Adds an outer border
    legend.margin = margin(0.6, 1, 0.6, 0.6),
    legend.spacing.y = unit(-0.0, "cm"),
    # legend.box = "horizontal"
  ) 

determ_9_p

ggsave("trajectories/determ_9_p.pdf", determ_9_p, w = width/3, h = height/4, units = "cm")
ggsave("trajectories/determ_9_p.jpg", determ_9_p, w = width/3, h = height/4, units = "cm")


## Deterministic 1 trajectory
determ_run_1 <- determ_polv(n = 6,
                            moi_mut_start = 100*0.0001, 
                            moi_wt_start = 100*(1-0.0001),
                            fit_func_in = fit_func_1$prob_surv,
                            v_prog = optim_params$optim_v_prog,
                            p2pfu = optim_params$optim_p2pfu
) %>%
  mutate(fit_type = "Full resistant fitness = 0.1")

#### Plotting 1 trajectory                                                  ####
determ_1_p <- ggplot(determ_run_1, aes(x = time, y = moi_type, color = type, linetype = fit_type)) +
  # geom_point(size = 3) + 
  geom_line(linewidth = 0.7) +
  scale_color_manual(labels = c("Res.", "Sus."),
                     values = color_in,
                     name = "Genotype") +
  theme_light() + 
  scale_linetype(guide = "none") +
  xlab("Passages") + 
  ylab("Viral density (genomes/cell)") + 
  scale_y_log10(breaks=10^(-6:10),
                labels = sapply(-6:10,function(i){parse(text = sprintf("10^%d",i))}),
                limits = c(10^-6.2, 10^2.4)) +
  scale_x_continuous(breaks = c(0:max(determ_run_9$time))) +
  # guides(
  #   color = "none"
  # ) +
  # ggtitle(paste0(unique(determ_run_1$fit_type))) +
  theme(
    ## Text size
    text = element_text(size= axis_text_size),
    axis.text = element_text(size= axis_text_size),
    ## Text color
    axis.text.x = element_text(colour="black"),
    axis.text.y = element_text(colour="black"),
    axis.title.y = element_text(),
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
    legend.position = "none",
    # legend.position = "inside",
    # legend.position.inside = c(0.66, 0.096),
    legend.text = element_text(size = legend_text_size),
    legend.title = element_blank(),
    # legend.key.width = unit(0, 'cm'),
    legend.key.size = unit(0.3, "cm"),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.2),  # Fills legend background & adds a border
    legend.box.background = element_rect(color = "black", linewidth = 0),  # Adds an outer border
    legend.margin = margin(0.6, 1, 0.6, 0.6),
    legend.spacing.y = unit(-0.0, "cm"),
    # legend.box = "horizontal"
  ) 

determ_1_p

ggsave("trajectories/determ_1_p.pdf", determ_1_p, w = width/3, h = height/4, units = "cm")
ggsave("trajectories/determ_1_p.jpg", determ_1_p, w = width/3, h = height/4, units = "cm")

## Deterministic 01 trajectory
determ_run_01 <- determ_polv(n = 6,
                            moi_mut_start = 100*0.0001, 
                            moi_wt_start = 100*(1-0.0001),
                            fit_func_in = fit_func_01$prob_surv,
                            v_prog = optim_params$optim_v_prog,
                            p2pfu = optim_params$optim_p2pfu
) %>%
  mutate(fit_type = "Full resistant fitness = 0.01")

#### Plotting 01 trajectory                                                  ####
determ_01_p <- ggplot(determ_run_01, aes(x = time, y = moi_type, color = type, linetype = fit_type)) +
  # geom_point(size = 3) + 
  geom_line(linewidth = 0.7) +
  scale_color_manual(labels = c("Res.", "Sus."),
                     values = color_in,
                     name = "Genotype") +
  theme_light() + 
  scale_linetype(guide = "none") +
  xlab("Passages") + 
  ylab("Viral density (genomes/cell)") + 
  scale_y_log10(breaks=10^(-6:10),
                labels = sapply(-6:10,function(i){parse(text = sprintf("10^%d",i))}),
                limits = c(10^-6.2, 10^2.4)) +
  scale_x_continuous(breaks = c(0:max(determ_run_9$time))) +
  # guides(
  #   color = "none"
  # ) +
  # ggtitle(paste0(unique(determ_run_01$fit_type))) +
  theme(
    ## Text size
    text = element_text(size= axis_text_size),
    axis.text = element_text(size= axis_text_size),
    ## Text color
    axis.text.x = element_text(colour="black"),
    axis.text.y = element_text(colour="black"),
    axis.title.y = element_text(),
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
    legend.position = "none",
    # legend.position = "inside",
    legend.position.inside = c(0.66, 0.096),
    legend.text = element_text(size = legend_text_size),
    legend.title = element_blank(),
    # legend.key.width = unit(0, 'cm'),
    legend.key.size = unit(0.3, "cm"),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.2),  # Fills legend background & adds a border
    legend.box.background = element_rect(color = "black", linewidth = 0),  # Adds an outer border
    legend.margin = margin(0.6, 1, 0.6, 0.6),
    legend.spacing.y = unit(-0.0, "cm"),
    # legend.box = "horizontal"
  ) 

determ_01_p


ggsave("trajectories/determ_01_p.pdf", determ_01_p, w = width/3, h = height/4, units = "cm")
ggsave("trajectories/determ_01_p.jpg", determ_01_p, w = width/3, h = height/4, units = "cm")


################################################################################
############################## Vector fields ###################################
################################################################################
suppressWarnings(dir.create(paste0("vector_fields/")))

mois <- 10^seq(from = -2, to = 2.35, length.out = 11)

# r_freq <- 10^seq(from = -4, to = 0, length.out = 11)
r_freq <- seq(from = 0, to = 1, length.out = 11)

## Set up parallel environment:
cl <- makeCluster((n_cores - 2))
registerDoParallel(cl)

## Run simulation
sims_9 <- foreach(i = mois, .combine = 'rbind', .packages = c('dplyr')) %:%
  foreach(j = r_freq, .combine = 'rbind', .packages = c('dplyr')) %dopar% {
    (determ_polv(n = 1,
                 moi_wt_start = (1-j)*i,
                 moi_mut_start = (j)*i,
                 fit_func_in = fit_func_9$prob_surv,
                 v_prog = optim_params$optim_v_prog,
                 p2pfu = optim_params$optim_p2pfu,
    ))
  }

stopImplicitCluster()
stopCluster(cl)

## Save the simulations
suppressWarnings(dir.create(paste0("sims")))
write.csv(sims_9, file = paste0("sims/phase_plain_sims_9.csv"))

sims_9 <- read.csv("sims/phase_plain_sims_9.csv")

## Clean simulation
sims_9_clean <- sims_9 %>%
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
sims_9_plot <- ggplot(sims_9_clean, aes(x = init_moi_tot, 
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
  # geom_segment(phase_plane_var, 
  #              mapping = aes(
  #                x = moi_tot, 
  #                y = pop_prop, 
  #                xend = (moi_tot+(delta_moi)), 
  #                yend = (pop_prop+(delta_p)),
  #                # color = factor(moi_tot, levels = moi_tot)
  #              ),
  #              arrow = arrow(length = unit(0.15, "cm")), 
  #              linewidth = 0.5, 
  #              color = "#717171") +
  # labs(color = "Total MOI") +
  # scale_colour_manual(values = c("#FFA61C", "#64A61C", rep("#000", 5))) +
  xlab("Total MOI") + 
  ylab(expression(f[Res])) + # paste0("\U0394 Freq. Resistance")) +
  # ggtitle(paste0(unique(determ_run_9$fit_type))) +
  theme(    
    ## Text size
    text = element_text(size= axis_text_size),
    axis.text = element_text(size= axis_text_size),
    ## Text color
    axis.text.x = element_text(colour="black"),
    axis.text.y = element_text(colour="black"),
    # axis.title.y = element_text(size=axis_text_size + axis_text_size * .4),
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
    legend.position = "none",
    # legend.position = "inside",
    # legend.position.inside = c(0.85,0.15),
    # legend.text = element_text(size = legend_text_size),
    legend.title = element_blank(),
    legend.background = element_rect(fill = "transparent"),
    legend.key = element_rect(fill = "transparent"),
    legend.key.width = unit(0, 'cm'),
    legend.key.size = unit(0.3, "cm")
  )

sims_9_plot

ggsave("vector_fields/sims_9_plot.pdf", sims_9_plot, w = width/3, h = height/4, units = "cm")
ggsave("vector_fields/sims_9_plot.jpg", sims_9_plot, w = width/3, h = height/4, units = "cm")



#### Res cost 01                                                            ####
## Set up parallel environment:
cl <- makeCluster((n_cores - 2))
registerDoParallel(cl)

## Run simulation
sims_1 <- foreach(i = mois, .combine = 'rbind', .packages = c('dplyr')) %:%
  foreach(j = r_freq, .combine = 'rbind', .packages = c('dplyr')) %dopar% {
    (determ_polv(n = 1,
                 moi_wt_start = (1-j)*i,
                 moi_mut_start = (j)*i,
                 fit_func_in = fit_func_1$prob_surv,
                 v_prog = optim_params$optim_v_prog,
                 p2pfu = optim_params$optim_p2pfu,
    ))
  }

stopImplicitCluster()
stopCluster(cl)

## Save the simulations
suppressWarnings(dir.create(paste0("sims")))
write.csv(sims_1, file = paste0("sims/phase_plain_sims_1.csv"))

sims_1 <- read.csv("sims/phase_plain_sims_1.csv")

## Clean simulation
sims_1_clean <- sims_1 %>%
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
sims_1_plot <- ggplot(sims_1_clean, aes(x = init_moi_tot, 
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
  # ggtitle(paste0(unique(determ_run_1$fit_type))) +
  # scale_y_continuous(trans="log10", 
  #                    breaks=10^(-4:0),
  #                    labels = sapply(c(-4:0),function(i){parse(text = sprintf("10^%d",i))}),
  #                    limits = c(10^(-4), 10^(0))
  # ) +
  # geom_segment(phase_plane_var, 
  #              mapping = aes(
  #                x = moi_tot, 
  #                y = pop_prop, 
  #                xend = (moi_tot+(delta_moi)), 
  #                yend = (pop_prop+(delta_p)),
  #                # color = factor(moi_tot, levels = moi_tot)
  #              ),
  #              arrow = arrow(length = unit(0.15, "cm")), 
  #              linewidth = 0.5, 
  #              color = "#717171") +
  # labs(color = "Total MOI") +
  # scale_colour_manual(values = c("#FFA61C", "#64A61C", rep("#000", 5))) +
  xlab("Total MOI") + 
  ylab(expression(f[Res])) + # paste0("\U0314 Freq. Resistance")) +
  theme(    
    ## Text size
    text = element_text(size= axis_text_size),
    axis.text = element_text(size= axis_text_size),
    ## Text color
    axis.text.x = element_text(colour="black"),
    axis.text.y = element_text(colour="black"),
    # axis.title.y = element_text(size=axis_text_size + axis_text_size * .4),
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
    legend.position = "none",
    # legend.position = "inside",
    # legend.position.inside = c(0.85,0.15),
    # legend.text = element_text(size = legend_text_size),
    legend.title = element_blank(),
    legend.background = element_rect(fill = "transparent"),
    legend.key = element_rect(fill = "transparent"),
    legend.key.width = unit(0, 'cm'),
    legend.key.size = unit(0.3, "cm")
  )

sims_1_plot

ggsave("vector_fields/sims_1_plot.pdf", sims_1_plot, w = width/3, h = height/4, units = "cm")
ggsave("vector_fields/sims_1_plot.jpg", sims_1_plot, w = width/3, h = height/4, units = "cm")

#### Res cost 01                                                            ####
## Set up parallel environment:
cl <- makeCluster((n_cores - 2))
registerDoParallel(cl)

## Run simulation
sims_01 <- foreach(i = mois, .combine = 'rbind', .packages = c('dplyr')) %:%
  foreach(j = r_freq, .combine = 'rbind', .packages = c('dplyr')) %dopar% {
    (determ_polv(n = 1,
                 moi_wt_start = (1-j)*i,
                 moi_mut_start = (j)*i,
                 fit_func_in = fit_func_01$prob_surv,
                 v_prog = optim_params$optim_v_prog,
                 p2pfu = optim_params$optim_p2pfu,
    ))
  }

stopImplicitCluster()
stopCluster(cl)

## Save the simulations
suppressWarnings(dir.create(paste0("sims")))
write.csv(sims_01, file = paste0("sims/phase_plain_sims_01.csv"))

sims_01 <- read.csv("sims/phase_plain_sims_01.csv")

## Clean simulation
sims_01_clean <- sims_01 %>%
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
sims_01_plot <- ggplot(sims_01_clean, aes(x = init_moi_tot, 
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
  # ggtitle(paste0(unique(determ_run_01$fit_type))) +
  # scale_y_continuous(trans="log10", 
  #                    breaks=10^(-4:0),
  #                    labels = sapply(c(-4:0),function(i){parse(text = sprintf("10^%d",i))}),
  #                    limits = c(10^(-4), 10^(0))
  # ) +
  # geom_segment(phase_plane_var, 
  #              mapping = aes(
  #                x = moi_tot, 
  #                y = pop_prop, 
  #                xend = (moi_tot+(delta_moi)), 
  #                yend = (pop_prop+(delta_p)),
  #                # color = factor(moi_tot, levels = moi_tot)
  #              ),
  #              arrow = arrow(length = unit(0.15, "cm")), 
  #              linewidth = 0.5, 
  #              color = "#717171") +
  # labs(color = "Total MOI") +
  # scale_colour_manual(values = c("#FFA61C", "#64A61C", rep("#000", 5))) +
  xlab("Total MOI") + 
  ylab(expression(f[Res])) + # paste0("\U0314 Freq. Resistance")) +
  theme(    
    ## Text size
    text = element_text(size= axis_text_size),
    axis.text = element_text(size= axis_text_size),
    ## Text color
    axis.text.x = element_text(colour="black"),
    axis.text.y = element_text(colour="black"),
    # axis.title.y = element_text(size=axis_text_size + axis_text_size * .4),
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
    legend.position = "none",
    # legend.position = "inside",
    # legend.position.inside = c(0.85,0.15),
    # legend.text = element_text(size = legend_text_size),
    legend.title = element_blank(),
    legend.background = element_rect(fill = "transparent"),
    legend.key = element_rect(fill = "transparent"),
    legend.key.width = unit(0, 'cm'),
    legend.key.size = unit(0.3, "cm")
  )

sims_01_plot

ggsave("vector_fields/sims_01_plot.pdf", sims_01_plot, w = width/3, h = height/4, units = "cm")
ggsave("vector_fields/sims_01_plot.jpg", sims_01_plot, w = width/3, h = height/4, units = "cm")


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
########################### Clinical trials ####################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

suppressWarnings(dir.create(paste0("clinical_trials/")))
## Range of resistance to run over
n_run <- 10
r_costs <- 10^seq(from = -0.0457, to = -2, length.out = n_run)

n_72 <- 70
n_24 <- 23

## Null dataframe to store everything
all_df <- NULL

cc <- 0

for(i in 1:n_run) {
  ## Generate the new fitness function
  temp_fit_func <- NULL
  
  res_c <- r_costs[i]
  m_rfun <- (res_c-1)/(60-0)
  rfun <- (m_rfun*sub_vec) + 1
  
  temp_fit_func <- data.frame(subunits = sub_vec, 
                           prob_surv = rfun*fit_func$prob_surv)
  
  cl <- makeCluster(n_cores - 1)
  registerDoParallel(cl)
  
  set.seed(1249128)
  trial_72_h <- foreach(k = 1:n_72, .combine = "rbind", 
                        .packages = c('dplyr')) %dopar% {
                          stoch_polv(n = 300, 
                                     moi_wt_start = 1, 
                                     moi_mut_start = 0, 
                                     t_pocap = 9, 
                                     id = k,
                                     imm_delay = clin_trial_params$optim_t_imm,
                                     c_pop = clin_trial_params$optim_c_pop, 
                                     imm_m = clin_trial_params$optim_imm_m, 
                                     imm_sd = clin_trial_params$optim_imm_sd,
                                     fit_func_in = temp_fit_func$prob_surv,
                                     v_prog = optim_params$optim_v_prog,
                                     p2pfu = optim_params$optim_p2pfu,
                                     seed_in = k*2
                          ) %>%
                            mutate(treatment = paste0("original"),
                                   res_cost = res_c) %>%
                            mutate(time = time/3) ## Changing from replications to days
                        }
  
  print(paste0("72h done"))
  
  set.seed(12491284)
  trial_24_h <- foreach(k = n_72:(n_72+n_24), .combine = "rbind") %dopar% {
    stoch_polv(n = 300, 
               moi_wt_start = 1, 
               moi_mut_start = 0, 
               t_pocap = 3, 
               id = k,
               imm_delay = clin_trial_params$optim_t_imm,
               c_pop = clin_trial_params$optim_c_pop, 
               imm_m = clin_trial_params$optim_imm_m, 
               imm_sd = clin_trial_params$optim_imm_sd,
               fit_func_in = temp_fit_func$prob_surv,
               v_prog = optim_params$optim_v_prog,
               p2pfu = optim_params$optim_p2pfu,
               seed_in = k*2
    ) %>% 
      mutate(treatment = paste0("original"),
             res_cost = res_c) %>%
      mutate(time = time/3) ## Changing from replications to days
  }
  
  print(paste0("24h done"))
  
  stopCluster(cl)
  stopImplicitCluster()
  
  trial <- rbind(trial_72_h, trial_24_h)
  
  write.csv(x = trial, file = paste0("clin_trial_res_cost_", round(res_c, 3), ".csv"))
  
  all_df <- rbind(all_df, trial)
  
  cc <- cc + 1
  print(cc)
  
}

## Save all_df for ease
write.csv(x = all_df, file = "clin_trial_res_cost_all.csv")


all_df <- read.csv("clin_trial_res_cost_all.csv")

sims_id <- all_df %>%
  group_by(wt_rep_abil, id, c_pop, res_cost) %>%
  filter(type == "resistant", time >= time_to_pocap/3) %>%
  mutate(fifty = pop_prop >= 0.5 & !duplicated(pop_prop >= 0.5),
         pop_prop = ifelse(is.na(pop_prop), 0, pop_prop)) %>% 
  reframe(max_time = max(time), 
          tot_v = sum(tot_pfu), 
          perc_res = sum(surv_pfu)/sum(as.numeric(tot_pfu)),
          max_v = max(tot_pfu),
          c_pop = c_pop,
          t_50 = ifelse(fifty, time, NA)) %>%
  # filter(!is.na(t_50)) %>%
  unique() %>%
  ungroup()

## Creating a stat summary for mean and standard dev
sims_stats_all <- sims_id %>%
  ungroup() %>%
  group_by(res_cost) %>%
  summarize(clear_mean = mean(max_time),
            clear_sd = sd(max_time),
            res_mean = mean(perc_res),
            res_sd = sd(perc_res),
            vl_mean = mean(tot_v),
            vl_sd = sd(tot_v)) %>%
  ungroup()

## Testing a polygon version of this
sim_sds <- sims_stats_all %>%
  mutate(
    clear_mean_plus_sd = clear_mean + clear_sd,
    clear_mean_minus_sd = clear_mean - clear_sd,
    res_mean_plus_sd = res_mean + res_sd,
    res_mean_minus_sd = res_mean - res_sd,
    vl_mean_plus_sd = vl_mean + vl_sd,
    vl_mean_minus_sd = vl_mean - vl_sd,
  )

## Changing the minimum values to 0, since I'm logging the y axis in the sum
## viral load, and I get an error
sim_sds$vl_mean_minus_sd <- ifelse(sim_sds$vl_mean_minus_sd < 0, 10^round(min(log10(sims_id$tot_v))), sim_sds$vl_mean_minus_sd)

#### Computing rolling average/sd                                           ####

# Define the window size (e.g., 100)
window_size <- 5

library(zoo)

## Compute rolling SD
sims_id <- sims_id %>%
  mutate(
    # Rolling standard deviation for max_time, perc_res, and tot_v
    clear_rolling_sd = rollapply(max_time, width = window_size, FUN = sd, fill = NA, align = "center"),
    res_rolling_sd = rollapply(perc_res, width = window_size, FUN = sd, fill = NA, align = "center"),
    vl_rolling_sd = rollapply(tot_v, width = window_size, FUN = sd, fill = NA, align = "center"),
    
    # Rolling variance (if needed)
    clear_rolling_mean = rollmean(max_time, k = window_size, fill = NA, align = "center"),
    res_rolling_mean = rollmean(perc_res, k = window_size, fill = NA, align = "center"),
    vl_rolling_mean = rollmean(tot_v, k = window_size, fill = NA, align = "center")
  ) %>%
  ungroup()

# Create a dataset for ribbon to show the +/- rolling standard deviations
sims_id <- sims_id %>%
  mutate(
    clear_mean_plus_sd = clear_rolling_mean + clear_rolling_sd,
    clear_mean_minus_sd = clear_rolling_mean - clear_rolling_sd,
    res_mean_plus_sd = res_rolling_mean + res_rolling_sd,
    res_mean_minus_sd = res_rolling_mean - res_rolling_sd,
    vl_mean_plus_sd = vl_rolling_mean + vl_rolling_sd,
    vl_mean_minus_sd = vl_rolling_mean - vl_rolling_sd
  )


#### Plot the clear times                                                   ####
# Calculate rolling average with a window size of 5
sims_id$rolling_avg_clear <- rollmean(sims_id$max_time, k = 5, fill = NA, align = "center")

## Calculate the fold difference from pocapavir
sims_id <- sims_id %>%
  mutate(fold_strin = sims_id$wt_rep_abil/min(sims_id$wt_rep_abil))



res_cost_clear_dates_p <- ggplot(sims_id, aes(x = res_cost, 
                               # group = wt_rep_abil,
                               # fill=factor(log(wt_rep_abil))
)) +
  # geom_boxplot(notch=FALSE, outlier.shape=NA,
  #              color = "black") +
  geom_jitter(  size = 0.6, 
                alpha = 0.4, 
                shape = 16,
                height  = 0,
                width = NULL,
                aes(y = as.numeric(max_time), color = perc_res)) +
  # geom_ribbon(
  #   # sim_sds,
  #   mapping = aes(ymin = clear_mean_minus_sd, ymax = clear_mean_plus_sd),
  #   alpha = 0.2) +
  # geom_line(sims_stats_all, mapping = aes(y = clear_mean), linewidth = 3) + ## Mean line
  # geom_line(aes(y = clear_rolling_mean), color = "black", size = 0.7) +
  scale_color_continuous(name = expression(italic(f) * "(Res.)"),
                         low = color_in[2],
                         high = color_in[1],
                         limits = c(0, 1),
  ) +
  geom_point(data = sims_stats_all,
             mapping = aes(y = clear_mean),
             color = "black",
             size = 1,
             shape = 3
  ) +
  geom_errorbar(data = sim_sds,
                aes(ymin=clear_mean_minus_sd, ymax=clear_mean_plus_sd),
                width=.08,
                # position=position_dodge(.05)
  ) +
  theme_light() +
  xlab("Fitness of a fully resistant capsid") + ylab("Clearance date (DPI)") +
  scale_x_log10(
    # breaks=10^(-5:0),
    # limits = c((min(sims_id$wt_rep_abil)-min(sims_id$wt_rep_abil)*0.1), 1.5),
    # labels = sapply(-5:0,function(j){parse(text = sprintf("10^%d",j))})
  ) +
  theme(
    ## Text size
    text = element_text(size= axis_text_size), 
    axis.text = element_text(size= axis_text_size),
    ## Text color
    axis.text.x = element_text(colour="black"),
    axis.text.y = element_text(colour="black"),
    axis.title.y = element_text(),
    axis.title.x = element_text(),
    ## Changing the lines in the plot 
    panel.grid.major.x = element_line(color = "lightgrey"), ## Major x lines
    # panel.grid.minor.x = element_blank(),                   ## Minor x lines
    panel.grid.major.y = element_line(color = "lightgrey"), ## Major y lines
    panel.grid.minor.y = element_line(color = "lightgrey"), ## Minor y lines
    ## For faceted plots: The strip is the top of the facet
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    ## Remove the background of the plot and panel
    panel.background = element_rect(fill='transparent'), 
    plot.background = element_rect(fill='transparent', color=NA),
    ## Legend stuff
    legend.position = "none",
    # legend.position = "inside",
    # legend.position.inside = c(0.35, 0.94), ## Upper left
    # legend.position.inside = c(0.65, 0.06), ## bottom right
    # legend.text = element_text(size = legend_text_size),
    # legend.title = element_blank(),
    # legend.background = element_rect(fill = "transparent"),
    # legend.key = element_rect(fill = "transparent"),
    # legend.key.width = unit(0, 'cm'),
    # legend.key.size = unit(0.3, "cm")
  )

res_cost_clear_dates_p

ggsave("clinical_trials/res_cost_clear_dates_p.pdf", res_cost_clear_dates_p, w = width/2, h = height/4, units = "cm")
ggsave("clinical_trials/res_cost_clear_dates_p.jpg", res_cost_clear_dates_p, w = width/2, h = height/4, units = "cm")


#### Plot the viral load                                                    ####
sims_id$vl_mean_minus_sd_adjusted <- ifelse(sims_id$vl_mean_minus_sd <= 0, 1000, sims_id$vl_mean_minus_sd)

res_cost_viral_load_p <- ggplot(sims_id, aes(x = res_cost, 
                               # group = wt_rep_abil,
                               # fill=factor(log(wt_rep_abil))
)) +
  # geom_boxplot(notch=FALSE, outlier.shape=NA,
  #              color = "black") +
  geom_jitter(mapping = aes(y = tot_v, 
                            color = perc_res
  ), 
  size = 0.6, 
  alpha = 0.4, 
  # shape = 16,
  height  = 0) +
  # geom_ribbon(sim_sds,
  # mapping = aes(ymin = vl_mean_minus_sd, ymax = vl_mean_plus_sd),
  # alpha = 0.2) +
  # geom_line(sims_stats_all, mapping = aes(y = vl_mean), linewidth = 3) + ## Mean line
  # geom_smooth(aes(y = tot_v)) +
  ## Rolling average
  # geom_ribbon(
  #   # data = sim_rolling_sds,
  #   mapping = aes(ymin = vl_mean_minus_sd_adjusted, ymax = vl_mean_plus_sd),
  #   alpha = 0.2) +
  
  # Add line for the mean values (from summary stats)
  # geom_line(data = sims_stats_rolling, mapping = aes(y = clear_mean), linewidth = 3) +
  
  # geom_line(aes(y = vl_rolling_mean), color = "black", size = 0.7) +
  # geom_line(aes(y = rolling_avg_vl), color = "black", size = 1.5) +
  # Add spline fit using GAM (Generalized Additive Model)
  # geom_smooth(mapping = aes(y = tot_v), 
  #             method = "gam", 
  #             formula = y ~ s(x),  # Use a smooth function for the data
  #             se = FALSE,  # Optionally, set to TRUE if you want a confidence interval
  #             color = "black", 
  #             size = 1.5) +
  geom_point(data = sims_stats_all,
             mapping = aes(y = vl_mean),
             color = "black",
             size = 1,
             shape = 3
  ) +
  geom_errorbar(data = sim_sds,
                aes(ymin=vl_mean_minus_sd, ymax=vl_mean_plus_sd),
                width=.08,
                # position=position_dodge(.05)
  ) +
  scale_color_continuous(name = expression(f[Res]),
                         low = color_in[2],
                         high = color_in[1],
                         limits = c(0, 1),
  ) +
  theme_light() +
  xlab("Fitness of a fully resistant capsid") + ylab("Sum total viral load") +
  scale_x_log10(
    # breaks=10^(-5:0),
    # limits = c((min(sims_id$wt_rep_abil)-min(sims_id$wt_rep_abil)*0.1), 1.5),
    # labels = sapply(-5:0,function(j){parse(text = sprintf("10^%d",j))})
  ) +
  scale_y_log10(breaks=10^(round(min(log10(sims_id$tot_v))):round(max(log10(sims_id$tot_v)))),
                # limits = 10^c(round(min(log10(sims_id$tot_v))), round(max(log10(sims_id$tot_v)))),
                labels = sapply(round(min(log10(sims_id$tot_v))):round(max(log10(sims_id$tot_v))),function(j){parse(text = sprintf("10^%d",j))}),
                expand = expansion(mult = c(0, 0.05))) +
  theme(
    ## Text size
    text = element_text(size= axis_text_size), 
    axis.text = element_text(size= axis_text_size),
    ## Text color
    axis.text.x = element_text(colour="black"),
    axis.text.y = element_text(colour="black"),
    axis.title.y = element_text(),
    axis.title.x = element_text(),
    ## Changing the lines in the plot 
    panel.grid.major.x = element_line(color = "lightgrey"), ## Major x lines
    # panel.grid.minor.x = element_blank(),                   ## Minor x lines
    panel.grid.major.y = element_line(color = "lightgrey"), ## Major y lines
    panel.grid.minor.y = element_line(color = "lightgrey"), ## Minor y lines
    ## For faceted plots: The strip is the top of the facet
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    ## Remove the background of the plot and panel
    panel.background = element_rect(fill='transparent'), 
    plot.background = element_rect(fill='transparent', color=NA),
    ## Legend stuff
    # legend.position = "none",
    # legend.position = "inside",
    # legend.position.inside = c(0.35, 0.94), ## Upper left
    legend.position = c(0.08, 0.715), ## bottom right
    legend.text = element_text(size = legend_text_size*.6),
    legend.title = element_text(size = legend_text_size*.3 + legend_text_size),
    # legend.background = element_rect(fill = "transparent"),
    legend.key = element_rect(fill = "transparent"),
    # legend.key.width = unit(0, 'cm'),
    legend.key.size = unit(0.2, "cm"),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.2),  # Fills legend background & adds a border
    legend.box.background = element_rect(color = "black", linewidth = 0),  # Adds an outer border
    legend.margin = margin(1, 1, 3, 1),
    legend.spacing.y = unit(-0.0, "cm"),
    legend.title.align = 0.5
  )

res_cost_viral_load_p

ggsave("clinical_trials/res_cost_viral_load_p.pdf", res_cost_viral_load_p, w = width/2, h = height/4, units = "cm")
ggsave("clinical_trials/res_cost_viral_load_p.jpg", res_cost_viral_load_p, w = width/2, h = height/4, units = "cm")


#### Plot everything together                                               ####

library(cowplot)

sfig_a <- fit_func_9_p
sfig_b <- fit_func_1_p
sfig_c <- fit_func_01_p

sfig_d <- determ_9_p
sfig_e <- determ_1_p
sfig_f <- determ_01_p

sfig_g <- sims_9_plot
sfig_h <- sims_1_plot
sfig_i <- sims_01_plot

sfig_j <- res_cost_clear_dates_p
sfig_k <- res_cost_viral_load_p



## Creating left top column
col1 <- plot_grid(sfig_a, sfig_d, sfig_g,
                  ncol = 1,
                  nrow = 3,
                  # rel_heights = c(0.25, 0.25, 0.5),
                  axis = "tblr", 
                  align = "h")

## Creating middle top column
col2 <- plot_grid(sfig_b, sfig_e, sfig_h,
                  ncol = 1,
                  nrow = 3,
                  # rel_heights = c(0.25, 0.25, 0.5),
                  axis = "tblr", 
                  align = "h")

## Creating middle top column
col3 <- plot_grid(sfig_c, sfig_f, sfig_i,
                  ncol = 1,
                  nrow = 3,
                  # rel_heights = c(0.25, 0.25, 0.5),
                  axis = "tblr", 
                  align = "h")

## Binding top section
top <- plot_grid(col3, col2, col1,
                 ncol = 3,
                 nrow = 1,
                 # rel_heights = c(0.25, 0.25, 0.5),
                 axis = "tblr", 
                 align = "h")

## Creating bottom row
row1 <- plot_grid(sfig_j, sfig_k,
                  ncol = 2,
                  nrow = 1,
                  # rel_heights = c(0.25, 0.25, 0.5),
                  axis = "tblr", 
                  align = "h")

## All together
sfig <- plot_grid(top, row1,
                  ncol = 1,
                  nrow = 2,
                  rel_heights = c(1, 1/3),
                  axis = "tblr", 
                  align = "h")

save_plot(plot = sfig, 
          filename = paste0("sfig_res_cost_wo_letters.pdf"), 
          base_height = height, 
          base_width = width, 
          units = "cm", 
          bg = "transparent", 
          dpi = 1200)

# 
# t24_1 <- readRDS("~/PhD/dominant_drug_targets/polio-res-eco-evo_update/250810_range_of_cpop/old_data_again/trial_24_h_new.RDS")
# t24_2 <- readRDS("~/PhD/dominant_drug_targets/polio-res-eco-evo_update/250810_range_of_cpop/old_data_again/trial_24_h.RDS")
# t24_3 <- readRDS("~/PhD/dominant_drug_targets/polio-res-eco-evo_update/250810_range_of_cpop/trial_24_h.RDS")
# flat_24_1 <- flatten_list(t24_1)
# flat_24_2 <- flatten_list(t24_2)
# flat_24_3 <- flatten_list(t24_3)
# 
# flat_24_1 <- bind_rows(flat_24_1)
# flat_24_2 <- bind_rows(flat_24_2)
# flat_24_3 <- bind_rows(flat_24_3)
# 
# flat_24_all <- rbind(flat_24_1, flat_24_2, flat_24_3)
# 
# 
# write.csv(file = "all_24.csv", x = flat_24_all)
# 
# 
# t72_1 <- readRDS("~/PhD/dominant_drug_targets/polio-res-eco-evo_update/250810_range_of_cpop/old_data_again/trial_72_h_new.RDS")
# t72_2 <- readRDS("~/PhD/dominant_drug_targets/polio-res-eco-evo_update/250810_range_of_cpop/old_data_again/trial_72_h.RDS")
# t72_3 <- readRDS("~/PhD/dominant_drug_targets/polio-res-eco-evo_update/250810_range_of_cpop/trial_72_h.RDS")
# flat_72_1 <- flatten_list(t72_1)
# flat_72_2 <- flatten_list(t72_2)
# flat_72_3 <- flatten_list(t72_3)
# 
# flat_72_1 <- bind_rows(flat_72_1)
# flat_72_2 <- bind_rows(flat_72_2)
# flat_72_3 <- bind_rows(flat_72_3)
# 
# flat_72_all <- rbind(flat_72_1, flat_72_2, flat_72_3)
# 
# 
# write.csv(file = "all_72.csv", x = flat_72_all)
