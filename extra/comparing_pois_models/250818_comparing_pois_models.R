## Writing code to compare the original Poisson entry model with the model that
## iteratively infects cells


## Setting up a remote connection for easy scripting


#remoter::client(port=54321,addr="localhost",prompt="r-on-cluster",timer=FALSE)

## Read in original functions
source("../scr/polv_DDT_functions.R")

## Read in new infection function
source("iterative_inf_model.R")

## Old model
source("original_inf_model.R")

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

## Create the filepath where things will be saved
filepath <- paste0("../runs/ddt_", init_fit_func, "_mu_", init_mu, "_", init_date, "/")

axis_text_size <- 6
legend_text_size <- 6

## Read in fitness function and parameters
fit_func <- read.csv(paste0(filepath, "dat_gen/params/logistic_fitness_function.csv"))
optim_params <- read.csv(paste0(filepath, "dat_gen/params/optim_params.csv"))
clin_trial_params <- read.csv(paste0(filepath, "dat_gen/params/optim_clin_trial_params.csv"))
# clin_trial_params <- NULL
# clin_trial_params$optim_t_imm <- 9.39
# clin_trial_params$optim_c_pop <- 51087
# clin_trial_params$optim_imm_m <- -1.32
# clin_trial_params$optim_imm_sd <- 0.51

#### Run clinical trials using either function                              ####
## Number of trials to run
n_trials <- 20

n_72 <- 70*n_trials

n_24 <- 23*n_trials

#### Original function
cl <- makeCluster(n_cores - 1)
registerDoParallel(cl)

set.seed(1249128)
trial_72_h <- foreach(k = 1:n_72, .combine = "rbind", 
                      .packages = c('dplyr')) %dopar% {
                        stoch_polv_og(n = 300, 
                                   moi_wt_start = 1, 
                                   moi_mut_start = 0, 
                                   t_pocap = 9, 
                                   id = k,
                                   imm_delay = clin_trial_params$optim_t_imm,
                                   c_pop = clin_trial_params$optim_c_pop, 
                                   imm_m = clin_trial_params$optim_imm_m, 
                                   imm_sd = clin_trial_params$optim_imm_sd,
                                   fit_func_in = fit_func$prob_surv,
                                   v_prog = optim_params$optim_v_prog,
                                   p2pfu = optim_params$optim_p2pfu,
                                   seed_in = k*2
                        ) %>%
                          mutate(treatment = paste0("original")) %>%
                          mutate(time = time/3) ## Changing from replications to days
                      }


set.seed(12491284)
trial_24_h <- foreach(k = n_72:(n_72+n_24), .combine = "rbind") %dopar% {
  stoch_polv_og(n = 300, 
             moi_wt_start = 1, 
             moi_mut_start = 0, 
             t_pocap = 3, 
             id = k,
             imm_delay = clin_trial_params$optim_t_imm,
             c_pop = clin_trial_params$optim_c_pop, 
             imm_m = clin_trial_params$optim_imm_m, 
             imm_sd = clin_trial_params$optim_imm_sd,
             fit_func_in = fit_func$prob_surv,
             v_prog = optim_params$optim_v_prog,
             p2pfu = optim_params$optim_p2pfu,
             seed_in = k*2
  ) %>% 
    mutate(treatment = paste0("original")) %>%
    mutate(time = time/3) ## Changing from replications to days
}

stopCluster(cl)
stopImplicitCluster()

trial_og <- rbind(trial_72_h, trial_24_h)

write.csv(x = trial_og, file = "original_function.csv")


## Run the alternative model
cl <- makeCluster(n_cores - 1)
registerDoParallel(cl)

set.seed(1249128)
trial_72_h_itr <- foreach(k = 1:n_72, .combine = "rbind", 
                      .packages = c('dplyr')) %dopar% {
                        stoch_polv_iter(n = 300, 
                                   moi_wt_start = 1, 
                                   moi_mut_start = 0, 
                                   t_pocap = 9, 
                                   id = k,
                                   imm_delay = clin_trial_params$optim_t_imm,
                                   c_pop = clin_trial_params$optim_c_pop, 
                                   imm_m = clin_trial_params$optim_imm_m, 
                                   imm_sd = clin_trial_params$optim_imm_sd,
                                   fit_func_in = fit_func$prob_surv,
                                   v_prog = optim_params$optim_v_prog,
                                   p2pfu = optim_params$optim_p2pfu,
                                   seed_in = k*2
                        ) %>%
                          mutate(treatment = paste0("alternative")) %>%
                          mutate(time = time/3) ## Changing from replications to days
                      }


set.seed(12491284)
trial_24_h_itr <- foreach(k = n_72:(n_72+n_24), .combine = "rbind") %dopar% {
  stoch_polv_iter(n = 300, 
             moi_wt_start = 1, 
             moi_mut_start = 0, 
             t_pocap = 3, 
             id = k,
             imm_delay = clin_trial_params$optim_t_imm,
             c_pop = clin_trial_params$optim_c_pop, 
             imm_m = clin_trial_params$optim_imm_m, 
             imm_sd = clin_trial_params$optim_imm_sd,
             fit_func_in = fit_func$prob_surv,
             v_prog = optim_params$optim_v_prog,
             p2pfu = optim_params$optim_p2pfu,
             seed_in = k*2
  ) %>% 
    mutate(treatment = paste0("alternative")) %>%
    mutate(time = time/3) ## Changing from replications to days
}

stopCluster(cl)
stopImplicitCluster()

trial_itr <- rbind(trial_72_h_itr, trial_24_h_itr)

write.csv(x = trial_itr, file = "alternative_function.csv")

#### Clean the data                                                         ####
trial_og <- read.csv("original_function.csv")
trial_itr <- read.csv("alternative_function.csv")

all_trials <- rbind(trial_itr, trial_og)

## Clean to find just the last days 
simsclean <- all_trials %>% group_by(treatment, id) %>% 
  filter(type == "resistant") %>%
  mutate(pop_prop = as.numeric(pop_prop)) %>%
  na.omit() %>%
  mutate(res = sum(surv_pfu)/sum(as.numeric(tot_pfu)) > 0.5,
         perc_res = sum(surv_pfu)/sum(as.numeric(tot_pfu)),
         log_perc_res = log(sum(surv_pfu))/log(sum(as.numeric(tot_pfu))),
         avg_pop_prop = mean(as.numeric(pop_prop)),
         tot_vl = sum(as.numeric(tot_pfu)),
  ) %>%
  select(time, res, perc_res, tot_vl, treatment, id, time_to_pocap, avg_pop_prop,
         log_perc_res, c_pop) %>% 
  filter(time == max(time)) %>%
  mutate(res = "Simulated") %>% ## Just coloring by if something was simulated or not
  # mutate(res = ifelse(res == TRUE, "Resistant", "Susceptible"),
  #        res_mean = ifelse(avg_pop_prop > 0.5, "Resistant", "Susceptible"),
  #        log_res = ifelse(log_perc_res > 0.5, "Resistant", "Susceptible")) %>%
  unique() %>%
  mutate(time = ceiling(time))

## Alt cleaning step
sims_id <- all_trials %>%
  group_by(wt_rep_abil, id, c_pop, treatment) %>%
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
  group_by(treatment) %>%
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




#### Plot clearance time comparison                                         ####
## Creating a dynamic size for the dot plot
dotsize <- (max(simsclean$time)-min(simsclean$time))/max(simsclean$time)-0.1


early_late_mapping <- c("#c94d4d", "#688cc8")


## Plotting
pois_clear_p <- 
  ggplot(simsclean, aes(x = treatment, group = treatment
  )) +
  # annotate("rect", ymin = 15.5, ymax = 17.5, xmin = 0.5, xmax = 2.5, color = "#f2f2f2", alpha=0.2) +
  # annotate("rect", ymin = 18.5, ymax = 21.5, xmin = 0.5, xmax = 2.5, color = "#f2f2f2", alpha=0.2) +
  # annotate("rect", ymin = 22.5, ymax = 28.5, xmin = 0.5, xmax = 2.5, color = "#f2f2f2", alpha=0.2) +
  # annotate("rect", ymin = 29.5, ymax = 42.5, xmin = 0.5, xmax = 2.5, color = "#f2f2f2", alpha=0.2) +
  # geom_boxplot(outliers = FALSE) +
  # geom_dotplot(binaxis = "y",
  #              stackdir = "center",
  #              # method="histodot",
  #              # stackdir = "centerwhole",
  #              binpositions = "all",
  #              stackgroups=TRUE,
  #              dotsize = dotsize,
  #              binwidth = 1,
  #              stroke = NA,
  #              # key_glyph = "point",
  #              # mapping = aes(fill = resistant)
  # ) +
  geom_violin(aes(y = time)) +
  geom_point(data = sims_stats_all,
             mapping = aes(y = clear_mean),
             color = "black",
             size = 2,
             shape = 3
  ) +
  geom_errorbar(data = sim_sds,
                mapping = aes(ymin=clear_mean_minus_sd, ymax=clear_mean_plus_sd),
                width=.08,
                # position=position_dodge(.05)
  ) +
  scale_fill_manual(values = early_late_mapping) +
  scale_y_continuous(limits = c(0, 43.5),
                     breaks = c(seq(from = 0, to = 15, by = 5), 18, 22, 29, 43),
                     expand = c(0,0)) +
  labs(x = "Treatment", 
       y = "Clearance time (DPI)") +
  theme_classic() +
  # geom_hline(yintercept = as.numeric(early_placebo) - 0.5, linetype = 2, linewidth = 0.3) +
  theme(
    text = element_text(size= axis_text_size),
    axis.text = element_text(size= axis_text_size),
    panel.grid.major.x = element_blank(), 
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(color = "lightgrey", linewidth = 0.2),
    panel.grid.minor.y = element_blank(), #element_line(color = "lightgrey"),
    axis.text.x = element_text(colour="black"),
    axis.text.y = element_text(colour="black"),
    axis.title.y = element_text(),
    axis.title.x = element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    panel.background = element_rect(fill='transparent'), 
    plot.background = element_rect(fill='transparent', color=NA),
    ## Legend stuff
    # legend.position = "inside",
    # legend.position.inside = c(0.18, 0.06),
    # legend.text = element_text(size = legend_text_size),
    # legend.title = element_blank(),
    # legend.background = element_rect(fill = "transparent"),
    # legend.key = element_rect(fill = "transparent"),
    # legend.key.width = unit(0, 'cm'),
    # legend.key.size = unit(0.3, "cm")
  )  

# pois_clear_p

ggsave("pois_clear_p.pdf", pois_clear_p, w = 8/2, h = 8/2, units = "cm")
ggsave("pois_clear_p.jpg", pois_clear_p, w = 8/2, h = 8/2, units = "cm")


#### Plot viral load and resistance frequency comparison                    ####
sim_sds$vl_mean_minus_sd_adjusted <- ifelse(sim_sds$vl_mean_minus_sd <= 0, 1000, sim_sds$vl_mean_minus_sd)

pois_viral_load_p <- ggplot(sims_id, aes(x = treatment, 
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
             size = 2,
             shape = 3
  ) +
  geom_errorbar(data = sim_sds,
                aes(ymin=vl_mean_minus_sd_adjusted, ymax=vl_mean_plus_sd),
                width=.08,
                # position=position_dodge(.05)
  ) +
  scale_color_continuous(name = expression(f[Res]),
                         low = color_in[2],
                         high = color_in[1],
                         limits = c(0, 1),
  ) +
  theme_light() +
  # xlab("Fitness of a fully resistant capsid") + 
  ylab("Sum Total Viral Load") +
  # scale_x_log10(
  #   # breaks=10^(-5:0),
  #   # limits = c((min(sims_id$wt_rep_abil)-min(sims_id$wt_rep_abil)*0.1), 1.5),
  #   # labels = sapply(-5:0,function(j){parse(text = sprintf("10^%d",j))})
  # ) +
  scale_y_log10(breaks=10^(round(min(log10(sims_id$tot_v))):round(max(log10(sims_id$tot_v)))),
                # limits = 10^c(round(min(log10(sims_id$tot_v))), round(max(log10(sims_id$tot_v)))),
                labels = sapply(round(min(log10(sims_id$tot_v))):round(max(log10(sims_id$tot_v))),
                                function(j){parse(text = sprintf("10^%d",j))}),
                expand = expansion(mult = c(0, 0.05))) +
  theme(
    ## Text size
    text = element_text(size= axis_text_size), 
    axis.text = element_text(size= axis_text_size),
    ## Text color
    axis.text.x = element_text(colour="black"),
    axis.text.y = element_text(colour="black"),
    axis.title.y = element_text(),
    axis.title.x = element_blank(),
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
    # legend.position = c(0.845, 0.24), ## bottom right
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

# pois_viral_load_p

ggsave("pois_viral_load_p.pdf", pois_viral_load_p, w = 9/2, h = 8/2, units = "cm")
ggsave("pois_viral_load_p.jpg", pois_viral_load_p, w = 9/2, h = 8/2, units = "cm")

