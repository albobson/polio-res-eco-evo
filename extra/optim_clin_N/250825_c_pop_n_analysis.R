### An analysis to run the clinical trials from the c_pop immune parameter
### fitting analysis.


## Read in files
init_c_pop <- c(30000/2, 30000*2)

## Date
init_date <- "2025-09-02"

print(paste0("init_c_pop: ", init_c_pop))

## Create the filepath where things will be saved
filepath <- paste0("init_c_pop_imm_train_results")

## Read in fitness function and parameters
fit_func <- read.csv(paste0("../runs/ddt_logistic_mu_2e-05_", paste0(init_date), "/dat_gen/params/logistic_fitness_function.csv"))
optim_params <- read.csv(paste0("../runs/ddt_logistic_mu_2e-05_", paste0(init_date), "/dat_gen/params/optim_params.csv"))

## Read in the fit data
clin_trial_params1 <- read.csv(paste0(filepath, "_params_c_pop_", init_c_pop[1], ".csv"))
clin_trial_params2 <- read.csv(paste0(filepath, "_params_c_pop_", init_c_pop[2], ".csv"))

## Read in Collett clinical trial data for comparison
coldat <- read.csv("../dat/collett_trial.csv")

## Source functions to run simulations
source("../scr/polv_DDT_functions.R")  ## Main functions


## Read in axis text size
axis_text_size <- 6

## Legend text size
legend_text_size <- 6

## This figure's size
width <-  5
height <- 3

## Number of cores to run the model
n_cores <- 3*2

#### Run the clinical trials using the fit data                             ####
print(paste0("First set of parameters"))

cl <- makeCluster(n_cores)
registerDoParallel(cl)

set.seed(1249128)
trial_72_h <- foreach(k = 1:70, .combine = "rbind", 
                      .packages = c('dplyr')) %dopar% {
                        stoch_polv(n = 300, 
                                   moi_wt_start = 1, 
                                   moi_mut_start = 0, 
                                   t_pocap = 9, 
                                   id = k,
                                   imm_delay = clin_trial_params1$optim_t_imm,
                                   c_pop = clin_trial_params1$optim_c_pop, 
                                   imm_m = clin_trial_params1$optim_imm_m, 
                                   imm_sd = clin_trial_params1$optim_imm_sd,
                                   fit_func_in = fit_func$prob_surv,
                                   v_prog = optim_params$optim_v_prog,
                                   p2pfu = optim_params$optim_p2pfu,
                                   seed_in = k*2
                        ) %>%
                          mutate(treatment = paste0("Simulated Pocapavir")) %>%
                          mutate(time = time/3) ## Changing from replications to days
                      }


set.seed(12491284)
trial_24_h <- foreach(k = 71:93, .combine = "rbind") %dopar% {
  stoch_polv(n = 300, 
             moi_wt_start = 1, 
             moi_mut_start = 0, 
             t_pocap = 3, 
             id = k,
             imm_delay = clin_trial_params1$optim_t_imm,
             c_pop = clin_trial_params1$optim_c_pop, #usually 8*10^3
             imm_m = clin_trial_params1$optim_imm_m, # Usually -1.6
             imm_sd = clin_trial_params1$optim_imm_sd, # Usually 0.5
             fit_func_in = fit_func$prob_surv,
             v_prog = optim_params$optim_v_prog,
             p2pfu = optim_params$optim_p2pfu,
             seed_in = k*2
  ) %>% 
    mutate(treatment = paste0("Simulated Pocapavir")) %>%
    mutate(time = time/3) ## Changing from replications to days
}

stopCluster(cl)
stopImplicitCluster()

trial1 <- rbind(trial_72_h, trial_24_h)
write.csv(x = trial1, file = "250825_first_clin_trial.csv")


## Second clinical trial
print(paste0("Second set of parameters"))

cl <- makeCluster(n_cores)
registerDoParallel(cl)

set.seed(1249128)
trial_72_h <- foreach(k = 1:70, .combine = "rbind", 
                      .packages = c('dplyr')) %dopar% {
                        stoch_polv(n = 300, 
                                   moi_wt_start = 1, 
                                   moi_mut_start = 0, 
                                   t_pocap = 9, 
                                   id = k,
                                   imm_delay = clin_trial_params2$optim_t_imm,
                                   c_pop = clin_trial_params2$optim_c_pop, 
                                   imm_m = clin_trial_params2$optim_imm_m, 
                                   imm_sd = clin_trial_params2$optim_imm_sd,
                                   fit_func_in = fit_func$prob_surv,
                                   v_prog = optim_params$optim_v_prog,
                                   p2pfu = optim_params$optim_p2pfu,
                                   seed_in = k*2
                        ) %>%
                          mutate(treatment = paste0("Simulated Pocapavir")) %>%
                          mutate(time = time/3) ## Changing from replications to days
                      }


set.seed(12491284)
trial_24_h <- foreach(k = 71:93, .combine = "rbind") %dopar% {
  stoch_polv(n = 300, 
             moi_wt_start = 1, 
             moi_mut_start = 0, 
             t_pocap = 3, 
             id = k,
             imm_delay = clin_trial_params2$optim_t_imm,
             c_pop = clin_trial_params2$optim_c_pop, #usually 8*10^3
             imm_m = clin_trial_params2$optim_imm_m, # Usually -1.6
             imm_sd = clin_trial_params2$optim_imm_sd, # Usually 0.5
             fit_func_in = fit_func$prob_surv,
             v_prog = optim_params$optim_v_prog,
             p2pfu = optim_params$optim_p2pfu,
             seed_in = k*2
  ) %>% 
    mutate(treatment = paste0("Simulated Pocapavir")) %>%
    mutate(time = time/3) ## Changing from replications to days
}

stopCluster(cl)
stopImplicitCluster()

trial2 <- rbind(trial_72_h, trial_24_h)
write.csv(x = trial2, file = "250825_second_clin_trial.csv")

#### Read in the original clinical trial                                    ####
og_trial <- read.csv(paste0("../runs/ddt_logistic_mu_2e-05_", init_date,"/dat_gen/sims/sim_clinical_trial.csv"))

# og_trial <- og_trial[,-1]

## Other trials
trial1 <- read.csv("250825_first_clin_trial.csv")
trial2 <- read.csv("250825_second_clin_trial.csv")

#### Plot the clinical trials                                               ####
## Clean the data
trial_full <- rbind(trial1, trial2, og_trial)

simsclean <- trial_full %>% group_by(treatment, id, c_pop) %>% 
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

## Fix to set the unsampled times
simsclean$time <- ifelse(simsclean$time %in%  16:17, 18,
                         ifelse(simsclean$time %in% 19:21, 22,
                                ifelse(simsclean$time %in% 23:28, 29,
                                       ifelse(simsclean$time %in% 30:42, 43,
                                              ifelse(simsclean$time > 43, 44,
                                                     simsclean$time)))))

## Early placebo 
coldat_p <- coldat %>% group_by(treatment) %>% 
  mutate(res = as.character(resistant), 
         clearange = as.numeric(clearange),
         source = "Clinical Trial") %>% 
  select(clearange, res, treatment, res) %>% 
  na.omit() %>% 
  # filter(treatment == "Pocapavir (n = 93)") %>%
  mutate(time = clearange, 
         treatment = ifelse(treatment == "Pocapavir (n = 93)", "Pocapavir", "Placebo"),
         # treatment = paste0(treatment, " ", res) ## if you want treatment to be stratified by resistance class
  ) %>%
  select(time, res, treatment) %>%
  mutate(shape = "circle") ## Adding a shape value for distinguising our own data

## Finding the earliest placebo clearer
early_placebo <- coldat_p %>%
  filter(treatment == "Placebo") %>%
  summarize(mint = min(time))

early_placebo <- early_placebo$mint[1]

coldat_p$earlylate <- ifelse(coldat_p$time >= early_placebo, "Late", "Early")

## Adding if they are early or late clearers
simsclean$earlylate <- ifelse(simsclean$time >= early_placebo, "Late", "Early")


## Creating a dynamic size for the dot plot
dotsize <- (max(simsclean$time)-min(simsclean$time))/max(simsclean$time)-0.1


## Plotting
tend_plot <- 
  ggplot(simsclean, aes(x = factor(round(c_pop, 0)), y = time,
  )) +
  annotate("rect", ymin = 15.5, ymax = 17.5, xmin = 0.5, xmax = 3.5, colour = "#f2f2f2", alpha=0.2) +
  annotate("rect", ymin = 18.5, ymax = 21.5, xmin = 0.5, xmax = 3.5, colour = "#f2f2f2", alpha=0.2) +
  annotate("rect", ymin = 22.5, ymax = 28.5, xmin = 0.5, xmax = 3.5, colour = "#f2f2f2", alpha=0.2) +
  annotate("rect", ymin = 29.5, ymax = 42.5, xmin = 0.5, xmax = 3.5, colour = "#f2f2f2", alpha=0.2) +
  # geom_boxplot(outliers = FALSE) +
  geom_dotplot(binaxis = "y", 
               stackdir = "center",
               # method="histodot",
               # stackdir = "centerwhole",
               binpositions = "all",
               stackgroups=TRUE,
               dotsize = dotsize,
               binwidth = 1,
               stroke = NA,
               # key_glyph = "point",
               # mapping = aes(fill = resistant)
  ) +
  # scale_fill_manual(values = early_late_mapping) +
  scale_y_continuous(limits = c(0, 43.5),
                     breaks = c(seq(from = 0, to = 15, by = 5), 18, 22, 29, 43),
                     expand = c(0,0)) +
  labs(x = "Treatment", 
       y = "Clearance time (DPI)") +
  theme_classic() +
  geom_hline(yintercept = early_placebo - 0.5, linetype = 2, linewidth = 0.3) +
  theme(text = element_text(size= axis_text_size), 
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

tend_plot


ggsave(tend_plot, file = paste0("clin_trial_comp.png"), h = height, w = width, units = "cm", bg = "transparent", dpi = 1200)
ggsave(tend_plot, file = paste0("clin_trial_comp.pdf"), h = height, w = width, units = "cm", bg = "transparent", dpi = 1200)

## Alternative bar plot of late vs early clearer
sims_sum <- simsclean %>%
  group_by(c_pop) %>%
  summarize(n_early = sum(earlylate == "Early"), n_late = sum(earlylate == "Late")) %>%
  mutate(c_pop = round(c_pop, 0))

sims_sum <- sims_sum %>%
  pivot_longer(
    cols = c(n_early, n_late), # The columns you want to pivot
    names_to = "earlylate",      # The new column for the original column names
    values_to = "count"       # The new column for the values
  )

sims_sum$perc <- sims_sum$count/93


bar_plot <- 
  ggplot(sims_sum, aes(x = factor(c_pop), y = perc, fill = factor(earlylate, levels = c("n_late", "n_early")))) +
  # geom_boxplot(outliers = FALSE) +
  geom_bar(position = "stack", stat = "identity") +
  # scale_fill_manual(values = early_late_mapping) +
  # scale_y_continuous(limits = c(0, 43.5),
  #                    breaks = c(seq(from = 0, to = 15, by = 5), 18, 22, 29, 43),
  #                    expand = c(0,0)) +
  labs(x = "Cell population size", 
       y = "Proportion") +
  theme_classic() +
  scale_fill_manual(
    name = "Clearers", # Sets the title of the legend
    values = c("n_late" = "black", # Assigns a specific color to a value
               "n_early" = "lightgrey"),
    labels = c("n_early" = "Early", # Sets the display text for each value
               "n_late" = "Late")
  ) +
  # geom_hline(yintercept = early_placebo - 0.5, linetype = 2, linewidth = 0.3) +
  theme(text = element_text(size= axis_text_size), 
        axis.text = element_text(size= axis_text_size),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "lightgrey", linewidth = 0.2),
        panel.grid.minor.y = element_blank(), #element_line(color = "lightgrey"),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        axis.title.y = element_text(),
        axis.title.x = element_text(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.background = element_rect(fill='transparent'), 
        plot.background = element_rect(fill='transparent', color=NA),
        ## Legend stuff
        # legend.position = "inside",
        # legend.position.inside = c(0.18, 0.06),
        legend.text = element_text(size = legend_text_size),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent"),
        legend.key.width = unit(0, 'cm'),
        legend.key.size = unit(0.3, "cm")
  )  

bar_plot


ggsave(bar_plot, file = paste0("clin_trial_barplot.png"), h = height, w = width, units = "cm", bg = "transparent", dpi = 1200)
ggsave(bar_plot, file = paste0("clin_trial_barplot.pdf"), h = height, w = width, units = "cm", bg = "transparent", dpi = 1200)
