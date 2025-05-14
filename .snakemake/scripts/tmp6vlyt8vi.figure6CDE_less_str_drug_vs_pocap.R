
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
    input = list('runs/ddt_logistic_mu_2e-05_2025-03-24/dat_gen/params/optim_params.csv', 'runs/ddt_logistic_mu_2e-05_2025-03-24/dat_gen/params/logistic_fitness_function.csv', 'runs/ddt_logistic_mu_2e-05_2025-03-24/dat_gen/sims/full_stringency_trials.csv'),
    output = list('runs/ddt_logistic_mu_2e-05_2025-03-24/res/fig6/C.rds', 'runs/ddt_logistic_mu_2e-05_2025-03-24/res/fig6/D.rds', 'runs/ddt_logistic_mu_2e-05_2025-03-24/res/fig6/E.rds'),
    params = list('2025-03-24', c('logistic'), c(2e-05), "date" = '2025-03-24', "fit_func" = c('logistic'), "mu" = c(2e-05)),
    wildcards = list(),
    threads = 1,
    log = list(),
    resources = list('mem_mb', 'disk_mb', 'tmpdir', 'log_loc', 'mfree', 'cluster_time', 'disk_free', 'runtime', 'n_cores', "mem_mb" = '3G', "disk_mb" = '4G', "tmpdir" = '/tmp/12679462.1.feder-short.q', "log_loc" = './log_scr/', "mfree" = '3G', "cluster_time" = '0:00:10:00', "disk_free" = '4G', "runtime" = 10, "n_cores" = '1'),
    config = list(),
    rule = 'fig6CDE',
    bench_iteration = as.numeric(NA),
    scriptdir = '/net/feder/vol1/home/alexrob/projects/dominant_drug_targets/pipeline/scr',
    source = function(...){
        wd <- getwd()
        setwd(snakemake@scriptdir)
        source(...)
        setwd(wd)
    }
)


######## snakemake preamble end #########
########################### Figure 6 C, D and E ################################

## Reason:

## This script will plot the average time to clearance, the average resistance
## frequency, and the average viral load of the pocapavir trial and the
## less stringent trial

#### Set Up                                                                 ####
## Set up where to save
## Date of the run
init_date <- as.character(snakemake@params[["date"]])

## Fitness function
init_fit_func <- as.character(snakemake@params[["fit_func"]])

## Mutation rate
init_mu <- as.numeric(snakemake@params[["mu"]])

## Create the filepath where things will be saved
filepath <- paste0("runs/ddt_", init_fit_func, "_mu_", init_mu, "_", init_date, "/")

## Read in fitness function and parameters
fit_func <- read.csv(paste0(filepath, "dat_gen/params/logistic_fitness_function.csv"))
optim_params <- read.csv(paste0(filepath, "dat_gen/params/optim_params.csv"))

## Source functions to run simulations
source("scr/polv_DDT_functions.R")  ## Main functions 
## Read in trial data
trial_sims <- read.csv(paste0(filepath, "dat_gen/sims/full_stringency_trials.csv"))

## Source functions to run simulations
source("scr/polv_DDT_functions.R")  ## Main functions 

## Create directory to save stuff (if not already created)
suppressWarnings(dir.create(paste0(filepath, "res/fig6")))

## Load the zoo package for rollapply()
library(zoo)

#### Clean data                                                             ####
## Group our data by id and drug stringency
## Converting to days
sims <- trial_sims %>%
  mutate(time = time/3,
         time_to_pocap = time_to_pocap/3) %>% ## Div by three for days
  filter(time >= time_to_pocap) ## removing any time pre pocap administration

## Finding the averages per each run (for box plots)
sims_id <- sims %>%
  group_by(wt_rep_abil, id, c_pop) %>%
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
  group_by(wt_rep_abil) %>%
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
window_size <- 80

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




#### Plotting                                                               ####

#### c. Clear time all ####
## Trying to create a rolling average
library(zoo)

# Calculate rolling average with a window size of 5
sims_id$rolling_avg_clear <- rollmean(sims_id$max_time, k = 100, fill = NA, align = "center")

plot6c <- ggplot(sims_id, aes(x = wt_rep_abil, 
                              # group = wt_rep_abil,
                              # fill=factor(log(wt_rep_abil))
)) +
  # geom_boxplot(notch=FALSE, outlier.shape=NA,
  #              color = "black") +
  geom_jitter(  size = 1.2, 
                alpha = 0.3, 
                shape = 16,
                height  = 0, aes(y = max_time, color = perc_res)) +
  # geom_ribbon(sim_sds, 
              # mapping = aes(ymin = clear_mean_minus_sd, ymax = clear_mean_plus_sd),
              # alpha = 0.2) +
  # geom_line(sims_stats_all, mapping = aes(y = clear_mean), linewidth = 3) + ## Mean line
  # Add ribbon for the rolling standard deviation
  geom_ribbon(
    # data = sim_rolling_sds,
              mapping = aes(ymin = clear_mean_minus_sd, ymax = clear_mean_plus_sd),
              alpha = 0.2) +
  
  # Add line for the mean values (from summary stats)
  # geom_line(data = sims_stats_rolling, mapping = aes(y = clear_mean), linewidth = 3) +
  
  geom_line(aes(y = clear_rolling_mean), color = "black", size = 1.5) +
  scale_color_continuous(name = expression(italic(f) * "(Res.)"),
                         low = "#688cc8",
                         high = "#c94d4d",
                         limits = c(0, 1),
  ) +
  theme_light() +
  xlab("P(WT Capsid Survival)") + ylab("Clearance Date (DPI)") +
  scale_x_log10(breaks=10^(-5:0),
                limits = c((min(sims_id$wt_rep_abil)-min(sims_id$wt_rep_abil)*0.1), 1.5),
                labels = sapply(-5:0,function(j){parse(text = sprintf("10^%d",j))})) +
  theme(text = element_text(size=18), 
        axis.text = element_text(size=14),
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "lightgrey"),
        panel.grid.minor.y = element_line(color = "lightgrey"),
        # legend.title = element_blank(),
        legend.position = "none",
        # title = element_blank(),
        # axis.text.x = element_blank(),
        axis.text.x = element_text(colour="black"),
        # axis.title.x=element_blank(),
        axis.text.y = element_text(colour="black"),
        strip.text = element_text(size = 22, color = "black"),
        strip.background =element_rect(fill="white"),
        # axis.title.y = element_blank(),
        # axis.title.x = element_blank(),
        ## Below is needed to fix some wonky ggbreak changes
        axis.text.x.top = element_blank(),
        axis.line.x.top = element_blank(),
        axis.ticks.x.top = element_blank())
plot6c

ggsave("res/fig6/C.png", plot6c, h = 8, w = 8, bg = "transparent")
ggsave("res/fig6/C.jpg", plot6c, h = 8, w = 8, bg = "transparent")
ggsave("res/fig6/C.svg", plot6c, h = 8, w = 8, bg = "transparent", device = 'svg')
saveRDS(plot6c, "res/fig6/C.rds")


#### D. Percent Resistant All ####
plot6d <- ggplot(sims_id, aes(x = wt_rep_abil,
                              # group = wt_rep_abil,
                              # fill=factor(log(wt_rep_abil))
)) +
  # geom_boxplot(notch=FALSE, outlier.shape=NA,
  #              color = "black") +
  geom_jitter(mapping = aes(y = perc_res, 
                            # color = perc_res
  ),
  size = 0.2, alpha = 0.3, height  = 0,
  ) +
  # geom_ribbon(sim_sds, 
              # mapping = aes(ymin = res_mean_minus_sd, ymax = res_mean_plus_sd),
              # alpha = 0.2) +
  # geom_line(sims_stats_all, mapping = aes(y = res_mean), linewidth = 3) + ## Mean line
  # scale_fill_grey(start = 0.3,
  #                 end = 1) +
  # scale_color_continuous(low = "#688cc8",
  #                        high = "#c94d4d",
  #                        limits = c(0, 1),
  # ) +
  geom_ribbon(
    # data = sim_rolling_sds,
    mapping = aes(ymin = res_mean_minus_sd, ymax = res_mean_plus_sd),
    alpha = 0.2) +
  
  # Add line for the mean values (from summary stats)
  # geom_line(data = sims_stats_rolling, mapping = aes(y = clear_mean), linewidth = 3) +
  
  geom_line(aes(y = res_rolling_mean), color = "black", size = 1.5) +
  theme_light() +
  xlab("P(WT Capsid Survival)") + ylab(expression(italic(f) * "(Resistant Genotype)")) +
  scale_x_log10(breaks=10^(-5:0),
                limits = c((min(sims_id$wt_rep_abil)-min(sims_id$wt_rep_abil)*0.1), 1.5),
                labels = sapply(-5:0,function(j){parse(text = sprintf("10^%d",j))})) +
  theme(text = element_text(size=18), 
        axis.text = element_text(size=14),
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "lightgrey"),
        panel.grid.minor.y = element_line(color = "lightgrey"),
        # legend.title = element_blank(),
        legend.position = "none",
        # title = element_blank(),
        # axis.text.x = element_blank(),
        # axis.title.x=element_blank(),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        strip.text = element_text(size = 22, color = "black"),
        strip.background =element_rect(fill="white"),
        # axis.title.y = element_blank(),
        # axis.title.x = element_blank(),
        ## Below is needed to fix some wonky ggbreak changes
        axis.text.x.top = element_blank(),
        axis.line.x.top = element_blank(),
        axis.ticks.x.top = element_blank())



ggsave(paste0(filepath, "res/fig6/D.png"), plot6d, h = 8, w = 8, bg = "transparent")
ggsave(paste0(filepath, "res/fig6/D.jpg"), plot6d, h = 8, w = 8, bg = "transparent")
ggsave(paste0(filepath, "res/fig6/D.svg"), plot6d, h = 8, w = 8, bg = "transparent", device = 'svg')
saveRDS(plot6d, paste0(filepath, "res/fig6/D.rds"))


#### E. Sum Total Viral Load All ####
## Trying to create a rolling average

## Replacing anything negative in the variance with 1000. I tried replacing with
## 1, but it wouldn't let me do that for some reason.
sims_id$vl_mean_minus_sd_adjusted <- ifelse(sims_id$vl_mean_minus_sd <= 0, 1000, sims_id$vl_mean_minus_sd)

plot6e <- ggplot(sims_id, aes(x = wt_rep_abil, 
                              # group = wt_rep_abil,
                              # fill=factor(log(wt_rep_abil))
)) +
  # geom_boxplot(notch=FALSE, outlier.shape=NA,
  #              color = "black") +
  geom_jitter(mapping = aes(y = tot_v, 
                            color = perc_res
  ), 
  size = 1.2, 
  alpha = 0.3, 
  # shape = 16,
  height  = 0) +
  # geom_ribbon(sim_sds,
              # mapping = aes(ymin = vl_mean_minus_sd, ymax = vl_mean_plus_sd),
              # alpha = 0.2) +
  # geom_line(sims_stats_all, mapping = aes(y = vl_mean), linewidth = 3) + ## Mean line
  # geom_smooth(aes(y = tot_v)) +
  ## Rolling average
  geom_ribbon(
    # data = sim_rolling_sds,
    mapping = aes(ymin = vl_mean_minus_sd_adjusted, ymax = vl_mean_plus_sd),
    alpha = 0.2) +
  
  # Add line for the mean values (from summary stats)
  # geom_line(data = sims_stats_rolling, mapping = aes(y = clear_mean), linewidth = 3) +
  
  geom_line(aes(y = vl_rolling_mean), color = "black", size = 1.5) +
  # geom_line(aes(y = rolling_avg_vl), color = "black", size = 1.5) +
  # Add spline fit using GAM (Generalized Additive Model)
  # geom_smooth(mapping = aes(y = tot_v), 
  #             method = "gam", 
  #             formula = y ~ s(x),  # Use a smooth function for the data
  #             se = FALSE,  # Optionally, set to TRUE if you want a confidence interval
  #             color = "black", 
  #             size = 1.5) +
  scale_color_continuous(name = expression(italic(f) * "(Res.)"),
                         low = "#688cc8",
                         high = "#c94d4d",
                         limits = c(0, 1),
  ) +
  theme_light() +
  xlab("P(WT Capsid Survival)") + ylab("Sum Total Viral Load") +
  scale_x_log10(breaks=10^(-5:0),
                limits = c((min(sims_id$wt_rep_abil)-min(sims_id$wt_rep_abil)*0.1), 1.5),
                labels = sapply(-5:0,function(j){parse(text = sprintf("10^%d",j))})) +
  scale_y_log10(breaks=10^(round(min(log10(sims_id$tot_v))):round(max(log10(sims_id$tot_v)))),
                limits = 10^c(round(min(log10(sims_id$tot_v))), round(max(log10(sims_id$tot_v)))),
                labels = sapply(round(min(log10(sims_id$tot_v))):round(max(log10(sims_id$tot_v))),function(j){parse(text = sprintf("10^%d",j))}),
                expand = expansion(mult = c(0, 0.05))) +
  theme(text = element_text(size=18), 
        axis.text = element_text(size=14),
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "lightgrey"),
        panel.grid.minor.y = element_line(color = "lightgrey"),
        # legend.title = element_blank(),
        # legend.position = "none",
        # title = element_blank(),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        strip.text = element_text(size = 22, color = "black"),
        strip.background =element_rect(fill="white"),
        # axis.title.y = element_blank(),
        # axis.title.x = element_blank(),
        ## Below is needed to fix some wonky ggbreak changes
        axis.text.x.top = element_blank(),
        axis.line.x.top = element_blank(),
        axis.ticks.x.top = element_blank())


ggsave(paste0(filepath, "res/fig6/E.png"), plot6e, h = 8, w = 8, bg = "transparent")
ggsave(paste0(filepath, "res/fig6/E.jpg"), plot6e, h = 8, w = 8, bg = "transparent")
ggsave(paste0(filepath, "res/fig6/E.svg"), plot6e, h = 8, w = 8, bg = "transparent", device = 'svg')
saveRDS(plot6e, paste0(filepath, "res/fig6/E.rds"))
