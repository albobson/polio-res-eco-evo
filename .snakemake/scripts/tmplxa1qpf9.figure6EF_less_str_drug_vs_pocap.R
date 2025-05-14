
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
    output = list('runs/ddt_logistic_mu_2e-05_2025-03-24/res/fig6/E.rds', 'runs/ddt_logistic_mu_2e-05_2025-03-24/res/fig6/F.rds'),
    params = list('2025-03-24', c('logistic'), c(2e-05), "date" = '2025-03-24', "fit_func" = c('logistic'), "mu" = c(2e-05)),
    wildcards = list(),
    threads = 1,
    log = list(),
    resources = list('mem_mb', 'disk_mb', 'tmpdir', 'log_loc', 'mfree', 'cluster_time', 'disk_free', 'runtime', 'n_cores', "mem_mb" = '3G', "disk_mb" = '4G', "tmpdir" = '/tmp/12723961.1.feder-short.q', "log_loc" = './log_scr/', "mfree" = '3G', "cluster_time" = '0:00:10:00', "disk_free" = '4G', "runtime" = 10, "n_cores" = '1'),
    config = list(),
    rule = 'fig6EF',
    bench_iteration = as.numeric(NA),
    scriptdir = '/net/feder/vol1/home/alexrob/projects/dominant_drug_targets/pipeline/scr/panels',
    source = function(...){
        wd <- getwd()
        setwd(snakemake@scriptdir)
        source(...)
        setwd(wd)
    }
)


######## snakemake preamble end #########
########################### Figure 6 E and F ###################################

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



## Making the color map for plotting
color_in <- c("#c94d4d", "#688cc8")

## geom_text() is in terms of mm, where element_text() is in terms of pt. I want
## my geom_text and axis text to be the same. Here, I am defining a conversion
## which will keep the text sizes the same between text types.
geom_text_conv = 0.3528

## I also want my axis text to be 12 pt font, so, setting that here
axis_text_size = 12

## I would like for the legend to be slightly smaller than the rest of the text.
## I'm setting that size here
legend_text_size = 8



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

## Calculate the fold difference from pocapavir
sims_id <- sims_id %>%
  mutate(fold_strin = sims_id$wt_rep_abil/min(sims_id$wt_rep_abil))
  


plot6e <- ggplot(sims_id, aes(x = fold_strin, 
                              # group = wt_rep_abil,
                              # fill=factor(log(wt_rep_abil))
)) +
  # geom_boxplot(notch=FALSE, outlier.shape=NA,
  #              color = "black") +
  geom_jitter(  size = 1, 
                alpha = 0.3, 
                shape = 16,
                height  = 0, 
                aes(y = max_time, color = perc_res)) +
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
  
  geom_line(aes(y = clear_rolling_mean), color = "black", size = 1) +
  scale_color_continuous(name = expression(italic(f) * "(Res.)"),
                         low = color_in[2],
                         high = color_in[1],
                         limits = c(0, 1),
  ) +
  theme_light() +
  xlab("Fold Less Stringent") + ylab("Clearance Date (DPI)") +
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
plot6e

suppressWarnings(dir.create(paste0(filepath, "res/")))
suppressWarnings(dir.create(paste0(filepath, "res/fig6/")))

ggsave(paste0(filepath, "res/fig6/E.png"), plot6e, h = 22/3, w = 18/3, unit = "cm", bg = "transparent", dpi = 300)
ggsave(paste0(filepath, "res/fig6/E.pdf"), plot6e, h = 22/3, w = 18/3, unit = "cm", bg = "transparent", dpi = 300)
ggsave(paste0(filepath, "res/fig6/E.svg"), plot6e, h = 22/3, w = 18/3, unit = "cm", bg = "transparent", dpi = 300)
saveRDS(plot6e, file = paste0(filepath, "res/fig6/E.rds"))


#### F. Sum Total Viral Load All                                            ####
## Trying to create a rolling average

## Replacing anything negative in the variance with 1000. I tried replacing with
## 1, but it wouldn't let me do that for some reason.
sims_id$vl_mean_minus_sd_adjusted <- ifelse(sims_id$vl_mean_minus_sd <= 0, 1000, sims_id$vl_mean_minus_sd)

plot6f <- ggplot(sims_id, aes(x = fold_strin, 
                              # group = wt_rep_abil,
                              # fill=factor(log(wt_rep_abil))
)) +
  # geom_boxplot(notch=FALSE, outlier.shape=NA,
  #              color = "black") +
  geom_jitter(mapping = aes(y = tot_v, 
                            color = perc_res
  ), 
  size = 1, 
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
  
  geom_line(aes(y = vl_rolling_mean), color = "black", size = 1) +
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
  xlab("Fold Less Stringent") + ylab("Sum Total Viral Load") +
  scale_x_log10(
    # breaks=10^(-5:0),
    # limits = c((min(sims_id$wt_rep_abil)-min(sims_id$wt_rep_abil)*0.1), 1.5),
    # labels = sapply(-5:0,function(j){parse(text = sprintf("10^%d",j))})
    ) +
  scale_y_log10(breaks=10^(round(min(log10(sims_id$tot_v))):round(max(log10(sims_id$tot_v)))),
                limits = 10^c(round(min(log10(sims_id$tot_v))), round(max(log10(sims_id$tot_v)))),
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
    legend.position = "inside",
    # legend.position.inside = c(0.35, 0.94), ## Upper left
    legend.position.inside = c(0.89, 0.15), ## bottom right
    legend.text = element_text(size = legend_text_size*.6),
    legend.title = element_text(size = legend_text_size),
    legend.background = element_rect(fill = "transparent"),
    legend.key = element_rect(fill = "transparent"),
    # legend.key.width = unit(0, 'cm'),
    legend.key.size = unit(0.2, "cm")
  ) 

plot6f

ggsave(paste0(filepath, "res/fig6/F.png"), plot6f, h = 22/3, w = 18/3, unit = "cm", bg = "transparent", dpi = 300)
ggsave(paste0(filepath, "res/fig6/F.pdf"), plot6f, h = 22/3, w = 18/3, unit = "cm", bg = "transparent", dpi = 300)
ggsave(paste0(filepath, "res/fig6/F.svg"), plot6f, h = 22/3, w = 18/3, unit = "cm", bg = "transparent", dpi = 300)
saveRDS(plot6f, file = paste0(filepath, "res/fig6/F.rds"))



# #### D. Percent Resistant All ####
# plot6d <- ggplot(sims_id, aes(x = wt_rep_abil,
#                               # group = wt_rep_abil,
#                               # fill=factor(log(wt_rep_abil))
# )) +
#   # geom_boxplot(notch=FALSE, outlier.shape=NA,
#   #              color = "black") +
#   geom_jitter(mapping = aes(y = perc_res, 
#                             # color = perc_res
#   ),
#   size = 0.2, alpha = 0.3, height  = 0,
#   ) +
#   # geom_ribbon(sim_sds, 
#               # mapping = aes(ymin = res_mean_minus_sd, ymax = res_mean_plus_sd),
#               # alpha = 0.2) +
#   # geom_line(sims_stats_all, mapping = aes(y = res_mean), linewidth = 3) + ## Mean line
#   # scale_fill_grey(start = 0.3,
#   #                 end = 1) +
#   # scale_color_continuous(low = "#688cc8",
#   #                        high = "#c94d4d",
#   #                        limits = c(0, 1),
#   # ) +
#   geom_ribbon(
#     # data = sim_rolling_sds,
#     mapping = aes(ymin = res_mean_minus_sd, ymax = res_mean_plus_sd),
#     alpha = 0.2) +
#   
#   # Add line for the mean values (from summary stats)
#   # geom_line(data = sims_stats_rolling, mapping = aes(y = clear_mean), linewidth = 3) +
#   
#   geom_line(aes(y = res_rolling_mean), color = "black", size = 1.5) +
#   theme_light() +
#   xlab("P(WT Capsid Survival)") + ylab(expression(italic(f) * "(Resistant Genotype)")) +
#   scale_x_log10(breaks=10^(-5:0),
#                 limits = c((min(sims_id$wt_rep_abil)-min(sims_id$wt_rep_abil)*0.1), 1.5),
#                 labels = sapply(-5:0,function(j){parse(text = sprintf("10^%d",j))})) +
#   theme(text = element_text(size=18), 
#         axis.text = element_text(size=14),
#         panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
#         panel.grid.major.y = element_line(color = "lightgrey"),
#         panel.grid.minor.y = element_line(color = "lightgrey"),
#         # legend.title = element_blank(),
#         legend.position = "none",
#         # title = element_blank(),
#         # axis.text.x = element_blank(),
#         # axis.title.x=element_blank(),
#         axis.text.x = element_text(colour="black"),
#         axis.text.y = element_text(colour="black"),
#         strip.text = element_text(size = 22, color = "black"),
#         strip.background =element_rect(fill="white"),
#         # axis.title.y = element_blank(),
#         # axis.title.x = element_blank(),
#         ## Below is needed to fix some wonky ggbreak changes
#         axis.text.x.top = element_blank(),
#         axis.line.x.top = element_blank(),
#         axis.ticks.x.top = element_blank())
# 
# 
# 
# ggsave(plot6d, file = paste0(filepath, "res/fig6/D.png"), h = 8, w = 8, bg = "transparent")
# ggsave(plot6d, file = paste0(filepath, "res/fig6/D.jpg"), h = 8, w = 8, bg = "transparent")
# ggsave(plot6d, file = paste0(filepath, "res/fig6/D.svg"), h = 8, w = 8, bg = "transparent", device = 'svg')
# saveRDS(plot6d, paste0(filepath, "res/fig6/D.rds"))