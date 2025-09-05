#### Testing the range of c_pop that gives us a similar outcome as to what we
#### saw in the clinical trial
print(paste0("Started"))

#### Set Up                                                                 ####
## Set up where to save
## Date of the run to pull parameters from
init_date <- "2025-09-02"

## Today's date
today_date <- "2025-09-02"

## Fitness function
init_fit_func <- "logistic"

## Mutation rate
init_mu <- 2e-5

## Plot dims
width <- 18
height <- 17

legend_text_size <- 6
axis_text_size <- 6

geom_text_conv = 0.3528

## Cores to use
cores_to_use <- 10*2


## Collett resistance frequency data Data
col_res_72 <- 25/(45+25)
col_res_24 <- 15/(15+8)


## Create the filepath where things will be saved
filepath <- paste0("runs/ddt_", init_fit_func, "_mu_", init_mu, "_", init_date, "/")


## Read in fitness function and parameters
fit_func <- read.csv(paste0(filepath, "dat_gen/params/logistic_fitness_function.csv"))
optim_params <- read.csv(paste0(filepath, "dat_gen/params/optim_params.csv"))


## Source functions to run simulations
source("scr/polv_DDT_functions.R")  ## Main functions 

## Read in the non-sharing function
# source("250810_range_of_cpop/polv_non_sharing.R")



## Range of c_pop to run across
c_pop_init <- c(round(10^seq(from = 1, to = 6, by = 0.25), 0))

## Range of v_tot to run over
v_init <- c(1, 2, 3, 4, 5, 6, 7, 10, 15, 100)

## Number of stochastic runs
n_72 <- 150
n_24 <- 150

print(paste0("Everything loaded"))
#### Running simulations                                                    ####
#### Run our simulations                                                    ####
# print(paste0("Cores to use: ", cores_to_use))
# 
# cl <- makeCluster(cores_to_use - 1)
# registerDoParallel(cl)
# 
# print(paste0("running 72 hour sims"))
# set.seed(1249128)
# trial_72_h <- foreach(k = 1:n_72,
#                       # .combine = "rbind",
#                       .packages = c('dplyr'),
#                       .errorhandling = "pass"
#                       ) %:%
#                 foreach(gamma = c_pop_init,
#                         # .combine = "rbind",
#                         .packages = c('dplyr'),
#                         .errorhandling = "pass"
#                         ) %:%
#                   foreach(v_in = v_init,
#                           # .combine = "rbind",
#                           .packages = c('dplyr'),
#                           .errorhandling = "pass"
#                           ) %dopar% {
#                             stoch_polv(n = (8*3),
#                                        moi_wt_start = v_in/gamma,
#                                        moi_mut_start = 0,
#                                        t_pocap = 9,
#                                        id = k,
#                                        c_pop = gamma,
#                                        imm_delay = 100,
#                                        imm_m = 1,
#                                        imm_sd = 1,
#                                        bg_mutat = 1/(2*10^-5),
#                                        fit_func_in = fit_func$prob_surv,
#                                        v_prog = optim_params$optim_v_prog,
#                                        p2pfu = 1,
#                                        seed_in = k*2
#                         ) %>%
#                           mutate(treatment = paste0("Simulated Pocapavir")) %>%
#                           mutate(time = time/3) ## Changing from replications to days
#                           }
# 
# 
# saveRDS(trial_72_h, paste0("250810_range_of_cpop/", today_date, "_trial_72_h.RDS"))
# 
# print(paste0("running 24 hour sims"))
# 
# set.seed(12491284)
# trial_24_h <- foreach(k = (n_72+1):(n_24+n_72),
#                       .packages = c('dplyr'),
#                       # .combine = "rbind",
#                       .errorhandling = "pass"
#                       ) %:%
#   foreach(gamma = c_pop_init,
#           # .combine = "rbind",
#           .packages = c('dplyr'),
#           .errorhandling = "pass"
#           ) %:%
#   foreach(v_in = v_init,
#           # .combine = "rbind",
#           .packages = c('dplyr'),
#           .errorhandling = "pass"
#           ) %dopar% {
#   stoch_polv(n = (8*3),
#              moi_wt_start = v_in/gamma,
#              moi_mut_start = 0,
#              t_pocap = 3,
#              id = k,
#              c_pop = gamma,
#              imm_delay = 100,
#              imm_m = 1,
#              imm_sd = 1,
#              bg_mutat = 1/(2*10^-5),
#              fit_func_in = fit_func$prob_surv,
#              v_prog = optim_params$optim_v_prog,
#              p2pfu = 1,
#              seed_in = k*2
#   ) %>%
#     mutate(treatment = paste0("Simulated Pocapavir")) %>%
#     mutate(time = time/3) ## Changing from replications to days
#           }
# 
# saveRDS(trial_24_h, paste0("250810_range_of_cpop/", today_date, "_trial_24_h.RDS"))
# 
# stopCluster(cl)
# stopImplicitCluster()
# 
# 
# #### Run the standard model                                                 ####
# cl <- makeCluster(cores_to_use - 1)
# registerDoParallel(cl)
# 
# print(paste0("running 72 hour sims -- Standard"))
# set.seed(1249128)
# trial_72_h_s <- foreach(k = 1:n_72,
#                       # .combine = "rbind",
#                       .packages = c('dplyr'),
#                       .errorhandling = "pass"
# ) %:%
#   foreach(gamma = c_pop_init,
#           # .combine = "rbind",
#           .packages = c('dplyr'),
#           .errorhandling = "pass"
#   ) %:%
#   foreach(v_in = v_init,
#           # .combine = "rbind",
#           .packages = c('dplyr'),
#           .errorhandling = "pass"
#   ) %dopar% {
#     stoch_polv_no_interaction(n = (8*3),
#                               moi_wt_start = v_in/gamma,
#                               moi_mut_start = 0,
#                               t_pocap = 9,
#                               id = k,
#                               c_pop = gamma,
#                               imm_delay = 100,
#                               imm_m = 1,
#                               imm_sd = 1,
#                               bg_mutat = 1/(2*10^-5),
#                               fit_func_in = fit_func$prob_surv,
#                               v_prog = optim_params$optim_v_prog,
#                               p2pfu = 1,
#                               seed_in = k*2
#     ) %>%
#       mutate(treatment = paste0("Simulated Pocapavir")) %>%
#       mutate(time = time/3) ## Changing from replications to days
#   }
# 
# ## Save
# saveRDS(trial_72_h_s, paste0("250810_range_of_cpop/", today_date, "_trial_72_h_s.RDS"))
# 
# print(paste0("running 24 hour sims -- Standard"))
# 
# set.seed(12491284)
# trial_24_h_s <- foreach(k = (n_72+1):(n_24+n_72),
#                       .packages = c('dplyr'),
#                       # .combine = "rbind",
#                       .errorhandling = "pass"
# ) %:%
#   foreach(gamma = c_pop_init,
#           # .combine = "rbind",
#           .packages = c('dplyr'),
#           .errorhandling = "pass"
#   ) %:%
#   foreach(v_in = v_init,
#           # .combine = "rbind",
#           .packages = c('dplyr'),
#           .errorhandling = "pass"
#   ) %dopar% {
#     stoch_polv_no_interaction(n = (8*3),
#                               moi_wt_start = v_in/gamma,
#                               moi_mut_start = 0,
#                               t_pocap = 3,
#                               id = k,
#                               c_pop = gamma,
#                               imm_delay = 100,
#                               imm_m = 1,
#                               imm_sd = 1,
#                               bg_mutat = 1/(2*10^-5),
#                               fit_func_in = fit_func$prob_surv,
#                               v_prog = optim_params$optim_v_prog,
#                               p2pfu = 1,
#                               seed_in = k*2
#     ) %>%
#       mutate(treatment = paste0("Simulated Pocapavir")) %>%
#       mutate(time = time/3) ## Changing from replications to days
#   }
# 
# saveRDS(trial_24_h_s, paste0("250810_range_of_cpop/", today_date, "_trial_24_h_s.RDS"))
# 
# stopCluster(cl)
# stopImplicitCluster()


## Create a function to flatten all of these into one dataframe
# flatten_list <- function(x, path = integer(), validate_fn = NULL) {
#   # helper to pretty-print path like [[1]][[3]][[2]]
#   path_to_label <- function(p) {
#     if (length(p) == 0) return("<root>")
#     paste0(paste0("[[", p, "]]"), collapse = "")
#   }
#   
#   # If it's a data.frame, optionally validate it and return as single-item list
#   if (is.data.frame(x)) {
#     if (!is.null(validate_fn) && !isTRUE(validate_fn(x))) {
#       stop(sprintf("Validation failed for data.frame at %s", path_to_label(path)))
#     }
#     return(list(x))
#   }
#   
#   # If it's a list, recurse element by element
#   if (is.list(x)) {
#     out <- list()
#     for (i in seq_along(x)) {
#       current_path <- c(path, i)
#       item <- x[[i]]
#       # recurse (the recursive call will stop with coordinates if needed)
#       out <- c(out, flatten_list(item, current_path, validate_fn))
#     }
#     return(out)
#   }
#   return(list())
# }
# 
# ## reread back in
# trial_24_h <- readRDS(paste0("250810_range_of_cpop/", today_date, "_trial_24_h.RDS"))
# trial_72_h <- readRDS(paste0("250810_range_of_cpop/", today_date, "_trial_72_h.RDS"))
# 
# ## Flatten the lists
# flat_72 <- flatten_list(trial_72_h)
# 
# ## bind rows
# df_72_h <- bind_rows(flat_72)
# 
# ## Flatten the lists
# flat_24 <- flatten_list(trial_24_h)
# 
# ## bind rows
# df_24_h <- bind_rows(flat_24)
# 
# trial <- rbind(df_72_h, df_24_h)
# 
# 
# ## save the data 
# write.csv(trial, paste0("250810_range_of_cpop/", today_date, "_c_pop_trial_data_interaction.csv"))
# 
# 
# ## Same thing for the standard model
# ## Read in the data
# trial_72_h_s <- readRDS(paste0("250810_range_of_cpop/", today_date, "_trial_72_h_s.RDS"))
# trial_24_h_s <- readRDS(paste0("250810_range_of_cpop/", today_date, "_trial_24_h_s.RDS"))
# 
# 
# ## Flatten the list
# flat_72_s <- flatten_list(trial_72_h_s)
# 
# ## bind rows
# df_72_h_s <- bind_rows(flat_72_s)
# 
# ## Flatten the list
# flat_24_s <- flatten_list(trial_24_h_s)
# 
# ## bind rows
# df_24_h_s <- bind_rows(flat_24_s)
# 
# standard_model <- rbind(df_72_h_s, df_24_h_s)
# 
# ## save the data 
# write.csv(standard_model, paste0("250810_range_of_cpop/", today_date, "_c_pop_trial_data_standard.csv"))
# 
# # write.csv(standard_model, "250810_range_of_cpop/c_pop_trial_data_first_run_standard.csv")
# print(paste0("Saved -- Standard"))
# 
# 
##### Cleaning the data                                                      ####
trial <- read.csv(paste0("250810_range_of_cpop/", today_date, "_c_pop_trial_data_interaction.csv"))

print(paste0("Cleaning"))

## Our model
df_24_h <- trial %>%
  filter(time_to_pocap == 3,
         type == "resistant")

df_72_h <- trial %>%
  filter(time_to_pocap == 9,
         type == "resistant")

trial_clean_24 <- df_24_h %>%
  dplyr::select(time, init_wt, c_pop, time_to_pocap, pop_prop, id) %>%
  group_by(id, init_wt, c_pop, time_to_pocap) %>%
  mutate(avg_r_freq_per_group = mean(as.numeric(pop_prop), na.rm = TRUE)) %>%
  filter(time == max(time)) %>%
  mutate(init_wt = init_wt*c_pop) %>%
  group_by(init_wt, c_pop, time_to_pocap) %>%
  summarize(n_late = sum(time >= 7), n_early = sum(time < 7),
            n_r_freq_ge_50 = sum(avg_r_freq_per_group >= 0.5),
            n_r_freq_le_50 = sum(avg_r_freq_per_group < 0.5)) %>%
  mutate(perc_late_clear = n_late/(n_early + n_late),
         perc_r_freq_ge_50 = n_r_freq_ge_50/(n_r_freq_ge_50 + n_r_freq_le_50))

trial_clean_72 <- df_72_h %>%
  dplyr::select(time, init_wt, c_pop, time_to_pocap, pop_prop, id) %>%
  group_by(id, init_wt, c_pop, time_to_pocap) %>%
  mutate(avg_r_freq_per_group = mean(as.numeric(pop_prop), na.rm = TRUE)) %>%
  filter(time == max(time)) %>%
  mutate(init_wt = init_wt*c_pop) %>%
  group_by(init_wt, c_pop, time_to_pocap) %>%
  summarize(n_late = sum(time >= 7), n_early = sum(time < 7),
            n_r_freq_ge_50 = sum(avg_r_freq_per_group >= 0.5),
            n_r_freq_le_50 = sum(avg_r_freq_per_group < 0.5)) %>%
  mutate(perc_late_clear = n_late/(n_early + n_late),
         perc_r_freq_ge_50 = n_r_freq_ge_50/(n_r_freq_ge_50 + n_r_freq_le_50))


## Combining the data frames
full_df <- rbind(trial_clean_24, trial_clean_72)

full_df$t_pocap <- ifelse(full_df$time_to_pocap == 3, "24h", "72h")

## Find the regions where 24h is preferred above 72h
filter_df <- full_df %>%
  select(init_wt, c_pop, time_to_pocap, t_pocap, perc_r_freq_ge_50)

## Pivoting
wide_df <- filter_df %>%
  group_by(init_wt, c_pop) %>%
  mutate(entry = row_number()) %>%
  ungroup() %>%
  pivot_wider(
    names_from = entry, # Use the new 'entry' column for new column names
    values_from = c(time_to_pocap, t_pocap, perc_r_freq_ge_50) # The columns to get values from
  )

wide_df <- wide_df %>%
  mutate(p72_m_p24 = perc_r_freq_ge_50_1 - perc_r_freq_ge_50_2,
         init_wt = round(init_wt, 0))

#### Standard model
## Reading in the standard model data:
standard_model <- read.csv(paste0("250810_range_of_cpop/", today_date, "_c_pop_trial_data_standard.csv"))


df_24_h_s <- standard_model %>%
  filter(time_to_pocap == 3,
         type == "resistant")

df_72_h_s <- standard_model %>%
  filter(time_to_pocap == 9,
         type == "resistant")

## Cleaning

print(paste0("Cleaning"))
trial_clean_24_s <- df_24_h_s %>%
  dplyr::select(time, init_wt, c_pop, time_to_pocap, pop_prop, id) %>%
  group_by(id, init_wt, c_pop, time_to_pocap) %>%
  mutate(avg_r_freq_per_group = mean(as.numeric(pop_prop), na.rm = TRUE)) %>%
  filter(time == max(time)) %>%
  mutate(init_wt = init_wt*c_pop) %>%
  group_by(init_wt, c_pop, time_to_pocap) %>%
  summarize(n_late = sum(time >= 7), n_early = sum(time < 7),
            n_r_freq_ge_50 = sum(avg_r_freq_per_group >= 0.5),
            n_r_freq_le_50 = sum(avg_r_freq_per_group < 0.5)) %>%
  mutate(perc_late_clear = n_late/(n_early + n_late),
         perc_r_freq_ge_50 = n_r_freq_ge_50/(n_r_freq_ge_50 + n_r_freq_le_50))

trial_clean_72_s <- df_72_h_s %>%
  dplyr::select(time, init_wt, c_pop, time_to_pocap, pop_prop, id) %>%
  group_by(id, init_wt, c_pop, time_to_pocap) %>%
  mutate(avg_r_freq_per_group = mean(as.numeric(pop_prop), na.rm = TRUE)) %>%
  filter(time == max(time)) %>%
  mutate(init_wt = init_wt*c_pop) %>%
  group_by(init_wt, c_pop, time_to_pocap) %>%
  summarize(n_late = sum(time >= 7), n_early = sum(time < 7),
            n_r_freq_ge_50 = sum(avg_r_freq_per_group >= 0.5),
            n_r_freq_le_50 = sum(avg_r_freq_per_group < 0.5)) %>%
  mutate(perc_late_clear = n_late/(n_early + n_late),
         perc_r_freq_ge_50 = n_r_freq_ge_50/(n_r_freq_ge_50 + n_r_freq_le_50))


## Combining the data frames
full_df_s <- rbind(trial_clean_24_s, trial_clean_72_s)
full_df_s$t_pocap <- ifelse(full_df_s$time_to_pocap == 3, "24h", "72h")

## First, find the regions where 24h is preferred above 72h
filter_df_s <- full_df_s %>%
  select(init_wt, c_pop, time_to_pocap, t_pocap, perc_r_freq_ge_50)

## Pivoting
wide_df_s <- filter_df_s %>%
  group_by(init_wt, c_pop) %>%
  mutate(entry = row_number()) %>%
  ungroup() %>%
  pivot_wider(
    names_from = entry, # Use the new 'entry' column for new column names
    values_from = c(time_to_pocap, t_pocap, perc_r_freq_ge_50) # The columns to get values from
  )

wide_df_s <- wide_df_s %>%
  mutate(p72_m_p24 = perc_r_freq_ge_50_1 - perc_r_freq_ge_50_2,
         init_wt = round(init_wt, 0))


### Plotting tileplots                                                      ####
## Defining the scale for fill values
global_lims <- range(c(wide_df$p72_m_p24, wide_df_s$p72_m_p24), na.rm = TRUE)

#### Interaction model
perc_res_tileplot <- ggplot(wide_df, aes(x = c_pop,
                                         y = factor(round(init_wt, 0)),
                                         fill = p72_m_p24)) +
  geom_tile() +
  scale_x_log10(
    breaks = c_pop_init[seq(1, length(c_pop_init), 2)],
    labels = sapply(
      round(log10(c_pop_init[seq(1, length(c_pop_init), 2)]), 1),
      function(i) { parse(text = sprintf("10^%f", i)) }
    ),
    # expand = c(0,0)
  ) +
  scale_fill_gradientn(
    colors = c("purple", "white", "darkgreen"),
    values = scales::rescale(c(global_lims[1], 0, global_lims[2]),
                             from = global_lims),
    limits = global_lims,
    name = "% R. inf. 24h - 72h"
  ) +
  ylab("Initial number of viruses") +
  xlab(expression("Cell population (" * gamma * ")")) +
  # theme_minimal() +
  theme(
    ## Text size
    legend.text = element_text(size = legend_text_size*.8),
    legend.title = element_text(size = legend_text_size),
    text = element_text(size= axis_text_size), 
    axis.text = element_text(size= axis_text_size),
    ## Text color
    axis.text.x = element_text(colour="black"),
    axis.text.y = element_text(colour="black"),
    axis.title.y = element_text(),
    axis.title.x = element_text(),
    
    ## Remove axis lines inside plot
    # axis.line = element_blank(),
    panel.border = element_blank(),
    
    ## For faceted plots: The strip is the top of the facet
    strip.text.x = element_text(size= axis_text_size),
    
    ## Remove the background of the plot and panel
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA),
    
    # panel.border = element_blank(),          # remove border around panel
    # plot.background = element_rect(fill = "transparent", colour = NA), 
    # panel.background = element_rect(fill = "transparent", colour = NA),
    # axis.line = element_blank(),             # remove axis lines
    # axis.ticks = element_blank(),            # remove tick marks if you want
    strip.background = element_blank(),       # remove facet borders
    
    ## Legend stuff
    legend.position = "none",
    legend.position.inside = c(0.87, 0.17),
    legend.key.width = unit(0.3, 'cm'),
    legend.key.size = unit(0.3, "cm"),
    
    ## Extra tweaks
    # plot.margin = margin(0, 3, 0, 3),          # no margin around plot
    panel.spacing = unit(0, "lines")           # no gap between panels
  ) +
  scale_y_discrete(expand = c(0, 0))   

ggsave(plot = perc_res_tileplot, 
       paste0("250810_range_of_cpop/", today_date, "percent_resistant_interaction.png"),
       h = height/2, w = width/2, unit = "cm")
ggsave(plot = perc_res_tileplot, 
       paste0("250810_range_of_cpop/", today_date, "percent_resistant_interaction.pdf"),
       h = height/2, w = width/2, unit = "cm")

#### Standard model
perc_res_tileplot_s <- ggplot(wide_df_s, aes(x = c_pop,
                                         y = factor(round(init_wt, 0)),
                                         fill = p72_m_p24)) +
  geom_tile() +
  scale_x_log10(
    breaks = c_pop_init[seq(1, length(c_pop_init), 2)],
    labels = sapply(
      round(log10(c_pop_init[seq(1, length(c_pop_init), 2)]), 1),
      function(i) { parse(text = sprintf("10^%f", i)) }
    ),
    # expand = c(0,0)
  ) +
  scale_fill_gradientn(
    colors = c("purple", "white", "darkgreen"),
    values = scales::rescale(c(global_lims[1], 0, global_lims[2]),
                             from = global_lims),
    limits = global_lims,
    name = paste0("% Res. inf. 24h - 72h")
  ) +
  ylab("Initial number of viruses") +
  xlab(expression("Cell population (" * gamma * ")")) +
  # theme_minimal() +
  theme(
    ## Text size
    legend.text = element_text(size = legend_text_size*.8),
    legend.title = element_text(size = legend_text_size),
    text = element_text(size= axis_text_size), 
    axis.text = element_text(size= axis_text_size),
    ## Text color
    axis.text.x = element_text(colour="black"),
    axis.text.y = element_text(colour="black"),
    axis.title.y = element_text(),
    axis.title.x = element_text(),
    
    ## Remove axis lines inside plot
    # axis.line = element_blank(),
    panel.border = element_blank(),
    
    ## For faceted plots: The strip is the top of the facet
    strip.text.x = element_text(size= axis_text_size),
    
    ## Remove the background of the plot and panel
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA),
    
    # panel.border = element_blank(),          # remove border around panel
    # plot.background = element_rect(fill = "transparent", colour = NA), 
    # panel.background = element_rect(fill = "transparent", colour = NA),
    # axis.line = element_blank(),             # remove axis lines
    # axis.ticks = element_blank(),            # remove tick marks if you want
    strip.background = element_blank(),       # remove facet borders
    
    ## Legend stuff
    legend.position = "inside",
    legend.position.inside = c(0.8, 0.4),
    legend.key.width = unit(0.3, 'cm'),
    legend.key.size = unit(0.3, "cm"),
    legend.title.align = 1,
    
    ## Extra tweaks
    # plot.margin = margin(0, 3, 0, 3),          # no margin around plot
    panel.spacing = unit(0, "lines")           # no gap between panels
  ) +
  scale_y_discrete(expand = c(0, 0))  

ggsave(plot = perc_res_tileplot_s, 
       paste0("250810_range_of_cpop/", today_date, "percent_resistant_standard.png"),
       h = height/2, w = width/2, unit = "cm")
ggsave(plot = perc_res_tileplot_s, 
       paste0("250810_range_of_cpop/", today_date, "percent_resistant_standard.pdf"),
       h = height/2, w = width/2, unit = "cm")


#### Pulling out the plots that correspond to important slices              ####
## Our model
slices <- full_df %>%
  filter(round(init_wt,0) == 1 |round(init_wt,0) == 3 | round(init_wt,0) == 100)

plot_out <- ggplot(slices, aes(x = c_pop, y = perc_r_freq_ge_50, color = t_pocap)) +
  geom_line(linewidth = 1) +
  labs(color='t Pocap.') +
  scale_x_log10(
    breaks = c_pop_init[seq(1, length(c_pop_init), 2)],
    labels = sapply(
      round(log10(c_pop_init[seq(1, length(c_pop_init), 2)]), 1),
      function(i) { parse(text = sprintf("10^%f", i)) }
    )
  ) +
  scale_color_manual(
    values = c(
      "darkgreen",
      "purple"
    ),
    name = expression(t[plain("poc")])
  ) +
  facet_wrap(~round(init_wt, 0), ncol = 1) +
  theme_minimal() +
  geom_hline(yintercept = col_res_24, linetype = 3) +
  annotate("text", x = 10^1.1, y = col_res_24 + 0.07, label = "24h %", 
           size = geom_text_conv*axis_text_size) +
  geom_hline(yintercept = col_res_72, linetype = 2) +
  annotate("text", x = 10^1.1, y = col_res_72 + 0.07, label = "72h %", 
           size = geom_text_conv*axis_text_size) +
  xlab(expression("Cell population (" * gamma * ")")) +
  ylab("Percent resistant infections"
  ) +
  facet_wrap(~factor(round(init_wt, 0), levels = c(100,3,1)), ncol = 1) +
  theme(
    ## Text size
    legend.text = element_text(size = legend_text_size*.8),
    legend.title = element_text(size = legend_text_size*1.5),
    text = element_text(size= axis_text_size), 
    axis.text = element_text(size= axis_text_size),
    ## Text color
    axis.text.x = element_text(colour="black"),
    axis.text.y = element_text(colour="black"),
    axis.title.y = element_text(),
    axis.title.x = element_text(),
    ## Changing the lines in the plot
    # panel.grid.major.x = element_line(color = "lightgrey"), ## Major x lines
    # panel.grid.minor.x = element_blank(),                   ## Minor x lines
    # panel.grid.major.y = element_line(color = "lightgrey"), ## Major y lines
    # panel.grid.minor.y = element_line(color = "lightgrey"), ## Minor y lines
    ## For faceted plots: The strip is the top of the facet
    # strip.background = element_blank(),
    strip.text.x = element_text(size= axis_text_size),
    ## Remove the background of the plot and panel
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA),
    ## Legend stuff
    # legend.position = "none",
    legend.position = "none",
    legend.position.inside = c(0.925, 0.1),
    # legend.text = element_text(size = legend_text_size),
    # legend.title = element_blank(),
    legend.background = element_rect(fill = "white"),
    # legend.key = element_rect(fill = "transparent"),
    legend.key.width = unit(0, 'cm'),
    legend.key.size = unit(0.3, "cm"),
    legend.title.align = 0.5
  )

ggsave(filename = paste0("250810_range_of_cpop/", today_date, "c_pop_by_init_wt_tall_slice_interaction.pdf"), 
       plot = plot_out,
       h = height/2, w = width/2, unit = "cm")


## Standard model
slices_s <- full_df_s %>%
  filter(round(init_wt,0) == 3 |round(init_wt,0) == 1 | round(init_wt,0) == 100)

plot_out_s <- ggplot(slices_s, aes(x = c_pop, y = perc_r_freq_ge_50, color = t_pocap)) +
  geom_line(linewidth = 1) +
  labs(color='t Pocap.') +
  scale_color_manual(
    values = c(
      "darkgreen",
      "purple"
    ),
    name = expression(t[plain("poc")])
  ) +
  scale_x_log10(
    breaks = c_pop_init[seq(1, length(c_pop_init), 2)],
    labels = sapply(
      round(log10(c_pop_init[seq(1, length(c_pop_init), 2)]), 1),
      function(i) { parse(text = sprintf("10^%f", i)) }
    )
  ) +
  theme_minimal() +
  geom_hline(yintercept = col_res_24, linetype = 3) +
  annotate("text", x = 10^1.1, y = col_res_24 + 0.07, label = "24h %", 
           size = geom_text_conv*axis_text_size) +
  geom_hline(yintercept = col_res_72, linetype = 2) +
  annotate("text", x = 10^1.1, y = col_res_72 + 0.07, label = "72h %", 
           size = geom_text_conv*axis_text_size) +
  xlab(expression("Cell population (" * gamma * ")")) +
  ylab("Percent resistant infections") +
  facet_wrap(~factor(round(init_wt, 0), levels = c(100,3,1)), ncol = 1) +
  theme(
    ## Text size
    legend.text = element_text(size = legend_text_size*.8),
    legend.title = element_text(size = legend_text_size*1.5),
    text = element_text(size= axis_text_size), 
    axis.text = element_text(size= axis_text_size),
    ## Text color
    axis.text.x = element_text(colour="black"),
    axis.text.y = element_text(colour="black"),
    axis.title.y = element_text(),
    axis.title.x = element_text(),
    ## Changing the lines in the plot
    # panel.grid.major.x = element_line(color = "lightgrey"), ## Major x lines
    # panel.grid.minor.x = element_blank(),                   ## Minor x lines
    # panel.grid.major.y = element_line(color = "lightgrey"), ## Major y lines
    # panel.grid.minor.y = element_line(color = "lightgrey"), ## Minor y lines
    ## For faceted plots: The strip is the top of the facet
    # strip.background = element_blank(),
    strip.text.x = element_text(size= axis_text_size),
    ## Remove the background of the plot and panel
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA),
    ## Legend stuff
    # legend.position = "none",
    legend.position = "inside",
    legend.position.inside = c(0.925, 0.1),
    # legend.text = element_text(size = legend_text_size),
    # legend.title = element_blank(),
    legend.background = element_rect(fill = "white"),
    # legend.key = element_rect(fill = "transparent"),
    legend.key.width = unit(0, 'cm'),
    legend.key.size = unit(0.3, "cm"),
    legend.title.align = 0.5
  )

ggsave(filename = paste0("250810_range_of_cpop/", today_date, "c_pop_by_init_wt_tall_slice_standard.pdf"), 
       plot = plot_out_s,
       h = height/2, w = width/2, unit = "cm")

#### Create the figure                                                      ####
library(cowplot)

## Defining the panels
sfig_a <- perc_res_tileplot

sfig_b <- perc_res_tileplot_s

sfig_c <- plot_out

sfig_d <- plot_out_s

row1 <- plot_grid(sfig_a, sfig_b,
                  ncol = 2,
                  nrow = 1,
                  # rel_heights = c(0.25, 0.25, 0.5),
                  axis = "tblr", 
                  align = "h")

row2 <- plot_grid(sfig_c, sfig_d,
                  ncol = 2,
                  nrow = 1,
                  # rel_heights = c(0.25, 0.25, 0.5),
                  axis = "tblr", 
                  align = "h")

sfig <- plot_grid(row1, row2,
                  ncol = 1,
                  nrow = 2,
                  rel_heights = c(0.5, 1),
                  axis = "tblr", 
                  align = "h")

save_plot(plot = sfig, 
          filename = paste0("250810_range_of_cpop/", today_date, "_full_sfig.pdf"), 
          base_height = height-3, 
          base_width = width-2, 
          units = "cm", 
          bg = "transparent", 
          dpi = 1200)

save_plot(plot = sfig, 
          filename = paste0("250810_range_of_cpop/", today_date, "_full_sfig.jpg"), 
          base_height = height, 
          base_width = width, 
          units = "cm", 
          bg = "transparent", 
          dpi = 1200)




#### All slices                                                             ####
## Interaction model
full_plot_out <- ggplot(full_df, aes(x = c_pop, y = perc_r_freq_ge_50, color = t_pocap)) +
  geom_line(linewidth = 1) +
  scale_x_log10(breaks = c_pop_init,
                labels = sapply(round(log10(c_pop_init),2),function(i){parse(text = sprintf("10^%f",i))})) +
  theme_minimal() +
  geom_hline(yintercept = col_res_24, linetype = 3) +
  annotate("text", x = 10^3.6, y = col_res_24 + 0.05, label = "24h %") +
  geom_hline(yintercept = col_res_72, linetype = 2) +
  annotate("text", x = 10^3.6, y = col_res_72 + 0.05, label = "72h %") +
  xlab(expression("Cell population (" * gamma * ")")) +
  ylab("Percent resistant infections") +
  facet_wrap(~round(init_wt, 0), ncol = 1) +
  theme(
    ## Text size
    legend.text = element_text(size = legend_text_size*.6),
    legend.title = element_text(size = legend_text_size*.3 + legend_text_size),
    text = element_text(size= axis_text_size), 
    axis.text = element_text(size= axis_text_size),
    ## Text color
    axis.text.x = element_text(colour="black"),
    axis.text.y = element_text(colour="black"),
    axis.title.y = element_text(),
    axis.title.x = element_text(),
    ## Changing the lines in the plot
    # panel.grid.major.x = element_line(color = "lightgrey"), ## Major x lines
    # panel.grid.minor.x = element_blank(),                   ## Minor x lines
    # panel.grid.major.y = element_line(color = "lightgrey"), ## Major y lines
    # panel.grid.minor.y = element_line(color = "lightgrey"), ## Minor y lines
    ## For faceted plots: The strip is the top of the facet
    # strip.background = element_blank(),
    strip.text.x = element_text(size= axis_text_size),
    ## Remove the background of the plot and panel
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA),
    ## Legend stuff
    # legend.position = "none",
    # legend.position = "inside",
    # legend.position.inside = c(0.18, 0.06),
    # legend.text = element_text(size = legend_text_size),
    # legend.title = element_blank(),
    legend.background = element_rect(fill = "transparent"),
    legend.key = element_rect(fill = "transparent"),
    legend.key.width = unit(0, 'cm'),
    legend.key.size = unit(0.3, "cm")
  )

ggsave(filename = paste0("250810_range_of_cpop/", today_date, "c_pop_by_init_wt_tall_interaction.pdf"), 
       plot = full_plot_out,
       h = 36, w =18, unit = "cm")

## Standard model
full_plot_out_s <- ggplot(full_df_s, aes(x = c_pop, y = perc_r_freq_ge_50, color = t_pocap)) +
  geom_line(linewidth = 1) +
  scale_x_log10(breaks = c_pop_init,
                labels = sapply(round(log10(c_pop_init),2),function(i){parse(text = sprintf("10^%f",i))})) +
  theme_minimal() +
  geom_hline(yintercept = col_res_24, linetype = 3) +
  annotate("text", x = 10^3.6, y = col_res_24 + 0.05, label = "24h %") +
  geom_hline(yintercept = col_res_72, linetype = 2) +
  annotate("text", x = 10^3.6, y = col_res_72 + 0.05, label = "72h %") +
  xlab(expression("Cell population (" * gamma * ")")) +
  ylab("Percent resistant infections") +
  facet_wrap(~round(init_wt, 0), ncol = 1) +
  theme(
    ## Text size
    legend.text = element_text(size = legend_text_size),
    legend.title = element_text(size = legend_text_size),
    text = element_text(size= axis_text_size), 
    axis.text = element_text(size= axis_text_size),
    ## Text color
    axis.text.x = element_text(colour="black"),
    axis.text.y = element_text(colour="black"),
    axis.title.y = element_text(),
    axis.title.x = element_text(),
    ## Changing the lines in the plot
    # panel.grid.major.x = element_line(color = "lightgrey"), ## Major x lines
    # panel.grid.minor.x = element_blank(),                   ## Minor x lines
    # panel.grid.major.y = element_line(color = "lightgrey"), ## Major y lines
    # panel.grid.minor.y = element_line(color = "lightgrey"), ## Minor y lines
    ## For faceted plots: The strip is the top of the facet
    # strip.background = element_blank(),
    strip.text.x = element_text(size= axis_text_size),
    ## Remove the background of the plot and panel
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA),
    ## Legend stuff
    # legend.position = "none",
    legend.position = "inside",
    legend.position.inside = c(0.9, 0.08),
    # legend.text = element_text(size = legend_text_size),
    # legend.title = element_blank(),
    # legend.background = element_rect(fill = "transparent"),
    # legend.key = element_rect(fill = "transparent"),
    legend.key.width = unit(0, 'cm'),
    legend.key.size = unit(0.3, "cm")
  )

ggsave(filename = paste0("250810_range_of_cpop/", today_date, "c_pop_by_init_wt_tall_standard.pdf"), 
       plot = full_plot_out_s,
       h = 36, w =18, unit = "cm")

