######################### Supplemental Figure 2 ################################

## Reason:

## The vector field gives an insight into the trends of the model but does not
## directly translate to how quickly those trends are realized. Here, I will
## plot steps of the simulations initated across the phase plane.

#### Set Up                                                                 ####

## Read in fitness function and parameters
fit_func <- read.csv("dat_gen/params/fitness_function.csv")
optim_params <- read.csv("dat_gen/params/optim_params.csv")

## Source functions to run simulations
source("scr/polv_DDT_functions.R")  ## Main functions 


#### Running simulations                                                    ####
## Initial MOIs
mois <- 10^seq(from = -2, to = 2.35, length.out = 11)

## Initial resistance frequencies
r_freq <- seq(from = 0, to = 1, length.out = 11)

## Max number of steps to run
max_time <- 50

## Set up parallel environment:
cl <- makeCluster(detectCores()-1)
registerDoParallel(cl)

## Run simulation
sims <- foreach(i = mois, .combine = 'rbind', .packages = c('dplyr')) %:%
  foreach(j = r_freq, .combine = 'rbind', .packages = c('dplyr')) %dopar% {
    (determ_polv(n = max_time,
                 moi_wt_start = (1-j)*i,
                 moi_mut_start = (j)*i,
                 fit_func_in = fit_func$prob_surv,
                 v_prog = optim_params$optim_v_prog,
                 p2pfu = optim_params$optim_p2pfu,
    ))
  }

stopImplicitCluster()
stopCluster(cl)

## Save for later just in case 
write.csv(sims, file = "dat_gen/sims/vector_field_by_time.csv")

## Clean simulation
sims_clean <- sims %>%
  select(type, time, moi_res, moi_wt, init_wt, init_mut, pop_prop) %>%
  filter(type == "resistant") %>%
  group_by(init_mut, init_wt, time) %>%
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
         moi_tot = (moi_res + moi_wt),
         dist_col = abs(delta_moi_tot) + abs(delta_p)) %>%
  mutate(x_end = (init_moi_tot+cos(atan(delta_p/delta_moi_tot))),
         y_end = (init_freq_mut+sin(atan(delta_p/delta_moi_tot)))) %>%
  ungroup()

## Generating a color indicator for each starting condition
sims_clean$unique_id <- as.numeric(interaction(sims_clean$init_moi_tot, 
                                               sims_clean$init_freq_mut))


#### Plot a few specific time points                                        ####
times_to_plot <- c(0, 1,3, 6, 12, 24, max_time)

for(t in 1:length(times_to_plot)) {
  temp_df <- temp_plot <- NULL
  
  temp_df <- sims_clean %>%
    filter(time == times_to_plot[t])
  
  temp_plot <- ggplot(temp_df, aes(x = moi_tot, y = pop_prop, color = unique_id)) +
    geom_point() +
    theme_light() +
    scale_color_viridis_c(option = "plasma") +
    ylim(c(0,1)) +
    scale_x_continuous(trans="log10", breaks=10^(-2:3),
                       labels = sapply(c(-2:3),function(i){parse(text = sprintf("10^%d",i))}),
                       limits = c(10^(-2.1), 10^(2.6))) +
    labs(color = "Total MOI") +
    xlab("Total MOI") + ylab("Resistance Frequency") +
    theme(text = element_text(size = 20), axis.text = element_text(size=20), 
          axis.title = element_text(size = 22),
          legend.text = element_text(size = 20), 
          legend.title = element_text(size = 20),
          legend.position = "none",
          axis.text.x = element_text(colour="black"),
          axis.text.y = element_text(colour="black"))
  
  suppressWarnings(dir.create("sup_res/fig2/"))
  
  ggsave(paste0("sup_res/fig2/", LETTERS[t], ".png"), temp_plot, h = 8, w = 8, bg = "transparent")
  ggsave(paste0("sup_res/fig2/", LETTERS[t], ".jpg"), temp_plot, h = 8, w = 8, bg = "transparent")
  ggsave(temp_plot, file = paste0("sup_res/fig2/", LETTERS[t], ".svg"), bg = "transparent", device = 'svg')
  saveRDS(temp_plot, paste0("sup_res/fig2/", LETTERS[t], ".rds"))
  
}