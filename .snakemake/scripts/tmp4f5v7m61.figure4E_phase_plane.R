
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
    input = list('runs/ddt_logistic_mu_2e-05_2025-03-24/dat_gen/params/optim_params.csv', 'runs/ddt_logistic_mu_2e-05_2025-03-24/dat_gen/params/logistic_fitness_function.csv', 'runs/ddt_logistic_mu_2e-05_2025-03-24/dat_gen/sims/traj_for_phase_plane.csv'),
    output = list('runs/ddt_logistic_mu_2e-05_2025-03-24/res/fig4/E.rds'),
    params = list('2025-03-24', c('logistic'), c(2e-05), "date" = '2025-03-24', "fit_func" = c('logistic'), "mu" = c(2e-05)),
    wildcards = list(),
    threads = 1,
    log = list(),
    resources = list('mem_mb', 'disk_mb', 'tmpdir', 'log_loc', 'mfree', 'cluster_time', 'disk_free', 'runtime', 'n_cores', "mem_mb" = '3G', "disk_mb" = '4G', "tmpdir" = '/tmp/12679870.1.feder-short.q', "log_loc" = './log_scr/', "mfree" = '3G', "cluster_time" = '0:00:45:00', "disk_free" = '4G', "runtime" = 45, "n_cores" = '10'),
    config = list(),
    rule = 'fig4E',
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


#### Running simulations                                                    ####
mois <- 10^seq(from = -2, to = 2.35, length.out = 11)

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




#### Plot                                                                   ####
plot4e <- ggplot(sims_clean, aes(x = init_moi_tot, y = init_freq_mut)) +
  geom_segment(aes(xend = (init_moi_tot+(delta_moi_tot/80)), yend = (init_freq_mut+(delta_p/80))),
               arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.6) +
  theme_light() +
  ylim(c(0,1)) +
  scale_x_continuous(trans="log10", breaks=10^(-2:3),
                     labels = sapply(c(-2:3),function(i){parse(text = sprintf("10^%d",i))}),
                     limits = c(10^(-2.1), 10^(2.6))) +
  geom_segment(phase_plane_var, mapping = aes(x = moi_tot, y = pop_prop, 
                                              xend = (moi_tot+(delta_moi)), 
                                              yend = (pop_prop+(delta_p)),
                                              color = factor(moi_tot,
                                                             levels = moi_tot)),
               arrow = arrow(length = unit(0.3, "cm")), linewidth = 1.4) +
  labs(color = "Total MOI") +
  scale_colour_manual(values = c("#FFA61C", "#64A61C", rep("#000", 5))) +
  xlab("Total MOI") + ylab("Resistance Frequency") +
  theme(text = element_text(size = 20), axis.text = element_text(size=20), 
        axis.title = element_text(size = 22),
        legend.text = element_text(size = 20), 
        legend.title = element_text(size = 20),
        legend.position = "none",
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"))

plot4e

ggsave(paste0(filepath, "res/fig4/E.png"), plot4e, h = 8, w = 8, bg = "transparent")
ggsave(paste0(filepath, "res/fig4/E.jpg"), plot4e, h = 8, w = 8, bg = "transparent")
ggsave(paste0(filepath, "res/fig4/E.svg"), plot4e, h = 8, w = 8, bg = "transparent", device = 'svg')
saveRDS(plot4e, paste0(filepath, "res/fig4/E.rds"))