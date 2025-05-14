
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
    input = list('runs/ddt_logistic_mu_2e-05_2025-03-24/dat_gen/params/optim_params.csv', 'runs/ddt_logistic_mu_2e-05_2025-03-24/dat_gen/params/logistic_fitness_function.csv'),
    output = list('runs/ddt_logistic_mu_2e-05_2025-03-24/res/fig4/A.rds'),
    params = list('2025-03-24', c('logistic'), c(2e-05), "date" = '2025-03-24', "fit_func" = c('logistic'), "mu" = c(2e-05)),
    wildcards = list(),
    threads = 1,
    log = list(),
    resources = list('mem_mb', 'disk_mb', 'tmpdir', 'log_loc', 'mfree', 'cluster_time', 'disk_free', 'runtime', 'n_cores', "mem_mb" = '3G', "disk_mb" = '4G', "tmpdir" = '/tmp/12668673.1.feder-short.q', "log_loc" = './log_scr/', "mfree" = '3G', "cluster_time" = '0:00:30:00', "disk_free" = '4G', "runtime" = 30, "n_cores" = '10'),
    config = list(),
    rule = 'fig4A',
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
############################### Figure 4 A #####################################

## Reason:

## Resistance suppression is an MOI dependent process.

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


#### Run simulations                                                        ####
range <- 0.0001
mois <- 10^seq(from = -2, to = 2, length.out = 50)

cl <- makeCluster(detectCores()-2)
registerDoParallel(cl)

low_freq_moi <- foreach(r = range, 
                        .combine = 'rbind', 
                        .packages = c('dplyr')) %:%
  foreach(m = mois, 
          .combine = 'rbind', 
          .packages = c('dplyr'))  %dopar% {
            (determ_polv(n = 1,
                         moi_mut_start = m*r, 
                         moi_wt_start = m*(1-r),
                         fit_func_in = fit_func$prob_surv,
                         v_prog = optim_params$optim_v_prog,
                         p2pfu = optim_params$optim_p2pfu
            ))
          }

stopCluster(cl)
stopImplicitCluster()


## Clean simulation
low_freq_clean <- low_freq_moi %>%
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
         dist_col = abs(delta_moi_tot) + abs(delta_p),
         moi_tot = (moi_wt + moi_res),
         last_p = lag(pop_prop, default = first(pop_prop))) %>%
  mutate(x_end = (init_moi_tot+cos(atan(delta_p/delta_moi_tot))),
         y_end = (init_freq_mut+sin(atan(delta_p/delta_moi_tot)))) %>%
  filter(time == 1) %>%
  ungroup()


plot4a <- ggplot(low_freq_clean, aes(x = init_moi_tot, y = delta_p)) +
  geom_line(linewidth = 1.3) +
  # scale_x_continuous(trans="log10", breaks=10^(-4:2),
  #                    labels = sapply(-4:2,function(i){parse(text = sprintf("10^%d",i))}),
  #                    limits = c(10^(-4), 10^(2))) +
  ylim(0, 0.15) +
  xlim(0, 102) +
  theme_minimal() +
  ## This adds the points to the plot
  # geom_point(data = df_points_clean, 
  #            mapping = aes(x = init_moi_tot, 
  #                          y = delta_p, 
  #                          color = factor(round(init_moi_tot, 2))),
  #            size = 5) +
  # scale_color_manual(values = c("#64A61C", "#FFA61C"),
  #                    name = "Total MOI") +
  xlab("Total MOI") + ylab(paste0("\U0394 Freq. Resistance")) +
  theme(text = element_text(size = 20), 
        axis.text = element_text(size=20), 
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        axis.title = element_text(size = 22),
        legend.text = element_text(size = 20), 
        legend.title = element_text(size = 20),
        legend.position = "none",
        plot.title = element_text(size=22, hjust = 0.4),
  ) 

ggsave(paste0(filepath, "res/fig4/A.png"), plot4a, h = 8, w = 8, bg = "transparent")
ggsave(paste0(filepath, "res/fig4/A.jpg"), plot4a, h = 8, w = 8, bg = "transparent")
ggsave(paste0(filepath, "res/fig4/A.svg"), plot4a, h = 8, w = 8, bg = "transparent", device = 'svg')
saveRDS(plot4a, paste0(filepath, "res/fig4/A.rds"))
