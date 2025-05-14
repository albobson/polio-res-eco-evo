
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
    output = list('runs/ddt_logistic_mu_2e-05_2025-03-24/res/fig4/B.rds', 'runs/ddt_logistic_mu_2e-05_2025-03-24/dat_gen/sims/traj_for_phase_plane.csv'),
    params = list('2025-03-24', c('logistic'), c(2e-05), "date" = '2025-03-24', "fit_func" = c('logistic'), "mu" = c(2e-05)),
    wildcards = list(),
    threads = 1,
    log = list(),
    resources = list('mem_mb', 'disk_mb', 'tmpdir', 'log_loc', 'mfree', 'cluster_time', 'disk_free', 'runtime', 'n_cores', "mem_mb" = '3G', "disk_mb" = '4G', "tmpdir" = '/tmp/12666916.1.feder-short.q', "log_loc" = './log_scr/', "mfree" = '3G', "cluster_time" = '0:00:10:00', "disk_free" = '4G', "runtime" = 10, "n_cores" = '1'),
    config = list(),
    rule = 'fig4B',
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
############################### Figure 4 B #####################################

## Reason:

## Over time, the viral population goes from a high MOI (where resistance is
## suppressed) to a low MOI (where resistance is expressed). At a low MOI,
## resistance can increase in frequency and become the dominant genotype.

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
curr_run <- determ_polv(n = 6,
                        moi_mut_start = 100*0.0001, 
                        moi_wt_start = 100*(1-0.0001),
                        fit_func_in = fit_func$prob_surv,
                        v_prog = optim_params$optim_v_prog,
                        p2pfu = optim_params$optim_p2pfu
) 

## Cleaning the data for plotting on phase plane
phase_plane_var <- curr_run %>%
  select(type, time, moi_res, moi_wt, pop_prop) %>%
  filter(type == "resistant") %>%
  mutate(moi_tot = (moi_res + moi_wt), pop_prop = as.numeric(pop_prop)) %>%
  arrange(time) %>%
  mutate(delta_p = (lead(pop_prop, default = last(pop_prop)) - pop_prop),
         delta_moi = (lead(moi_tot, default = last(moi_tot)) - moi_tot),
         dist_col = abs((delta_p)+(delta_moi)/100),
         moi_color = factor(round(moi_tot,5), 
                            levels = unique(as.numeric(round(moi_tot,5))))) %>%
  ungroup()

## Saving this data to be used for the phase plane
suppressWarnings(dir.create("dat_gen/sims"))
write.csv(phase_plane_var, file = "dat_gen/sims/traj_for_phase_plane.csv")


#### Plotting trajectory                                                    ####
plot4b <- ggplot(curr_run, aes(x = time, y = moi_type, color = type)) +
  geom_point(size = 3) + 
  geom_line(linewidth = 1) +
  scale_color_manual(labels = c("Resistant", "Susceptible"),
                     values = c("#c94d4d", "#688cc8"),
                     name = "Genotype") +
  theme_light() + 
  xlab("Passages") + ylab("MOI") + 
  scale_y_log10(breaks=10^(-6:10),
                labels = sapply(-6:10,function(i){parse(text = sprintf("10^%d",i))})) +
  scale_x_continuous(breaks = c(0:max(curr_run$time))) +
  theme(text = element_text(size = 20), axis.text = element_text(size=20), 
        axis.title = element_text(size = 22),
        legend.text = element_text(size = 20), 
        legend.title = element_blank(),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        legend.position=c(0.85,0.2))

plot4b

ggsave(paste0(filepath, "res/fig4/B.png"), plot4b, h = 8, w = 8, bg = "transparent")
ggsave(paste0(filepath, "res/fig4/B.jpg"), plot4b, h = 8, w = 8, bg = "transparent")
ggsave(paste0(filepath, "res/fig4/B.svg"), plot4b, h = 8, w = 8, bg = "transparent", device = 'svg')
saveRDS(plot4b, paste0(filepath, "res/fig4/B.rds"))