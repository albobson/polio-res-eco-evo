
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
    input = list('runs/ddt_logistic_mu_2e-05_2025-03-24/dat_gen/params/logistic_fitness_function.csv', 'runs/ddt_logistic_mu_2e-05_2025-03-24/dat_gen/params/tan_fit_df.csv'),
    output = list('runs/ddt_logistic_mu_2e-05_2025-03-24/res/fig3/B.rds'),
    params = list('2025-03-24', c('logistic'), c(2e-05), "date" = '2025-03-24', "fit_func" = c('logistic'), "mu" = c(2e-05)),
    wildcards = list(),
    threads = 1,
    log = list(),
    resources = list('mem_mb', 'disk_mb', 'tmpdir', 'log_loc', 'mfree', 'cluster_time', 'disk_free', 'runtime', 'n_cores', "mem_mb" = '3G', "disk_mb" = '4G', "tmpdir" = '/tmp/12668491.1.feder-short.q', "log_loc" = './log_scr/', "mfree" = '3G', "cluster_time" = '0:00:10:00', "disk_free" = '4G', "runtime" = 10, "n_cores" = '1'),
    config = list(),
    rule = 'fig3B',
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
############################### Figure 3 B #####################################

## Reason:

## This figure pannel shows what the fitness function looks like for our
## analysis. For a comparison of other fitness functions, please look at
## supplemental figure 1.

#### Set Up                                                                 ####

## libraries
require(dplyr)
require(ggplot2)

## Set up where to save
## Date of the run
init_date <- as.character(snakemake@params[["date"]])

## Fitness function
init_fit_func <- as.character(snakemake@params[["fit_func"]])

## Mutation rate
init_mu <- as.numeric(snakemake@params[["mu"]])

## Create the filepath where things will be saved
filepath <- paste0("runs/ddt_", init_fit_func, "_mu_", init_mu, "_", init_date, "/")

print(filepath)

#### Read in the data                                                       ####
## Tanner data points
fit_df <- read.csv(file = paste0(filepath, "dat_gen/params/tan_fit_df.csv"))

print(fit_df)

## Logistic fit
fix_logistic <- read.csv(file = paste0(filepath, "dat_gen/params/logistic_fitness_function.csv"))

print(fix_logistic)

#### Visualize                                                              ####
func_plot <- ggplot(fit_df, aes(x = subunits, y = prob_surv)) +
  geom_line(data = fix_logistic, 
            aes(group = fit_type),
            linewidth = 1) +
  geom_point(color = "black") +
  xlab("Number of Resistant Subunits") +
  ylab("Probability of Survival") +
  scale_y_continuous(trans="log10",
                     limits = c(0.0001,1),
                     breaks=10^(-4:0),
                     labels = c(sapply(-4:-1,function(i){parse(text = sprintf("10^%d",i))}), "1 ")) +
  theme_light() +
  theme(text = element_text(size = 20), axis.text = element_text(size=20), 
        axis.title = element_text(size = 22),
        legend.text = element_text(size = 20), legend.title = element_text(size = 20),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        legend.position = "none",
  )

print("Plotted")

dir.create(paste0(filepath, "res/"))
dir.create(paste0(filepath, "res/fig3/"))

print("Dirs created")

ggsave(paste0(filepath, "res/fig3/logistic_B.png"), func_plot, h = 8, w = 8, bg = "transparent")
ggsave(paste0(filepath, "res/fig3/logistic_B.jpg"), func_plot, h = 8, w = 8, bg = "transparent")
ggsave(paste0(filepath, "res/fig3/logistic_B.svg"), func_plot, h = 8, w = 8, bg = "transparent", device = 'svg')
saveRDS(paste0(filepath, func_plot, "res/fig3/logistic_B.rds"))

print("Figs saved")