
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
    input = list('dat/collett_trial.csv', 'runs/ddt_logistic_mu_2e-05_2025-03-24/dat_gen/params/optim_params.csv', 'runs/ddt_logistic_mu_2e-05_2025-03-24/dat_gen/params/logistic_fitness_function.csv', 'runs/ddt_logistic_mu_2e-05_2025-03-24/dat_gen/sims/sim_clinical_trial.csv'),
    output = list('runs/ddt_logistic_mu_2e-05_2025-03-24/res/fig5/C.rds'),
    params = list('2025-03-24', c('logistic'), c(2e-05), "date" = '2025-03-24', "fit_func" = c('logistic'), "mu" = c(2e-05)),
    wildcards = list(),
    threads = 1,
    log = list(),
    resources = list('mem_mb', 'disk_mb', 'tmpdir', 'log_loc', 'mfree', 'cluster_time', 'disk_free', 'runtime', 'n_cores', "mem_mb" = '3G', "disk_mb" = '4G', "tmpdir" = '/tmp/12680208.1.feder-short.q', "log_loc" = './log_scr/', "mfree" = '3G', "cluster_time" = '0:00:10:00', "disk_free" = '4G', "runtime" = 10, "n_cores" = '1'),
    config = list(),
    rule = 'fig5C',
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
############################### Figure 5 C #####################################

## Reason:

## We want to plot an early clearer from our simulated clinical trial

#### Set Up                                                                 ####

## Source functions to run simulations
source("scr/polv_DDT_functions.R")  ## Main functions 

## Read in Collett data
coldat <- read.csv("dat/collett_trial.csv")

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

## Read in clincial trial simulations
trial <- read.csv(paste0(filepath, "dat_gen/sims/sim_clinical_trial.csv"))

## Create directory to save stuff (if not already created)
suppressWarnings(dir.create(paste0(filepath, "res/fig5")))

## An individual early clearer who cleared susceptible
early <- trial %>%
  group_by(id) %>%
  mutate(res = sum(surv_pfu)/sum(as.numeric(tot_pfu)) > 0.5,
         max_date = max(time)) %>%
  ungroup()

## Finding the ranges of dates
date_range <- sort(unique(early$max_date))

## Save max date for our scaling
max_day <- ceiling(last(date_range))

## Find the list of the people who cleared 3rd earliest (to make the plot more interesting)
early_id <- early %>%
  filter(max_date == date_range[3]) %>%
  select(id) %>% 
  unique()

## We'll just select the first one
f_early_id <- as.vector(unlist(early_id[1,]))

early_sim <- early %>%
  filter(id == f_early_id)

## Now going to plot this
early_sim_p <- ggplot(early_sim, aes(x = time, y = moi_type, color = type)) +
  # geom_point(size = ) + 
  geom_line(linewidth = 1) +
  scale_color_manual(values = c("#c94d4d", "#688cc8"),
                     name = "Genotype",
                     labels = c("Resistant", "Susceptible"),) +
  theme_light() + 
  xlab("Number of Replication Cycles") + ylab("MOI") + 
  scale_x_continuous(limits = c(0, max_day)) +
  scale_y_log10(breaks=10^(-6:10),
                labels = sapply(-6:10,function(i){parse(text = sprintf("10^%d",i))})) +
  theme(text = element_text(size = 20), axis.text = element_text(size=20), 
        axis.title = element_text(size = 22),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        rect = element_rect(fill = "transparent"),
        legend.text = element_text(size = 14), 
        legend.title = element_text(size = 16),
        legend.position=c(0.85,0.85)
        )

# early_sim_p

ggsave(early_sim_p, file = paste0(filepath, "res/fig5/C.png"), bg = "transparent")
ggsave(early_sim_p, file = paste0(filepath, "res/fig5/C.jpg"), bg = "transparent")
ggsave(early_sim_p, file = paste0(filepath, "res/fig5/C.svg"), bg = "transparent", device = 'svg')
saveRDS(early_sim_p, paste0(filepath, "res/fig5/C.rds"))
