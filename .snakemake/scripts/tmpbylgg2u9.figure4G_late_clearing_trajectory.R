
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
    output = list('runs/ddt_logistic_mu_2e-05_2025-03-24/res/fig4/G.rds'),
    params = list('2025-03-24', c('logistic'), c(2e-05), "date" = '2025-03-24', "fit_func" = c('logistic'), "mu" = c(2e-05)),
    wildcards = list(),
    threads = 1,
    log = list(),
    resources = list('mem_mb', 'disk_mb', 'tmpdir', 'log_loc', 'mfree', 'cluster_time', 'disk_free', 'runtime', 'n_cores', "mem_mb" = '3G', "disk_mb" = '4G', "tmpdir" = '/tmp/12724533.1.feder-short.q', "log_loc" = './log_scr/', "mfree" = '3G', "cluster_time" = '0:00:10:00', "disk_free" = '4G', "runtime" = 10, "n_cores" = '1'),
    config = list(),
    rule = 'fig4G',
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
############################### Figure 4 G #####################################

## Reason:

## We want to plot a late clearer from our simulated clinical trial

#### Set Up                                                                 ####

## Source functions to run simulations
source("scr/polv_DDT_functions.R")  ## Main functions 

## Library
library(cowplot)

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



#### Data Cleaning                                                          ####
early <- trial %>%
  group_by(id) %>%
  mutate(res = sum(surv_pfu)/sum(as.numeric(tot_pfu)) > 0.5,
         max_date = max(time)) %>%
  ungroup()

## Finding the ranges of dates
date_range <- sort(unique(early$max_date))

## Find the list of the people who cleared 3rd earliest (to make the plot more interesting)
early_id <- early %>%
  filter(max_date == date_range[3]) %>%
  select(id) %>% 
  unique()



## An individual late clearer who cleared resistant
late_id <- early %>%
  filter(max_date == max(date_range)) %>%
  select(id) %>% 
  unique()

## We'll just select the first one
f_late_id <- as.vector(unlist(late_id[1,]))

late_sim <- early %>%
  filter(id == f_late_id)

## A function to simplify the x axis:
evens <- function(x) subset(x, x %% 2 == 0)

#### Plotting                                                               ####
late_sim_p <- ggplot(late_sim, aes(x = time, y = moi_type, color = type)) +
  # geom_point(size = ) + 
  geom_line(linewidth = 1) +
  scale_color_manual(labels = c("Resistant", "Susceptible"),
                     values = color_in,
                     name = "Genotype") +
  theme_light() + 
  xlab("Replications") + 
  ylab("MOI") + 
  scale_x_continuous() +
  scale_y_log10(breaks=10^(evens(-6:10)),
                labels = sapply(evens(-6:10), function(i){parse(text = sprintf("10^%d",i))})) +
  theme(
    ## Text size
    text = element_text(size= axis_text_size), 
    axis.text = element_text(size= axis_text_size),
    ## Text color
    axis.text.x = element_text(colour="black"),
    axis.text.y = element_text(colour="black"),
    axis.title.y = element_text(),
    axis.title.x = element_blank(),
    ## Changing the lines in the plot 
    panel.grid.major.x = element_line(color = "lightgrey"), ## Major x lines
    panel.grid.minor.x = element_blank(),                   ## Minor x lines
    panel.grid.major.y = element_line(color = "lightgrey"), ## Major y lines
    panel.grid.minor.y = element_blank(), ## Minor y lines
    ## For faceted plots: The strip is the top of the facet
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    ## Remove the background of the plot and panel
    panel.background = element_rect(fill='transparent'), 
    plot.background = element_rect(fill='transparent', color=NA),
    ## Legend stuff
    # legend.position = "none",
    legend.position = "inside",
    legend.position.inside = c(0.85,0.83),
    legend.text = element_text(size = legend_text_size),
    legend.title = element_blank(),
    legend.background = element_rect(fill = "transparent"),
    legend.key = element_rect(fill = "transparent"),
    legend.key.width = unit(0, 'cm'),
    legend.key.size = unit(0.3, "cm")
  ) 

late_sim_p

## text plot for "Late clearer" 
late_clear <- data.frame(cond = 0, name = "Late Clearer")

## Plotting just pocapavir
late_clear_p <- ggplot(late_clear, aes(cond, name, label = name)) +
  geom_text(size = axis_text_size * geom_text_conv) +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.text.x = element_blank(),
        text = element_text(size= axis_text_size),
        strip.text = element_blank(),
        axis.text.y=element_blank(),
        strip.background = element_blank())

late_clear_p

## Full plot
late_p <- plot_grid(late_clear_p, NULL, late_sim_p,
                      rel_heights = c(0.3, -0.1, 1),
                      ncol = 1, axis = "l", align = "v")
late_p

late_list <- list(late_clear_p, late_sim_p)


## Save
ggsave(late_p, file = paste0(filepath, "res/fig4/G.png"), h = 22/3/2, w = 18/2, unit = "cm", bg = "transparent", dpi = 300)
ggsave(late_p, file = paste0(filepath, "res/fig4/G.pdf"), h = 22/3/2, w = 18/2, unit = "cm", bg = "transparent", dpi = 300)
ggsave(late_p, file = paste0(filepath, "res/fig4/G.svg"), h = 22/3/2, w = 18/2, unit = "cm", bg = "transparent", dpi = 300)
saveRDS(late_list, paste0(filepath, "res/fig4/G.rds"))

