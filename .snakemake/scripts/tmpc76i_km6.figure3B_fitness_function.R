
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
    resources = list('mem_mb', 'disk_mb', 'tmpdir', 'log_loc', 'mfree', 'cluster_time', 'disk_free', 'runtime', 'n_cores', "mem_mb" = '3G', "disk_mb" = '4G', "tmpdir" = '/tmp/12722196.1.feder-short.q', "log_loc" = './log_scr/', "mfree" = '3G', "cluster_time" = '0:00:10:00', "disk_free" = '4G', "runtime" = 10, "n_cores" = '1'),
    config = list(),
    rule = 'fig3B',
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
############################### Figure 3 B #####################################

## Reason:

## This figure pannel shows what the fitness function looks like for our
## analysis. For a comparison of other fitness functions, please look at
## supplemental figure 1.

#### Set Up                                                                 ####

## libraries
require(dplyr)
require(ggplot2)
library(cowplot)

## Set up where to save
## Date of the run
init_date <- as.character(snakemake@params[["date"]])

## Fitness function
init_fit_func <- as.character(snakemake@params[["fit_func"]])

## Mutation rate
init_mu <- as.numeric(snakemake@params[["mu"]])

## Create the filepath where things will be saved
filepath <- paste0("runs/ddt_", init_fit_func, "_mu_", init_mu, "_", init_date, "/")

## Making the color map for plotting
color_in <- c("#c94d4d", "#688cc8")

## geom_text() is in terms of mm, where element_text() is in terms of pt. I want
## my geom_text and axis text to be the same. Here, I am defining a conversion
## which will keep the text sizes the same between text types.
geom_text_conv = 0.3528

## I also want my axis text to be 12 pt font, so, setting that here
axis_text_size = 14

## I would like for the legend to be slightly smaller than the rest of the text.
## I'm setting that size here
legend_text_size = 10

#### Read in the data                                                       ####
## Tanner data points
fit_df <- read.csv(file = paste0(filepath, "dat_gen/params/tan_fit_df.csv"))

## Logistic fit
fix_logistic <- read.csv(file = paste0(filepath, "dat_gen/params/logistic_fitness_function.csv"))

#### Visualize                                                              ####
func_plot <- ggplot(fit_df, aes(x = subunits, y = prob_surv)) +
  geom_line(data = fix_logistic, 
            aes(group = fit_type),
            linewidth = 1) +
  geom_point(color = "black") +
  xlab("Resistant Subunits") +
  ylab("Probability of Survival") +
  scale_y_continuous(trans="log10",
                     limits = c(0.0001,1),
                     breaks=10^(-4:0),
                     labels = c(sapply(-4:-1,function(i){parse(text = sprintf("10^%d",i))}), "1 ")) +
  theme_light() +
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
    # legend.position.inside = c(0.18, 0.06),
    # legend.text = element_text(size = legend_text_size),
    # legend.title = element_blank(),
    # legend.background = element_rect(fill = "transparent"),
    # legend.key = element_rect(fill = "transparent"),
    # legend.key.width = unit(0, 'cm'),
    # legend.key.size = unit(0.3, "cm")
  ) 

func_plot

suppressWarnings(dir.create(paste0(filepath, "res/")))
suppressWarnings(dir.create(paste0(filepath, "res/fig3/")))

ggsave(paste0(filepath, "res/fig3/B.png"), func_plot, h = 3.5*0.72, w = 3, unit = "in", bg = "transparent")
ggsave(paste0(filepath, "res/fig3/B.pdf"), func_plot, h = 3.5*0.72, w = 3, unit = "in", bg = "transparent")
ggsave(paste0(filepath, "res/fig3/B.svg"), func_plot, h = 3.5*0.72, w = 3, unit = "in", bg = "transparent", device = 'svg')
saveRDS(func_plot, file = paste0(filepath, "res/fig3/B.rds"))
