
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
    input = list('dat/collett_trial.csv'),
    output = list('runs/ddt_logistic_mu_2e-05_2025-03-24/res/fig2/B.rds'),
    params = list('2025-03-24', c('logistic'), c(2e-05), "date" = '2025-03-24', "fit_func" = c('logistic'), "mu" = c(2e-05)),
    wildcards = list(),
    threads = 1,
    log = list(),
    resources = list('mem_mb', 'disk_mb', 'tmpdir', 'log_loc', 'mfree', 'cluster_time', 'disk_free', 'runtime', 'n_cores', "mem_mb" = '3G', "disk_mb" = '4G', "tmpdir" = '/tmp/12666795.1.feder-short.q', "log_loc" = './log_scr/', "mfree" = '3G', "cluster_time" = '0:00:10:00', "disk_free" = '4G', "runtime" = 10, "n_cores" = '1'),
    config = list(),
    rule = 'fig2B',
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
################################## Figure 2 B ##################################

## Reason:
## This is data from Collett et al (2017). It shows the results of the pocapavir
## clinical trial.

#### Set Up                                                                 ####

## Libraries
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

## Create this filepath
suppressWarnings(dir.create(paste0(filepath)))

## Read in data
coldat <- read.csv("dat/collett_trial.csv")

## Making the color map for plotting
color_in <- c("#c94d4d", "#688cc8")

names(color_in) <- c("Resistant", "Susceptible")

#### . . . . . . Analysis                                                   ####
## Calculating the number of resistant infections
colstat <- coldat %>%
  group_by(treatment, resistant) %>%
  summarize(prop = length(resistant))

#### . . . . . . Plotting as dot plots                                      ####

box_color <- "#b5b5b5"

tend_plot <- 
  ggplot(coldat, aes(
    x = treatment, 
    y = clearange,
  )) +
  annotate("rect", ymin = 15.5, ymax = 17.5, xmin = 0.5, xmax = 2.5, colour = "white", fill = box_color, alpha=0.2) +
  annotate("rect", ymin = 18.5, ymax = 21.5, xmin = 0.5, xmax = 2.5, colour = "white", fill = box_color, alpha=0.2) +
  annotate("rect", ymin = 22.5, ymax = 28.5, xmin = 0.5, xmax = 2.5, colour = "white", fill = box_color, alpha=0.2) +
  annotate("rect", ymin = 29.5, ymax = 42.5, xmin = 0.5, xmax = 2.5, colour = "white", fill = box_color, alpha=0.2) +
  # geom_boxplot(outliers = FALSE) +
  geom_dotplot(binaxis = "y", 
               stackdir = "center",
               # method="histodot",
               # stackdir = "centerwhole",
               binpositions = "all",
               stackgroups=TRUE,
               dotsize = 0.6,
               binwidth = 1,
               stroke = 0,
               mapping = aes(fill = resistant)
  ) +
  scale_fill_manual(values = c("#c94d4d", "#688cc8"),
  ) +
  scale_y_continuous(limits = c(0, 43.5),
                     breaks = c(seq(from = 0, to = 15, by = 5), 18, 22, 29, 43),
                     expand = c(0,0)
                     ) +
  scale_x_discrete(expand = c(0,0)) +
  labs(x = "Treatment", y = "Clearance Date (DPI)", fill = "Infection") +
  theme_classic() +
  theme(text = element_text(size=20), axis.text = element_text(size=18),
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "lightgrey"),
        panel.grid.minor.y = element_line(color = "lightgrey"),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        strip.text = element_text(size = 22, color = "black"),
        strip.background =element_rect(fill="white"),
        ## Below is needed to fix some wonky ggbreak changes
        axis.text.x.top = element_blank(),
        axis.line.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        ## For some reason, I can't change my legend position. GGbreak problem?
        legend.position = "inside",
        legend.position.inside = c(0.15, 0.1),
        legend.text = element_text(size=14),
        legend.title = element_text(size = 16),
        legend.background = element_rect(fill = "transparent")
        # legend.box.background = element_rect(fill = "transparent")
  ) 

tend_plot

ggsave(paste0(filepath, "res/fig2/B.png"), tend_plot, h = 8, w = 8, bg = "transparent")
ggsave(paste0(filepath, "res/fig2/B.jpg"), tend_plot, h = 8, w = 8, bg = "transparent")
ggsave(paste0(filepath, "res/fig2/B.svg"), tend_plot, h = 6.8, w = 6, bg = "transparent", device = 'svg')
saveRDS(tend_plot, paste0(filepath, "res/fig2/B.rds"))
