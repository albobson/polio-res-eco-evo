
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
    output = list('runs/ddt_logistic_mu_2e-05_2025-03-24/res/fig5/D.rds'),
    params = list('2025-03-24', c('logistic'), c(2e-05), "date" = '2025-03-24', "fit_func" = c('logistic'), "mu" = c(2e-05)),
    wildcards = list(),
    threads = 1,
    log = list(),
    resources = list('mem_mb', 'disk_mb', 'tmpdir', 'log_loc', 'mfree', 'cluster_time', 'disk_free', 'runtime', 'n_cores', "mem_mb" = '3G', "disk_mb" = '4G', "tmpdir" = '/tmp/12724629.1.feder-short.q', "log_loc" = './log_scr/', "mfree" = '3G', "cluster_time" = '0:00:10:00', "disk_free" = '4G', "runtime" = 10, "n_cores" = '1'),
    config = list(),
    rule = 'fig5D',
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
############################## Figure 5 D ######################################

## Reason:

## This will plot the range of dominance that will be assessed.

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
fit_func$fit_type <- "Pocapavir"

optim_params <- read.csv(paste0(filepath, "dat_gen/params/optim_params.csv"))

## Source functions to run simulations
source("scr/polv_DDT_functions.R")  ## Main functions 


## Making the color map for plotting
color_in <- c("#c94d4d", "#688cc8")

## Create colors for the funcs
custom_colors <- c("Pocapavir" = "black", # black
                   "Full Dominant" = "darkgrey") # Saddle Brown

## geom_text() is in terms of mm, where element_text() is in terms of pt. I want
## my geom_text and axis text to be the same. Here, I am defining a conversion
## which will keep the text sizes the same between text types.
geom_text_conv = 0.3528

## I also want my axis text to be 12 pt font, so, setting that here
axis_text_size = 12

## I would like for the legend to be slightly smaller than the rest of the text.
## I'm setting that size here
legend_text_size = 8

#### Set up less stringent drugs                                            ####
## Set the range of dominance
## Define the range for which we will scale the function
min_scale <- 0
max_scale <- 60

scale_range <- seq(min_scale, max_scale, length.out = max_scale+1)

## Storage df
scale_store <- NULL

## generate functions across scale range
for(i in 1:length(scale_range)) {
  temp_scale <- NULL
  temp_df <- NULL
  
  temp_scale <- fit_func_logistic(k = 100, 
                                  min_val = min(fit_func$prob_surv),
                                  subunits = 60, mid = i)
  
  ## Set our minimum and maximum values to the values that we want
  temp_scale[1] <- min(fit_func$prob_surv)
  temp_scale[length(temp_scale)] <- max(fit_func$prob_surv)
  
  temp_df <- data.frame(subunits = 0:60, prob_surv = temp_scale, fit_type = scale_range[i])
  
  scale_store <- rbind(scale_store, temp_df)
}

scale_store


#### Plot                                                                   ####
## Define color range
num_lines <- length(scale_range)

gray_range <- gray(seq(0, 0.8, length.out = num_lines))
gray_range


#### Visualize                                                              ####
func_plot <- ggplot() +
  geom_line(data = scale_store,
            mapping = aes(group = fit_type,
                          x = subunits, 
                          y = prob_surv, 
                          color = factor(fit_type)
            ),
            linetype = "solid",
            linewidth = 1
  ) +
  scale_linetype(guide = "none") +
  xlab("Resistant Subunits") +
  ylab("Probability of Survival") +
  scale_y_continuous(trans="log10",
                     # limits = c(0.0001,1),
                     breaks=10^(-4:0),
                     labels = c(sapply(-4:-1,function(i){parse(text = sprintf("10^%d",i))}), "1 ")) +
  theme_light() +
  scale_color_manual(
    values = gray_range) +
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
    # legend.position.inside = c(0.35, 0.94), ## Upper left
    # legend.position.inside = c(0.65, 0.06), ## bottom right
    # legend.text = element_text(size = legend_text_size),
    # legend.title = element_blank(),
    # legend.background = element_rect(fill = "transparent"),
    # legend.key = element_rect(fill = "transparent"),
    # legend.key.width = unit(0, 'cm'),
    # legend.key.size = unit(0.3, "cm")
  ) 

func_plot

suppressWarnings(dir.create(paste0(filepath, "res/")))
suppressWarnings(dir.create(paste0(filepath, "res/fig5/")))

ggsave(paste0(filepath, "res/fig5/D.png"), func_plot, h = 22/3, w = 18/3, unit = "cm", bg = "transparent", dpi = 300)
ggsave(paste0(filepath, "res/fig5/D.pdf"), func_plot, h = 22/3, w = 18/3, unit = "cm", bg = "transparent", dpi = 300)
ggsave(paste0(filepath, "res/fig5/D.svg"), func_plot, h = 22/3, w = 18/3, unit = "cm", bg = "transparent", dpi = 300)
saveRDS(func_plot, file = paste0(filepath, "res/fig5/D.rds"))



