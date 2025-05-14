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

## Color scheme for resistant and susceptible
color_in <- snakemake@params[["rs_colors"]]

## Read in axis text size
axis_text_size <- as.numeric(snakemake@params[["axis_text_size"]])

## Legend text size
legend_text_size <- as.numeric(snakemake@params[["lege_text_size"]])

## This figure's size
width <-  as.numeric(snakemake@params[["fig3_dim"]])[1]
height <- as.numeric(snakemake@params[["fig3_dim"]])[2]

## Create the filepath where things will be saved
filepath <- paste0("runs/ddt_", init_fit_func, "_mu_", init_mu, "_", init_date, "/")


## geom_text() is in terms of mm, where element_text() is in terms of pt. I want
## my geom_text and axis text to be the same. Here, I am defining a conversion
## which will keep the text sizes the same between text types.
geom_text_conv = 0.3528

#### Read in the data                                                       ####
## Tanner data points
fit_df <- read.csv(file = paste0(filepath, "dat_gen/params/tan_fit_df.csv"))

## Logistic fit
fix_logistic <- read.csv(file = paste0(filepath, "dat_gen/params/logistic_fitness_function.csv"))

#### Visualize                                                              ####
func_plot <- ggplot(fit_df, aes(x = subunits, y = prob_surv)) +
  geom_line(data = fix_logistic, 
            aes(group = fit_type),
            linewidth = 0.7) +
  geom_point(color = "black",
             size = 0.8) +
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

suppressWarnings(dir.create(paste0(filepath, "res/")))
suppressWarnings(dir.create(paste0(filepath, "res/fig3/")))

ggsave(paste0(filepath, "res/fig3/B.png"), func_plot, h = height/2, w = width/2, units = "cm", bg = "transparent", dpi = 1200)
ggsave(paste0(filepath, "res/fig3/B.pdf"), func_plot, h = height/2, w = width/2, units = "cm", bg = "transparent", dpi = 1200)
ggsave(paste0(filepath, "res/fig3/B.svg"), func_plot, h = height/2, w = width/2, units = "cm", bg = "transparent", dpi = 1200, device = 'svg')
saveRDS(func_plot, file = paste0(filepath, "res/fig3/B.rds"))
