############################### Figure 4 H #####################################

## Reason:

## We want to plot an early clearer from our simulated clinical trial

#### Set Up                                                                 ####
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
width <-  as.numeric(snakemake@params[["fig4_dim"]])[1]
height <- as.numeric(snakemake@params[["fig4_dim"]])[2]

## Create the filepath where things will be saved
filepath <- paste0("runs/ddt_", init_fit_func, "_mu_", init_mu, "_", init_date, "/")

## Source functions to run simulations
source("scr/polv_DDT_functions.R")  ## Main functions 

## Read in Collett data
coldat <- read.csv("dat/collett_trial.csv")

## Read in fitness function and parameters
fit_func <- read.csv(paste0(filepath, "dat_gen/params/logistic_fitness_function.csv"))
optim_params <- read.csv(paste0(filepath, "dat_gen/params/optim_params.csv"))

## Read in clincial trial simulations
trial <- read.csv(paste0(filepath, "dat_gen/sims/sim_clinical_trial.csv"))

## A function which will be used that only selects even numbers
evens <- function(x) subset(x, x %% 2 == 0)

## geom_text() is in terms of mm, where element_text() is in terms of pt. I want
## my geom_text and axis text to be the same. Here, I am defining a conversion
## which will keep the text sizes the same between text types.
geom_text_conv = 0.3528

#### Data Cleaning                                                          ####
## Create directory to save stuff (if not already created)
suppressWarnings(dir.create(paste0(filepath, "res/fig4")))

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
  geom_line(linewidth = 0.7) +
  scale_color_manual(labels = c("Res.", "Sus."),
                     values = color_in,
                     name = "Genotype") +
  theme_light() + 
  xlab("Replications") + 
  ylab("Genomes/cell") + 
  scale_x_continuous(limits = c(0, max(early$max_date))) +
  scale_y_log10(breaks=10^(evens(-6:10)),
                labels = sapply(evens(-6:10), function(i){parse(text = sprintf("10^%d",i))}),
                limits = c(10^-4, 10^2.5)) +
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
    legend.position.inside = c(0.91,0.8),
    legend.text = element_text(size = legend_text_size),
    legend.title = element_blank(),
    # legend.background = element_rect(fill = "transparent"),
    # legend.key = element_rect(fill = "transparent"),
    legend.key.width = unit(0, 'cm'),
    legend.key.size = unit(0.0, "cm"),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.2),  # Fills legend background & adds a border
    legend.box.background = element_rect(color = "black", linewidth = 0),  # Adds an outer border
    legend.margin = margin(0.6, 0.6, 0.6, 2),
  ) 

## text plot for "Early clearer" 
early_clear <- data.frame(cond = 0, name = "Early Clearer")

## Plotting just "Early Clearer"
early_clear_p <- ggplot(early_clear, aes(cond, name, label = name)) +
  geom_text(size = axis_text_size * geom_text_conv) +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.text.x = element_blank(),
        text = element_text(size= axis_text_size),
        strip.text = element_blank(),
        axis.text.y=element_blank(),
        strip.background = element_blank())

## Full plot
early_p <- plot_grid(early_clear_p, NULL, early_sim_p,
                     rel_heights = c(0.4, -0.2, 1),
                    ncol = 1, axis = "l", align = "v")

## text plot for x axis labels
reps_df <- data.frame(cond = 0, name = "Replications")

## Plotting just "Early Clearer"
reps_p <- ggplot(reps_df, aes(cond, name, label = name)) +
  geom_text(size = axis_text_size * geom_text_conv) +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.text.x = element_blank(),
        text = element_text(size= axis_text_size),
        strip.text = element_blank(),
        axis.text.y=element_blank(),
        strip.background = element_blank())

early_list <- list(early_clear_p, early_sim_p, reps_p)


## Save
ggsave(early_p, file = paste0(filepath, "res/fig4/H.png"), h = height/2/2, w = width/3, units = "cm", bg = "transparent", dpi = 1200)
ggsave(early_p, file = paste0(filepath, "res/fig4/H.pdf"), h = height/2/2, w = width/3, units = "cm", bg = "transparent", dpi = 1200)
ggsave(early_p, file = paste0(filepath, "res/fig4/H.svg"), h = height/2/2, w = width/3, units = "cm", bg = "transparent", dpi = 1200, device = 'svg')
saveRDS(early_list, paste0(filepath, "res/fig4/H.rds"))

