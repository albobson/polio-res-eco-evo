############################## Figure 5 D ######################################

## Reason:

## This will plot the range of worse drugs that will be assessed.

#### Set Up                                                                 ####

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
width <-  as.numeric(snakemake@params[["fig5_dim"]])[1]
height <- as.numeric(snakemake@params[["fig5_dim"]])[2]

## Find how many variations we want to run 
n_trials <- as.numeric(snakemake@params[["n_trials"]])

## Create the filepath where things will be saved
filepath <- paste0("runs/ddt_", init_fit_func, "_mu_", init_mu, "_", init_date, "/")

## Read in fitness function and parameters
fit_func <- read.csv(paste0(filepath, "dat_gen/params/logistic_fitness_function.csv"))
optim_params <- read.csv(paste0(filepath, "dat_gen/params/optim_params.csv"))

## Source functions to run simulations
source("scr/polv_DDT_functions.R")  ## Main functions 

## geom_text() is in terms of mm, where element_text() is in terms of pt. I want
## my geom_text and axis text to be the same. Here, I am defining a conversion
## which will keep the text sizes the same between text types.
geom_text_conv = 0.3528

#### Creating other functions                                               ####
## Generate a function that is less stringent than pocapavir
scaled_vector <- function(og_probs, scaler) {
  ## New vector for storing
  new_vec <- NULL
  
  new_vec <- ((1- scaler * og_probs[1])/(1-og_probs[1])) * (og_probs - og_probs[1]) + scaler * og_probs[1]
  
  return(new_vec)
}

# ## Define the range for which we will scale the function
# min_scale <- 1
# max_scale <- 1000
# 
# scale_range <- 10^log10(seq(min_scale, max_scale, length.out = n_trials))

min_scale <- 0
max_scale <- n_trials-1

scale_range <- 10^c(min_scale:max_scale)

## Storage df
scale_store <- NULL
## generate functions across scale range
for(i in 1:length(scale_range)) {
  temp_scale <- NULL
  temp_df <- NULL
  
  temp_scale <- scaled_vector(og_probs = fit_func$prob_surv, scaler = scale_range[i])
  
  temp_df <- data.frame(subunits = 0:60, prob_surv = temp_scale, fit_type = scale_range[i])
  
  scale_store <- rbind(scale_store, temp_df)
}

#### Plot                                                                   ####
## Define color range
num_lines <- length(scale_range)

gray_range <- gray(seq(0, 0.8, length.out = num_lines))


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
  scale_color_manual(name = "Fold Weaker",
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
    # legend.position = "none",
    legend.position = "inside",
    # legend.position.inside = c(0.35, 0.94), ## Upper left
    legend.position.inside = c(0.73, 0.19), ## bottom right
    legend.text = element_text(size = legend_text_size),
    # legend.title = element_blank(),
    # legend.background = element_rect(fill = "transparent"),
    legend.key = element_rect(fill = "transparent"),
    legend.key.width = unit(0, 'cm'),
    legend.key.height = unit(0, "cm"),
    legend.key.size = unit(0.0, "cm"),
    legend.background = element_rect(fill = "white", color = "black"),  # Fills legend background & adds a border
    legend.box.background = element_rect(color = "black", linewidth = 0),  # Adds an outer border
    legend.margin = margin(1, 1, 1, 2),
  ) 

suppressWarnings(dir.create(paste0(filepath, "res/")))
suppressWarnings(dir.create(paste0(filepath, "res/fig5/")))

ggsave(paste0(filepath, "res/fig5/D.png"), func_plot, h = height/2, w = width/3, units = "cm", bg = "transparent", dpi = 1200)
ggsave(paste0(filepath, "res/fig5/D.pdf"), func_plot, h = height/2, w = width/3, units = "cm", bg = "transparent", dpi = 1200)
ggsave(paste0(filepath, "res/fig5/D.svg"), func_plot, h = height/2, w = width/3, units = "cm", bg = "transparent", dpi = 1200, device = 'svg')
saveRDS(func_plot, file = paste0(filepath, "res/fig5/D.rds"))



