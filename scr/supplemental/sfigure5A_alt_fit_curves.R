############################## sFigure 5 A ######################################

## Reason:

## This will plot the more dominant alongside pocapavir

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
width <-  as.numeric(snakemake@params[["sfig5_dim"]])[1]
height <- as.numeric(snakemake@params[["sfig5_dim"]])[2]

## Create the filepath where things will be saved
filepath <- paste0("runs/ddt_", init_fit_func, "_mu_", init_mu, "_", init_date, "/")

## Read in fitness function and parameters
fit_func <- read.csv(paste0(filepath, "dat_gen/params/logistic_fitness_function.csv"))
fit_func$fit_type <- "Pocapavir"

optim_params <- read.csv(paste0(filepath, "dat_gen/params/optim_params.csv"))

## Source functions to run simulations
source("scr/polv_DDT_functions.R")  ## Main functions 


## geom_text() is in terms of mm, where element_text() is in terms of pt. I want
## my geom_text and axis text to be the same. Here, I am defining a conversion
## which will keep the text sizes the same between text types.
geom_text_conv = 0.3528

#### Creating other functions                                               ####
## Create a fully dominant fitness curve
high_dom <- data.frame(subunits = 0:60,
                       prob_surv = fit_func_logistic(k = 100, 
                                                     min_val = min(fit_func$prob_surv),
                                                     subunits = 60, mid = 59.5),
                       fit_type = "Full Dom."
) 

all_funcs <- rbind(fit_func, high_dom)

## Create colors for them
custom_colors <- c("Pocapavir" = "black", # black
                   "Full Dom." = "darkgrey") # Saddle Brown

#### Visualize                                                              ####
func_plot <- ggplot() +
  geom_line(data = fit_func,
            mapping = aes(group = fit_type,
                          x = subunits, 
                          y = prob_surv, 
                          color = fit_type
            ),
            linetype = "solid",
            linewidth = 0.7
  ) +
  geom_line(data = high_dom,
            mapping = aes(group = fit_type,
                          x = subunits, 
                          y = prob_surv, 
                          color = fit_type, 
            ),
            linetype = "longdash",
            linewidth = 0.7
  ) +
  xlab("Res. Subunit Composition") +
  ylab("Probability of Survival") +
  scale_y_continuous(trans="log10",
                     limits = c(0.0003,1),
                     breaks=10^(-4:0),
                     labels = c(sapply(-4:-1,function(i){parse(text = sprintf("10^%d",i))}), "1 ")) +
  theme_light() +
  scale_color_manual(values = custom_colors) +
  # labs(color = "") +
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
    legend.position.inside = c(0.27, 0.92), ## Upper left
    # legend.position.inside = c(0.65, 0.06), ## bottom right
    legend.text = element_text(size = legend_text_size),
    legend.title = element_blank(),
    legend.key.width = unit(0, 'cm'),
    legend.key.size = unit(0.0, "cm"),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.2),  # Fills legend background & adds a border
    legend.box.background = element_rect(color = "black", linewidth = 0),  # Adds an outer border
    legend.margin = margin(0.6, 0.6, 0.6, 2),
  ) 

suppressWarnings(dir.create(paste0(filepath, "res/")))
suppressWarnings(dir.create(paste0(filepath, "res/sup/")))
suppressWarnings(dir.create(paste0(filepath, "res/sup/sfig5/")))

ggsave(paste0(filepath, "res/sup/sfig5/A.png"), func_plot, h = height/2, w = width/3, units = "cm", bg = "transparent", dpi = 1200)
ggsave(paste0(filepath, "res/sup/sfig5/A.pdf"), func_plot, h = height/2, w = width/3, units = "cm", bg = "transparent", dpi = 1200)
ggsave(paste0(filepath, "res/sup/sfig5/A.svg"), func_plot, h = height/2, w = width/3, units = "cm", bg = "transparent", dpi = 1200, device = 'svg')
saveRDS(func_plot, file = paste0(filepath, "res/sup/sfig5/A.rds"))

