############################### Sup. Fig 4 #####################################

## Reason:

## This will plot the full trajectories of the clinical trial runs. There will
## be two plots: early clearers and late clearers

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
width <-  as.numeric(snakemake@params[["sfig4_dim"]])[1]
height <- as.numeric(snakemake@params[["sfig4_dim"]])[2]

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

## Read in Collett data
coldat <- read.csv("dat/collett_trial.csv")

## Read in clincial trial simulations
trial <- read.csv(paste0(filepath, "dat_gen/sims/sim_clinical_trial.csv"))


#### Data cleaning                                                          ####
## Find the earliest placebo from the clinical trial
coldat_p <- coldat %>% group_by(treatment) %>% 
  mutate(res = as.character(resistant), 
         clearange = as.numeric(clearange),
         source = "Clinical Trial") %>% 
  select(clearange, res, treatment, res) %>% 
  na.omit() %>% 
  # filter(treatment == "Pocapavir (n = 93)") %>%
  mutate(time = clearange, 
         treatment = ifelse(treatment == "Pocapavir (n = 93)", "Pocapavir", "Placebo"),
         # treatment = paste0(treatment, " ", res) ## if you want treatment to be stratified by resistance class
  ) %>%
  select(time, res, treatment) %>%
  mutate(shape = "circle") ## Adding a shape value for distinguising our own data

early_placebo <- coldat_p %>%
  filter(treatment == "Placebo") %>%
  summarize(mint = min(time))

early_placebo <- early_placebo$mint[1]


## Separate simulations into early vs late clearers
simsclean <- trial %>% 
  group_by(id) %>% 
  mutate(earlylate = ifelse(max(time) >= early_placebo, "Late", "Early"))

## A function to simplify the y axis:
evens <- function(x) subset(x, x %% 2 == 0)


sup_plot4 <- ggplot(simsclean, aes(x = time, y = moi_type, color = type, group = interaction(id, type))) +
  geom_point(size = 0.7, alpha = 0.4) + 
  geom_line(linewidth = 0.7, alpha = 0.4) +
  scale_color_manual(labels = c("Res.", "Sus."),
                    values = color_in,
                    name = "Genotype") +
  theme_light() + 
  xlab("Days post infection (DPI)") + 
  ylab("Viral Density (Genomes/cell)") + 
  xlim(c(0, round(max(simsclean$time), 0))) +
  scale_y_log10(breaks=10^(evens(-6:10)),
                labels = sapply(evens(-6:10), function(i){parse(text = sprintf("10^%d",i))})) +
  # scale_x_continuous(breaks = c(0:max(simsclean$time))) +
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
    panel.grid.minor.x = element_blank(),                   ## Minor x lines
    panel.grid.major.y = element_line(color = "lightgrey"), ## Major y lines
    panel.grid.minor.y = element_blank(), ## Minor y lines
    ## For faceted plots: The strip is the top of the facet
    strip.background = element_blank(),
    strip.text = element_text(color = "black"),
    ## Remove the background of the plot and panel
    panel.background = element_rect(fill='transparent'), 
    plot.background = element_rect(fill='transparent', color=NA),
    ## Legend stuff
    # legend.position = "none",
    legend.position = "inside",
    legend.position.inside = c(0.93, 0.945),
    legend.text = element_text(size = legend_text_size),
    legend.title = element_blank(),
    # legend.key.width = unit(0, 'cm'),
    legend.key.size = unit(0.3, "cm"),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.2),  # Fills legend background & adds a border
    legend.box.background = element_rect(color = "black", linewidth = 0),  # Adds an outer border
    legend.margin = margin(0.6, 0.6, 0.6, 0.6),
    legend.spacing.y = unit(-0.0, "cm"),
    # legend.box = "horizontal"
  )  +
  facet_grid(rows = vars(earlylate))

suppressWarnings(dir.create(paste0(filepath, "res/")))
suppressWarnings(dir.create(paste0(filepath, "res/sup/")))
suppressWarnings(dir.create(paste0(filepath, "res/sup/sfig4")))

## Save
ggsave(paste0(filepath, "res/sup/sfig4/sfig4.png"), sup_plot4, h = height, w = width, units = "cm", bg = "transparent", dpi = 1200)
ggsave(paste0(filepath, "res/sup/sfig4/sfig4.pdf"), sup_plot4, h = height, w = width, units = "cm", bg = "transparent", dpi = 1200)
saveRDS(sup_plot4, file = paste0(filepath, "res/sup/sfig4/sfig4.rds"))

