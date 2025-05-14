############################### SFigure 3 ######################################

## Reason:

## This will show the distribution of progeny phenotypes after each round of rep.

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
# axis_text_size <- as.numeric(snakemake@params[["axis_text_size"]])
axis_text_size <- 6

## Legend text size
# legend_text_size <- as.numeric(snakemake@params[["lege_text_size"]])
legend_text_size = 6

## This figure's size
width <-  as.numeric(snakemake@params[["sfig3_dim"]])[1]
height <- as.numeric(snakemake@params[["sfig3_dim"]])[2]

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



#### Running simulation                                                     ####
## Res freq. same as seen in unexposed cultures
moi100_df0001 <- determ_polv(n = 6, 
                             moi_mut_start = 100*0.0001, 
                             moi_wt_start = 100*(1-0.0001),
                             fit_func_in = fit_func$prob_surv,
                             v_prog = optim_params$optim_v_prog,
                             p2pfu = optim_params$optim_p2pfu,
                             report_subunit_dist = TRUE) 

moi100_df0001 <- moi100_df0001 %>%
  filter(type != "total")

moi100_df0001$type <- factor(moi100_df0001$type, levels = c("susceptible", "resistant"))

moi100_df0001 <- moi100_df0001 %>%
  filter(time != 0) %>%
  mutate(time = time - 1)

## Plotting
moi100 <- ggplot(moi100_df0001, aes(x = subs, y = dens_b4_pocap, fill = type)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(rows = vars(time), labeller = labeller(time = function(x) paste0("t = ", x))) +
  theme_light() +
  xlab("Res. Capsid Subunits") +
  ylab("Density") +
  scale_fill_manual(labels = c("Sus.", "Res."),
                     values = rev(color_in),
                     name = "Genotype") +
  scale_y_continuous(breaks = c(0, 0.5, 1), 
                     limits = c(-0.1, 1.1),
                     expand = c(0,0)
                     ) +
  xlim(c(-1, 61)) +
  # ggtitle("Capsid phenotype distribution by genotype") +
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
    strip.text = element_text(color = "black", size= axis_text_size),
    ## Remove the background of the plot and panel
    panel.background = element_rect(fill='transparent'), 
    plot.background = element_rect(fill='transparent', color=NA),
    ## Legend stuff
    # legend.position = "none",
    legend.position = "inside",
    legend.position.inside = c(0.1, 0.05),
    legend.text = element_text(size = legend_text_size),
    legend.title = element_blank(),
    # legend.key.width = unit(0, 'cm'),
    legend.key.size = unit(0.3, "cm"),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.2),  # Fills legend background & adds a border
    legend.box.background = element_rect(color = "black", linewidth = 0),  # Adds an outer border
    legend.margin = margin(0.6, 0.6, 0.6, 0.6),
    legend.spacing.y = unit(-0.0, "cm"),
    # legend.box = "horizontal"
  ) 


suppressWarnings(dir.create(paste0(filepath, "res/")))
suppressWarnings(dir.create(paste0(filepath, "res/sup/")))
suppressWarnings(dir.create(paste0(filepath, "res/sup/sfig3")))

ggsave(paste0(filepath, "res/sup/sfig3/sfig3.png"), moi100, h = height, w = width, units = "cm", bg = "transparent", dpi = 1200)
ggsave(paste0(filepath, "res/sup/sfig3/sfig3.pdf"), moi100, h = height, w = width, units = "cm", bg = "transparent", dpi = 1200)
saveRDS(moi100, file = paste0(filepath, "res/sup/sfig3/sfig3.rds"))

