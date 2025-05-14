############################### SFigure 1 ######################################

## Reason:

## This will show the distribution of co-infected cells as a function of the MOI

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
width <-  as.numeric(snakemake@params[["sfig1_dim"]])[1]
height <- as.numeric(snakemake@params[["sfig1_dim"]])[2]

## Create the filepath where things will be saved
filepath <- paste0("runs/ddt_", init_fit_func, "_mu_", init_mu, "_", init_date, "/")

## geom_text() is in terms of mm, where element_text() is in terms of pt. I want
## my geom_text and axis text to be the same. Here, I am defining a conversion
## which will keep the text sizes the same between text types.
geom_text_conv = 0.3528

## Load library
library(ggplot2)

#### Generating data frames                                                 ####
## Generating x values
x <- seq(from = 0, to = 150, by = 1)
lam01 <- dpois(x = x, lambda = 1)
lam10 <- dpois(x = x, lambda = 10)
lam50 <- dpois(x = x, lambda = 50)
lam100 <- dpois(x = x, lambda = 100)
df01 <- data.frame(x, dens = lam01, MOI = "1")
df10 <- data.frame(x, dens = lam10, MOI = "10")
df50 <- data.frame(x, dens = lam50, MOI = "50")
df100 <- data.frame(x, dens = lam100, MOI = "100")
df_moi <- rbind(df01, df10, df50, df100)

df_moi$MOI <- factor(df_moi$MOI, levels = c(1, 10, 50, 100))


#### Plot                                                                   ####
moi_plot <- ggplot(df_moi, aes(x = x, y = dens, color = MOI)) + 
  geom_line(linewidth = 0.7) +
  xlab("Number of viruses per cell") +
  ylab("Probability") + labs(title = "") +
  theme_light() +
  scale_color_grey(start = 0.1, end = 0.9) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  # geom_vline(xintercept = 1, linetype = "dashed", color = "grey") +
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
    legend.position.inside = c(0.945, 0.87), ## Upper right
    # legend.position.inside = c(0.655, 0.07), ## bottom right
    legend.text = element_text(size = legend_text_size),
    # legend.title = element_blank(),
    legend.key.width = unit(0, 'cm'),
    legend.key.size = unit(0.0, "cm"),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.2),  # Fills legend background & adds a border
    legend.box.background = element_rect(color = "black", linewidth = 0),  # Adds an outer border
    legend.margin = margin(1, 1, 1, 2),
  ) 

suppressWarnings(dir.create(paste0(filepath, "res/")))
suppressWarnings(dir.create(paste0(filepath, "res/sup/")))
suppressWarnings(dir.create(paste0(filepath, "res/sup/sfig1")))

ggsave(paste0(filepath, "res/sup/sfig1/sfig1.png"), moi_plot, h = height, w = width, units = "cm", bg = "transparent", dpi = 1200)
ggsave(paste0(filepath, "res/sup/sfig1/sfig1.pdf"), moi_plot, h = height, w = width, units = "cm", bg = "transparent", dpi = 1200)
saveRDS(moi_plot, file = paste0(filepath, "res/sup/sfig1/sfig1.rds"))

