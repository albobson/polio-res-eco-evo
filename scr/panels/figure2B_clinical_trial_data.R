################################## Figure 2 B ##################################

## Reason:
## This is data from Collett et al (2017). It shows the results of the pocapavir
## clinical trial.

#### Set Up                                                                 ####

## Libraries
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
width <-  as.numeric(snakemake@params[["fig2_dim"]])[1]
height <- as.numeric(snakemake@params[["fig2_dim"]])[2]

## Create the filepath where things will be saved
filepath <- paste0("runs/ddt_", init_fit_func, "_mu_", init_mu, "_", init_date, "/")

## Create this filepath
suppressWarnings(dir.create(paste0(filepath)))

## Read in data
coldat <- read.csv("dat/collett_trial.csv")

names(color_in) <- c("Res.", "Sus.")

## Figure 2A is put together kind of strangely. For that reason, we have to
## slightly rescale this plot in order to get its bounds to match up to the
## figure 2A's bounds:

## geom_text() is in terms of mm, where element_text() is in terms of pt. I want
## my geom_text and axis text to be the same. Here, I am defining a conversion
## which will keep the text sizes the same between text types.
geom_text_conv = 0.3528

#### . . . . . . Analysis                                                   ####
## Calculating the number of resistant infections
colstat <- coldat %>%
  group_by(treatment, resistant) %>%
  summarize(prop = length(resistant)) 

coldat$resistant <- ifelse(coldat$resistant == "Resistant", "Res.", "Sus.")
#### . . . . . . Plotting as dot plots                                      ####

box_color <- "#b5b5b5"

## Creating a dynamic size for the dot plot
dotsize = (max(coldat$clearange)-min(coldat$clearange))/max(coldat$clearange)

## Test with formatting similar to plot 2A
tend_plot_only <- 
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
               dotsize = dotsize,
               binwidth = 1,
               stroke = NA,
               key_glyph = "point",
               mapping = aes(fill = resistant)
  ) +
  scale_fill_manual(values = color_in,
  ) +
  guides(fill = guide_legend(override.aes = list(color = color_in,
                                                 size = 1))) +
  scale_y_continuous(limits = c(0, 43.5),
                     breaks = c(seq(from = 0, to = 15, by = 5), 18, 22, 29, 43),
                     expand = c(0,0)
  ) +
  scale_x_discrete(expand = c(0,0)) +
  labs(x = "Treatment", y = "Measured Clearance Time (DPI)", fill = "Infection") +
  theme_classic() +
  theme(text = element_text(size= axis_text_size), 
        axis.text = element_text(size= axis_text_size),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "lightgrey", linewidth = 0.2),
        panel.grid.minor.y = element_blank(), #element_line(color = "lightgrey"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(colour="black"),
        axis.title.y = element_text(),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.background = element_rect(fill='transparent'), 
        plot.background = element_rect(fill='transparent', color=NA),
        ## Legend stuff
        legend.position = "inside",
        legend.position.inside = c(0.84, 0.83),
        legend.text = element_text(size = legend_text_size),
        legend.title = element_blank(),
        # legend.background = element_rect(fill = "transparent"),
        legend.background = element_rect(fill = "white", color = "black", linewidth = 0.2),  # Fills legend background & adds a border
        legend.box.background = element_rect(color = "white", linewidth = 0.2),  # Adds an outer border
        legend.key = element_rect(fill = "transparent"),
        legend.margin = margin(0.6, 0.6, 0.6, 0.6),
        legend.key.width = unit(-00, 'cm'),
        legend.key.size = unit(0.2, "cm")
  ) 

## Plotting clinical trial title
clin_group_df <- coldat %>%
  select(treatment) %>%
  mutate(cond = 0) %>%
  unique()

clin_group_df <- data.frame(cond = c(0,0), treatment = c("Placebo\n(n = 48)", "Pocapavir\n(n = 93)"))

clin_group_p <- ggplot(clin_group_df, aes(cond, cond, label = treatment)) +
  geom_text(size = axis_text_size * geom_text_conv) +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  facet_grid(~treatment) +
  theme(panel.grid = element_blank(), axis.text.x = element_blank(),
        text = element_text(size= axis_text_size * geom_text_conv),
        strip.text = element_blank(),
        axis.text.y=element_blank(),
        strip.background = element_blank())


full_clin <- plot_grid(NULL, NULL, tend_plot_only, NULL, clin_group_p, ## Adding NULL to move plots closer together
                       rel_heights = c(0.3, -0.15, 1.2, -0.1, 0.3),
                            ncol = 1, axis = "l", align = "v")


full_clin_list <- list(tend_plot_only, clin_group_p)


suppressWarnings(dir.create(paste0(filepath, "res/")))

ggsave(paste0(filepath, "res/fig2/B.png"), full_clin, h = height, w = width/2, units = "cm", bg = "transparent", dpi = 1200)
ggsave(paste0(filepath, "res/fig2/B.pdf"), full_clin, h = height, w = width/2, units = "cm", bg = "transparent", dpi = 1200)
ggsave(paste0(filepath, "res/fig2/B.svg"), full_clin, h = height, w = width/2, units = "cm", bg = "transparent", device = 'svg')
saveRDS(full_clin_list, paste0(filepath, "res/fig2/B.rds"))


# 
# clin_group_df <- data.frame(x = c(0,1), y = c(1,1), treatment = c("Placebo", "Pocapavir"))
# 
# clin_num_df <- data.frame(x = c(0,1), y = c(0,0), treatment = c("(n = 48)", "(n = 93)"))
# 
# clin_group_p <- ggplot(mapping = aes(x = x, y = y)) +
#   geom_text(data = clin_group_df,
#             size = axis_text_size * geom_text_conv,
#             mapping = aes(label = treatment)
#   ) +
#   geom_text(data = clin_num_df,
#             size = axis_text_size * geom_text_conv,
#             mapping = aes(label = treatment)
#   ) +
#   labs(x = NULL, y = NULL) +
#   theme_minimal() +
#   theme(panel.grid = element_blank(), axis.text.y = element_blank(),
#         axis.text.x = element_blank(),
#         text = element_text(size=axis_text_size),
#         legend.position = "none",
#         strip.text = element_blank(),
#         strip.background = element_blank())
# 
# clin_group_p
# 
# full_clin <- plot_grid(NULL, NULL, tend_plot_only, NULL, clin_group, ## Adding NULL to move plots closer together
#                        rel_heights = c(0.3, -0.1, 1.2, -0.05, 0.3),
#                        ncol = 1, axis = "l", align = "v")
# 
# 
# full_clin_list <- list(tend_plot_only, clin_group_p)