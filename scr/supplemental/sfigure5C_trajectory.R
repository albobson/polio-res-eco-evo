############################## Figure 5 C ######################################

## Reason:

## This will showcase how resistance can emerge *due* to pocapavir's high
## effectiveness in a single trajectory plot

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

## Create colors for the funcs
custom_colors <- c("Pocapavir" = "black", # black
                   "Full Dom." = "darkgrey") # Saddle Brown

## geom_text() is in terms of mm, where element_text() is in terms of pt. I want
## my geom_text and axis text to be the same. Here, I am defining a conversion
## which will keep the text sizes the same between text types.
geom_text_conv = 0.3528

#### Set up less stringent drugs                                            ####
high_dom <- data.frame(subunits = 0:60,
                       prob_surv = fit_func_logistic(k = 100, 
                                                     min_val = min(fit_func$prob_surv),
                                                     subunits = 60, mid = 59.5),
                       fit_type = "Full Dom."
) 


#### Run simulations                                                        ####
pocap_run <- determ_polv(n = 6,
                         moi_mut_start = 100*0.0001, 
                         moi_wt_start = 100*(1-0.0001),
                         fit_func_in = fit_func$prob_surv,
                         v_prog = optim_params$optim_v_prog,
                         p2pfu = optim_params$optim_p2pfu
) %>%
  mutate(fit_type = "Pocapavir")

high_dom_run <- determ_polv(n = 6,
                            moi_mut_start = 100*0.0001, 
                            moi_wt_start = 100*(1-0.0001),
                            fit_func_in = high_dom$prob_surv,
                            v_prog = optim_params$optim_v_prog,
                            p2pfu = optim_params$optim_p2pfu
) %>%
  mutate(fit_type = "Full Dom.")

tot_df <- rbind(pocap_run, high_dom_run)

tot_df$fit_type <- factor(tot_df$fit_type, levels = c("Pocapavir", "Full Dom."))

#### Plotting trajectory                                                    ####
plot5c <- ggplot(tot_df, aes(x = time, y = moi_type, color = type, linetype = fit_type)) +
  # geom_point(size = 3) + 
  geom_line(linewidth = 0.7) +
  scale_color_manual(labels = c("Res.", "Sus."),
                     values = c("#c94d4d", "#688cc8"),
                     name = "Genotype") +
  theme_light() + 
  # scale_linetype(guide = "none") +
  xlab("Passages") + 
  ylab("Viral Density (Genomes/cell)") + 
  scale_y_log10(breaks=10^(-6:10),
                labels = sapply(-6:10,function(i){parse(text = sprintf("10^%d",i))}),
                limits = c(10^-5, 10^2.4)) +
  scale_x_continuous(breaks = c(0:max(tot_df$time))) +
  guides(
    color = "none"
  ) +
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
    strip.text.x = element_blank(),
    ## Remove the background of the plot and panel
    panel.background = element_rect(fill='transparent'), 
    plot.background = element_rect(fill='transparent', color=NA),
    ## Legend stuff
    # legend.position = "none",
    legend.position = "inside",
    legend.position.inside = c(0.68, 0.1),
    legend.text = element_text(size = legend_text_size),
    legend.title = element_blank(),
    # legend.key.width = unit(0, 'cm'),
    legend.key.size = unit(0.3, "cm"),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.2),  # Fills legend background & adds a border
    legend.box.background = element_rect(color = "black", linewidth = 0),  # Adds an outer border
    legend.margin = margin(0.6, 0.6, 0.6, 2),
    legend.spacing.y = unit(-0.0, "cm"),
    # legend.box = "horizontal"
  ) 

suppressWarnings(dir.create(paste0(filepath, "res/")))
suppressWarnings(dir.create(paste0(filepath, "res/sup/")))
suppressWarnings(dir.create(paste0(filepath, "res/sup/sfig5/")))

ggsave(paste0(filepath, "res/sup/sfig5/C.png"), plot5c, h = height/2, w = width/3, units = "cm", bg = "transparent", dpi = 1200)
ggsave(paste0(filepath, "res/sup/sfig5/C.pdf"), plot5c, h = height/2, w = width/3, units = "cm", bg = "transparent", dpi = 1200)
ggsave(paste0(filepath, "res/sup/sfig5/C.svg"), plot5c, h = height/2, w = width/3, units = "cm", bg = "transparent", dpi = 1200, device = "svg")
saveRDS(plot5c, paste0(filepath, "res/sup/sfig5/C.rds"))

