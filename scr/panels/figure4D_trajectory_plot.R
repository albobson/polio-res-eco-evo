############################### Figure 4 D #####################################

## Reason:

## Over time, the viral population goes from a high MOI (where resistance is
## suppressed) to a low MOI (where resistance is expressed). At a low MOI,
## resistance can increase in frequency and become the dominant genotype.

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
width <-  as.numeric(snakemake@params[["fig4_dim"]])[1]
height <- as.numeric(snakemake@params[["fig4_dim"]])[2]

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

#### Run simulations                                                        ####
curr_run <- determ_polv(n = 6,
                        moi_mut_start = 100*0.0001, 
                        moi_wt_start = 100*(1-0.0001),
                        fit_func_in = fit_func$prob_surv,
                        v_prog = optim_params$optim_v_prog,
                        p2pfu = optim_params$optim_p2pfu
) 

## Cleaning the data for plotting on phase plane
phase_plane_var <- curr_run %>%
  select(type, time, moi_res, moi_wt, pop_prop) %>%
  filter(type == "resistant") %>%
  mutate(moi_tot = (moi_res + moi_wt), pop_prop = as.numeric(pop_prop)) %>%
  arrange(time) %>%
  mutate(delta_p = (lead(pop_prop, default = last(pop_prop)) - pop_prop),
         delta_moi = (lead(moi_tot, default = last(moi_tot)) - moi_tot),
         dist_col = abs((delta_p)+(delta_moi)/100),
         moi_color = factor(round(moi_tot,5), 
                            levels = unique(as.numeric(round(moi_tot,5))))) %>%
  ungroup()

## Saving this data to be used for the phase plane
suppressWarnings(dir.create(paste0(filepath, "dat_gen/sims")))
write.csv(phase_plane_var, file = paste0(filepath, "dat_gen/sims/traj_for_phase_plane.csv"))


#### Plotting trajectory                                                    ####

## A function to simplify the x axis:
evens <- function(x) subset(x, x %% 2 == 0)


plot4d <- ggplot(curr_run, aes(x = time, y = moi_type, color = type)) +
  geom_point(size = 0.7) + 
  geom_line(linewidth = 0.7) +
  scale_color_manual(labels = c("Res.", "Sus."),
                     values = color_in,
                     name = "Genotype") +
  # theme_light() + 
  xlab("Replications") + 
  ylab("Viral Density (Genomes/cell)") + 
  scale_y_log10(breaks=10^(evens(-6:10)),
                labels = sapply(evens(-6:10), function(i){parse(text = sprintf("10^%d",i))})) +
  scale_x_continuous(breaks = c(0:max(curr_run$time))) +
  guides(fill = guide_legend(override.aes = list(size = 0.3))) +
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
    legend.position.inside = c(0.865,0.14),
    legend.text = element_text(size = legend_text_size),
    legend.title = element_blank(),
    # legend.background = element_rect(fill = "transparent"),
    legend.key = element_rect(fill = "transparent"),
    legend.key.width = unit(0, 'cm'),
    legend.key.size = unit(0.0, "cm"),
    legend.background = element_rect(fill = "white", color = "black",linewidth = 0.2),  # Fills legend background & adds a border
    legend.box.background = element_rect(color = "black", linewidth = 0),  # Adds an outer border
    legend.margin = margin(0.6, 0.6, 0.6, 2),
  ) 
        
suppressWarnings(dir.create(paste0(filepath, "res/fig4/")))

ggsave(paste0(filepath, "res/fig4/D.png"), plot4d, h = height/2, w = width/3, units = "cm", bg = "transparent", dpi = 1200)
ggsave(paste0(filepath, "res/fig4/D.pdf"), plot4d, h = height/2, w = width/3, units = "cm", bg = "transparent", dpi = 1200)
ggsave(paste0(filepath, "res/fig4/D.svg"), plot4d, h = height/2, w = width/3, units = "cm", bg = "transparent", dpi = 1200, device = "svg")
saveRDS(plot4d, paste0(filepath, "res/fig4/D.rds"))

