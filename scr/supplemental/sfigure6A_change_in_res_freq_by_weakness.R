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

## Cores
init_cores <- as.numeric(snakemake@resources[["n_cores"]])

## This figure's size
width <-  as.numeric(snakemake@params[["sfig6_dim"]])[1]
height <- as.numeric(snakemake@params[["sfig6_dim"]])[2]

## Find how many variations we want to run 
n_trials <- as.numeric(snakemake@params[["n_trials"]])

print(paste0("n_trials: ", n_trials))

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
## Generate a function that is less stringent than pocapavir
scaled_vector <- function(og_probs, scaler) {
  ## New vector for storing
  new_vec <- NULL
  
  new_vec <- ((1- scaler * og_probs[1])/(1-og_probs[1])) * (og_probs - og_probs[1]) + scaler * og_probs[1]
  
  return(new_vec)
}

# ## Define the range for which we will scale the function
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

############################### Run sims #######################################
range <- 0.0001
mois <- 10^seq(from = -2, to = 2.5, length.out = 50)

cl <- makeCluster(init_cores-1)
registerDoParallel(cl = cl)

low_freq_moi <- foreach(w = scale_range, 
                        .combine = 'rbind', 
                        .packages = c('dplyr')) %:%
  foreach(m = mois, 
          .combine = 'rbind', 
          .packages = c('dplyr'))  %dopar% {
            ## Filter to what we want
            fit_func <- scale_store[which(scale_store$fit_type == w),]
            ## Run suims
            determ_polv(n = 1,
                        moi_mut_start = m*range, 
                        moi_wt_start = m*(1-range),
                        fit_func_in = fit_func$prob_surv,
                        v_prog = optim_params$optim_v_prog,
                        p2pfu = optim_params$optim_p2pfu
            ) %>%
              mutate(weakness = w)
          }

stopImplicitCluster()
stopCluster(cl)

## Clean simulation
low_freq_clean <- low_freq_moi %>%
  select(type, time, moi_res, moi_wt, init_wt, init_mut, pop_prop, weakness) %>%
  filter(type == "resistant") %>%
  group_by(init_mut, init_wt, weakness) %>%
  arrange(time) %>%
  mutate(pop_prop = as.numeric(pop_prop),
         init_mut = as.numeric(init_mut),
         init_wt = as.numeric(init_wt),
         moi_res = as.numeric(moi_res),
         moi_wt = as.numeric(moi_wt)) %>%
  mutate(delta_p = (pop_prop - lag(pop_prop, default = first(pop_prop))),
         delta_moi_tot = ((moi_res+moi_wt) - lag((moi_res+moi_wt), default = first((moi_res+moi_wt)))),
         init_moi_tot = (init_wt+init_mut),
         init_freq_mut = init_mut/(init_wt+init_mut),
         dist_col = abs(delta_moi_tot) + abs(delta_p),
         moi_tot = (moi_wt + moi_res),
         last_p = lag(pop_prop, default = first(pop_prop))) %>%
  mutate(x_end = (init_moi_tot+cos(atan(delta_p/delta_moi_tot))),
         y_end = (init_freq_mut+sin(atan(delta_p/delta_moi_tot)))) %>%
  filter(time == 1) %>%
  ungroup()



################################ Plot ##########################################
## Define color range
num_lines <- length(scale_range)

# gray_range <- gray(seq(0, 0.9, length.out = num_lines))
gray_range <- c("#000000", "#696969", "#A9A9A9", "#D3D3D3")

splot <- ggplot(low_freq_clean, aes(x = init_moi_tot, y = delta_p, color = factor(weakness))) +
  geom_line(linewidth = 0.7) +
  scale_x_continuous(trans="log10", breaks=10^(-2:3),
                     labels = sapply(c(-2:3),function(i){parse(text = sprintf("10^%d",i))}),
                     limits = c(10^(-2.1), 10^(2.6))) +
  scale_y_continuous(trans="log10", 
                     breaks=10^(-4:0),
                     labels = sapply(c(-4:0),function(i){parse(text = sprintf("10^%d",i))}),
                     limits = c(10^(-4), 10^(0))
  ) +
  scale_color_manual(name = expression("\u00d7-Weaker"),
                     values = gray_range) +
  # ylim(0, round(max(low_freq_clean$delta_p), 1)) +
  # xlim(0, 102) +
  # theme_minimal() +
  ## This adds the points to the plot
  # geom_point(data = df_points_clean, 
  #            mapping = aes(x = init_moi_tot, 
  #                          y = delta_p, 
  #                          color = factor(round(init_moi_tot, 2))),
  #            size = 5) +
  # scale_color_manual(values = c("#64A61C", "#FFA61C"),
  #                    name = "Total MOI") +
  xlab("Total MOI") + 
  ylab(expression(Delta * " " * italic(f)[Res])) +# paste0("\U0394 Freq. Resistance")) +
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
    # legend.position.inside = c(0.35, 0.94), ## Upper left
    legend.position.inside = c(0.828, 0.83), ## bottom right
    legend.text = element_text(size = legend_text_size),
    # legend.title = element_blank(),
    # legend.background = element_rect(fill = "transparent"),
    legend.key = element_rect(fill = "transparent"),
    legend.key.width = unit(0, 'cm'),
    legend.key.height = unit(0, "cm"),
    legend.key.size = unit(0.0, "cm"),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.2),  # Fills legend background & adds a border
    legend.box.background = element_rect(color = "black", linewidth = 0),  # Adds an outer border
    legend.margin = margin(1, 1, 1, 2),
  ) 

suppressWarnings(dir.create(paste0(filepath, "res/")))
suppressWarnings(dir.create(paste0(filepath, "res/sup/")))
suppressWarnings(dir.create(paste0(filepath, "res/sup/sfig6")))

ggsave(paste0(filepath, "res/sup/sfig6/A.png"), splot, h = height, w = width/3, units = "cm", bg = "transparent", dpi = 1200)
ggsave(paste0(filepath, "res/sup/sfig6/A.pdf"), splot, h = height, w = width/3, units = "cm", bg = "transparent", dpi = 1200)
ggsave(paste0(filepath, "res/sup/sfig6/A.svg"), splot, h = height, w = width/3, units = "cm", bg = "transparent", dpi = 1200, device = 'svg')
saveRDS(splot, file = paste0(filepath, "res/sup/sfig6/A.rds"))

