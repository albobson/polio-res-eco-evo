############################### Figure 5 B #####################################

## Reason:

## Less stringent single step resistance change

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

#### Set up less stringent drugs                                            ####
scaled_vector <- function(og_probs, scaler) {
  ## New vector for storing
  new_vec <- NULL
  
  new_vec <- ((1- scaler * og_probs[1])/(1-og_probs[1])) * (og_probs - og_probs[1]) + scaler * og_probs[1]
  
  return(new_vec)
}

scaled_fit_func <- fit_func

## Define the scaler
scale <- 100

scaled_fit_func$prob_surv <- scaled_vector(og_probs = fit_func$prob_surv, scaler = scale)

scaled_fit_func$fit_type <- "Less Stringent Drug"

#### Run simulations                                                        ####
range <- 0.0001
mois <- 10^seq(from = -2, to = 2.3, length.out = 50)

cl <- makeCluster(detectCores()-2)
registerDoParallel(cl)

low_freq_moi <- foreach(r = range, 
                        .combine = 'rbind', 
                        .packages = c('dplyr')) %:%
  foreach(m = mois, 
          .combine = 'rbind', 
          .packages = c('dplyr'))  %dopar% {
            (determ_polv(n = 1,
                         moi_mut_start = m*r, 
                         moi_wt_start = m*(1-r),
                         fit_func_in = fit_func$prob_surv,
                         v_prog = optim_params$optim_v_prog,
                         p2pfu = optim_params$optim_p2pfu
            ))
          }

less_str <- foreach(r = range, 
                    .combine = 'rbind', 
                    .packages = c('dplyr')) %:%
  foreach(m = mois, 
          .combine = 'rbind', 
          .packages = c('dplyr'))  %dopar% {
            (determ_polv(n = 1,
                         moi_mut_start = m*r, 
                         moi_wt_start = m*(1-r),
                         fit_func_in = scaled_fit_func$prob_surv,
                         v_prog = optim_params$optim_v_prog,
                         p2pfu = optim_params$optim_p2pfu
            ))
          }

stopCluster(cl)
stopImplicitCluster()


## Clean simulations
low_freq_clean <- low_freq_moi %>%
  select(type, time, moi_res, moi_wt, init_wt, init_mut, pop_prop) %>%
  filter(type == "resistant") %>%
  group_by(init_mut, init_wt) %>%
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
  ungroup() %>%
  mutate(fit_type = "Pocapavir")

less_str_clean <- less_str %>%
  select(type, time, moi_res, moi_wt, init_wt, init_mut, pop_prop) %>%
  filter(type == "resistant") %>%
  group_by(init_mut, init_wt) %>%
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
  ungroup() %>%
  mutate(fit_type = "100x Weaker")

tot_df <- rbind(less_str_clean, low_freq_clean)


plot6b <- ggplot() +
  geom_line(data = low_freq_clean, mapping = aes(x = init_moi_tot, y = delta_p, color = fit_type),
            linewidth = 0.7, 
            # color = "black",
            linetype = "solid") +
  geom_line(data = less_str_clean, mapping = aes(x = init_moi_tot, y = delta_p, color = fit_type),
            linewidth = 0.7, 
            # color = "darkgrey",
            linetype = "dashed") +
  scale_y_continuous(trans="log10", 
                     breaks=10^(-4:0),
                     labels = sapply(c(-4:0),function(i){parse(text = sprintf("10^%d",i))}),
                     limits = c(10^(-4), 10^(0))
  ) +
  scale_x_continuous(trans="log10", breaks=10^(-2:3),
                     labels = sapply(c(-2:3),function(i){parse(text = sprintf("10^%d",i))}),
                     limits = c(10^(-2.1), 10^(2.6))) +
  xlab("Total MOI") + 
  ylab(expression(Delta * " f(Res.)")) +# paste0("\U0394 Freq. Resistance")) +
  # ylim(0, round(max(low_freq_clean$delta_p), 1)) +
  scale_color_manual(labels = c(paste0(scale, "x Weaker"), "Pocapavir"),
                     values = c("darkgrey", "black")) +
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
  scale_linetype(guide = "none") +
  xlab("Total MOI") + 
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
    legend.position.inside = c(0.7, 0.92), ## Upper left
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
suppressWarnings(dir.create(paste0(filepath, "res/fig5/")))

ggsave(paste0(filepath, "res/fig5/B.png"), plot6b, h = height/2, w = width/3, units = "cm", bg = "transparent", dpi = 1200)
ggsave(paste0(filepath, "res/fig5/B.pdf"), plot6b, h = height/2, w = width/3, units = "cm", bg = "transparent", dpi = 1200)
ggsave(paste0(filepath, "res/fig5/B.svg"), plot6b, h = height/2, w = width/3, units = "cm", bg = "transparent", dpi = 1200, device = 'svg')
saveRDS(plot6b, file = paste0(filepath, "res/fig5/B.rds"))

