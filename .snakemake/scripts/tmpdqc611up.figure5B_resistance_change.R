
######## snakemake preamble start (automatically inserted, do not edit) ########
library(methods)
Snakemake <- setClass(
    "Snakemake",
    slots = c(
        input = "list",
        output = "list",
        params = "list",
        wildcards = "list",
        threads = "numeric",
        log = "list",
        resources = "list",
        config = "list",
        rule = "character",
        bench_iteration = "numeric",
        scriptdir = "character",
        source = "function"
    )
)
snakemake <- Snakemake(
    input = list('runs/ddt_logistic_mu_2e-05_2025-03-24/dat_gen/params/optim_params.csv', 'runs/ddt_logistic_mu_2e-05_2025-03-24/dat_gen/params/logistic_fitness_function.csv'),
    output = list('runs/ddt_logistic_mu_2e-05_2025-03-24/res/fig5/B.rds'),
    params = list('2025-03-24', c('logistic'), c(2e-05), "date" = '2025-03-24', "fit_func" = c('logistic'), "mu" = c(2e-05)),
    wildcards = list(),
    threads = 1,
    log = list(),
    resources = list('mem_mb', 'disk_mb', 'tmpdir', 'log_loc', 'mfree', 'cluster_time', 'disk_free', 'runtime', 'n_cores', "mem_mb" = '3G', "disk_mb" = '4G', "tmpdir" = '/tmp/12722635.1.feder-short.q', "log_loc" = './log_scr/', "mfree" = '3G', "cluster_time" = '0:00:10:00', "disk_free" = '4G', "runtime" = 10, "n_cores" = '10'),
    config = list(),
    rule = 'fig5B',
    bench_iteration = as.numeric(NA),
    scriptdir = '/net/feder/vol1/home/alexrob/projects/dominant_drug_targets/pipeline/scr/panels',
    source = function(...){
        wd <- getwd()
        setwd(snakemake@scriptdir)
        source(...)
        setwd(wd)
    }
)


######## snakemake preamble end #########
############################### Figure 5 B #####################################

## Reason:

## More dominant drug single step resistance change

#### Set Up                                                                 ####

## Set up where to save
## Date of the run
init_date <- as.character(snakemake@params[["date"]])

## Fitness function
init_fit_func <- as.character(snakemake@params[["fit_func"]])

## Mutation rate
init_mu <- as.numeric(snakemake@params[["mu"]])

## Create the filepath where things will be saved
filepath <- paste0("runs/ddt_", init_fit_func, "_mu_", init_mu, "_", init_date, "/")

## Read in fitness function and parameters
fit_func <- read.csv(paste0(filepath, "dat_gen/params/logistic_fitness_function.csv"))
fit_func$fit_type <- "Pocapavir"

optim_params <- read.csv(paste0(filepath, "dat_gen/params/optim_params.csv"))

## Source functions to run simulations
source("scr/polv_DDT_functions.R")  ## Main functions 


## Making the color map for plotting
color_in <- c("#c94d4d", "#688cc8")

## Create colors for the funcs
custom_colors <- c("Pocapavir" = "black", # black
                   "Full Dominant" = "darkgrey") # Saddle Brown

## geom_text() is in terms of mm, where element_text() is in terms of pt. I want
## my geom_text and axis text to be the same. Here, I am defining a conversion
## which will keep the text sizes the same between text types.
geom_text_conv = 0.3528

## I also want my axis text to be 12 pt font, so, setting that here
axis_text_size = 12

## I would like for the legend to be slightly smaller than the rest of the text.
## I'm setting that size here
legend_text_size = 8

print("initialized")

#### Set up less stringent drugs                                            ####
high_dom <- data.frame(subunits = 0:60,
                       prob_surv = fit_func_logistic(k = 100, 
                                                     min_val = min(fit_func$prob_surv),
                                                     subunits = 60, mid = 59.5),
                       fit_type = "Full Dominant"
) 

print("high_dom created")

library(doParallel)

#### Run simulations                                                        ####
range <- 0.0001
mois <- 10^seq(from = -2, to = 2, length.out = 50)

cl <- makeCluster(detectCores()-2)
registerDoParallel(cl)

print("cluster created")

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

print("low_freq_moi created")

dom_freq_moi <- foreach(r = range, 
                    .combine = 'rbind', 
                    .packages = c('dplyr')) %:%
  foreach(m = mois, 
          .combine = 'rbind', 
          .packages = c('dplyr'))  %dopar% {
            (determ_polv(n = 1,
                         moi_mut_start = m*r, 
                         moi_wt_start = m*(1-r),
                         fit_func_in = high_dom$prob_surv,
                         v_prog = optim_params$optim_v_prog,
                         p2pfu = optim_params$optim_p2pfu
            ))
          }


print("dom_freq_moi created")

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

dom_freq_clean <- dom_freq_moi %>%
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
  mutate(fit_type = "100x Less Stringent")

print("dfs cleaned")

tot_df <- rbind(dom_freq_clean, low_freq_clean)


plot5b <- ggplot() +
  geom_line(data = low_freq_clean, mapping = aes(x = init_moi_tot, y = delta_p, color = fit_type),
            linewidth = 1, 
            # color = "black",
            linetype = "solid") +
  geom_line(data = dom_freq_clean, mapping = aes(x = init_moi_tot, y = delta_p, color = fit_type),
            linewidth = 1, 
            # color = "darkgrey",
            linetype = "longdash") +
  # scale_x_continuous(trans="log10", breaks=10^(-4:2),
  #                    labels = sapply(-4:2,function(i){parse(text = sprintf("10^%d",i))}),
  #                    limits = c(10^(-4), 10^(2))) +
  ylim(0, round(max(low_freq_clean$delta_p), 1)) +
  scale_color_manual(labels = c("Fully Dominant", "Pocapavir"),
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
  xlab("Total MOI") + ylab(expression(Delta * " Freq. Resistance")) +# paste0("\U0394 Freq. Resistance")) +
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
    legend.position.inside = c(.8, 0.95), ## bottom right (order: H, V)
    legend.text = element_text(size = legend_text_size),
    legend.title = element_blank(),
    legend.background = element_rect(fill = "transparent"),
    legend.key = element_rect(fill = "transparent"),
    legend.key.width = unit(0, 'cm'),
    legend.key.size = unit(0.3, "cm")
  ) 

print("plotted")


suppressWarnings(dir.create(paste0(filepath, "res/")))
suppressWarnings(dir.create(paste0(filepath, "res/fig5/")))


ggsave(paste0(filepath, "res/fig5/B.png"), plot5b, h = 22/3, w = 18/2, unit = "cm", bg = "transparent", dpi = 300)
ggsave(paste0(filepath, "res/fig5/B.pdf"), plot5b, h = 22/3, w = 18/2, unit = "cm", bg = "transparent", dpi = 300)
ggsave(paste0(filepath, "res/fig5/B.svg"), plot5b, h = 22/3, w = 18/2, unit = "cm", bg = "transparent", device = 'svg', dpi = 300)
saveRDS(plot5b, file = paste0(filepath, "res/fig5/B.rds"))

print("saved")