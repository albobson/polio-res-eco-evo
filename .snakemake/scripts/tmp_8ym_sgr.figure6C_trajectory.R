
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
    output = list('runs/ddt_logistic_mu_2e-05_2025-03-24/res/fig6/C.rds'),
    params = list('2025-03-24', c('logistic'), c(2e-05), "date" = '2025-03-24', "fit_func" = c('logistic'), "mu" = c(2e-05)),
    wildcards = list(),
    threads = 1,
    log = list(),
    resources = list('mem_mb', 'disk_mb', 'tmpdir', 'log_loc', 'mfree', 'cluster_time', 'disk_free', 'runtime', 'n_cores', "mem_mb" = '3G', "disk_mb" = '4G', "tmpdir" = '/tmp/12723950.1.feder-short.q', "log_loc" = './log_scr/', "mfree" = '3G', "cluster_time" = '0:00:10:00', "disk_free" = '4G', "runtime" = 10, "n_cores" = '1'),
    config = list(),
    rule = 'fig6C',
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
############################## Figure 6 C ######################################

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

## I also want my axis text to be 12 pt font, so, setting that here
axis_text_size = 12

## I would like for the legend to be slightly smaller than the rest of the text.
## I'm setting that size here
legend_text_size = 6

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
pocap_run <- determ_polv(n = 6,
                         moi_mut_start = 100*0.0001, 
                         moi_wt_start = 100*(1-0.0001),
                        fit_func_in = fit_func$prob_surv,
                        v_prog = optim_params$optim_v_prog,
                        p2pfu = optim_params$optim_p2pfu
) %>%
  mutate(fit_type = "Pocapavir")

less_str_run <- determ_polv(n = 6,
                            moi_mut_start = 100*0.0001, 
                            moi_wt_start = 100*(1-0.0001),
                             fit_func_in = scaled_fit_func$prob_surv,
                             v_prog = optim_params$optim_v_prog,
                             p2pfu = optim_params$optim_p2pfu
) %>%
  mutate(fit_type = "100x Less Stringent")

tot_df <- rbind(pocap_run, less_str_run)

tot_df$fit_type <- factor(tot_df$fit_type, levels = c("Pocapavir", "100x Less Stringent"))

#### Plotting trajectory                                                    ####
plot6c <- ggplot(tot_df, aes(x = time, y = moi_type, color = type, linetype = fit_type)) +
  # geom_point(size = 3) + 
  geom_line(linewidth = 1) +
  scale_color_manual(labels = c("Res.", "Sus."),
                     values = c("#c94d4d", "#688cc8"),
                     name = "Genotype") +
  theme_light() + 
  # scale_linetype(guide = "none") +
  xlab("Passages") + ylab("MOI") + 
  scale_y_log10(breaks=10^(-6:10),
                labels = sapply(-6:10,function(i){parse(text = sprintf("10^%d",i))}),
                limits = c(10^-5, 10^2.4)) +
  scale_x_continuous(breaks = c(0:max(tot_df$time))) +
  guides(
    color = guide_legend(order = 1),
    linetype = guide_legend(order = 2)
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
    legend.position.inside = c(0.5,0.05),
    legend.text = element_text(size = legend_text_size),
    legend.title = element_blank(),
    legend.background = element_rect(fill = "transparent"),
    legend.key = element_rect(fill = "transparent"),
    # legend.key.width = unit(0, 'cm'),
    legend.key.size = unit(0.3, "cm"),
    legend.box = "horizontal"
  ) 


plot6c

suppressWarnings(dir.create(paste0(filepath, "res/fig6/")))

ggsave(paste0(filepath, "res/fig6/C.png"), plot6c, h = 22/3, w = 18/3, unit = "cm", bg = "transparent", dpi = 300)
ggsave(paste0(filepath, "res/fig6/C.pdf"), plot6c, h = 22/3, w = 18/3, unit = "cm", bg = "transparent", dpi = 300)
ggsave(paste0(filepath, "res/fig6/C.svg"), plot6c,h = 22/3, w = 18/3, unit = "cm", bg = "transparent", dpi = 300, device = "svg")
saveRDS(plot6c, paste0(filepath, "res/fig6/C.rds"))

