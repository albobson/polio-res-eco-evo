
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
    output = list('runs/ddt_logistic_mu_2e-05_2025-03-24/res/fig3/D.rds'),
    params = list('2025-03-24', c('logistic'), c(2e-05), "date" = '2025-03-24', "fit_func" = c('logistic'), "mu" = c(2e-05)),
    wildcards = list(),
    threads = 1,
    log = list(),
    resources = list('mem_mb', 'disk_mb', 'tmpdir', 'log_loc', 'mfree', 'cluster_time', 'disk_free', 'runtime', 'n_cores', "mem_mb" = '3G', "disk_mb" = '4G', "tmpdir" = '/tmp/12722233.1.feder-short.q', "log_loc" = './log_scr/', "mfree" = '3G', "cluster_time" = '0:00:10:00', "disk_free" = '4G', "runtime" = 10, "n_cores" = '7'),
    config = list(),
    rule = 'fig3D',
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
############################### Figure 3 D #####################################

## Reason:

## This shows the model's ability to capture the DDT phenominon observed in cell
## culture. It uses the fitness function that was generated in figure 3B and
## requires the respective parameters that were generated through the
## parameter_optimization.R script.

#### Set Up                                                                 ####

## libraries
require(dplyr)
require(tidyr)
require(ggplot2)
require(cowplot)

## Set up where to save
## Date of the run
init_date <- as.character(snakemake@params[["date"]])

## Fitness function
init_fit_func <- as.character(snakemake@params[["fit_func"]])

## Mutation rate
init_mu <- as.numeric(snakemake@params[["mu"]])

## Create the filepath where things will be saved
filepath <- paste0("runs/ddt_", init_fit_func, "_mu_", init_mu, "_", init_date, "/")

## geom_text() is in terms of mm, where element_text() is in terms of pt. I want
## my geom_text and axis text to be the same. Here, I am defining a conversion
## which will keep the text sizes the same between text types.
geom_text_conv = 0.3528

## I also want my axis text to be 12 pt font, so, setting that here
axis_text_size = 14

## I would like for the legend to be slightly smaller than the rest of the text.
## I'm setting that size here
legend_text_size = 10



## Read in fitness function and parameters
fit_func <- read.csv(paste0(filepath, "dat_gen/params/logistic_fitness_function.csv"))
optim_params <- read.csv(paste0(filepath, "dat_gen/params/optim_params.csv"))


## Source functions to run simulations
source("scr/polv_DDT_functions.R")  ## Main functions 
source("scr/model_optimization_funcs.R")  ## Optimization functions

## Read in Tanner cell culture data
tan_pure <- read.csv("dat/tanner_pure_culture.csv") ## Mono culture data
tan_mix <- read.csv("dat/tanner_mix_culture.csv")   ## Co-culture data

#### . . . . . . My Data                                                    ####
## Generating the data
cl <- makeCluster((detectCores() - 1))
registerDoParallel(cl)

my_data_pocap <- foreach(wt = c(0, 5, 15, 50, 100), 
                         .combine = 'rbind', 
                         .packages = c('dplyr')) %dopar% {
                           (determ_polv(n = 1,
                                        moi_mut_start = 10, 
                                        moi_wt_start = wt,
                                        v_prog = optim_params$optim_v_prog,
                                        p2pfu = optim_params$optim_p2pfu,
                                        max_vpc = optim_params$optim_max_vpc,
                                        fit_func_in = fit_func$prob_surv
                           ))
                         }

stopCluster(cl)
stopImplicitCluster()

## Cleaning the data
my_data_clean <- my_data_pocap %>%
  select(type, time, init_wt, init_mut, wt_rep_abil, surv_pfu, tot_pfu, pop_prop) %>%
  filter(time == 1, type == "resistant") %>%
  group_by(init_wt, wt_rep_abil) %>%
  summarize(sum_pfu = sum(surv_pfu)) %>%
  # filter(time == 1, type == "susceptible") %>%
  mutate(moi_wt = factor(init_wt)) %>%
  mutate(pocap = ifelse(wt_rep_abil == 1, "-", "+")) %>%
  mutate(pocap = factor(pocap, levels = c("+", "-")))

my_data_clean$sim <- "Simulated"

my_data_clean$sim <- factor(my_data_clean$sim, levels = c("Simulated", "Cell Culture"))

## Plotting my data
my_p <- ggplot(my_data_clean, aes(x = moi_wt, y = sum_pfu, fill = sim)) +
  geom_bar(stat="identity", position = "dodge", fill = "black") +
  theme_classic() + 
  labs(fill = "Type") +
  xlab("WT MOI") + ylab("Res. Viral Yield (PFU/mL)") + 
  scale_y_continuous(trans="log10",
                     breaks=10^(0:10),
                     # expand = c(0.01, 0),
                     labels = sapply(0:10,function(i){parse(text = sprintf("10^%d",i))})) +
  coord_cartesian(ylim=c(10^4,10^9)) +
  theme(text = element_text(size= axis_text_size), 
        axis.text = element_text(size= axis_text_size),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "lightgrey"),
        panel.grid.minor.y = element_line(color = "lightgrey"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(colour="black"),
        axis.title.y = element_text(),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.title = element_blank(), 
        legend.position = "none",
        panel.background = element_rect(fill='transparent'), 
        plot.background = element_rect(fill='transparent', color=NA))

#### Creating X axis labels                                                 ####
## Cleaning Data
tan <- tan_mix %>%
  mutate(moi_wt = factor(moi_wt)) %>%
  mutate(pocap = ifelse(pocap == "yes", "+", "-")) %>%
  mutate(pocap = factor(pocap, levels = c("+", "-"))) %>% ## To fix the ordering 
  filter(mutation == "VP3-A24V", strain == "Mahoney", pocap == "+")

## Adding the appropriate labels
longer_tan <- data.frame(lapply(tan, as.character)) ## Turning everything to character

## Renaming the columns according to how I want them plotted
longer_tan <- longer_tan %>%
  mutate("Susceptible MOI" = moi_wt,
         "Resistant MOI" = moi_mut,
         "Pocapavir" = pocap) %>%
  select(-moi_mut, -pfu) %>%
  mutate(cond = moi_wt) %>% ## Adding a condition column
  pivot_longer(-c(cond, pocap)) %>%## Pivoting longer by condition
  mutate(pocap = factor(pocap, levels = c("+", "-"))) 

## Making a dataframe for pocapavir
tan_pocap_df <- longer_tan %>%
  filter(name == "Pocapavir", cond == "0") %>%
  mutate(name = "Simulated") ## Changing this to observed for the figure

## Plotting just pocapavir
tan_pocap <- ggplot(tan_pocap_df, aes(cond, name, label = name)) +
  geom_text(size = axis_text_size * geom_text_conv) +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  facet_grid(~pocap) +
  theme(panel.grid = element_blank(), axis.text.x = element_blank(),
        text = element_text(size= axis_text_size),
        strip.text = element_blank(),
        axis.text.y=element_blank(),
        strip.background = element_blank())

tan_pocap

## Making a dataframe for the MOIs
tan_moi_df_res <- longer_tan %>%
  filter(name != "Pocapavir", name != "moi_wt", name == "Resistant MOI") %>%
  mutate(value = factor(value, levels = c("0", "5", "10", "15", "50", "100")),
         cond = factor(cond, levels = c("0", "5", "10", "15", "50", "100")))

tan_moi_df_sus <- longer_tan %>%
  filter(name != "Pocapavir", name != "moi_wt", name == "Susceptible MOI") %>%
  mutate(value = factor(value, levels = c("0", "5", "10", "15", "50", "100")),
         cond = factor(cond, levels = c("0", "5", "10", "15", "50", "100")))

## Plotting just pocapavir
tan_mois <- ggplot() +
  geom_text(tan_moi_df_res, 
            mapping = aes(cond, name, label = value, color = name), 
            size = axis_text_size * geom_text_conv) +
  geom_text(tan_moi_df_sus, 
            mapping = aes(cond, name, label = value, color = name), 
            size = axis_text_size * geom_text_conv) +
  labs(x = NULL, y = NULL) +
  scale_color_manual(values = c("#c94d4d", "#688cc8")) +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.text.y = element_text(colour = color_in),
        axis.text.x = element_blank(),
        text = element_text(size=axis_text_size),
        legend.position = "none",
        strip.text = element_blank(),
        strip.background = element_blank())

## The full plot together
## Note, I had to make some "null" plots in order to manually adjust spacing
full_my <- plot_grid(tan_pocap, NULL, my_p, NULL, tan_mois, ## Adding NULL to move plots closer together
                      rel_heights = c(0.3, -0.1, 1.2, -0.05, 0.3),
                      ncol = 1, axis = "l", align = "v")
full_my

full_my_list <- list(tan_pocap, my_p, tan_mois)

ggsave(paste0(filepath,"res/fig3/D.png"), full_my, h = 3.5, w = 3.5, units = "in", bg = "transparent", dpi = 300)
ggsave(paste0(filepath,"res/fig3/D.pdf"), full_my, h = 3.5, w = 3.5, units = "in", bg = "transparent", dpi = 300)
ggsave(paste0(filepath,"res/fig3/D.svg"), full_my, h = 3.5, w = 3.5, units = "in", bg = "transparent", device = 'svg')
saveRDS(full_my_list, paste0(filepath, "res/fig3/D.rds"))
