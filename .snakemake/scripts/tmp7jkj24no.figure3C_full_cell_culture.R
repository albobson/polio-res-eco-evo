
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
    input = list('dat/tanner_pure_culture.csv', 'dat/tanner_mix_culture.csv'),
    output = list('runs/ddt_logistic_mu_2e-05_2025-03-24/res/fig3/C.rds'),
    params = list('2025-03-24', c('logistic'), c(2e-05), "date" = '2025-03-24', "fit_func" = c('logistic'), "mu" = c(2e-05)),
    wildcards = list(),
    threads = 1,
    log = list(),
    resources = list('mem_mb', 'disk_mb', 'tmpdir', 'log_loc', 'mfree', 'cluster_time', 'disk_free', 'runtime', 'n_cores', "mem_mb" = '3G', "disk_mb" = '4G', "tmpdir" = '/tmp/12680028.1.feder-short.q', "log_loc" = './log_scr/', "mfree" = '3G', "cluster_time" = '0:00:10:00', "disk_free" = '4G', "runtime" = 10, "n_cores" = '1'),
    config = list(),
    rule = 'fig3C',
    bench_iteration = as.numeric(NA),
    scriptdir = '/net/feder/vol1/home/alexrob/projects/dominant_drug_targets/pipeline/scr',
    source = function(...){
        wd <- getwd()
        setwd(snakemake@scriptdir)
        source(...)
        setwd(wd)
    }
)


######## snakemake preamble end #########
############################### Figure 3 C #####################################

## Reason:

## This is plotting the full cell co-culture data for comparison with my model

#### Set Up                                                                 ####

## libraries
require(dplyr)
require(tidyr)
require(ggplot2)
require(cowplot)

## Read in Tanner cell culture data
tan_pure <- read.csv("dat/tanner_pure_culture.csv") ## Mono culture data
tan_mix <- read.csv("dat/tanner_mix_culture.csv")   ## Co-culture data

## Set up where to save
## Date of the run
init_date <- as.character(snakemake@params[["date"]])

## Fitness function
init_fit_func <- as.character(snakemake@params[["fit_func"]])

## Mutation rate
init_mu <- as.numeric(snakemake@params[["mu"]])

## Create the filepath where things will be saved
filepath <- paste0("runs/ddt_", init_fit_func, "_mu_", init_mu, "_", init_date, "/")

## Create this filepath
suppressWarnings(dir.create(paste0(filepath)))

#### Cleaning Data                                                          ####
tan <- tan_mix %>%
  mutate(moi_wt = factor(moi_wt)) %>%
  mutate(pocap = ifelse(pocap == "yes", "+", "-")) %>%
  mutate(pocap = factor(pocap, levels = c("+", "-"))) %>% ## To fix the ordering 
  filter(mutation == "VP3-A24V", strain == "Mahoney", pocap == "+")

## Plotting tanner data
tan_p <- ggplot(tan, aes(x = moi_wt, y = pfu, fill = "black")) +
  geom_bar(stat="identity", position = "dodge", fill = "black") +
  theme_classic() + 
  xlab("WT MOI") + ylab("Viral Yield (PFU/mL)") + 
  scale_y_continuous(trans="log10",
                     breaks=10^(0:10),
                     # expand = c(0.01, 0),
                     labels = sapply(0:10,function(i){parse(text = sprintf("10^%d",i))})) +
  coord_cartesian(ylim=c(10^4,10^9)) +
  theme(text = element_text(size=20), 
        axis.text = element_text(size=18),
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
        panel.background = element_blank()) 

## Adding the appropreate labels
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
  filter(name == "Pocapavir", cond == "0")

## Plotting just pocapavir
tan_pocap <- ggplot(tan_pocap_df, aes(cond, name, label = value)) +
  geom_text(size = 10) +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  facet_grid(~pocap) +
  theme(panel.grid = element_blank(), axis.text.x = element_blank(),
        text = element_text(size=22),
        strip.text = element_blank(),
        strip.background = element_blank())

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
            size = 9) +
  geom_text(tan_moi_df_sus, 
            mapping = aes(cond, name, label = value, color = name), 
            size = 9) +
  labs(x = NULL, y = NULL) +
  scale_color_manual(values = c("#c94d4d", "#688cc8")) +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.text.x = element_blank(),
        text = element_text(size=20),
        legend.position = "none") +
  facet_grid(~pocap)+
  theme(panel.grid = element_blank(), axis.text.x = element_blank(),
        text = element_text(size=22),
        strip.text = element_blank(),
        strip.background = element_blank())

## The full plot together
## Note, I had to make some "null" plots in order to manually adjust spacing
full_tan <- plot_grid(tan_p, NULL, tan_mois,
                      rel_heights = c(1.2, -0.05, 0.3),
                      ncol = 1, axis = "l", align = "v")

ggsave(paste0(filepath, "res/fig3/C.png"), full_tan, h = 8, w = 8, bg = "transparent")
ggsave(paste0(filepath, "res/fig3/C.jpg"), full_tan, h = 8, w = 8, bg = "transparent")
ggsave(paste0(filepath, "res/fig3/C.svg"), full_tan, h = 8, w = 8, bg = "transparent", device = 'svg')
saveRDS(full_tan, paste0(filepath, "res/fig3/C.rds"))