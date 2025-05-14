
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
    output = list('runs/ddt_logistic_mu_2e-05_2025-03-24/res/fig2/A.rds'),
    params = list('2025-03-24', c('logistic'), c(2e-05), "date" = '2025-03-24', "fit_func" = c('logistic'), "mu" = c(2e-05)),
    wildcards = list(),
    threads = 1,
    log = list(),
    resources = list('mem_mb', 'disk_mb', 'tmpdir', 'log_loc', 'mfree', 'cluster_time', 'disk_free', 'runtime', 'n_cores', "mem_mb" = '3G', "disk_mb" = '4G', "tmpdir" = '/tmp/12666807.1.feder-short.q', "log_loc" = './log_scr/', "mfree" = '3G', "cluster_time" = '0:00:10:00', "disk_free" = '4G', "runtime" = 10, "n_cores" = '1'),
    config = list(),
    rule = 'fig2A',
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
################################## Figure 2 A ##################################

## Reason:
## This is data from Tanner et al (2014). It is intended to show the susceptible
## interference phenominon exhibited by pocapavir

#### Set Up                                                                 ####

## Libraries 
library(dplyr)
library(ggplot2)
library(cowplot)
library(tidyr)

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

## Read in Tanner cell culture data
tan_pure <- read.csv("dat/tanner_pure_culture.csv") ## Mono culture data
tan_mix <- read.csv("dat/tanner_mix_culture.csv")   ## Co-culture data

## Making the color map for plotting
color_in <- c("#c94d4d", "#688cc8")

#### Cleaning data                                                          ####
## Separating Tanner data into pocap/no pocap
tan_trunc <- tan_mix %>%
  filter(mutation == "VP3-A24V", 
         strain == "Mahoney") %>%
  mutate(moi_wt = factor(moi_wt)) %>%
  filter(moi_wt == 0 | moi_wt == 100) %>%
  mutate(pocap = ifelse(pocap == "yes", "+", "-")) %>%
  mutate(pocap = factor(pocap, levels = c("-", "+"))) ## To fix the ordering

#### . . . . . . Plotting                                                   ####
## Plotting tanner data
tan_trunc_p <- ggplot(tan_trunc, aes(x = moi_wt, y = pfu, fill = "black")) +
  geom_bar(stat="identity", position = "dodge", fill = "black") +
  theme_classic() + 
  xlab("WT MOI") + ylab("Viral Yield (PFU/mL)") + 
  scale_y_continuous(trans="log10",
                     breaks=10^(0:10),
                     # expand = c(0.01, 0),
                     labels = sapply(0:10,function(i){parse(text = sprintf("10^%d",i))})) +
  coord_cartesian(ylim=c(10^4,10^9)) +
  theme(text = element_text(size= 22), 
        axis.text = element_text(size= 22),
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
        plot.background = element_rect(fill='transparent', color=NA)) +
  facet_grid(~pocap)

## Adding labels
longer_tan <- data.frame(lapply(tan_trunc, as.character)) ## Turning everything to character

## Renaming the columns according to how I want them ploted
longer_tan <- longer_tan %>%
  mutate("Susceptible MOI" = moi_wt,
         "Resistant MOI" = moi_mut,
         "Pocapavir" = pocap) %>%
  select(-moi_mut, -pfu) %>%
  mutate(cond = moi_wt) %>% ## Adding a condition column
  pivot_longer(-c(cond, pocap)) %>%
  mutate(pocap = factor(pocap, levels = c("-", "+"))) ## Pivoting longer by condition

## Making a dataframe for pocapavir
tan_pocap_df <- longer_tan %>%
  filter(name == "Pocapavir", cond == "0")

## Plotting pocapavir title
tan_pocap <- ggplot(tan_pocap_df, aes(cond, name, label = value)) +
  geom_text(size = 16) +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  facet_grid(~pocap) +
  theme(panel.grid = element_blank(), axis.text.x = element_blank(),
        text = element_text(size= 26),
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

## Plotting MOI titles
tan_mois <- ggplot() +
  geom_text(tan_moi_df_res, 
            mapping = aes(cond, name, label = value, color = name), 
            size = 10) +
  geom_text(tan_moi_df_sus, 
            mapping = aes(cond, name, label = value, color = name), 
            size = 10) +
  labs(x = NULL, y = NULL) +
  scale_color_manual(values = color_in) +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.text.y = element_text(colour = color_in),
        axis.text.x = element_blank(),
        text = element_text(size=26),
        legend.position = "none",
        strip.text = element_blank(),
        strip.background = element_blank()) +
  facet_grid(~pocap)

## The full plot together
## Note, I had to make some "null" plots in order to manually adjust spacing
full_tan_trunc <- plot_grid(tan_pocap, NULL, tan_trunc_p, NULL, tan_mois, ## Adding NULL to move plots closer together
                            rel_heights = c(0.3, -0.1, 1.2, -0.05, 0.3),
                            ncol = 1, axis = "l", align = "v")

full_tan_list <- list(tan_pocap, tan_trunc_p, tan_mois)


suppressWarnings(dir.create(paste0(filepath, "res/")))

ggsave(paste0(filepath,"res/fig2/A.png"), full_tan_trunc, h = 8, w = 8, bg = "transparent")
ggsave(paste0(filepath,"res/fig2/A.jpg"), full_tan_trunc, h = 8, w = 8, bg = "transparent")
ggsave(paste0(filepath,"res/fig2/A.svg"), full_tan_trunc, h = 8, w = 8, bg = "transparent", device = 'svg')
saveRDS(full_tan_list, paste0(filepath, "res/fig2/A.rds"))
