
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
    input = list('dat/collett_trial.csv', 'runs/ddt_logistic_mu_2e-05_2025-03-24/dat_gen/params/optim_params.csv', 'runs/ddt_logistic_mu_2e-05_2025-03-24/dat_gen/params/logistic_fitness_function.csv', 'runs/ddt_logistic_mu_2e-05_2025-03-24/dat_gen/sims/sim_clinical_trial.csv'),
    output = list('runs/ddt_logistic_mu_2e-05_2025-03-24/res/fig4/F.rds'),
    params = list('2025-03-24', c('logistic'), c(2e-05), "date" = '2025-03-24', "fit_func" = c('logistic'), "mu" = c(2e-05)),
    wildcards = list(),
    threads = 1,
    log = list(),
    resources = list('mem_mb', 'disk_mb', 'tmpdir', 'log_loc', 'mfree', 'cluster_time', 'disk_free', 'runtime', 'n_cores', "mem_mb" = '3G', "disk_mb" = '4G', "tmpdir" = '/tmp/12723959.1.feder-short.q', "log_loc" = './log_scr/', "mfree" = '3G', "cluster_time" = '0:00:10:00', "disk_free" = '4G', "runtime" = 10, "n_cores" = '1'),
    config = list(),
    rule = 'fig4F',
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
############################### Figure 5 A #####################################

## Reason:

## By incorporating stochastic immune clearance, we can run simulations of the
## pocapavir clinical trial. We observe a similar bifurcation in outcomes.

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

## Making the color map for plotting
color_in <- c("#c94d4d", "#688cc8")

## geom_text() is in terms of mm, where element_text() is in terms of pt. I want
## my geom_text and axis text to be the same. Here, I am defining a conversion
## which will keep the text sizes the same between text types.
geom_text_conv = 0.3528

## I also want my axis text to be 12 pt font, so, setting that here
axis_text_size = 12

## I would like for the legend to be slightly smaller than the rest of the text.
## I'm setting that size here
legend_text_size = 8


## Source functions to run simulations
source("scr/polv_DDT_functions.R")  ## Main functions 

## Read in Collett data
coldat <- read.csv("dat/collett_trial.csv")

## Read in clincial trial simulations
trial <- read.csv(paste0(filepath, "dat_gen/sims/sim_clinical_trial.csv"))


#### Cleaning Collett data                                                  ####
coldat_p <- coldat %>% group_by(treatment) %>% 
  mutate(res = as.character(resistant), 
         clearange = as.numeric(clearange),
         source = "Clinical Trial") %>% 
  select(clearange, res, treatment, res) %>% 
  na.omit() %>% 
  # filter(treatment == "Pocapavir (n = 93)") %>%
  mutate(time = clearange, 
         treatment = ifelse(treatment == "Pocapavir (n = 93)", "Pocapavir", "Placebo"),
         # treatment = paste0(treatment, " ", res) ## if you want treatment to be stratified by resistance class
  ) %>%
  select(time, res, treatment) %>%
  mutate(shape = "circle") ## Adding a shape value for distinguising our own data

## Finding the earliest placebo clearer
early_placebo <- coldat_p %>%
  filter(treatment == "Placebo") %>%
  summarize(mint = min(time))

early_placebo <- early_placebo$mint[1]

coldat_p$earlylate <- ifelse(coldat_p$time >= early_placebo, "Late", "Early")


#### Cleaning simulated trial                                               ####
simsclean <- trial %>% group_by(treatment, id) %>% 
  filter(type == "resistant") %>%
  mutate(pop_prop = as.numeric(pop_prop)) %>%
  na.omit() %>%
  mutate(res = sum(surv_pfu)/sum(as.numeric(tot_pfu)) > 0.5,
         perc_res = sum(surv_pfu)/sum(as.numeric(tot_pfu)),
         log_perc_res = log(sum(surv_pfu))/log(sum(as.numeric(tot_pfu))),
         avg_pop_prop = mean(as.numeric(pop_prop)),
         tot_vl = sum(as.numeric(tot_pfu)),
  ) %>%
  select(time, res, perc_res, tot_vl, treatment, id, time_to_pocap, avg_pop_prop,
         log_perc_res) %>% 
  filter(time == max(time)) %>%
  mutate(res = "Simulated") %>% ## Just coloring by if something was simulated or not
  # mutate(res = ifelse(res == TRUE, "Resistant", "Susceptible"),
  #        res_mean = ifelse(avg_pop_prop > 0.5, "Resistant", "Susceptible"),
  #        log_res = ifelse(log_perc_res > 0.5, "Resistant", "Susceptible")) %>%
  unique() %>%
  mutate(time = ceiling(time))

## Fix the fact that not every day was sampled
simsclean$time <- ifelse(simsclean$time %in%  16:17, 18,
                         ifelse(simsclean$time %in% 19:21, 22,
                                ifelse(simsclean$time %in% 23:28, 29,
                                       ifelse(simsclean$time %in% 30:42, 43,
                                              ifelse(simsclean$time > 43, 44,
                                                     simsclean$time)))))


## Adding if they are early or late clearers
simsclean$earlylate <- ifelse(simsclean$time >= early_placebo, "Late", "Early")


#### . . . Plotting                                                         ####
# library(RColorBrewer)
# early_late_pal <- data.frame(earlylate = c("Late", "Early"),
#                              color = c("#9467bd","#2ca02c"))

# Create a named vector of colors
# early_late_mapping <- setNames(early_late_pal$color, early_late_pal$earlylate)

## Binding
tot_df_w_placebo <- rbind(coldat_p, simsclean)

## Removing placebo 
tot_df_og <- tot_df_w_placebo %>%
  filter(treatment != "Placebo")

## Changing the names from pocapavir to observed and from simulated pocapavir to simulated
tot_df <- tot_df_og
tot_df$treatment <- ifelse(tot_df_og$treatment == "Pocapavir", "Observed", "Simulated")

## Creating a dynamic size for the dot plot
dotsize <- (max(tot_df$time)-min(tot_df$time))/max(tot_df$time)-0.1


## Plotting
tend_plot <- 
  ggplot(tot_df, aes(x = treatment, y = time,
  )) +
  annotate("rect", ymin = 15.5, ymax = 17.5, xmin = 0.5, xmax = 2.5, colour = "#f2f2f2", alpha=0.2) +
  annotate("rect", ymin = 18.5, ymax = 21.5, xmin = 0.5, xmax = 2.5, colour = "#f2f2f2", alpha=0.2) +
  annotate("rect", ymin = 22.5, ymax = 28.5, xmin = 0.5, xmax = 2.5, colour = "#f2f2f2", alpha=0.2) +
  annotate("rect", ymin = 29.5, ymax = 42.5, xmin = 0.5, xmax = 2.5, colour = "#f2f2f2", alpha=0.2) +
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
               # key_glyph = "point",
               # mapping = aes(fill = resistant)
  ) +
  # scale_fill_manual(values = early_late_mapping) +
  scale_y_continuous(limits = c(0, 43.5),
                     breaks = c(seq(from = 0, to = 15, by = 5), 18, 22, 29, 43),
                     expand = c(0,0)) +
  labs(x = "Treatment", 
       y = "Clearance Date (DPI)") +
  theme_classic() +
  geom_hline(yintercept = early_placebo - 0.5, linetype = 2, linewidth = 0.5) +
  theme(text = element_text(size= axis_text_size), 
        axis.text = element_text(size= axis_text_size),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "lightgrey"),
        panel.grid.minor.y = element_line(color = "lightgrey"),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        axis.title.y = element_text(),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.background = element_rect(fill='transparent'), 
        plot.background = element_rect(fill='transparent', color=NA),
        ## Legend stuff
        legend.position = "inside",
        legend.position.inside = c(0.18, 0.06),
        legend.text = element_text(size = legend_text_size),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent"),
        legend.key.width = unit(0, 'cm'),
        legend.key.size = unit(0.3, "cm")
  )  

tend_plot


suppressWarnings(dir.create(paste0(filepath, "res/fig4")))

## Save
ggsave(tend_plot, file = paste0(filepath, "res/fig4/F.png"), h = 22/3, w = 18/2, unit = "cm", bg = "transparent", dpi = 300)
ggsave(tend_plot, file = paste0(filepath, "res/fig4/F.pdf"), h = 22/3, w = 18/2, unit = "cm", bg = "transparent", dpi = 300)
ggsave(tend_plot, file = paste0(filepath, "res/fig4/F.svg"), h = 22/3, w = 18/2, unit = "cm", bg = "transparent", dpi = 300, device = 'svg')
saveRDS(tend_plot, paste0(filepath, "res/fig4/F.rds"))
