############################## Figure 6 D ######################################

## Reason:

## This will plot the simulated pocapavir next to the simulated less stringent
## drug

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

## Source functions to run simulations
source("scr/polv_DDT_functions.R")  ## Main functions

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


#### Cleaning Collett data                                                  ####
## Read in Collett data
coldat <- read.csv("dat/collett_trial.csv")

## Cleaning
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


#### Cleaning simulated data                                                ####
## Read in the clinical trial and the less stringent clinical trial
pocap_trial <- read.csv(paste0(filepath, "dat_gen/sims/sim_clinical_trial.csv"))
less_trial <- read.csv(paste0(filepath, "dat_gen/sims/less_str_clinical_trial.csv"))

## Clean
pocap_clean <- pocap_trial %>% group_by(treatment, id) %>% 
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

## Clean
less_clean <- less_trial %>% group_by(treatment, id) %>% 
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



#### . . . Plotting                                                         ####
# library(RColorBrewer)
# early_late_pal <- data.frame(earlylate = c("Late", "Early"),
#                              color = c("#9467bd","#2ca02c"))

# Create a named vector of colors
# early_late_mapping <- setNames(early_late_pal$color, early_late_pal$earlylate)

## Binding
tot_trial_df <- rbind(pocap_clean, less_clean)

## Changing their names to be shorter
tot_trial_df$treatment <- ifelse(tot_trial_df$treatment == "Simulated Pocapavir", "Pocapavir", "100x Less Stringent")

## Setting it as a factor to change the ordering in the plot
tot_trial_df$treatment <- factor(tot_trial_df$treatment, c("Pocapavir", "100x Less Stringent"))


## Creating a dynamic size for the dot plot
dotsize <- (max(tot_trial_df$time)-min(tot_trial_df$time))/max(tot_trial_df$time)-0.1


## Plotting
tend_plot <- 
  ggplot(tot_trial_df, aes(x = treatment, y = time, color = avg_pop_prop,
  )) +
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
               # mapping = aes(fill = avg_pop_prop)
  ) +
  scale_fill_gradient(low = color_in[2], high = color_in[1]) +
  scale_y_continuous(limits = c(0, 43.5),
                     breaks = c(seq(from = 0, to = 15, by = 5), 18, 22, 29, 43),
                     expand = c(0,0)) +
  labs(x = "Treatment", 
       y = "Clearance Date (DPI)") +
  theme_classic() +
  geom_hline(aes(yintercept = early_placebo - 0.5), linetype = 2, linewidth = 0.5) +
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
        # legend.position = "inside",
        # legend.position.inside = c(0.18, 0.06),
        # legend.text = element_text(size = legend_text_size),
        # legend.title = element_blank(),
        # legend.background = element_rect(fill = "transparent"),
        # legend.key = element_rect(fill = "transparent"),
        # legend.key.width = unit(0, 'cm'),
        # legend.key.size = unit(0.3, "cm")
  )  

tend_plot


suppressWarnings(dir.create(paste0(filepath, "res/fig6")))

## Save
ggsave(tend_plot, file = paste0(filepath, "res/fig6/D.png"), h = 22/3, w = 18/3, unit = "cm", bg = "transparent", dpi = 300)
ggsave(tend_plot, file = paste0(filepath, "res/fig6/D.pdf"), h = 22/3, w = 18/3, unit = "cm", bg = "transparent", dpi = 300)
ggsave(tend_plot, file = paste0(filepath, "res/fig6/D.svg"), h = 22/3, w = 18/3, unit = "cm", bg = "transparent", dpi = 300, device = 'svg')
saveRDS(tend_plot, paste0(filepath, "res/fig6/D.rds"))
