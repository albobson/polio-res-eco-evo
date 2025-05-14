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

## Color scheme for resistant and susceptible
color_in <- unlist(snakemake@params[["rs_colors"]])

## Read in axis text size
axis_text_size <- as.numeric(snakemake@params[["axis_text_size"]])

## Legend text size
legend_text_size <- as.numeric(snakemake@params[["lege_text_size"]])

## This figure's size
width <-  unlist(as.numeric(snakemake@params[["fig2_dim"]]))[1]
height <- unlist(as.numeric(snakemake@params[["fig2_dim"]]))[2]

## Create the filepath where things will be saved
filepath <- paste0("runs/ddt_", init_fit_func, "_mu_", init_mu, "_", init_date, "/")

## Create this filepath
suppressWarnings(dir.create(paste0(filepath)))

## Read in Tanner cell culture data
tan_pure <- read.csv("dat/tanner_pure_culture.csv") ## Mono culture data
tan_mix <- read.csv("dat/tanner_mix_culture.csv")   ## Co-culture data

## geom_text() is in terms of mm, where element_text() is in terms of pt. I want
## my geom_text and axis text to be the same. Here, I am defining a conversion
## which will keep the text sizes the same between text types.
geom_text_conv = 0.3528



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
  xlab("WT MOI") + ylab("Res. Viral Yield (PFU/mL)") + 
  scale_y_continuous(trans="log10",
                     breaks=10^(0:10),
                     # expand = c(0.01, 0),
                     labels = sapply(0:10,function(i){parse(text = sprintf("10^%d",i))})) +
  coord_cartesian(ylim=c(10^6,10^9)) +
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
        plot.background = element_rect(fill='transparent', color=NA)) +
  facet_grid(~pocap)

## Adding labels
longer_tan <- data.frame(lapply(tan_trunc, as.character)) ## Turning everything to character

## Renaming the columns according to how I want them ploted
longer_tan <- longer_tan %>%
  mutate("Sus. MOI" = moi_wt,
         "Res. MOI" = moi_mut,
         "Pocapavir" = pocap) %>%
  select(-moi_mut, -pfu) %>%
  mutate(cond = moi_wt) %>% ## Adding a condition column
  pivot_longer(-c(cond, pocap)) %>%
  mutate(pocap = factor(pocap, levels = c("-", "+"))) ## Pivoting longer by condition

## Making a dataframe for pocapavir
tan_pocap_df <- longer_tan %>%
  filter(name == "Pocapavir", cond == "0") %>%
  mutate(pocap_value = paste0(value, " ", name))

## Plotting pocapavir title
tan_pocap <- ggplot(tan_pocap_df, aes(cond, name, label = pocap_value)) +
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
  filter(name != "Pocapavir", name != "moi_wt", name == "Res. MOI") %>%
  mutate(value = factor(value, levels = c("0", "5", "10", "15", "50", "100")),
         cond = factor(cond, levels = c("0", "5", "10", "15", "50", "100")))

tan_moi_df_sus <- longer_tan %>%
  filter(name != "Pocapavir", name != "moi_wt", name == "Sus. MOI") %>%
  mutate(value = factor(value, levels = c("0", "5", "10", "15", "50", "100")),
         cond = factor(cond, levels = c("0", "5", "10", "15", "50", "100")))

## Plotting MOI titles
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
        strip.background = element_blank()) +
  facet_grid(~pocap)

## The full plot together
## Note, I had to make some "null" plots in order to manually adjust spacing
full_tan_trunc <- plot_grid(tan_pocap, NULL, tan_trunc_p, NULL, tan_mois, ## Adding NULL to move plots closer together
                            rel_heights = c(0.3, -0.15, 1.2, -0.1, 0.3),
                            ncol = 1, axis = "l", align = "v")

full_tan_list <- list(tan_pocap, tan_trunc_p, tan_mois)


suppressWarnings(dir.create(paste0(filepath, "res/")))

ggsave(paste0(filepath,"res/fig2/A.png"), full_tan_trunc, h = height, w = width/2, units = "cm", bg = "transparent", dpi = 1200)
ggsave(paste0(filepath,"res/fig2/A.pdf"), full_tan_trunc, h = height, w = width/2, units = "cm", bg = "transparent", dpi = 1200)
ggsave(paste0(filepath,"res/fig2/A.svg"), full_tan_trunc, h = height, w = width/2, units = "cm", bg = "transparent", device = 'svg')
saveRDS(full_tan_list, paste0(filepath, "res/fig2/A.rds"))


# ggsave(paste0(filepath,"res/fig2/A.png"), full_tan_trunc, h = 5, w = 4, units = "cm", bg = "transparent", dpi = 300)
# ggsave(paste0(filepath,"res/fig2/A.pdf"), full_tan_trunc, h = 5, w = 4, units = "cm", bg = "transparent", dpi = 300)
# ggsave(paste0(filepath,"res/fig2/A.svg"), full_tan_trunc, h = 5, w = 4, units = "cm", bg = "transparent", device = 'svg')