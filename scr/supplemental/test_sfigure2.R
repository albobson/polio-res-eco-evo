############################### SFigure 2 ######################################

## Reason:

## This will show the distribution of progeny phenotypes as a function of
## resistance frequency and MOI

#### Set Up                                                                 ####

library(ggridges)

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
width <-  as.numeric(snakemake@params[["sfig2_dim"]])[1]
height <- as.numeric(snakemake@params[["sfig2_dim"]])[2]

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

print("stuff read in")
#### Run the simulations                                                    ####
range <- 10^seq(from = -4, to = -1, by = 1)
range <- c(range, 0.5)
mois <- c(1, 10, 100)

cl <- makeCluster(detectCores() - 2)
registerDoParallel(cl)
print("Cluster created")

df <- foreach(r = range, 
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
                         p2pfu = optim_params$optim_p2pfu,
                         report_subunit_dist = TRUE) %>%
               mutate(tot_init_moi = m))
          }


print("Sims run")
stopImplicitCluster()
stopCluster(cl)
print("Cluster stopped")

## Clean the data
df_clean <- df %>%
  group_by(tot_init_moi, type) %>%
  filter(time == 1, type != "total") %>%
  mutate(prop_res = init_moi_res/(tot_init_moi)) %>%
  pivot_longer(cols = c(dens, dens_b4_pocap), names_to = "dens_at", values_to = "dens_val")

df_clean$dens_at <- factor(df_clean$dens_at,
                           levels = c("dens_b4_pocap", "dens"),
                           labels = c("Before pocapavir", "After pocapavir"))

df_clean <- df_clean %>%
  filter(dens_at == "Before pocapavir")
print("Data cleaned")


df_clean$type <- factor(df_clean$type, levels = c("susceptible", "resistant"))

df_clean$prop_res <- as.character(df_clean$prop_res)

## Plot the data
sub_dens_b4_plot <- ggplot(df_clean, aes(x = subs, y = dens_val, fill = type)) +
  geom_bar(stat = "identity", position = "dodge") +
  # geom_density_ridges(stat = "binline", bins = 61, scale = 0.95, draw_baseline = FALSE, aes(height = stat(density))) +
  # scale_y_log10(breaks=10^(-4:0),
  #               labels = sapply(-4:0,function(i){parse(text = sprintf("10^%d",i))})) +
  # facet_grid(cols = vars(tot_init_moi), labeller = labeller(tot_init_moi = function(x) paste0("MOI = ", x)),
  #            rows = vars(dens_at)) +
  facet_grid(cols = vars(tot_init_moi), labeller = labeller(tot_init_moi = function(x) paste0("MOI = ", x), 
                                                            prop_res = as_labeller(function(x) {
                                                              lapply(x, function(val) bquote(f[Res] == .(as.numeric(val))))
                                                            }, default = label_parsed)
                                                            ),
             rows = vars(prop_res)) +
  scale_fill_manual(labels = c("Sus.", "Res."),
                    values = rev(color_in),
                    name = "Genotype") +
  scale_alpha_continuous(range = c(0, 1), 
                         name = "Progeny \n Density",
                         limits = c(0, 1)
  ) +
  # ylab(expression("Initial " * italic(f)[Res])) +
  ylab("Density") +
  xlab("Res. Capsid Subunits") +
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
    panel.grid.major.x = element_line(color = "#F4F4F3"), ## Major x lines
    panel.grid.minor.x = element_blank(),                   ## Minor x lines
    panel.grid.major.y = element_line(color = "#F4F4F3"), ## Major y lines
    panel.grid.minor.y = element_blank(), ## Minor y lines
    ## For faceted plots: The strip is the top of the facet
    strip.background = element_blank(),
    strip.text = element_text(color = "black"),
    ## Remove the background of the plot and panel
    panel.background = element_rect(fill='transparent'), 
    plot.background = element_rect(fill='transparent', color=NA),
    ## Legend stuff
    # legend.position = "none",
    # legend.position = "inside",
    # legend.position.inside = c(0.1, 0.06),
    legend.text = element_text(size = legend_text_size),
    # legend.title = element_blank(),
    # legend.key.width = unit(0, 'cm'),
    legend.key.size = unit(0.3, "cm"),
    # legend.background = element_rect(fill = "white", color = "black"),  # Fills legend background & adds a border
    # legend.box.background = element_rect(color = "black", linewidth = 0),  # Adds an outer border
    legend.margin = margin(0.6, 0.6, 0.6, 0.6),
    # legend.spacing.y = unit(-0.0, "cm"),
    # legend.box = "horizontal"
  )

print("Data Plotted")

suppressWarnings(dir.create(paste0(filepath, "res/")))
suppressWarnings(dir.create(paste0(filepath, "res/sup/")))
suppressWarnings(dir.create(paste0(filepath, "res/sup/sfig2")))

ggsave(paste0(filepath, "res/sup/sfig2/sfig2.png"), sub_dens_b4_plot, h = height, w = width, units = "cm", bg = "transparent", dpi = 1200)
ggsave(paste0(filepath, "res/sup/sfig2/sfig2.pdf"), sub_dens_b4_plot, h = height, w = width, units = "cm", bg = "transparent", dpi = 1200)
saveRDS(sub_dens_b4_plot, file = paste0(filepath, "res/sup/sfig2/sfig2.rds"))

print("Data Saved")
