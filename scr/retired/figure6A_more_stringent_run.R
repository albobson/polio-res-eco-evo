############################## Figure 6 A ######################################

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


#### Run simulations                                                        ####
curr_run_a <- determ_polv(n = 6,
                        moi_mut_start = 0, 
                        moi_wt_start = 100,
                        fit_func_in = fit_func$prob_surv,
                        v_prog = optim_params$optim_v_prog,
                        p2pfu = optim_params$optim_p2pfu
) 


#### Plotting trajectory                                                    ####
plot6a <- ggplot(curr_run_a, aes(x = time, y = moi_type, color = type)) +
  geom_point(size = 3) + 
  geom_line(linewidth = 1) +
  scale_color_manual(labels = c("Resistant", "Susceptible"),
                     values = c("#c94d4d", "#688cc8"),
                     name = "Genotype") +
  theme_light() + 
  xlab("Passages") + ylab("MOI") + 
  scale_y_log10(breaks=10^(-6:10),
                labels = sapply(-6:10,function(i){parse(text = sprintf("10^%d",i))}),
                limits = c(10^-5, 10^2.4)) +
  scale_x_continuous(breaks = c(0:max(curr_run_a$time))) +
  theme(text = element_text(size = 20), axis.text = element_text(size=20), 
        axis.title = element_text(size = 22),
        legend.text = element_text(size = 20), 
        legend.title = element_blank(),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        legend.position=c(0.85,0.2))

suppressWarnings(dir.create(paste0(filepath, "res/fig6")))

ggsave(paste0(filepath, "res/fig6/A.png"), plot6a, h = 8, w = 8, bg = "transparent")
ggsave(paste0(filepath, "res/fig6/A.jpg"), plot6a, h = 8, w = 8, bg = "transparent")
ggsave(paste0(filepath, "res/fig6/A.svg"), plot6a, h = 8, w = 8, bg = "transparent", device = 'svg')
saveRDS(plot6a, paste0(filepath, "res/fig6/A.rds"))