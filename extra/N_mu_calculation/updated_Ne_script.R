#### Script to run the Ne extinction calculation                            ####


#### Set-Up                                                                 ####
suppressWarnings(suppressMessages({
  library(dplyr)
  library(ggplot2)
}))

setwd("~/PhD/dominant_drug_targets/book/250512_Ne/")


axis_text_size <- 6


width = 7
height = 5






## Read in the files
data <- read.csv(paste0("all_stoch-v2_var.csv"))

## Create a new results/ directory
dir.create(file.path("results/"), showWarnings = FALSE)

## Making sure everything that is a number is being read in as a number
numerics <- c("moi_type","moi_res","moi_wt","time",
              "init_mut","init_wt","pop_prop", "nu")

data[numerics] <- lapply(data[numerics], as.numeric)

## Find the extinction probabilities for each N*u
data_nu <- data %>%
  ## Since we're looking at _total_ extinction probability of the population I
  ## think we should remove one of the virus types since we already have a
  ## tot_pfu column.
  filter(time != 0, type == "resistant") %>%
  select(nu, id, tot_pfu) %>%
  group_by(nu) %>%
  reframe(
    n_u = nu,
    total_rows = n(),
    zero_count = sum(tot_pfu == 0),
    proportion_zero = zero_count / total_rows
  )


write.csv(data_nu, "results/extinction_probabilities.csv")

breaks10 <- 10^seq(min(log10(data_nu$n_u)), max(log10(data_nu$n_u)), by = 1)

#### Plotting ####
# p <- ggplot(data_nu, aes(x = n_u, y = proportion_zero)) +
#   geom_line() +
#   theme_light() +
#   xlab("N * u") + ylab("Extinction Probability") +
#   # scale_x_log10(breaks = breaks10) +
#   scale_x_continuous(trans="log10",breaks=breaks10,
#                      labels = sapply(seq(min(log10(data_nu$n_u)), max(log10(data_nu$n_u)), by = 1),
#                                      function(i){parse(text = sprintf("10^%d",i))})) +
#   coord_cartesian(ylim=c(0,1)) +
#   labs(title = paste0("Stochastic (V2) Extinction Probability")) +
#   theme(axis.text = element_text(size = 13))
# 
# print(p)
# 
# ggsave(filename = "results/result_stoch-v2_n-u.png",
#        width = 16, height = 9, units = "in")

# #### MOI Specifically ####
# ## Find the extinction probabilities for each N*u
# data_m <- data %>%
#   ## Since we're looking at _total_ extinction probability of the population I
#   ## think we should remove one of the virus types since we already have a
#   ## tot_pfu column.
#   filter(time != 0, type == "susceptible") %>%
#   select(nu, id, tot_pfu, moi_type, init_wt) %>%
#   group_by(nu) %>%
#   reframe(
#     moi = init_wt,
#     total_rows = n(),
#     zero_count = sum(tot_pfu == 0),
#     proportion_zero = zero_count / total_rows
#   )
# 
# breaks10 <- 10^seq(min(log10(data_m$moi)), max(log10(data_m$moi)), by = 1)
# 
# #### Plotting ####
# p <- ggplot(data_m, aes(x = moi, y = proportion_zero)) +
#   geom_line() +
#   theme_light() +
#   xlab("Initial Susceptible MOI") + ylab("Extinction Probability") +
#   geom_vline(aes(xintercept = 0.000153678, color = "red"), linetype = "dashed") +
#   annotate(geom="text", x=0.000953678, y=.15, label="WT 'Replication Ability'",
#            color="red") +
#   # scale_x_log10(breaks = breaks10) +
#   scale_x_continuous(trans="log10",breaks=breaks10,
#                      labels = sapply(seq(min(log10(data_m$moi)), max(log10(data_m$moi)), by = 1),
#                                      function(i){parse(text = sprintf("10^%f",i))})) +
#   coord_cartesian(ylim=c(0,1)) +
#   labs(title = paste0("Stochastic (V2) Extinction Probability - MOI")) +
#   theme(axis.text = element_text(size = 13))
# 
# print(p)
# 
# ggsave(filename = "results/result_stoch-v2_moi.png",
#        width = 16, height = 9, units = "in")


## For loop for each time step

max_run = max(data$id)

numerics <- c("moi_type","moi_res","moi_wt", "time", "id",
              "init_mut","init_wt","pop_prop", "nu")

data[numerics] <- lapply(data[numerics], as.numeric)

# data$time <- as.character(data$time)
# data$id <- as.character(data$id)

## In order to account for time points which were not included in the data
## (since the simulations ended early). For this reason, I am generating a
## dummy data frame with all of the possible time points and then will left
## join my data frame to incorporate the missing rows into my calculation.
time_dum <- expand.grid(time = 0:max(data$time), id = 0:max(data$id), nu = unique(data$nu))

## Add the initial MOI data to the time_dum df
moi_nu <- data %>%
  select(nu, init_wt) %>%
  distinct()

time_dum <- left_join(time_dum, moi_nu, by = "nu")

for(i in 0:max(data$time)) {
  ## Filtering the dummy df
  time_dum_temp <- time_dum %>% filter(time == i)
  
  ## Find the extinction probabilities for each N*u
  temp <- data %>%
    ## Since we're looking at _total_ extinction probability of the population I
    ## think we should remove one of the virus types since we already have a
    ## tot_pfu column. Then, filtering out anything that is 0 or NA.
    filter(time == i, type == "susceptible", tot_pfu > 0) %>%
    select(time, nu, id, tot_pfu)
  ## Now joining together the temp df and the time_dummy_df
  result <- left_join(time_dum_temp, temp, by = c("time", "id", "nu")) %>%
    group_by(nu) %>%
    reframe(
      moi = init_wt,
      total_still_alive = sum(!is.na(tot_pfu)),
      proportion_zero = 1 - (total_still_alive / max_run)
    ) %>%
    distinct() %>%
    ungroup()
  
  write.csv(temp, paste0("results/extinction_probabilities_t", i, ".csv"))
  
  breaks10nu <- 10^seq(min(log10(data$nu)), max(log10(data$nu)), by = 1)
  
  #### Plotting ####
  nup <- ggplot(result, aes(x = nu, y = proportion_zero)) +
    geom_line() +
    theme_light() +
    xlab(expression("Host cell population size" ~ "\u00d7" ~ "mutation rate (" * gamma %*% mu * ")")) + ylab("Extinction Probability") +
    scale_x_continuous(trans="log10",breaks=breaks10nu,
                       labels = sapply(seq(min(log10(data$nu)), max(log10(data$nu)), by = 1),
                                       function(i){parse(text = sprintf("10^%d",i))})) +
    coord_cartesian(ylim=c(0,1)) +
    # labs(title = paste0(i, " Replications")) +
    theme_light() +
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
      # panel.grid.minor.x = element_blank(),                   ## Minor x lines
      panel.grid.major.y = element_line(color = "lightgrey"), ## Major y lines
      panel.grid.minor.y = element_line(color = "lightgrey"), ## Minor y lines
      ## For faceted plots: The strip is the top of the facet
      strip.background = element_blank(),
      strip.text.x = element_blank(),
      ## Remove the background of the plot and panel
      panel.background = element_rect(fill='transparent'), 
      plot.background = element_rect(fill='transparent', color=NA),
      ## Legend stuff
      legend.position = "none",
      # legend.position = "inside",
      # legend.position.inside = c(0.18, 0.06),
      # legend.text = element_text(size = legend_text_size),
      # legend.title = element_blank(),
      # legend.background = element_rect(fill = "transparent"),
      # legend.key = element_rect(fill = "transparent"),
      # legend.key.width = unit(0, 'cm'),
      # legend.key.size = unit(0.3, "cm")
    ) 
  
  ggsave(plot = nup,
         filename = paste0("results/result_t", i,"_n-u.pdf"),
         w = width, h = height, units = "cm", bg = "transparent", dpi = 1200)
  
  # breaks10moi <- 10^seq(min(log10(result$moi)), max(log10(result$moi)), by = 1)
  # 
  # moip <- ggplot(result, aes(x = moi, y = proportion_zero)) +
  #   geom_line() +
  #   theme_light() +
  #   xlab("Initial Susceptible MOI") + ylab("Extinction Probability") +
  #   # geom_vline(aes(xintercept = 0.000153678, color = "red"), linetype = "dashed") +
  #   # annotate(geom="text", x=0.000953678, y=.15, label="WT 'Replication Ability'",
  #   #          color="red") +
  #   scale_x_continuous(trans="log10",breaks=breaks10moi,
  #                      labels = sapply(seq(min(log10(result$moi)), max(log10(result$moi)), by = 1),
  #                                      function(i){parse(text = sprintf("10^%f",i))})) +
  #   coord_cartesian(ylim=c(0,1)) +
  #   labs(title = paste0("Stochastic (V2) Extinction Probability, ", i, " Replications - MOI")) +
  #   theme(axis.text = element_text(size = 13))
  # 
  # ggsave(plot = moip,
  #        filename = paste0("results/result_stoch-v2_t", i,"_moi.png"),
  #        width = 16, height = 9, units = "in")
  
}

# suppressWarnings(dir.create(paste0(filepath, "res/")))
# suppressWarnings(dir.create(paste0(filepath, "res/fig3/")))
# 
# ggsave(paste0(filepath, "res/fig3/B.png"), func_plot, h = height/2, w = width/2, units = "cm", bg = "transparent", dpi = 1200)
# ggsave(paste0(filepath, "res/fig3/B.pdf"), func_plot, h = height/2, w = width/2, units = "cm", bg = "transparent", dpi = 1200)
# ggsave(paste0(filepath, "res/fig3/B.svg"), func_plot, h = height/2, w = width/2, units = "cm", bg = "transparent", dpi = 1200, device = 'svg')
# saveRDS(func_plot, file = paste0(filepath, "res/fig3/B.rds"))


