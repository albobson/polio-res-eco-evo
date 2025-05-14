#### Function to optimize (Note, I'm going to use foreach to parallelize the sims) ####
func_to_optim <- function(data, par, fit_func_in) {
  burst <- exp(par[1])                   ## Params[1] is the log(burst size)
  ppfu <- exp(par[2])                    ## Params[2] is the log(p2pfu)
  max_v_p_c <- par[3]                    ## I've removed this for now (it doesn't do anything inside of the function)
  
  ## Generate a dataframe to store data and a df to bind with
  df_test <- NULL
  df_it <- NULL
  
  ## Run simulations based on this fitness function
  df_test <- foreach(sus_moi = c(0, 5, 15, 50, 100),
                     .combine = 'rbind', 
                     .export = 'determ_polv',
                     .packages = c('dplyr')) %dopar% {
                       ## Run simulation
                       determ_polv(
                         n = 1,
                         moi_mut_start = 10,
                         moi_wt_start = sus_moi,
                         fit_func_in = fit_func_in,
                         v_prog = burst,
                         p2pfu = ppfu,
                         max_vpc = max_v_p_c
                       )
                     }
  
  ## Now we need to calculate the differences. First need to clean sim data
  df_test <- df_test %>%
    filter(time == 1, type == "resistant") %>%
    mutate(moi_mut = init_mut, 
           moi_wt = init_wt, 
           sim_pfu = surv_pfu, 
           sim_tot_pfu = tot_pfu) %>%
    dplyr::select(moi_mut, moi_wt, sim_pfu, sim_tot_pfu)
  
  
  df_test
  
  ## Now compare to the "real" data. Compare to resistant PFU first
  full_df <- merge(df_test, data)
  
  full_df$diff <- (log(full_df$pfu) - log(full_df$sim_pfu))^2
  
  sum_sqr_dif_log <- sum(full_df$diff)

  ## Saving the parameters and optimization variable to visualize later
  ## First updating counter
  cc <<- cc+1
  
  ## Now recoding the values
  vals[[cc]] <<- c(cc, par, sum_sqr_dif_log)
  
  print(paste0("iteration: ", cc, ", Val: ", sum_sqr_dif_log))
  return(sum_sqr_dif_log)
}

#### Function to plot bar plots ####
sim_bar_plot <- function(data, sim_name, plot_res_only = TRUE) {
  
  ## Clean up the data
  df_optim_clean <- data %>%
    select(type, time, init_wt, init_mut, wt_rep_abil, surv_pfu, tot_pfu, pop_prop) %>%
    filter(time == 1) %>%
    group_by(init_wt, wt_rep_abil) %>%
    filter(type == "resistant") %>%
    # mutate(sum_pfu = sum(surv_pfu)) %>%
    # filter(time == 1, type == "susceptible") %>%
    mutate(moi_wt = factor(init_wt), tot_pfu = as.numeric(tot_pfu)) %>%
    mutate(pocap = ifelse(wt_rep_abil == 1, "-", "+")) %>%
    mutate(pocap = factor(pocap, levels = c("+", "-"))) %>%
    mutate(sim = sim_name)
  
  
  ## If plot_res_only == TRUE, only plot the resistant virus
  
  ## Plot data
  optim_plot <- ggplot(df_optim_clean, aes(x = moi_wt, fill = sim)) +
    {if(plot_res_only)geom_bar(stat="identity", position = "dodge", mapping = aes(y = surv_pfu))}+
    {if(!plot_res_only)geom_bar(stat="identity", position = "dodge", mapping = aes(y = tot_pfu))}+
    scale_fill_manual(values = c("optim_mah_a24v" = "darkgreen",
                                 "optim_his_a24v" = "brown",
                                 "optim_his_i94f" = "purple",
                                 "Cell Culture" = "black")
    ) +
    theme_classic() + 
    labs(fill = "Type") +
    xlab("WT MOI") + ylab("Resistant Viral Yield (PFU/mL)") + 
    scale_y_continuous(trans="log10",
                       breaks=10^(0:10),
                       # expand = c(0.01, 0),
                       labels = sapply(0:10,function(i){parse(text = sprintf("10^%d",i))})) +
    coord_cartesian(ylim=c(10^4,10^log10(max(df_optim_clean$tot_pfu)))) +
    {if(plot_res_only)ggtitle(paste0("Resistant Virus - ", sim_name))}+
    {if(!plot_res_only)ggtitle(paste0("Total Virus - ", sim_name))}+
    theme(text = element_text(size=14), 
          axis.text = element_text(size=12),
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
          plot.title = element_text(hjust = 0.5)
          # panel.background = element_blank()
    ) 
  
  return(optim_plot)
  
}


#### Function to plot cell culture plots ####
cell_bar_plot <- function(data, title) {
  cell_df <- data %>%
    mutate(moi_wt = factor(moi_wt),
           surv_pfu = as.numeric(pfu)) %>%
    mutate(sim = "Cell Culture")
  
  ## Plot data
  cell_plot <- ggplot(cell_df, aes(x = moi_wt, fill = sim)) +
    geom_bar(stat="identity", position = "dodge", mapping = aes(y = surv_pfu)) +
    scale_fill_manual(values = c("optim_mah_a24v" = "darkgreen",
                                 "optim_his_a24v" = "brown",
                                 "optim_his_i94f" = "purple",
                                 "Cell Culture" = "black")
    ) +
    theme_classic() + 
    labs(fill = "Type") +
    xlab("WT MOI") + ylab("Resistant Viral Yield (PFU/mL)") + 
    scale_y_continuous(trans="log10",
                       breaks=10^(0:10),
                       # expand = c(0.01, 0),
                       labels = sapply(0:10,function(i){parse(text = sprintf("10^%d",i))})) +
    coord_cartesian(ylim=c(10^4,10^9)) +
    ggtitle(paste0("Cell Culture - ", title)) +
    theme(text = element_text(size=14), 
          axis.text = element_text(size=12),
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
          plot.title = element_text(hjust = 0.5)
          # panel.background = element_blank()
    ) 
  
  return(cell_plot)
  
}

#### Function to plot fitness functions ####
plot_single_fit_func <- function(fit_func_df, fit_df, title, line_color) {
  ggplot(fit_df, aes(x = subunits, y = prob_surv)) +
  geom_line(data = fit_func_df,
            aes(color = fit_type,
                group = fit_type),
            linewidth = 1) +
    scale_color_manual(values = line_color,
                       name = "Fit Type") +
  geom_point(color = "black") +
  xlab("Number of Resistant Subunits") +
  ylab("Probability of Survival") +
  scale_y_continuous(trans="log10",
                     limits = c(0.0001,1),
                     breaks=10^(-4:0),
                     labels = c(sapply(-4:-1,function(i){parse(text = sprintf("10^%d",i))}), "1 ")) +
    ggtitle(title) +
  theme_light() +
  theme(text = element_text(size = 20), axis.text = element_text(size=20), 
        axis.title = element_text(size = 22),
        legend.text = element_text(size = 20), 
        legend.title = element_text(size = 20),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"))
}
