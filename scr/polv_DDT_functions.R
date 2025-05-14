#### Document Timeline ####

## 230823 - Retired old document, which used nested for loops and arrays to
##          generate probabilities. I figured out a way to make things run
##          faster by NOT using for loops, so I'm updating all of my code
##          accordingly. I am also going to clean up some of the functions
##          that I've wanted to clean up for some time now.

## 240418 - As the project has changed over time, I have honed in on specific
##          models and functions. This script will clean up all of my various
##          models that I've created over the last few months and will be,
##          hopefully, much easier to read.
##          In particular, I am simplifying the core function to perform all of
##          the various functions that I had spread across multiple functions. 

## 240925 - Keeping only necessary functions.


#### TO DO ####

    ## Pull out the fitness vector from the model? Or perhaps, make it optional,
    ## where a different fitness vector could be given instead?
        ## 240925 - Done

#### Libraries ####
## Should I load each of these within each function? I can't decide.
## General use
library(dplyr)
library(tidyr)

## Parallelization
library(foreach)
library(doParallel)
library(parallel)

## Plotting
library(ggplot2)
library(cowplot)

#### Logistic fitness function ####
fit_func_logistic <- function(min_val, k, subunits, mid) {
    ((1-min_val)/(1+exp(-k*((0:subunits)-mid))))+min_val
}


#### Deterministic core model ####
## A function which calculates the change in population of each virus type over 
## time. 

## Where:
## n = the number of generations
## ld = The # of susceptible capsid subunits needed to neutralize the virus
## moi_wt_start = A vector of initial  wild type MOIs to test
## moi_mut_start = A vector of initial mutant MOIs to test
## c_pop = The cell population (refreshed at each time step)
## v_prog = The number of viral particles produced with each infection
## p2pfu = The particle to PFU ratio of the virus
## subunits = The number of capsid subunits of the virus
## bg_mutat = The mutation rate of the virus between wt and mut and vice versa
## max_vpc = The maximum number of viruses that can enter a cell.
## fit_func_in = Vector with survival probabilities associated with capsid
##               subunit composition
## wt_rep_abil = The replication ability of the wild type virus in the presence 
##               of pocapavir. The first of the fitness vector
## report_subunit_dist = TRUE returns a breakdown of viral population stratified
##                       by subunit composition

determ_polv <-
  function(n = 1,
           moi_wt_start = 1,
           moi_mut_start = 1,
           c_pop = 2.9 * 10 ^ 6,
           v_prog = 5243.65,
           p2pfu = 24.46,
           subunits = 60,
           bg_mutat = 10000,
           max_vpc = 500,
           wt_rep_abil,
           fit_func_in,
           report_subunit_dist = FALSE ## TRUE returns a breakdown of viral population stratified by subunit composition
           ) { 
    
    #### Initial set up ####
    ## Condition counter for progress bar
    cc <- 1
    ## progress bar to report progress
    pb <- txtProgressBar(min = 0, max = n, initial = 0, style = 3)
    df_dens <- NULL
    
    ## Loop through up to n generations
    for (t in 1:n) {
      ## If n=1, then use the starting MOI, otherwise, determine MOI based
      ## on the previous run
      if (t == 1) {
        moi_wt <- moi_wt_start
        moi_mut <- moi_mut_start
      } else {
        moi_wt <- new_wt_v / c_pop
        moi_mut <- new_mut_v / c_pop
      }
      
      ## Population of the susceptible and mutant virus
      v_pop_wt <- c_pop * moi_wt
      v_pop_mut <- c_pop * moi_mut
      v_pop_tot <- v_pop_wt + v_pop_mut
      
      ## Calculating the total MOI
      moi_tot <- v_pop_tot / c_pop
      
      ## Calculating the proportion of susceptible to resistant virus
      v_prop <- v_pop_mut / v_pop_tot
      
      ## The range of the "true MOI" that we're interested in (viruses per
      ## cell)
      ## A function of MOI tot
      vpc <- qpois(0.99999999, (moi_tot + 1))
      if (vpc > max_vpc) {
        vpc <- max_vpc
      }
      
      #### Data Frames ####
      ## Create a df with all of the combinations of viruses/focal genomes
      df_res <- expand.grid(viruses_per_cell = 1:vpc, 
                            res_per_cell = 0:vpc, 
                            subs = 0:subunits) %>%
        filter(viruses_per_cell >= res_per_cell) %>%
        mutate(type = "res") %>%
        arrange(viruses_per_cell, res_per_cell, subs)
      
      df_sus <- expand.grid(viruses_per_cell = 1:vpc, 
                            sus_per_cell = 0:vpc, 
                            subs = 0:subunits) %>%
        filter(viruses_per_cell >= sus_per_cell) %>%
        mutate(type = "sus") %>%
        arrange(viruses_per_cell, sus_per_cell, subs)
      
      ## Fill out the DF with the initial infection probabilities
      
      ###### Resistant
      ## dpois(df_res$viruses_per_cell, moi_tot) = the probability of a cell
      ## being infected with `viruses_per_cell` number of total viruses
      df_res$prob_inf_w_n_v <- dpois(df_res$viruses_per_cell, moi_tot)
      
      ## dbinom(df_res$res_per_cell, df_res$viruses_per_cell, v_prop) = The
      ## probability of a cell being infected with `res_per_cell` mutant
      ## genomes, given it was infected with `viruses_per_cell` total viruses
      df_res$prob_inf_w_n_r <- 
        dbinom(df_res$res_per_cell, df_res$viruses_per_cell, v_prop)
      
      ## dbinom(df_res$subs, subunits,
      ## (df_res$res_per_cell/df_res$viruses_per_cell)) = the probability of
      ## generating a virus with `subs` subunits, given it was infected with
      ## `viruses_per_cell` viruses, `res_per_cell` of them being resistant
      df_res$prob_prog_w_n_res_subs <- 
        dbinom(df_res$subs, subunits, 
               (df_res$res_per_cell/df_res$viruses_per_cell))
      
      ## (df_res$res_per_cell/df_res$viruses_per_cell) + ((1 -
      ## (df_res$res_per_cell / df_res$viruses_per_cell)) / bg_mutat) -
      ## ((df_res$res_per_cell / df_res$viruses_per_cell) / bg_mutat) = the
      ## probability of a resistant genome existing within a cell, given that
      ## there were `res_per_cell` number of resistant genomes,
      ## `viruses_per_cell` total viruses, and that the frequency of conversion
      ## between the two genotypes is `bg_mutat`
      df_res$prob_res_genome <- 
        ## Probability w/o mutation
        (df_res$res_per_cell/df_res$viruses_per_cell) +
        ## Plus the frequency of susceptible mutating to resistant
        ((1 - (df_res$res_per_cell / df_res$viruses_per_cell)) / bg_mutat) -
        ## Minus the frequency of resistance mutating back to susceptible
        ((df_res$res_per_cell / df_res$viruses_per_cell) / bg_mutat)
      
      ## The total probability of a virus with a resistant genome being produced
      ## by a cell with the given infection composition, wrapped with the
      ## specific number of resistant subunits is the intersection of all of the
      ## prior probabilities
      df_res$tot_prob <- df_res$prob_inf_w_n_v * df_res$prob_inf_w_n_r *
        df_res$prob_prog_w_n_res_subs * df_res$prob_res_genome
      
      #### Susceptible
      ## Note: Since I documented the version above so well, I am going to skip
      ## documenting here. Basically everything that I wrote before holds true
      ## here unless noted.
      df_sus$prob_inf_w_n_v <- dpois(df_sus$viruses_per_cell, moi_tot)
      
      df_sus$prob_inf_w_n_r <- 
        dbinom(df_sus$sus_per_cell, df_sus$viruses_per_cell, (1 - v_prop))
      
      df_sus$prob_prog_w_n_sus_subs <- 
        dbinom(df_sus$subs, subunits, 
               (1-(df_sus$sus_per_cell/df_sus$viruses_per_cell)))
      
      df_sus$prob_sus_genome <- 
        (df_sus$sus_per_cell/df_sus$viruses_per_cell) +
        ((1 - (df_sus$sus_per_cell / df_sus$viruses_per_cell)) / bg_mutat) -
        ((df_sus$sus_per_cell / df_sus$viruses_per_cell) / bg_mutat)
      
      ## Total probability of a susceptible genome wrapped in a capsid with some
      ## number of resistant subunits, given the initial infection composition
      ## and the system's MOI
      df_sus$tot_prob <- df_sus$prob_inf_w_n_v * df_sus$prob_inf_w_n_r *
        df_sus$prob_prog_w_n_sus_subs * df_sus$prob_sus_genome
      
      ## Finding the number of uninfected cells
      dist <- dpois(0:vpc, moi_tot)
      cell_dist <- dist * c_pop
      
      #### Total number of viruses #### 
      ## Now, multiply the probabilities by the total number of cells, the burst
      ## size per cell, and divide by the particle to plaque forming unit ratio.
      
      ## NOTE: The sum of the probabilities of both data frames will be close
      ## to, but not equal to 1. This is because there are cells that will be
      ## infected by 0 viruses. Because the probability of being infected by 0
      ## viruses has been removed from the data frames, we can multiply by the
      ## total cell population and burst size.
      
      df_res$tot_virus <- df_res$tot_prob * (c_pop * (v_prog/p2pfu))
      df_sus$tot_virus <- df_sus$tot_prob * (c_pop * (v_prog/p2pfu))
      
      #### Drug action ####
      ## fit_func_in is our fitness probability associated with 
      
      ## Now pairing these values up with their associated subunits
      subsun <- 0:subunits
      pocap_vec <- data.frame(fit_func = fit_func_in, subs = subsun)
      
      ## Now multiply the corresponding fitness value for each subunit by the
      ## number of viruses with those subunits. Left join combines the data
      ## frames by the subunits columns. Then I create the multiplied vector,
      ## drop everything else, and assign it as the new column.
      df_res$tot_virus_after_pocap <- left_join(df_res, pocap_vec,
                                                by = join_by(subs)) %>%
        mutate(tot_virus_after_pocap = tot_virus*fit_func) %>%
        select(tot_virus_after_pocap)
      
      df_sus$tot_virus_after_pocap <- left_join(df_sus, pocap_vec,
                                                by = join_by(subs)) %>%
        mutate(tot_virus_after_pocap = tot_virus*fit_func) %>%
        select(tot_virus_after_pocap)
      
      ## Then, assuming that the mutation has no effect on the binding
      ## affinity or the ability to infect the next cell, we can assume
      ## that the sum of these arrays, is all of the viruses that are able
      ## to escape.
      
      ## Total virus before pocap applied
      new_mut_v_b4_pocap <- sum(df_res$tot_virus)
      new_wt_v_b4_pocap <- sum(df_sus$tot_virus)
      
      ## Virus after pocap applied
      new_mut_v <- sum(df_res$tot_virus_after_pocap)
      new_wt_v <- sum(df_sus$tot_virus_after_pocap)
      
      ## The first place in our cell_infection distribution is the number
      ## of cells infected by 0 viruses. These are the "surviving" cells,
      ## assuming all cells die once infected.
      cell_surv <- cell_dist[1]
      
      #### Checking work ####
      ## Double checking that all of the probabilities add up to 1.
      sumsmut <- sum(df_res$tot_prob)
      sumswt <- sum(df_sus$tot_prob)
      
      #### Tracking capsid subunit compositions ####
      ## resistant virus
      res_sub_df <- df_res %>%
        select(subs, tot_virus, tot_virus_after_pocap) %>%
        group_by(subs) %>%
        summarize(dens_b4_pocap = sum(tot_virus)/new_mut_v_b4_pocap,
                dens = sum(tot_virus_after_pocap)/new_mut_v) %>%
        mutate(time = t, type = "resistant", 
               type_moi = moi_mut, tot_moi = moi_tot,
               init_moi_res = moi_mut_start,
               init_moi_sus = moi_wt_start)
        
      ## Susceptible virus
      sus_sub_df <- df_sus %>%
        select(subs, tot_virus, tot_virus_after_pocap) %>%
        group_by(subs) %>%
        summarize(dens_b4_pocap = sum(tot_virus)/new_wt_v_b4_pocap,
                  dens = sum(tot_virus_after_pocap)/new_wt_v) %>%
        mutate(time = t, type = "susceptible", 
               type_moi = moi_wt, tot_moi = moi_tot,
               init_moi_res = moi_mut_start,
               init_moi_sus = moi_wt_start)
      
      sub_dens_df <- rbind(res_sub_df, sus_sub_df)
      
      ## Total virus
      tot_sub_df <- sub_dens_df %>%
        group_by(subs) %>%
        summarize(dens_b4_pocap = sum(dens_b4_pocap),
                  dens = sum(dens)) %>%
        mutate(time = t, 
               type = "total",
               type_moi = moi_tot, tot_moi = moi_tot,
               init_moi_res = moi_mut_start,
               init_moi_sus = moi_wt_start) %>%
        ungroup() %>%
        mutate(dens_b4_pocap = dens_b4_pocap/sum(dens_b4_pocap),
               dens = dens/sum(dens))
      
      ## All together
      sub_dens_df <- rbind(sub_dens_df, tot_sub_df)
      
      #### Generating dfs ####
      
      ## Generating a long dataframe to hold all of our information
      ## if the time == 1, generate a new dataframe
      if(t==1) {
        df_mut <- data.frame(type = "resistant",          ## resistant virus
                             time = 0,                    ## time
                             moi_type = c(moi_mut_start),
                             moi_res = c(moi_mut_start),
                             moi_wt = c(moi_wt_start),
                             init_wt = moi_wt_start,
                             init_mut = moi_mut_start,
                             c_pop = c(c_pop),
                             subunits = subunits,
                             v_prog = c(v_prog),
                             p2pfu = c(p2pfu),
                             bg_mut_rate = bg_mutat,
                             wt_rep_abil = fit_func_in[1],
                             uninf_cells = c(cell_surv),
                             surv_pfu = moi_mut_start*c_pop,
                             tot_pfu = (moi_wt_start*c_pop+moi_mut_start*c_pop),
                             pop_prop = moi_mut_start/(moi_mut_start+moi_wt_start),
                             check_sums = (sum(sumsmut)+sum(sumswt)+dist[1]))
        
        df_wt <- data.frame(type = "susceptible",
                            time = 0,
                            moi_type = c(moi_wt_start),
                            moi_res = c(moi_mut_start),
                            moi_wt = c(moi_wt_start),
                            init_wt = moi_wt_start,
                            init_mut = moi_mut_start,
                            c_pop = c(c_pop),
                            subunits = subunits,
                            v_prog = c(v_prog),
                            p2pfu = c(p2pfu),
                            bg_mut_rate = bg_mutat,
                            wt_rep_abil = fit_func_in[1],
                            uninf_cells = c(cell_surv),
                            surv_pfu = moi_wt_start*c_pop,
                            tot_pfu = (moi_wt_start*c_pop+moi_mut_start*c_pop),
                            pop_prop = moi_wt_start/(moi_mut_start+moi_wt_start),
                            check_sums = (sum(sumsmut)+sum(sumswt)+dist[1]))
        
        df_long <- rbind(df_mut, df_wt)
        
        ## Subunit densities
        ## resistant virus
        res_sub_df0 <- data.frame(time = 0,
                                 type = "resistant",
                                 subs = 0:subunits, 
                                 ## Density of the population
                                 dens = c(rep(0, 60), 1),
                                 dens_b4_pocap = c(rep(0, 60), 1),
                                 type_moi = moi_mut_start, 
                                 tot_moi = moi_mut_start+moi_wt_start,
                                 init_moi_res = moi_mut_start,
                                 init_moi_sus = moi_wt_start)
            ## If there is 0 to start, account for that
        if(moi_mut_start == 0){
          res_sub_df0 <- data.frame(time = 0,
                                    type = "resistant",
                                    subs = 0:subunits, 
                                    ## Density of the population
                                    dens = c(rep(0, 60), 0),
                                    dens_b4_pocap = c(rep(0, 60), 0),
                                    type_moi = moi_mut_start, 
                                    tot_moi = moi_mut_start+moi_wt_start,
                                    init_moi_res = moi_mut_start,
                                    init_moi_sus = moi_wt_start)
        }
        
        sus_sub_df0 <- data.frame(time = 0, 
                                 type = "susceptible",
                                 subs = 0:subunits, 
                                 ## Density of the population
                                 dens = c(1, rep(0, 60)),
                                 dens_b4_pocap = c(1, rep(0, 60)),
                                 type_moi = moi_wt_start, tot_moi = moi_mut_start+moi_wt_start,
                                 init_moi_res = moi_mut_start,
                                 init_moi_sus = moi_wt_start)
        ## If there is 0 to start, account for that
        if(moi_wt_start == 0){
          sus_sub_df0 <- data.frame(time = 0,
                                    type = "susceptible",
                                    subs = 0:subunits, 
                                    ## Density of the population
                                    dens = c(rep(0, 60), 0),
                                    dens_b4_pocap = c(rep(0, 60), 0),
                                    type_moi = moi_wt_start, tot_moi = moi_mut_start+moi_wt_start,
                                    init_moi_res = moi_mut_start,
                                    init_moi_sus = moi_wt_start)
        }
        
        tot_sub_df0 <- data.frame(time = 0, 
                                 type = "total",
                                 subs = 0:subunits, 
                                 ## Density of the population
                                 dens =  c(1*(1-moi_wt_start/(moi_wt_start+moi_mut_start)), 
                                           rep(0, 59), 1*(moi_wt_start/(moi_wt_start+moi_mut_start))),
                                 dens_b4_pocap = c(1*(1-moi_wt_start/(moi_wt_start+moi_mut_start)), 
                                                   rep(0, 59), 1*(moi_wt_start/(moi_wt_start+moi_mut_start))),
                                 type_moi = moi_mut_start+moi_wt_start, tot_moi = moi_mut_start+moi_wt_start,
                                 init_moi_res = moi_mut_start,
                                 init_moi_sus = moi_wt_start)
        ## All together
        sub_dens_df0 <- rbind(res_sub_df0, sus_sub_df0, tot_sub_df0)
        sub_dens_df <- rbind(sub_dens_df0, sub_dens_df)
      }
      
      ## Create a new row to the current dataframe
      new_row_mut <- c("resistant",
                       t,
                       (new_mut_v/c_pop),
                       (new_mut_v/c_pop),
                       (new_wt_v/c_pop),
                       moi_wt_start,
                       moi_mut_start,
                       c_pop,
                       subunits,
                       v_prog,
                       p2pfu,
                       bg_mutat,
                       fit_func_in[1],
                       cell_surv,
                       new_mut_v,
                       new_mut_v+new_wt_v,
                       new_mut_v/(new_wt_v+new_mut_v),
                       (sum(sumsmut)+sum(sumswt)+dist[1]))
      
      new_row_wt <- c("susceptible",
                      t,
                      (new_wt_v/c_pop),
                      (new_mut_v/c_pop),
                      (new_wt_v/c_pop),
                      moi_wt_start,
                      moi_mut_start,
                      c_pop,
                      subunits,
                      v_prog,
                      p2pfu,
                      bg_mutat,
                      fit_func_in[1],
                      cell_surv,
                      new_wt_v,
                      new_mut_v+new_wt_v,
                      new_wt_v/(new_wt_v+new_mut_v),
                      (sum(sumsmut)+sum(sumswt)+dist[1]))
      
      ##bind to new dataframe
      df_long <- rbind(df_long, new_row_mut, new_row_wt)
      df_dens <- rbind(df_dens, sub_dens_df)
      
      ## Update progress bar
      setTxtProgressBar(pb,cc)
      
      ## Increase the condition counter by 1 since condition has been finished
      cc <- cc+1
    }
    
    close(pb)
    
    ## Fixing variables ##
    numericsl <- c("moi_type","moi_res","moi_wt","c_pop","wt_rep_abil",
                   "v_prog","p2pfu","uninf_cells",
                   "check_sums","time", "bg_mut_rate",
                   "init_mut","init_wt","surv_pfu")
    df_long[numericsl] <- lapply(df_long[numericsl], as.numeric)
    
    ## If you just want trajectories, just return trajectories
    if(report_subunit_dist == FALSE){
       return(df_long)     
    }
    
    ## If you want trajectories and subunit distributions, return both as a list
    if(report_subunit_dist == TRUE){
      return(df_dens)
    }
}


#### Stochastic core model ####
stoch_polv <- function(n = 1,
                       moi_wt_start = 1,
                       moi_mut_start = 1,
                       t_pocap = 0,
                       imm_delay = 9,
                       imm_m = -1.6,
                       imm_sd = 0.5,
                       c_pop = 2.9 * 10 ^ 6,
                       v_prog = 5244,
                       p2pfu = 25,
                       subunits = 60,
                       bg_mutat = 10000,
                       max_vpc = 1000, ## Set at a number higher than ever experienced in a simulation for now. This is a factor we decided not to incorporate into this model due to a lack of empirical data
                       id = 1,
                       fit_func_in,
                       seed_in = round(runif(1, 0, 10000))) {
  #### Initial set up ####
  ## Set the seed:
  set.seed(seed_in)
  
  ## Condition counter (only usable when running the function in sequence)
  cc <- 1
  ## progress bar to report progress
  pb <- txtProgressBar(min = 0, max = n, initial = 0, style = 3)
  
  ## Immune sensitivity
  immune_sensitivity <- rlnorm(1, imm_m, imm_sd) 
  
  ## Loop through up to n generations
  for (t in 1:n) {
    ## If n=1, then use the starting MOI, otherwise, determine MOI based
    ## on the previous run
    if (t == 1) {
      moi_wt <- moi_wt_start
      moi_mut <- moi_mut_start
    } else {
      moi_wt <- new_wt_v / c_pop
      moi_mut <- new_mut_v / c_pop
    }
    
    ## Reset the fitness function
    curr_fit_func <- fit_func_in
    
    ## Population of the susceptible and mutant virus
    v_pop_wt <- c_pop * moi_wt
    v_pop_mut <- c_pop * moi_mut
    v_pop_tot <- v_pop_wt + v_pop_mut
    
    ## Calculating the total MOI
    moi_tot <- v_pop_tot / c_pop
    
    ## Calculating the proportion of susceptible to resistant virus
    v_prop <- v_pop_mut / v_pop_tot
    
    ## The range of the "true MOI" that we're interested in (viruses per
    ## cell)
    ## A function of MOI tot
    vpc <- qpois(0.99999999, (moi_tot + 1))
    if (length(vpc) > max_vpc) {
      vpc <- 0:max_vpc
    }
    
    #### Initial data frames ####
    df <- data.frame(viruses_per_cell = 0:vpc)
    
    ## Initial infection probabilities ##
    ## dpois(df_res$viruses_per_cell, moi_tot) = the probability of a cell
    ## being infected with `viruses_per_cell` number of total viruses
    df$prob_inf_w_n_v <- dpois(df$viruses_per_cell, moi_tot)
    
    ## Random pulls across the c_population given the infection probs
    df$inf_w_n_v = rmultinom(n = 1, size = c_pop, prob = df$prob_inf_w_n_v)
    
    ## We can now remove the 0s from our dataframe, since they have no viruses
    df_0s <- df[1,] ## Saving for later, just in case desired.
    
    df <- df[-1,]
    
    ## Creating vectors with the # of viruses per cell
    df$res_per_cell <- lapply(df$viruses_per_cell, function(x) 0:x)
    
    ## Generate the probability of having the exact number of res viruses in
    ## a cell given some number of total viruses.
    f_gen_vf_prob <- function(x,y) c(dbinom(x, y, v_prop))
    
    df$prob_inf_w_n_r <- mapply(f_gen_vf_prob, df$res_per_cell, df$viruses_per_cell)
    
    has.na <- sapply(df$prob_inf_w_n_r, function(x) any(is.na(x)))
    
    if(any(has.na) == TRUE) break ## Check that we have viruses
    
    ## Given the probability of having the exact number of resistant viruses
    ## in each cell with a given number of total viruses, pull the number of
    ## viruses across a multinomial distribution, where the size is the total
    ## number of cells in that res/total category
    f_gen_vf_tot <- function(x, y) rmultinom(n = 1, size = x, prob = y)
    
    df$tot_inf_w_n_r <- mapply(f_gen_vf_tot, df$inf_w_n_v, df$prob_inf_w_n_r)
    
    ## Now that we have these, we want to unnest the list of viral populations
    ## in order to calculate our next probabilities easier
    df <- df %>% tidyr::unnest(cols = c(res_per_cell, prob_inf_w_n_r, tot_inf_w_n_r))
    
    ## Finding how many total viruses in the cell by generating a random burst
    ## size per cell and then summing all of the bursts (for each category)
    f_bursts <- function(x) sum(rpois(x,v_prog))
    
    df$tot_v_in_cells <- mapply(f_bursts, df$tot_inf_w_n_r)
    
    ## Assuming that without mutation, both viruses replicate at the same
    ## rate. Using a binomial distribution, with the probability being the
    ## ratio of resistant to susceptible in the cell to start. Then,
    ## subtracting the resistant viruses from the total, giving the
    ## susceptible viruses
    f_gen_res <- function(x, y, z) rbinom(n = 1, size = x, prob = (y/z))
    
    df$tot_r_b4_mut <- mapply(f_gen_res, x = df$tot_v_in_cells,
                              y = df$res_per_cell,
                              z = df$viruses_per_cell)
    
    df$tot_s_b4_mut <- df$tot_v_in_cells - df$tot_r_b4_mut
    
    ## Doing mutation to find how many susceptible and how many resistant
    ## viruses there are in the cell total
    ## Subtracting the total from the probability that it STAYS the same. This
    ## way, I can use a really small bg_mutat without any problem
    f_prob_mutat <- function(x) rbinom(1, size = x, prob = (1/bg_mutat))
    
    ## Finding the number of mutants for each virus type in each category
    df$n_res_2_sus <- unlist(lapply(df$tot_r_b4_mut, f_prob_mutat))
    
    df$n_sus_2_res <- unlist(lapply(df$tot_s_b4_mut, f_prob_mutat))
    
    ## Adding/subtracting mutants from total
    df$tot_r_b4_faulty <- df$tot_r_b4_mut - df$n_res_2_sus + df$n_sus_2_res
    
    df$tot_s_b4_faulty <- df$tot_s_b4_mut - df$n_sus_2_res + df$n_res_2_sus
    
    ## Finding the number of faulty viruses (P:PFU ratio) from the total
    f_faulty <- function(x) x - rbinom(1, size = x, prob = 1/p2pfu)
    df$n_res_faulty <- unlist(lapply(df$tot_r_b4_faulty, f_faulty))
    df$n_sus_faulty <- unlist(lapply(df$tot_s_b4_faulty, f_faulty))
    
    ## Subtracting out faulty viruses
    df$tot_r <- df$tot_r_b4_faulty - df$n_res_faulty
    df$tot_s <- df$tot_s_b4_faulty - df$n_sus_faulty
    
    ## Make sure there are still viruses in the system
    if(sum(df$tot_r + df$tot_s)==0) break
    
    ## Function to generate probability vector for subunit composition
    f_prob_gen_r_subs_4_r <- function(x, y) list(dbinom(0:subunits, subunits, (x/y)))
    
    ## The probability of generating a capsid with n resistant subunits
    df$prob_gen_res_subs_r <- mapply(f_prob_gen_r_subs_4_r, df$res_per_cell, df$viruses_per_cell)
    df$prob_gen_res_subs_s <- mapply(f_prob_gen_r_subs_4_r, df$res_per_cell, df$viruses_per_cell)
    
    ## Function to pull viruses into subunit categories
    f_gen_v_subs_large <- function(x, y, max_chunk) {
      ## Number of chunks to run
      chunks <- x / max_chunk
      
      ## If the number of chunks is less than 1, just run as normal
      if(chunks <= 1) {
        ## Multinomial pull from these probs
        tot_prog <- rmultinom(1, x, y)
      } 
      ## If there are more pulls than the computer can handle, run in batches
      if(chunks > 1) {
        ## floor() rounds down to the nearest integer
        new_temp <- rmultinom(n = floor(chunks), size = max_chunk, prob = y)
        new_viruses <- rowSums(new_temp)
        ## Now add the remainder to the new_viruses vector
        remainder <- chunks - floor(chunks)
        ## Rounding just in case
        new_remand <- rmultinom(n = 1, size = round(max_chunk*remainder, 0), prob = y)
        tot_prog <- new_remand+new_viruses
      }
      return(list(tot_prog))
    }
    
    df$tot_gen_res_subs_r <- mapply(f_gen_v_subs_large, 
                                    x = df$tot_r, 
                                    y = df$prob_gen_res_subs_r,
                                    max_chunk = 2147483647)
    
    df$tot_gen_res_subs_s <- mapply(f_gen_v_subs_large,
                                    x = df$tot_s, 
                                    y = df$prob_gen_res_subs_r,
                                    max_chunk = 2147483647)
    
    ## Turn them into doubles so that R can add them together
    df$tot_gen_res_subs_r <- lapply(df$tot_gen_res_subs_r, as.double)
    
    df$tot_gen_res_subs_s <- lapply(df$tot_gen_res_subs_s, as.double)
    
    ## Add all of the corresponding indices together 
    res_vec <- Reduce(function(x,y) x + y, df$tot_gen_res_subs_r)
    
    sus_vec <- Reduce(function(x,y) x + y, df$tot_gen_res_subs_s)
    
    #### Drug action ####
    ## fit_func_in has the survival probabilities associated with capsid subunit
    
    ## If the t_pocap variable is greater than the current generation, reassign
    ## fit_func_in to be a vector of 1s rather than the normal fitness vector
    if(t_pocap > t) {
      curr_fit_func <- rep(1, subunits+1)
    }
    
    ## Now treating these as the probability of survival pairing these values
    ## up with their associated subunits
    res_post_pocap <- mapply(rbinom, n = 1, size = res_vec, prob = curr_fit_func)
    sus_post_pocap <- mapply(rbinom, n = 1, size = sus_vec, prob = curr_fit_func)
    
    
    ## Then, assuming that the mutation has no effect on the binding
    ## affinity or the ability to infect the next cell, we can assume
    ## that the sum of these arrays, is all of the viruses that are able
    ## to escape.
    
    new_mut_v <- sum(res_post_pocap)
    new_wt_v <- sum(sus_post_pocap)
    
    #### Immune system ####
    
    ## imm_delay allows for the modulation of the immune system delay. If the
    ## generation is greater than the immune system delay, then there is an
    ## additional binomial sampling step which escalates over time and is
    ## proportional to the number of infections that occurred in the previous
    ## generation.
    if(t > imm_delay){
      imm_delay <- as.numeric(imm_delay)
      new_mut_v <- rbinom(1, new_mut_v, prob = exp(-(t - imm_delay) * immune_sensitivity))
      new_wt_v <- rbinom(1, new_wt_v, prob = exp(-(t - imm_delay) * immune_sensitivity))
    }
    
    ## Recording the number of cells that did not get infected
    cell_surv <- as.vector(df_0s$inf_w_n_v)
    
    #### Generating dfs ####
    ## Generating a long dataframe to hold all of our information
    ## if it is the first run (cc==1) and the time == the first time
    ## input, generate a new dataframe
    if(t==1) {
      df_mut <- data.frame(type = "resistant",          ## resistant virus
                           time = 0,                    ## time
                           moi_type = c(moi_mut_start),
                           moi_res = c(moi_mut_start),
                           moi_wt = c(moi_wt_start),
                           init_wt = moi_wt_start,
                           init_mut = moi_mut_start,
                           c_pop = c(c_pop),
                           v_prog = c(v_prog),
                           p2pfu = c(p2pfu),
                           wt_rep_abil = fit_func_in[1],
                           time_to_pocap = t_pocap,
                           bg_mut_rate = bg_mutat,
                           uninf_cells = c(cell_surv),
                           surv_pfu = moi_mut_start*c_pop,
                           tot_pfu = (moi_wt_start*c_pop+moi_mut_start*c_pop),
                           pop_prop = moi_mut_start/(moi_mut_start+moi_wt_start),
                           imm_delay = imm_delay,
                           imm_m = imm_m,
                           imm_sd = imm_sd,
                           id = id,
                           seed = seed_in)
      
      df_wt <- data.frame(type = "susceptible",
                          time = 0,
                          moi_type = c(moi_wt_start),
                          moi_res = c(moi_mut_start),
                          moi_wt = c(moi_wt_start),
                          init_wt = moi_wt_start,
                          init_mut = moi_mut_start,
                          c_pop = c(c_pop),
                          time_to_pocap = t_pocap,
                          v_prog = c(v_prog),
                          p2pfu = c(p2pfu),
                          wt_rep_abil = fit_func_in[1],
                          bg_mut_rate = bg_mutat,
                          uninf_cells = c(cell_surv),
                          surv_pfu = moi_wt_start*c_pop,
                          tot_pfu = (moi_wt_start*c_pop+moi_mut_start*c_pop),
                          pop_prop = moi_wt_start/(moi_mut_start+moi_wt_start),
                          imm_delay = imm_delay,
                          imm_m = imm_m,
                          imm_sd = imm_sd,
                          id = id,
                          seed = seed_in)
      
      df_long <- rbind(df_mut, df_wt)
    }
    ## Create a new row to the current dataframe
    new_row_mut <- c("resistant",
                     t,
                     new_mut_v/c_pop,
                     new_mut_v/c_pop,
                     new_wt_v/c_pop,
                     moi_wt_start,
                     moi_mut_start,
                     c_pop,
                     v_prog,
                     p2pfu,
                     fit_func_in[1],
                     t_pocap,
                     bg_mutat,
                     cell_surv,
                     new_mut_v,
                     new_mut_v+new_wt_v,
                     new_mut_v/(new_wt_v+new_mut_v),
                     imm_delay,
                     imm_m,
                     imm_sd,
                     id,
                     seed_in)
    
    new_row_wt <- c("susceptible",
                    t,
                    new_wt_v/c_pop,
                    new_mut_v/c_pop,
                    new_wt_v/c_pop,
                    moi_wt_start,
                    moi_mut_start,
                    c_pop,
                    v_prog,
                    p2pfu,
                    fit_func_in[1],
                    t_pocap,
                    bg_mutat,
                    cell_surv,
                    new_wt_v,
                    new_mut_v+new_wt_v,
                    new_wt_v/(new_wt_v+new_mut_v),
                    imm_delay,
                    imm_m,
                    imm_sd,
                    id,
                    seed_in)
    
    ##bind to new dataframe
    df_long <- rbind(df_long, new_row_mut, new_row_wt)
    
    ## Update progress bar
    setTxtProgressBar(pb,cc)
    
    ## Increase the condition counter by 1 since condition has been finished
    cc <- cc+1
  }
  
  close(pb)
  
  ## Fixing variables ##
  numericsl <- c("moi_type","moi_res","moi_wt","c_pop",
                 "v_prog","p2pfu","uninf_cells","time", 
                 "bg_mut_rate","init_mut","init_wt","surv_pfu","wt_rep_abil",
                 "imm_sd", "imm_m", "imm_delay")
  df_long[numericsl] <- lapply(df_long[numericsl], as.numeric)
  return(df_long)
}

#### Main plotting function ####
## df = An output simulation dataframe from one of the previous functions
## stoch = Only use if you want to plot multiple stochastic functions per plot
plot_traj_polv <- function(df, stoch = FALSE, plot_title = TRUE, title_text) {
  
  if(stoch == FALSE) {
    s <- ggplot(df, aes(x = time, y = surv_pfu, color = type)) +
      geom_point(size = 6) + 
      geom_line() +
      scale_color_manual(values = c("#c94d4d", "#688cc8")) +
      theme_light() + 
      xlab("Number of Replication Cycles") + ylab("Viral yield (PFU/mg)") + 
      if(plot_title) {ggtitle(title)} +
      scale_y_log10(breaks=10^(0:10),
                         labels = sapply(0:10,function(i){parse(text = sprintf("10^%d",i))})) +
      theme(text = element_text(size=18), axis.text = element_text(size=18),
            panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
            legend.title = element_blank(),
            plot.title = element_text(size=20)) 
  }
  
  if(stoch == TRUE) {
    s <- ggplot(df, aes(x = time, y = surv_pfu, color = type, 
                            group = interaction(type, id))) +
        geom_point(size = 6, alpha = 0.5) + 
        geom_line() +
        scale_color_manual(values = c("#c94d4d", "#688cc8")) +
        theme_light() + 
        xlab("Number of Replication Cycles") + ylab("Viral yield (PFU/mg)") + 
        if(plot_title) {ggtitle(title)} +
        scale_y_continuous(trans="log10",breaks=10^(0:10),
                           labels = sapply(0:10,function(i){parse(text = sprintf("10^%d",i))})) +
        theme(text = element_text(size=18), axis.text = element_text(size=18),
              panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
              legend.title = element_blank(),
              plot.title = element_text(size=20)) 
  }
  
  return(s)
}

