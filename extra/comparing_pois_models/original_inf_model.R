stoch_polv_og <- function(n = 1,
                       moi_wt_start = 1,
                       moi_mut_start = 1,
                       t_pocap = 0,
                       imm_delay = 9,
                       imm_m = -1.6,
                       imm_sd = 0.5,
                       c_pop = 2.9 * 10 ^ 6,
                       v_prog = 203,
                       p2pfu = 1,
                       subunits = 60,
                       bg_mutat = 1/(2*10^-5),
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
    # browser()
    
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
    vpc <- qpois(0.999999999999, (moi_tot + 1))
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
    
    # 
    # #### Initial data frames ####
    # df <- data.frame(viruses_per_cell = 0:vpc)
    # 
    # ## Probabiilty of being infected with this number of genomes (after below 
    # ## update, now just used for tracking)
    # df$prob_inf_w_n_v <- dpois(df$viruses_per_cell, moi_tot)
    # 
    # ## Create a vector for the actual number of cells in each category to add to
    # df$inf_w_n_v <-  rep(0, (vpc+1))
    # 
    # 
    # ## Initial infections probabilities ##
    # ## Updated after the first round of reviewer comments (7/31/2025)
    # ## "Fill up" the cells, one at a time, with viruses, based on a pois. pull.
    # ## Take viruses and cells out of the population as they are pulled.
    # ## Here, we assume indipendence between pulls
    # 
    # ## Create temporary variables to iterate over
    # temp_v_tot <- v_pop_tot
    # temp_c_pop <- c_pop
    # temp_moi_tot <- moi_tot
    # # v_tot_tracker <- NULL
    # # moi_tracker <- NULL
    # 
    # ## Check that we have viruses
    # 
    # if(v_pop_tot <= 0) break ## Check that we have viruses
    # if(moi_tot <= 0) break
    # 
    # ## While there are cells available:
    # while(temp_c_pop > 0) {
    #   ## Have to break the loop once the MOI is 0 or less (temp_moi_tot <= 0)
    #   if(temp_moi_tot <= 0) {
    #     ## Make sure to add back any cells that were not infected
    #     df$inf_w_n_v[1] <- df$inf_w_n_v[1] + temp_c_pop
    #     break
    #   }
    #   # moi_tracker <- c(moi_tracker, temp_moi_tot)
    #   # v_tot_tracker <- c(v_tot_tracker, temp_v_tot)
    #   ## Create a variable for the new viruses
    #   new_v <- NULL
    #   
    #   ## If this is the last cell, assign it all of the remaining viruses
    #   if (temp_c_pop == 1) {
    #     new_v <- temp_v_tot
    #   }
    #   
    #   ## If there are more than one cells left, pull the number from pois
    #   if (temp_c_pop > 1) {
    #     new_v <- rpois(1, temp_moi_tot)
    #     ## If the number pulled is larger than the total number of viruses, just
    #     ## assign the total number of viruses left.
    #     if (new_v > temp_v_tot) {
    #       new_v <- temp_v_tot
    #     }
    #   }
    #   
    #   ## Track the assigned cell (add 1 to account for 0 viruses)
    #   df$inf_w_n_v[new_v+1] <- df$inf_w_n_v[new_v+1] + 1 
    #   ## Subtract the new_v from total virus count
    #   temp_v_tot <- temp_v_tot - new_v 
    #   ## Subtract from total cell population count
    #   temp_c_pop <-  temp_c_pop - 1 
    #   ## Calculate new MOI
    #   temp_moi_tot <- temp_v_tot/temp_c_pop
    #   # cc <- cc+1
    # }
    # 
    # ## We can now remove the 0s from our dataframe, since they have no viruses
    # df_0s <- df[1,] ## Saving for later, just in case desired.
    # 
    # df <- df[-1,]
    # 
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