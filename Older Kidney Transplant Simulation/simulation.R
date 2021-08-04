`" 
*-------------------------------------------------------------------------------
* Older Kidney Transplant Simulation Model (OKTSM)
* Simulation of Older Renal Transplantation (SORT)
* simulation.R
* Date Created: 04/19/2021
* Last Updated: 05/04/2021
* Matt Kaufmann, Stanford University
*-------------------------------------------------------------------------------
"`

##### Load Packages #####
library(haven)
library(dplyr)
library(survival)
library(survminer)

##### Load data #####
in_path <- "/Volumes/Nephrology/Projects/UNOS/KaufmannM/data/"
out_path <- "/Volumes/Nephrology/Projects/UNOS/KaufmannM/results/"
data = as.matrix(read.csv(file = paste0(in_path, "simulation.csv"), sep=","))

##### Functions for Simulation #####
build.data <- function(seed=11111, n) {
  ### Function Inputs
  # seed: used to replicate results
  # n: sample size used for simulation
  ### Create bootstrapped data set
  set.seed(seed)
  num_patients <- n
  
  boot_data <<- data[sample(nrow(data), num_patients, replace=T),]
  sim_data <- boot_data[,-1]
  rownames(sim_data) <- boot_data[,1]
  
  graft_loss <<- boot_data[,20]
  months_to_gl <<- boot_data[,21]
  gs_death <<- boot_data[,22]
  months_to_gs_mort <<- boot_data[,23]
  gl_death <<- boot_data[,24]
  months_to_gl_mort <<- boot_data[,25]
  
  boot_data <<- boot_data[,-c(20:25)]
  
  data_func_graft <<- boot_data
  data_graft_loss <<- cbind(boot_data, rep(0, num_patients))
  obs_df <<- as.data.frame(data)
  obs_df <<- obs_df %>%
    mutate(race = ifelse(black==0 & hispanic==0 & other_race==0, 0, 
                         ifelse(black==1 & hispanic==0 & other_race==0, 1, 
                                ifelse(black==0 & hispanic==1 & other_race==0, 2, 3)))) %>%
    mutate(race = factor(race, levels = c(0, 1, 2, 3),
                         labels = c("White", "Black", 
                                    "Hispanic", "Other")))
}

run.simulation <- function(seed=11111,
                           n, 
                           t, 
                           data_gs = data_func_graft,
                           data_gl = data_graft_loss) {
  ### Function Inputs
  # seed: used to replicate results
  # n: sample size used for simulation
  # t: number of months to run simulation over
  ### Create bootstrapped data set
  set.seed(seed)
  num_patients <- n
  months <- t
  
  build.data(n=num_patients)
  ### Graft loss (AFT Weibull)
  gl_coef <- cbind(c(-0.0088464, -0.1638935, -0.1649998, 0.1786882, 0.3689721, 
                     -0.0484216, -0.0790271, -0.0205087, 0.2752222, 0.1820622,
                     0.0042934, -0.0056192, 0.0864764, -0.1689237, -0.4211846,
                     -0.2048911, -0.121687, -0.7027984, 8.441661))
  gl_p <- 0.8409265
  gl_lambda <- exp(-gl_p * data_gs %*% gl_coef)
  
  ### Death w/ functioning graft (PH Gompertz)
  gs_mort_coef <- cbind(c(0.0527071, 0.1176476, -0.2356672, -0.4074649, 
                          -0.4080931, 0.0513649, 0.169329, 0.0050798, 
                          -0.0136294, -0.1041224, -0.0490057, 0.0063064, 
                          0.0548066, 0.4064592, 0.3953489, 0.2723943, 0.0319272,
                          0.222242, -7.238323))
  gs_mort_gamma <- 0.0174279
  gs_mort_lambda <- exp(data_gs %*% gs_mort_coef)
  
  ### Death after graft loss (put months2gl variable last)
  dial_mort_coef <- cbind(c(-0.016629, 0.2210384, 0.8659112, 0.2514827, 
                            0.2345926, -0.0204829, 0.8456908, 0.0053321, 
                            0.143826, -0.5176801, 0.0307328, 0.0192943, 
                            -0.0075228, -0.2775279, 0.2060207, -0.3552075, 
                            -0.1486936, -0.0311216, 1.660609, 0.005687))
  dial_mort_p <- 0.7257093
  
  ### Same day graft loss - mortality 
  gl_mort_coef <- cbind(c(0.0021596, -0.1084306, 0.0030603, -0.1479822, 
                          -0.3039797, -0.0390912, 0.1028756, -0.0059598,
                          -0.1124294, 0.2996683, 0.0733479, -0.0032874, 
                          0.0414696, -0.016484, -0.3080031, 0.1165943, 
                          -0.1381915, 0.3061243, -0.5784331, 0.0098482))

  ### Initialize vectors and matrices
  gs_mort_results <- rep(0, num_patients)
  gs_mort_trace <- matrix(nrow = num_patients, ncol=months)
  
  gl_results <- rep(0, num_patients)
  gl_trace <- matrix(nrow = num_patients, ncol=months)
  
  gl_mort_results <- rep(0, num_patients)
  gl_mort_trace <- matrix(nrow = num_patients, ncol=months)
  
  m2gl <- rep(0, num_patients)
  m2gs_mort <- rep(0, num_patients)
  m2gl_mort <- rep(0, num_patients)
  
  ### Run Simulation
  sim_time <- proc.time()  # start timer
  
  ### Start first 30 days
  
  for (i in 1:months) {
    ### First we will see who has died w/ functioning graft
    gs_mort_haz <- gs_mort_lambda * exp(gs_mort_gamma*i) # vector of hazards
    gs_mort_prob <- 1 - exp(-gs_mort_haz)  # Graft loss probability vector
    rand_gs_mort_surv <- runif(num_patients, min = 0, max = 1)  # random number vector
    x <- as.integer(rand_gs_mort_surv < gs_mort_prob)  # vector of possible events
    new_gs_mort <- as.integer(x==1 & pmax(x, gs_mort_results) > gs_mort_results)
    gs_mort_results <- pmax(x, gs_mort_results)   # create vector that combines those with previous deaths to new deaths
    match <- as.integer(gs_mort_results == gl_results & gl_results==1)
    gs_mort_results <- gs_mort_results - match
    new_gs_mort <- new_gs_mort - match
    m2gs_mort <- m2gs_mort + new_gs_mort*i
    gs_mort_trace[,i] <- gs_mort_results  
    
    ### Next we see who died after graft failure
    data_gl[,21] <- m2gl # Update data set
    dial_mort_lambda <- exp(-dial_mort_p * data_gl %*% dial_mort_coef) # Update lambda vector
    dial_mort_haz <- dial_mort_p * dial_mort_lambda * i^(dial_mort_p - 1)  # Graft loss hazard vector
    dial_mort_prob <- 1 - exp(-dial_mort_haz)  # Graft loss probability vector
    rand_dial_surv <- runif(num_patients, min = 0, max = 1)  # Vector of random probabilities from uniform
    z <- as.integer(rand_dial_surv < dial_mort_prob)  # vector of those who had event in period i
    new_gl_mort <- as.integer(z==1 & pmax(z, gl_mort_results) > gl_mort_results)
    gl_mort_results <- pmax(z, gl_mort_results)  # create vector that combines those with previous graft loss to new graft loss
    gl_mort_results <- gl_results * gl_mort_results
    new_gl_mort <- gl_results * new_gl_mort
    m2gl_mort <- m2gl_mort + new_gl_mort*i 
    gl_mort_trace[,i] <- gl_mort_results
    
    ### Next we will see who experiencing graft loss
    gl_haz <- gl_p * gl_lambda * i^(gl_p - 1)  # Graft loss hazard vector
    gl_prob <- 1 - exp(-gl_haz)  # Graft loss probability vector
    rand_gl_surv <- runif(num_patients, min = 0, max = 1)  # Vector of random probabilities from uniform
    y <- as.integer(rand_gl_surv < gl_prob)  # vector of those who had event in period i
    new_gl <- as.integer(y==1 & pmax(y, gl_results) > gl_results)
    gl_results <- pmax(y, gl_results)  # create vector that combines those with previous graft loss to new graft loss
    match <- as.integer(gs_mort_results == gl_results & gs_mort_results==1)
    gl_results <- gl_results - match
    new_gl <- new_gl - match
    m2gl <- m2gl + new_gl*i
    gl_trace[,i] <- gl_results
    
    ### Day of graft loss death
    gl_mort_prob <- exp(data_gl %*% gl_mort_coef) / 
      (1 + exp(data_gl %*% gl_mort_coef))
    rand_gl_mort_surv <- runif(num_patients, min = 0, max = 1)  # Vector of random probabilities from uniform
    a <- as.integer(rand_gl_mort_surv < gl_mort_prob)  # vector of those who had event in period i
    new_gl_mort <- as.integer(a==1 & pmax((new_gl *a), gl_mort_results) > gl_mort_results)
    gl_mort_results <- pmax(new_gl * a, gl_mort_results)
    m2gl_mort <- m2gl_mort + new_gl_mort*i 
    gl_mort_trace[,i] <- gl_mort_results
  }
  
  ### Build dataset for analysis
  results <- as.data.frame(cbind(boot_data, gl_trace[,months], 
                                  gs_mort_trace[,months], gl_mort_trace[,months],
                                  m2gl, m2gs_mort, m2gl_mort))
  results <- results %>%
    rename(
      sim_gl = V21,
      sim_gs_mort = V22,
      sim_gl_mort = V23,
      sim_gl_time = m2gl,
      sim_gs_mort_time = m2gs_mort,
      sim_gl_mort_time = m2gl_mort
    ) %>%
    mutate(sim_gl_time = ifelse(sim_gl==0 & sim_gs_mort==0, 
                                months, sim_gl_time)) %>%
    mutate(sim_gl_time = ifelse(sim_gl==0 & sim_gs_mort==1, 
                                sim_gs_mort_time, sim_gl_time)) %>%
    mutate(sim_gs_mort_time = ifelse(sim_gs_mort==0 & sim_gl==0, 
                                     months, sim_gs_mort_time)) %>%
    mutate(sim_gs_mort_time = ifelse(sim_gs_mort==0 & sim_gl==1, 
                                     sim_gl_time, sim_gs_mort_time)) %>%  
    mutate(sim_gl_mort_time = ifelse(sim_gl_mort==0, months, sim_gl_mort_time))
  
  df_test <<- cbind(results, graft_loss, months_to_gl, gs_death, months_to_gs_mort)
  df_test <<- df_test %>%
    mutate(graft_loss = ifelse(graft_loss==1 & months_to_gl > months, 0, graft_loss)) %>%
    mutate(months_to_gl = ifelse(months_to_gl>months, months, months_to_gl)) %>%
    mutate(gs_death = ifelse(gs_death==1 & months_to_gs_mort > months, 0, gs_death)) %>%
    mutate(months_to_gs_mort = ifelse(months_to_gs_mort>months, months, months_to_gs_mort)) %>%
    mutate(race = ifelse(black==0 & hispanic==0 & other_race==0, 0, 
                         ifelse(black==1 & hispanic==0 & other_race==0, 1, 
                                ifelse(black==0 & hispanic==1 & other_race==0, 2, 3))))
  
  results_sim <- as.data.frame(cbind(boot_data, gl_trace[,months], 
                                 gs_mort_trace[,months], gl_mort_trace[,months],
                                 m2gl, m2gs_mort, m2gl_mort))
  results_sim <- results_sim %>%
    rename(
      graft_loss = V21,
      gs_death = V22,
      gl_death = V23,
      gl_time = m2gl,
      gs_death_time = m2gs_mort,
      gl_death_time = m2gl_mort
    ) %>%
    mutate(gl_time = ifelse(graft_loss==0 & gs_death==0, 
                                months, gl_time)) %>%
    mutate(gl_time = ifelse(graft_loss==0 & gs_death==1, 
                                gs_death_time, gl_time)) %>%
    mutate(gs_death_time = ifelse(gs_death==0 & graft_loss==0, 
                                     months, gs_death_time)) %>%
    mutate(gs_death_time = ifelse(gs_death==0 & graft_loss==1, 
                                  gl_time, gs_death_time)) %>%
    mutate(gl_death_time = gl_death_time - gl_time) %>%
    mutate(gl_death_time = ifelse(gl_death==0 & graft_loss==1,
                                  -gl_time+months, gl_death_time)) %>%
    mutate(gl_death_time = ifelse(gl_death_time < 0, NA, gl_death_time)) %>%
    mutate(race = ifelse(black==0 & hispanic==0 & other_race==0, 0, 
                         ifelse(black==1 & hispanic==0 & other_race==0, 1, 
                                ifelse(black==0 & hispanic==1 & other_race==0, 2, 3)))) %>%
    mutate(sim = rep(1, num_patients)) %>%
    mutate(race = factor(race, levels = c(0, 1, 2, 3),
                         labels = c("White", "Black", 
                                           "Hispanic", "Other"))) %>%
    mutate(sim = factor(sim, levels = c(0, 1), 
                        labels = c("Actual", "Simulated")))
  # results_sim$sim <- rep(1, num_patients)
  
  results_real <- as.data.frame(cbind(boot_data, graft_loss, gs_death, gl_death,
                                      months_to_gl, months_to_gs_mort,
                                      months_to_gl_mort))
  results_real <- results_real %>%
    rename(
      gl_time = months_to_gl,
      gs_death_time = months_to_gs_mort,
      gl_death_time = months_to_gl_mort
    ) %>%
    mutate(graft_loss = ifelse(graft_loss==1 & gl_time > months, 0, graft_loss)) %>%
    mutate(gl_time = ifelse(gl_time>months, months, gl_time)) %>%
    mutate(gs_death = ifelse(gs_death==1 & gs_death_time > months, 0, gs_death)) %>%
    mutate(gs_death_time = ifelse(gs_death_time>months, months, gs_death_time)) %>%
#    mutate(gl_death = ifelse(gl_death==1 & gl_death_time > (-gl_time+months), 
#                             0, gl_death)) %>%
#    mutate(gl_death_time = ifelse(gl_death_time > -gl_time+months, 
#                                  -gl_time+months, gl_death_time)) %>%
    mutate(race = ifelse(black==0 & hispanic==0 & other_race==0, 0, 
                         ifelse(black==1 & hispanic==0 & other_race==0, 1, 
                                ifelse(black==0 & hispanic==1 & other_race==0, 2, 3)))) %>%
    mutate(race = ifelse(black==0 & hispanic==0 & other_race==0, 0, 
                         ifelse(black==1 & hispanic==0 & other_race==0, 1, 
                                ifelse(black==0 & hispanic==1 & other_race==0, 2, 3)))) %>%
    mutate(sim = rep(0, num_patients)) %>%
    mutate(race = factor(race, levels = c(0, 1, 2, 3),
                         labels = c("White", "Black", 
                                    "Hispanic", "Other"))) %>%
    mutate(sim = factor(sim, levels = c(0, 1), 
                        labels = c("Actual", "Simulated")))
  # results_real$sim <- rep(0, num_patients)
  df_test_long <<- rbind(results_sim, results_real)
  
  sim_df <<- results_sim
}

##### Run Simulation #####
run.simulation(n=100000, t=120)

##### SCRAPS #####
'" 
### Wait list mortality
wait_mort_coef <- cbind(c())
wait_mort_p <- 
dial_mort_lambda <- exp(-wait_mort_p * boot_data %*% wait_mort_coef)

### First month after Tx outcomes
beta_gl <- cbind(c(0.0251736, -0.0105749, 0.261029, 0.0002929, -0.0145875, 
                   0.0718415, -0.0190684, 0.0129367, -0.0780347, 0.0753779,
                   -0.1140512, 0.0332946, -0.051162, 0.0901054, -0.430521,
                   0.2728332, 0.0041982, 1.067956, -0.0458052, -5.550331))
beta_dgf <- cbind(c(-0.0037742, 0.0017812, -0.0071346, 0.0269121, -0.0032198,
                    -0.0171439, -0.0008169, -0.0004284, -0.0090288, -0.020568,
                    0.0016661, -0.0054831, 0.0066592, -0.0131414, 0.0264853,
                    -0.1003867, -0.0058592, ))
                    
risk.table = TRUE,                       # Add risk table
tables.theme = theme_cleantable(),  # Clean risk table
palette = "jco"

(1-0.8984) - (sum(gl_trace[,60])/num_patients)
(1-0.7698) - (sum(gs_mort_trace[,60])/num_patients)
 
(1-0.9708) - (sum(gl_trace[,12])/num_patients)
(1-0.9539) - (sum(gs_mort_trace[,12])/num_patients)
"'