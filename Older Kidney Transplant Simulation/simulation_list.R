`" 
*-------------------------------------------------------------------------------
* Simulation of Waitlist Outcomes
* simulation_list.R
* Date Created: 07/02/2021
* Last Updated: 07/28/2021
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
data = as.matrix(read.csv(file = paste0(in_path, "simulation_list.csv"), sep=","))

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
  
  event           <<- boot_data[,16]
  months_to_event <<- boot_data[,17]
  event_cd        <<- boot_data[,18]
  
  boot_data <<- boot_data[,-c(16:18)]
  
  data_list <<- boot_data
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
                           data_list = boot_data) {
  ### Function Inputs
  # seed: used to replicate results
  # n: sample size used for simulation
  # t: number of months to run simulation over
  ### Create bootstrapped data set
  set.seed(seed)
  num_patients <- n
  months <- t
  
  build.data(n=num_patients)
  
  ### Time to any event (AFT Weibull)
  event_coef <- cbind(c(-0.0248207, 0.0164917, 0.0496562, 0.019184,
                        0.2171164, 0.0697977, -0.0145885, 0.1982676,
                        -0.2521326, 0.1572624, 0.1823128, -0.0344955, 
                        -0.0363011, -0.2121378, 0.0265579, 3.204003))
  event_p <- 1.032089
  event_lambda <- exp(-event_p * data_list %*% event_coef)
  
  ### Which event? (Multinomial Logit)
  # Base Outcome: Other Removal
  ddtx_beta <- cbind(c(-0.0690856, -0.0275894, 0.0575006, 0.0042321,
                       0.065416, -0.0229673, 0.1318429, -0.0589578, 
                       0.6172187, -0.4337746, -0.4624144, -0.4878962, 
                       -0.3108714, 0.4528148, -0.0772435, 0.2973837))
  ldtx_beta <- cbind(c(-0.0529223, 0.1351953, -0.7818026, -0.2728606,
                       -0.5658424, -0.4933644, 0.1441382, -0.9363874,
                       0.0917832, -0.1546752, -0.2690706, -0.5287634, 
                       -0.3341261, 0.6440933, -0.1102809, 0.2476341))
  mort_beta <- cbind(c(-0.0699549, 0.0898672, -0.1353044, -0.1416642,
                       -0.0220675, -0.0118303, -0.0248519, -0.2473104,
                       0.1463788, 0.0534331, 0.0805427, 0.2356591, 
                       0.1181792, 0.1204255, -0.0250232, -0.6208005))
  ddtx_exb <- exp(data_list %*% ddtx_beta)  
  ldtx_exb <- exp(data_list %*% ldtx_beta)
  mort_exb <- exp(data_list %*% mort_beta)
  
  ddtx_prob   <- ddtx_exb / (1 + ddtx_exb + ldtx_exb + mort_exb)
  ldtx_prob   <- ldtx_exb / (1 + ddtx_exb + ldtx_exb + mort_exb)
  mort_prob   <- mort_exb / (1 + ddtx_exb + ldtx_exb + mort_exb)
  remove_prob <- 1 / (1 + ddtx_exb + ldtx_exb + mort_exb)
  
  ddtx_range  <- ddtx_prob
  ldtx_range  <- ddtx_range + ldtx_prob
  mort_range  <- ldtx_range + mort_prob
  
  ### Initialize vectors and matrices
  event_results   <- rep(0, num_patients)
  event_trace     <- matrix(nrow = num_patients, ncol=months)
  
  ddtx_results    <- rep(0, num_patients)
  ddtx_trace      <- matrix(nrow = num_patients, ncol=months)
  
  ldtx_results    <- rep(0, num_patients)
  ldtx_trace      <- matrix(nrow = num_patients, ncol=months)
  
  mort_results    <- rep(0, num_patients)
  mort_trace      <- matrix(nrow = num_patients, ncol=months)
  
  remove_results  <- rep(0, num_patients)
  remove_trace    <- matrix(nrow = num_patients, ncol=months)
  
  m2event         <- rep(0, num_patients)
  
  ### Run Simulation
  sim_time <- proc.time()  # start timer
  i <- 0
  while (all(event == 1) == FALSE) {
    i = i + 1
    ### Any event?
    event_haz <- event_p * event_lambda * i^(event_p - 1)  # Any event hazard vector
    event_prob <- 1 - exp(-event_haz)  # Any event probability vector
    rand_event_surv <- runif(num_patients, min = 0, max = 1)  # Vector of random probabilities from uniform
    y <- as.integer(rand_event_surv < event_prob)  # vector of those who had event in period i
    new_event <- as.integer(y==1 & pmax(y, event_results) > event_results)
    event_results <- pmax(y, event_results)  # create vector that combines those with previous events and new events
    event_trace[,i] <- event_results
    m2event <- m2event + new_event*i
    
    event <- event_results
    ### Which Event
    rand_num_v  <- runif(num_patients, min = 0, max = 1)  # Vector of random probabilities from uniform
    
    a <- as.integer(rand_num_v <= ddtx_range)
    b <- as.integer((rand_num_v > ddtx_range) & (rand_num_v <= ldtx_range))
    c <- as.integer((rand_num_v > ldtx_range) & (rand_num_v <= mort_range))
    d <- as.integer(rand_num_v > mort_range)
    
    new_ddtx   <- as.integer(a==1 & 
                               pmax(a, ddtx_results) > ddtx_results &
                               new_event==1)
    new_ldtx   <- as.integer(b==1 & 
                               pmax(b, ldtx_results) > ldtx_results &
                               new_event==1)
    new_mort   <- as.integer(c==1 & 
                               pmax(c, mort_results) > mort_results &
                               new_event==1)
    new_remove <- as.integer(d==1 & 
                               pmax(d, remove_results) > remove_results &
                               new_event==1)
    
    ddtx_results   <- pmax(new_ddtx, ddtx_results)
    ldtx_results   <- pmax(new_ldtx, ldtx_results)
    mort_results   <- pmax(new_mort, mort_results)
    remove_results <- pmax(new_remove, remove_results)
    
    ddtx_trace[,i]   <- ddtx_results
    ldtx_trace[,i]   <- ldtx_results
    mort_trace[,i]   <- mort_results
    remove_trace[,i] <- remove_results
  }
  
  ### Build dataset for analysis
  ddtx_trace   <- ddtx_trace[, 1:i]
  ldtx_trace   <- ldtx_trace[, 1:i]
  mort_trace   <- mort_trace[, 1:i]
  remove_trace <- remove_trace[, 1:i]
  event_trace  <- event_trace[, 1:i]
  
  results_sim <<- as.data.frame(cbind(data_list, event_trace[,i],
                                     ddtx_trace[,i], ldtx_trace[,i],
                                     mort_trace[,i], remove_trace[,i],
                                     m2event))
  results_sim <- results_sim %>%
    rename(
      event = V17,
      ddtx = V18,
      ldtx = V19,
      mort = V20,
      remove = V21,
      months_to_event = m2event
    ) %>%
    mutate(months_to_event = ifelse(event==0, i, months_to_event)) %>%
    mutate(race = ifelse(black==0 & hispanic==0 & other_race==0, 0,
                         ifelse(black==1 & hispanic==0 & other_race==0, 1,
                                ifelse(black==0 & hispanic==1 & other_race==0, 2, 3)))) %>%
    mutate(sim = rep(1, num_patients)) %>%
    mutate(race = factor(race, levels = c(0, 1, 2, 3),
                         labels = c("White", "Black",
                                    "Hispanic", "Other"))) %>%
    mutate(sim = factor(sim, levels = c(0, 1),
                        labels = c("Actual", "Simulated")))


  results_real <- as.data.frame(cbind(data_list, event, months_to_event,
                                      event_cd))
  results_real <- results_real %>%
    mutate(race = ifelse(black==0 & hispanic==0 & other_race==0, 0,
                         ifelse(black==1 & hispanic==0 & other_race==0, 1,
                                ifelse(black==0 & hispanic==1 & other_race==0, 2, 3)))) %>%
    mutate(ddtx = ifelse(event_cd==0, 1, 0)) %>%
    mutate(ldtx = ifelse(event_cd==1, 1, 0)) %>%
    mutate(mort = ifelse(event_cd==2, 1, 0)) %>%
    mutate(remove = ifelse(event_cd==3, 1, 0)) %>%
    mutate(sim = rep(0, num_patients)) %>%
    mutate(race = factor(race, levels = c(0, 1, 2, 3),
                         labels = c("White", "Black",
                                    "Hispanic", "Other"))) %>%
    mutate(sim = factor(sim, levels = c(0, 1),
                        labels = c("Actual", "Simulated"))) %>%
    select(!event_cd)

  df_test_long <<- rbind(results_sim, results_real)

  sim_df <<- results_sim
}

run.simulation2 <- function(seed=11111,
                            n, 
                            t, 
                            data_list = boot_data) {
  ### Function Inputs
  # seed: used to replicate results
  # n: sample size used for simulation
  # t: number of months to run simulation over
  ### Create bootstrapped data set
  set.seed(seed)
  num_patients <- n
  months <- t
  
  build.data(n=num_patients)
  
  ### Deceased Donor Tx (AFT Weibull)
  ddtx_coef   <- cbind(c(0.0050035, 0.0906279, -0.1136044, -0.0524093,
                         0.0958358, 0.0208221, -0.0831058, 0.1089517,
                         -0.6740915, 0.5164808, 0.5505225, 0.3014152,
                         0.1673722, -0.4502303, 0.0888333, 4.497449))
  ddtx_p      <- 0.8057152
  ddtx_lambda <- exp(-ddtx_p * data_list %*% ddtx_coef)
  
  ### Living Donor Tx (PH Gompertz)
  ldtx_coef   <- cbind(c(0.0066903, 0.0868512, -0.7872054, -0.1730859,
                         -0.6899881, -0.4436386, 0.0113922, -0.9649756, 
                         -0.0205408, -0.0856782, -0.1890228, -0.247052,
                         -0.2171749, 0.4456786, -0.0662714, -3.37321))
  ldtx_gamma  <- -0.0558557
  ldtx_lambda <- exp(data_list %*% ldtx_coef)
  
  ### Wait list mortality (AFT Log-logistic)
  mort_coef   <- cbind(c(0.0029614, -0.0433458, 0.0557588, 0.0691281, 
                         0.1577234, 0.0602654, 0.0454398, 0.2553751, 
                         -0.0968927, -0.0228495, -0.0404601, -0.3906527,
                         -0.2343233, -0.0794903, 0.0107482, 4.586921))
  mort_gamma  <- 0.7179732
  mort_lambda <- exp(-data_list %*% mort_coef)  
  
  ### Other wait list removal (AFT Weibull)
  remove_coef   <- cbind(c(-0.04617, 0.0149515, -0.0269867, -0.0258846,
                           0.1306117, 0.0357642, -0.0082171, 0.0360195, 
                           -0.0060275, 0.0330261, 0.0396322, -0.1706369, 
                           -0.0911433, -0.0243874, -0.015688, 4.286393))
  remove_p      <- 1.685768
  remove_lambda <- exp(-remove_p * data_list %*% remove_coef)
  
  ### Initialize vectors and matrices
  event_results   <- rep(0, num_patients)
  event_trace     <- matrix(nrow = num_patients, ncol=months)
  
  ddtx_results    <- rep(0, num_patients)
  ddtx_trace      <- matrix(nrow = num_patients, ncol=months)
  
  ldtx_results    <- rep(0, num_patients)
  ldtx_trace      <- matrix(nrow = num_patients, ncol=months)
  
  mort_results    <- rep(0, num_patients)
  mort_trace      <- matrix(nrow = num_patients, ncol=months)
  
  remove_results  <- rep(0, num_patients)
  remove_trace    <- matrix(nrow = num_patients, ncol=months)
  
  m2event         <- rep(0, num_patients)
  events          <- rep(0, num_patients)
  
  ### Run Simulation
  sim_time <- proc.time()  # start timer
  for (i in 1:months) {
    ### Deceased Donor Tx?
    # Generate vector of ddtx probabilities
    ddtx_prob <- 1 - exp(-(ddtx_p * ddtx_lambda * i^(ddtx_p - 1)))
    
    # Determine who has event by using random uniform draws and determine who
    # had a new event this period
    rand_ddtx_surv <- runif(num_patients, min = 0, max = 1)  # Vector of random probabilities from uniform
    y              <- as.integer(rand_ddtx_surv < ddtx_prob)  # vector of those who had event in period i
    new_ddtx       <- as.integer(y==1 & pmax(y, ddtx_results) > ddtx_results)
    
    # Check check that they are still at risk for event
    ddtx_results <- pmax(y, ddtx_results)
    match        <- as.integer(ddtx_results == pmax(ldtx_results, mort_results, 
                                                    remove_results) & 
                                 pmax(ldtx_results, mort_results, remove_results)==1)
    ddtx_results <- ddtx_results - match
    new_ddtx     <- new_ddtx - match
    
    # Create vector the tracks cumulative ddtx events, create patient trace, and 
    # record time to event
    ddtx_trace[,i] <- ddtx_results
    
    
    ### Living Donor Tx?
    # Generate vector of ldtx probabilities
    ldtx_prob <- 1 - exp(-(ldtx_lambda * exp(ldtx_gamma*i)))  # Graft loss probability vector
    
    # Determine who has event by using random uniform draws and determine who
    # had a new event this period
    rand_ldtx_surv <- runif(num_patients, min = 0, max = 1)  # random number vector
    x              <- as.integer(rand_ldtx_surv < ldtx_prob)  # vector of possible events
    new_ldtx       <- as.integer(x==1 & pmax(x, ldtx_results) > ldtx_results)
    
    # Check check that they are still at risk for event
    ldtx_results <- pmax(x, ldtx_results)
    match        <- as.integer(ldtx_results == pmax(ddtx_results, mort_results, 
                                                    remove_results) & 
                                 pmax(ddtx_results, mort_results, remove_results)==1)
    ldtx_results <- ldtx_results - match
    new_ldtx     <- new_ldtx - match    
    
    # Create vector the tracks cumulative ldtx events, create patient trace
    ldtx_trace[,i] <- ldtx_results
    
    
    ### Wait list mortality?
    # Generate vector of mortality probabilities
    mort_survival <- (1 + (mort_lambda*i)^(1/mort_gamma))^(-1)    # S(t)
    mort_density  <- ((mort_lambda^(1/mort_gamma))*(i^(1/mort_gamma - 1))) /
      (mort_gamma * (1 + (mort_lambda*i)^(1/mort_gamma))^2) # f(t)
    mort_prob     <- 1 - exp(-(mort_density / mort_survival))
    
    # Determine who has event by using random uniform draws and determine who
    # had a new event this period
    rand_mort_surv <- runif(num_patients, min = 0, max = 1)  # random number vector
    x              <- as.integer(rand_mort_surv < mort_prob)  # vector of possible events
    new_mort       <- as.integer(x==1 & pmax(x, mort_results) > mort_results)
    
    # Check check that they are still at risk for event
    mort_results <- pmax(x, mort_results)
    match        <- as.integer(mort_results == pmax(ddtx_results, ldtx_results, 
                                                    remove_results) & 
                                 pmax(ddtx_results, ldtx_results, remove_results)==1)
    mort_results <- mort_results - match
    new_mort     <- new_mort - match    
    
    # Create vector the tracks cumulative ddtx events, create patient trace
    mort_trace[,i] <- mort_results    
    
    
    ### Other List Removal?
    # Generate vector of list removal probabilities
    remove_prob <- 1 - exp(-(remove_p * remove_lambda * i^(remove_p - 1)))
    
    # Determine who has event by using random uniform draws and determine who
    # had a new event this period
    rand_remove_surv <- runif(num_patients, min = 0, max = 1)  # Vector of random probabilities from uniform
    y                <- as.integer(rand_remove_surv < remove_prob)  # vector of those who had event in period i
    new_remove       <- as.integer(y==1 & pmax(y, remove_results) > remove_results)
    
    # Check check that they are still at risk for event
    remove_results <- pmax(y, remove_results)
    match          <- as.integer(remove_results == pmax(ddtx_results, 
                                                        ldtx_results, mort_results) 
                                 & pmax(ddtx_results, ldtx_results, mort_results)==1)
    remove_results <- remove_results - match
    new_remove     <- new_remove - match 
    
    # Create vector the tracks cumulative removal events, create patient trace
    remove_trace[,i] <- remove_results
    
    ### Time to any event
    event_trace[,i] <- pmax(ddtx_results, ldtx_results, 
                            mort_results, remove_results)
    m2event         <- m2event + pmax(new_ddtx, new_ldtx, 
                                      new_mort, new_remove)*i 
  }
  
  ### Build dataset for analysis
  results_sim <- as.data.frame(cbind(data_list, event_trace[,months],
                                     ddtx_trace[,months], ldtx_trace[,months],
                                     mort_trace[,months], remove_trace[,months],
                                     m2event))
  results_sim <- results_sim %>%
    rename(
      event = V17,
      ddtx = V18,
      ldtx = V19,
      mort = V20,
      remove = V21,
      months_to_event = m2event
    ) %>%
    mutate(months_to_event = ifelse(event==0, months, months_to_event)) %>%
    mutate(race = ifelse(black==0 & hispanic==0 & other_race==0, 0, 
                         ifelse(black==1 & hispanic==0 & other_race==0, 1, 
                                ifelse(black==0 & hispanic==1 & other_race==0, 2, 3)))) %>%
    mutate(sim = rep(1, num_patients)) %>%
    mutate(race = factor(race, levels = c(0, 1, 2, 3),
                         labels = c("White", "Black", 
                                    "Hispanic", "Other"))) %>%
    mutate(sim = factor(sim, levels = c(0, 1), 
                        labels = c("Actual", "Simulated")))
  
  
  results_real <- as.data.frame(cbind(data_list, event, months_to_event, 
                                      event_cd))
  results_real <- results_real %>%
    mutate(race = ifelse(black==0 & hispanic==0 & other_race==0, 0, 
                         ifelse(black==1 & hispanic==0 & other_race==0, 1, 
                                ifelse(black==0 & hispanic==1 & other_race==0, 2, 3)))) %>%
    mutate(ddtx = ifelse(event_cd==0, 1, 0)) %>%
    mutate(ldtx = ifelse(event_cd==1, 1, 0)) %>%
    mutate(mort = ifelse(event_cd==2, 1, 0)) %>%
    mutate(remove = ifelse(event_cd==3, 1, 0)) %>%
    mutate(sim = rep(0, num_patients)) %>%
    mutate(race = factor(race, levels = c(0, 1, 2, 3),
                         labels = c("White", "Black", 
                                    "Hispanic", "Other"))) %>%
    mutate(sim = factor(sim, levels = c(0, 1), 
                        labels = c("Actual", "Simulated"))) %>%
    select(!event_cd)
  
  df_test_long <<- rbind(results_sim, results_real)
  
  sim_df <<- results_sim
}

##### Run Simulation #####
run.simulation(n=100000, t=600)
run.simulation2(n=100000, t=360)
