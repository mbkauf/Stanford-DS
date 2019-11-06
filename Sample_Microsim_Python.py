################################################################################
################# Microsimulation Model in Python ##############################
################################################################################
# This code creates a simple microsimulation using python
#
# Author: Matt Kaufmann
# Date Created: July 30th, 2019
# Last Updated: August 10th, 2019
# Version used: python3
#
############################### Install Packages ###############################
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

############################### Start Code #####################################
# Model inputs
n_i = 1000                            # number of simulated individuals
n_t = 10                               # time horizon, 30 cycles
v_n = ["H", "S1", "S2", "D"]           # the model states: Healthy (H), Sick (S1), Sicker (S2), Dead (D)
n_s = len(v_n)                         # the number of health states
v_M_1 = ["H"] * n_i                    # everyone begins in the healthy state
d_c = 0.03                             # discounting of costs is 3%
d_e = 0.03                             # discounting of QALYs is 3%
v_Trt = ["No Treatment", "Treatment"]  # store the strategy names
seed = 1                               # manually change the seed that is used

# Transition probabilities (per cycle)
p_HD = 0.005                 # probability to die when healthy
p_HS1 = 0.15                 # probability to become sick when healthy
p_S1H = 0.5                  # probability to become healthy when sick
p_S1S2 = 0.105               # probability to become sicker when sick
rr_S1 = 3                    # rate ratio of death when sick vs healthy
rr_S2 = 10                   # rate ratio of death when sicker vs healthy
r_HD = -np.log(1 - p_HD)     # rate of death when healthy
r_S1D = rr_S1 * r_HD         # rate of death when sick
r_S2D = rr_S2 * r_HD         # rate of death when sicker
p_S1D = 1 - np.exp(- r_S1D)  # probability to die when sick
p_S2D = 1 - np.exp(- r_S2D)  # probability to die when sicker

# Cost and utility inputs
c_H = 2000     # cost of remaining one cycle healthy
c_S1 = 4000    # cost of remaining one cycle sick
c_S2 = 15000   # cost of remaining one cycle sicker
c_Trt = 12000  # cost of treatment (per cycle)

u_H = 1        # utility when healthy
u_S1 = 0.75    # utility when sick
u_S2 = 0.5     # utility when sicker
u_Trt = 0.95   # utility when sick(er) and being treated

# PSA vars
num_iter = 100   # number of iterations for PSA

##################################### Functions ###########################################

# The MicroSim function for the simple microsimulation of the 'Sick-Sicker' model keeps track of what happens to each
# individual during each cycle.

def microsim(v_M_1, n_i, n_t, v_n, d_c, d_e, TR_out=True, TS_out=True, Trt=False, psa_seed=0):
  # Arguments:
  # v_M_1:   vector of initial states for individuals
  # n_i:     number of individuals
  # n_t:     total number of cycles to run the model
  # v_n:     vector of health state names
  # d_c:     discount rate for costs
  # d_e:     discount rate for health outcome (QALYs)
  # TR_out:  should the output include a microsimulation trace? (default is TRUE)
  # TS_out:  should the output include a matrix of transitions between states? (default is TRUE)
  # Trt:     are the n.i individuals receiving treatment? (scalar with a Boolean value, default is FALSE)
  # seed:    starting seed number for random number generator (default is 1)
  # Makes use of:
  # probs:   function for the estimation of transition probabilities
  # costs:   function for the estimation of cost state values
  # effs:    function for the estimation of state specific health outcomes (QALYs)

  # following code is used to create discount rate vectors for costs and qalys
  v_dwc = np.array([])
  v_dwe = np.array([])
  v_cycle = np.array([])
  for i in range(0, n_t + 1):
    result_c = 1 / ((1 + d_c) ** (np.maximum(0, i - 1)))  # calculate the cost discount weights
    result_e = 1 / ((1 + d_e) ** (np.maximum(0, i - 1)))  # calculate the QALY discount weights
    v_dwc = np.append(v_dwc, result_c)  # appends discount rate to vector v_dwc
    v_dwe = np.append(v_dwe, result_e)  # appends discount rate to vector v_dwe
    v_cycle = np.append(v_cycle, "Cycle" + str(i))

  # create the matrix capturing the state name/costs/health outcomes for all individuals at each time point
  m_M = pd.DataFrame("NA", index=(range(1, n_i + 1)), columns=v_cycle)
  m_C = pd.DataFrame(0, index=(range(1, n_i + 1)), columns=v_cycle)
  m_E = pd.DataFrame(0.00, index=(range(1, n_i + 1)), columns=v_cycle)

  m_M['Cycle0'] = "H"  # indicate the initial health state

  for i in range(1, n_i+1):
    np.random.seed(seed + i + psa_seed)  # set the seed for every individual for the random number generator
    m_C.iat[i - 1, 0] = costs(m_M.iloc[i - 1, 0], Trt)  # estimate costs
    m_E.iat[i - 1, 0] = effs(m_M.iloc[i - 1, 0], Trt)   # estimate QALYs

    for t in range(1, n_t + 1):
      v_p = probs(m_M.iloc[i - 1, t - 1])  # calculate the transition probabilities at cycle t

      m_M.iat[i - 1, t] = np.random.choice(v_n, p=v_p)    # sample the next health state and store that state
      m_C.iat[i - 1, t] = costs(m_M.iloc[i - 1, t], Trt)  # estimate costs per individual
      m_E.iat[i - 1, t] = effs(m_M.iloc[i - 1, t], Trt)   # estimate QALYs per individual

    # close the loop for the time points 
    print(np.round(i / n_i * 100, 2), "% done")  # display the progress of the simulation

  # close the loop for the individuals 

  tc = m_C.dot(v_dwc)  # total (discounted) cost per individual
  te = m_E.dot(v_dwe)  # total (discounted) QALYs per individual

  tc_hat = np.mean(tc)  # average (discounted) cost
  te_hat = np.mean(te)  # average (discounted) QALYs

  if TS_out:  # create a  matrix of transitions across states
    ts = pd.DataFrame("NA", index=(range(1, n_i + 1)), columns=range(0, n_t + 1))
    ts.iloc[:, 0] = m_M.iloc[:, 0]
    for i in range(1, n_t + 1):
      ts.iloc[:, i] = m_M.iloc[:, i-1] + "->" + m_M.iloc[:, i]  # transitions from one state to the other

  else:
    ts = None

  if TR_out:  # create a trace from the individual trajectories
    tr = pd.DataFrame(0, index=v_cycle, columns=v_n)
    for state in v_n:
      for cycle in v_cycle:
        tr.at[cycle, state] = (m_M[cycle] == state).sum()
    tr = tr / n_i
  else:
    tr = None

  return m_M, m_C, m_E, tc, te, tc_hat, te_hat, ts, tr  # return results of simulation


def probs(M_it):  # TODO: Can this be done in a matrix?
  v_p_it = [0] * n_s  # create vector of state transition probabilities

  # update v_p_it with the appropriate probabilities
  if M_it == "H":
    v_p_it = [1 - p_HS1 - p_HD, p_HS1, 0, p_HD]  # transition probabilities when healthy
  if M_it == "S1":
    v_p_it = [p_S1H, 1 - p_S1H - p_S1S2 - p_S1D, p_S1S2, p_S1D]  # transition probabilities when sick
  if M_it == "S2":
    v_p_it = [0, 0, 1 - p_S2D, p_S2D]  # transition probabilities when sicker
  if M_it == "D":
    v_p_it = [0, 0, 0, 1]  # transition probabilities when dead

  return v_p_it


### Costs function
# The Costs function estimates the costs at every cycle.

def costs(M_it, Trt=False):  # TODO: Can this be done in a vectorized approach?
  # M_it: health state occupied by individual i at cycle t (character variable)
  # Trt:  is the individual being treated? (default is FALSE)

  c_it = 0  # by default the cost for everyone is zero
  if M_it == "H":
    c_it = c_H  # update the cost if healthy
  if M_it == "S1":
    c_it = c_S1 + c_Trt * Trt  # update the cost if sick conditional on treatment
  if M_it == "S2":
    c_it = c_S2 + c_Trt * Trt  # update the cost if sicker conditional on treatment
  return c_it  # return the costs


### Health outcome function
# The Effs function to update the utilities at every cycle.

def effs(M_it, Trt=False, cl=1):  # TODO: Can this be done in a vectorized approach?

  # M_it: health state occupied by individual i at cycle t (character variable)
  # Trt:  is the individual treated? (default is FALSE)
  # cl:   cycle length (default is 1)

  u_it = 0  # by default the utility for everyone is zero
  if M_it == "H":
    u_it = u_H  # update the utility if healthy
  if M_it == "S1":
    u_it = Trt * u_Trt + (1 - Trt) * u_S1  # update the utility if sick conditional on treatment
  if M_it == "S2":
    u_it = u_S2  # update the utility if sicker
  qalys = u_it * cl  # calculate the QALYs during cycle t
  return qalys  # return the QALYs


def psa():
  v_results = {}  # empty vector of results
  for i in range(num_iter):
    psa_no_trt = microsim(v_M_1, n_i, n_t, v_n, d_c, d_e, Trt=False, psa_seed=i)
    psa_trt = microsim(v_M_1, n_i, n_t, v_n, d_c, d_e, Trt=True, psa_seed=i)
    v_C = [psa_no_trt[5], psa_trt[5]]
    v_E = [psa_no_trt[6], psa_trt[6]]
    delta_C = v_C[1] - v_C[0]  # calculate incremental costs
    delta_E = v_E[1] - v_E[0]  # calculate incremental QALYs
    v_results.update({i: {'Inc. Cost': delta_C, 'Inc. QALYs': delta_E}})
  psa_results = pd.DataFrame(v_results).T
  fig = psa_results.plot(kind='scatter', x='Inc. QALYs', y='Inc. Cost')
  return psa_results, fig


##################################### Run the simulation ##################################
sim_no_trt = microsim(v_M_1, n_i, n_t, v_n, d_c, d_e, Trt=False)  # run for no treatment
sim_trt = microsim(v_M_1, n_i, n_t, v_n, d_c, d_e, Trt=True)  # run for treatment
################################# Cost-effectiveness analysis #############################

#
# # store the mean costs (and the MCSE) of each strategy in a new variable v.C (vector costs)
v_C = [sim_no_trt[5], sim_trt[5]]
se_C = [np.std(sim_no_trt[3]), np.std(sim_trt[3])] / np.sqrt(n_i)
# # store the mean QALYs (and the MCSE) of each strategy in a new variable v.E (vector health outcomes)
v_E = [sim_no_trt[6], sim_trt[6]]
se_E = [np.std(sim_no_trt[4]), np.std(sim_trt[4])] / np.sqrt(n_i)
#
delta_C = v_C[1] - v_C[0]  # calculate incremental costs
delta_E = v_E[1] - v_E[0]  # calculate incremental QALYs
se_delta_E = np.std(sim_trt[4] - sim_no_trt[4]) / np.sqrt(n_i)  # Monte Carlo squared error (MCSE) of incremental costs
se_delta_C = np.std(sim_trt[3] - sim_no_trt[3]) / np.sqrt(n_i)  # Monte Carlo squared error (MCSE) of incremental QALYs
icer = delta_C / delta_E  # calculate the ICER
results = [delta_C, delta_E, icer]  # store the values in a new variable

#
# # Create full incremental cost-effectiveness analysis table
dict_results = {'Costs': {'No Treatment': np.round(v_C[0]), 'Treatment': np.round(v_C[1])},
                'QALYs': {'No Treatment': np.round(v_E[0], 3), 'Treatment': np.round(v_E[1], 3)},
                'Inc. Costs': {'No Treatment': '', 'Treatment': np.round(delta_C)},
                'Inc. QALYs': {'No Treatment': '', 'Treatment': np.round(delta_E, 3)},
                'ICER': {'No Treatment': '', 'Treatment': np.round(icer)}}
table_micro = pd.DataFrame(dict_results)
print(table_micro)

psa()
