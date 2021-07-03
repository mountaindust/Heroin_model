import pandas as pd
import numpy as np
from SPAHR_class import SPAHR_ABM

# Loading parameters into model
params = {}
params['turtle-count']  = 2000
params['num-ticks']     = 1000
params['delta_t']       = 0.01
params['mu']            = 0.0071

params['m_tilde']       = -0.0056
params['b_tilde']       = 0.27 
params['c_tilde']       = -0.027
params['beta_A']        = 0.000878
params['beta_P']        = 0.0000654
params['theta_1']       = 0.222

params['gamma']         = 0.00505
params['epsilon']       = 2.53
params['theta_2']       = 0.236

params['zeta']          = 0.198
params['theta_3']       = 19.7
params['d_tilde']       = 0.000977
params['e_tilde']       = 0.00883

params['nu']            = 0.000531
params['mu_H']          = 0.0466

params['sigma']         = 0.102

params['P0']  = 0.095
params['A0']  = 0.00710
params['H0']  = 0.000465
params['R0']  = 0.00507
params['S0']  = 1 - (params['P0'] + params['A0'] + params['H0'] + params['R0'])

# These need only to be input if alpha and mu_a are constant in the model
#params['alpha'] = 0.27
#params['mu_a']  = 0.00883

model = SPAHR_ABM(parameters = params, constants = [0, 0])
# constants = [1,1] tells us that [alpha, mu_a] are constant

model.run_model(n_runs = 10, relapse_mem = False)

model.data_to_csv(file_name = "6-27_timedep_SPAHR_ABM_faster.csv")
