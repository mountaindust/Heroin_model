import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode

#initial population values
S_0 = 0.87
P_0 = 0.1
A_0 = 0.02
R_0 = 0.01

#temporal info
tstart = 0
tstop = 10000

#parameters
params = {}
params['alpha'] = 0.2                   #S->P: prescription rate
params['beta'] = 0.006                  #total S->A addiction rate
params['delta'] = 0.09                  #R->S: finish recovery
params['epsilon'] = 0.76                #P->S rate
params['gamma'] = 1 - params['epsilon'] #P->A
params['xi'] = 0.505                    #fraction of beta due to P
params['zeta'] = 0.25                   #A->R rate of starting treatment
params['nu'] = 1-params['delta']/2      #R->A treatment relapse due to A
params['mu'] = 0.008237                 #nomral death rate
params['mu_star'] = 0.00834             #addiction death rate
params['sigma'] = 1-params['delta']-params['mu']    #R->A natural treatment relapse

def update_params(new_params):
    '''Update the params dict with new values contained in new_params'''
    global params
    for key,val in new_params.items():
        if key in params.keys():
            params[key] = val
        else:
            print('Could not find parameter {}.'.format(key))



def opioid_odes(t, X, params):
    '''Specifies the opioid model as a system of ODEs.'''
    pass



def solve_odes(S0,P0,A0,R0,tstart,tstop,params):
    '''Solve opioid_odes with given initial conditions, time, and params.'''
    pass
    