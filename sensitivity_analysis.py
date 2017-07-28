'''Use SALib to perform Sobol sensitivity analysis on the model parameters'''

from SALib.sample import saltelli
from SALib.analyze import sobol
import numpy as np

def main():
    '''Uses the ODE version of the heroin model and conducts parameter sensitivity'''

    #Define the parameter space

    problem = {
        'num_vars': 13, #number of parameters
        'names': ['alpha', 'beta', 'delta1', 'delta2',
                  'mu1', 'mu2', 'irate', 'hrate', 'epsilon',
                  'xi', 'sigma', 'gamma_0', 'dr_add'],
        'bounds': [[0.01, 0.5], [0.01, 0.5], [0, 0.2], [0, 0.8],
                   [0, 0.15], [0, 0.3], [0, 1], [0, 1], [0.01, 0.8],
                   [0, 1], [0, 0.95], [0.1, 10], [0.00001, 0.95]]
    }