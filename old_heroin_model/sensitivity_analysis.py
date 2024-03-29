'''Use SALib to perform Sobol sensitivity analysis on the model parameters'''

import sys, os
import pickle
import argparse
from multiprocessing import Pool
from SALib.sample import saltelli
from SALib.analyze import sobol
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import heroin_model

default_N = os.cpu_count()
parser = argparse.ArgumentParser()
parser.add_argument("-N", type=int, default=1000,
                    help="obtain N*(2D+2) samples from parameter space")
parser.add_argument("-n", "--ncores", type=int,
                    help="number of cores, defaults to {}".format(default_N))
parser.add_argument("-o", "--filename", type=str, 
                    help="filename to write output to, no extension",
                    default='analysis')

### Define a model wrapper based on the parameter space in main() ###
def run_model(alpha, beta, delta1, delta2, mu1, mu2, irate, hrate,
                epsilon, xi, sigma, gamma_0, r, g):
    # Length to run each model realization
    endt = 1000
    # Copy default parameter dict
    params = dict(heroin_model.params)
    # Replace parameter values
    params['alpha'] = alpha
    params['beta'] = beta
    params['delta1'] = delta1
    params['delta2'] = delta2
    params['mu1'] = mu1
    params['mu2'] = mu2
    params['irate'] = irate
    params['hrate'] = hrate
    params['epsilon'] = epsilon
    params['xi'] = xi
    params['sigma'] = sigma
    params['gamma_0'] = gamma_0
    params['r'] = r
    params['g'] = g
    # Get initial conditions
    s0 = heroin_model.s0
    e0 = heroin_model.e0
    i0 = heroin_model.i0
    h0 = heroin_model.h0
    r0 = heroin_model.r0
    # Run model
    try:
        result = heroin_model.heroin_solve_odes(s0,e0,i0,h0,r0,
                                                version=1,tsteps=endt,params=params)
    except:
        return (sys.exc_info()[1], 
                None, None, None, None, None, None, None, None, None)
    # Return just the end value and std of each variable
    if result[4] is not None:
        return (np.mean(result[0][-50:]), np.mean(result[1][-50:]), 
                np.mean(result[2][-50:]), np.mean(result[3][-50:]), 
                np.mean(result[4][-50:]), np.std(result[0][-50:]),
                np.std(result[1][-50:]), np.std(result[2][-50:]),
                np.std(result[3][-50:]), np.std(result[4][-50:]))
    else:
        return (np.mean(result[0][-50:]), np.mean(result[1][-50:]), 
                np.mean(result[2][-50:]), np.mean(result[3][-50:]), 
                0, np.std(result[0][-50:]),
                np.std(result[1][-50:]), np.std(result[2][-50:]),
                np.std(result[3][-50:]), 0)



def main(N, filename, pool=None):
    '''Uses the ODE version of the heroin model and conducts parameter sensitivity'''

    ### Define the parameter space ###

    problem = {
        'num_vars': 14, #number of parameters
        'names': ['alpha', 'beta', 'delta1', 'delta2',
                  'mu1', 'mu2', 'irate', 'hrate', 'epsilon',
                  'xi', 'sigma', 'gamma_0', 'r', 'g'],
        'bounds': [[0.01, 0.5], [0.01, 0.5], [0, 0.2], [0, 0.8],
                   [0, 0.15], [0, 0.3], [0, 1], [0, 1], [0.01, 0.8],
                   [0, 1], [0, 0.95], [0.1, 10], [0.00001, 0.95], [0, 0.5]]
    }

    ### Create an N by num_var matrix of parameter values ###
    param_values = saltelli.sample(problem, N, calc_second_order=True)

    ### Run model ###
    print('Examining the parameter space.')
    poolsize = os.cpu_count()
    chunksize = param_values.shape[0]//poolsize
    output = pool.starmap(run_model, param_values, chunksize=chunksize)

    ### Parse and save the output ###
    print('Saving and reviewing the results...')
    param_values = pd.DataFrame(param_values, columns=problem['names'])
    # write data to temporary location in case of errors
    with open("raw_result_data.pickle", "wb") as f:
        result = {'output':output, 'param_values':param_values}
        pickle.dump(result, f)
    error_num = 0
    error_places = []
    for n, result in enumerate(output):
        if result[1] is None:
            error_num += 1
            error_places.append(n)
    if error_num > 0:
        print("Errors discovered in output.")
        print("Parameter locations: {}".format(error_places))
        print("Please review pickled output.")
        return
    print('Parsing the results...')
    output = np.array(output)
    # Resave as dataframe in hdf5
    store = pd.HDFStore(filename+'.h5')
    store['param_values'] = param_values
    store['raw_output'] = pd.DataFrame(output, columns=['S', 'E', 'I', 'H', 'R',
                                       'std_S', 'std_E', 'std_I', 'std_H', 'std_R'])
    os.remove('raw_result_data.pickle')
    

    ### Analyze the results and view using Pandas ###
    # Conduct the sobol analysis and pop out the S2 results to a dict
    std_sum = output[:,5:].sum(1)
    S2 = {}
    S_sens = sobol.analyze(problem, output[:,0], calc_second_order=True)
    S2['S'] = pd.DataFrame(S_sens.pop('S2'), index=problem['names'],
                           columns=problem['names'])
    S2['S_conf'] = pd.DataFrame(S_sens.pop('S2_conf'), index=problem['names'],
                           columns=problem['names'])
    E_sens = sobol.analyze(problem, output[:,1], calc_second_order=True)
    S2['E'] = pd.DataFrame(E_sens.pop('S2'), index=problem['names'],
                           columns=problem['names'])
    S2['E_conf'] = pd.DataFrame(E_sens.pop('S2_conf'), index=problem['names'],
                           columns=problem['names'])
    I_sens = sobol.analyze(problem, output[:,2], calc_second_order=True)
    S2['I'] = pd.DataFrame(I_sens.pop('S2'), index=problem['names'],
                           columns=problem['names'])
    S2['I_conf'] = pd.DataFrame(I_sens.pop('S2_conf'), index=problem['names'],
                           columns=problem['names'])
    H_sens = sobol.analyze(problem, output[:,3], calc_second_order=True)
    S2['H'] = pd.DataFrame(H_sens.pop('S2'), index=problem['names'],
                           columns=problem['names'])
    S2['H_conf'] = pd.DataFrame(H_sens.pop('S2_conf'), index=problem['names'],
                           columns=problem['names'])
    R_sens = sobol.analyze(problem, output[:,4], calc_second_order=True)
    S2['R'] = pd.DataFrame(R_sens.pop('S2'), index=problem['names'],
                           columns=problem['names'])
    S2['R_conf'] = pd.DataFrame(R_sens.pop('S2_conf'), index=problem['names'],
                           columns=problem['names'])
    var_sens = sobol.analyze(problem, std_sum, calc_second_order=True)
    S2['var'] = pd.DataFrame(var_sens.pop('S2'), index=problem['names'],
                             columns=problem['names'])
    S2['var_conf'] = pd.DataFrame(var_sens.pop('S2_conf'), index=problem['names'],
                                  columns=problem['names'])
    # Convert the rest to a pandas dataframe
    S_sens = pd.DataFrame(S_sens,index=problem['names'])
    E_sens = pd.DataFrame(E_sens,index=problem['names'])
    I_sens = pd.DataFrame(I_sens,index=problem['names'])
    H_sens = pd.DataFrame(H_sens,index=problem['names'])
    R_sens = pd.DataFrame(R_sens,index=problem['names'])
    var_sens = pd.DataFrame(var_sens,index=problem['names'])
    # Plot
    plot_S1_ST(S_sens, E_sens, I_sens, H_sens, R_sens, False)
    plot_var(var_sens, False)
    plt.show()
    
    ### Save the analysis ###
    print('Saving...')
    store['S_sens'] = S_sens
    store['E_sens'] = E_sens
    store['I_sens'] = I_sens
    store['H_sens'] = H_sens
    store['R_sens'] = R_sens
    store['var_sens'] = var_sens
    for key in S2.keys():
        store['S2/'+key] = S2[key]
    store.close()



def load_data(filename):
    '''Load analysis data from previous run and return for examination'''

    return pd.HDFStore(filename)



def plot_S1_ST(S_sens, E_sens, I_sens, H_sens, R_sens, show=True):
    # Gather the S1 and ST results
    S1 = pd.concat({'S':S_sens['S1'], 'E':E_sens['S1'], 'I':I_sens['S1'], 
                   'H':H_sens['S1'], 'R':R_sens['S1']}, axis=1)
    ST = pd.concat({'S':S_sens['ST'], 'E':E_sens['ST'], 'I':I_sens['ST'], 
                   'H':H_sens['ST'], 'R':R_sens['ST']}, axis=1)
    # Plot
    fig, axes = plt.subplots(ncols=2, figsize=(8, 4.5))
    S1.plot.bar(stacked=True, ax=axes[0], title='First-order indices')
    ST.plot.bar(stacked=True, ax=axes[1], title='Total-order indices')
    plt.tight_layout()
    if show:
        plt.show()
    return (fig, axes)



def plot_var(var_sens, show=True):
    fig, axes = plt.subplots(ncols=2, figsize=(8, 4.5))
    var_sens['S1'].plot.bar(ax=axes[0], title='First-order variance')
    var_sens['ST'].plot.bar(ax=axes[1], title='Total-order variance')
    plt.tight_layout()
    if show:
        plt.show()
    return (fig, axes)



if __name__ == "__main__":
    args = parser.parse_args()
    if args.ncores is None:
        with Pool() as pool:
            main(args.N, filename=args.filename, pool=pool)
    else:
        with Pool(args.ncores) as pool:
            main(args.N, filename=args.filename, pool=pool)