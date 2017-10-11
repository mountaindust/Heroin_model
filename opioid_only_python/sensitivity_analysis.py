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
import opioid_model

default_N = os.cpu_count()
parser = argparse.ArgumentParser()
parser.add_argument("-N", type=int, default=1000,
                    help="obtain N*(2D+2) samples from parameter space")
parser.add_argument("-n", "--ncores", type=int,
                    help="number of cores, defaults to {}".format(default_N))
parser.add_argument("-o", "--filename", type=str, 
                    help="filename to write output to, no extension",
                    default='analysis')


def run_reduced_model(alpha,beta,delta,epsilon,zeta,nu,mu,mu_star,sigma):
    '''Defines a model wrapper based on the parameter space in main()'''
    # Length to run each model
    tstart = 0
    tstop =10000
    # Copy default parameter dict
    params = dict(opioid_model.params)
    # Set gamma and xi equal to zero
    params['gamma'] = 0
    params['xi'] = 0
    # Replace other parameter values
    params['alpha'] = alpha
    params['beta'] = beta
    params['delta'] = delta
    params['epsilon'] = epsilon
    params['zeta'] = zeta
    params['nu'] = nu
    params['mu'] = mu
    params['mu_star'] = mu_star
    params['sigma'] = sigma
    # Get initial conditions
    S_0 = opioid_model.S_0
    P_0 = opioid_model.P_0
    A_0 = opioid_model.A_0
    R_0 = opioid_model.R_0
    # Run model
    try:
        result = opioid_model.solve_odes(S_0,P_0,A_0,R_0,tstart,tstop,params)
    except:
        return (sys.exc_info()[1],None,None,None)
    # Return just the mean end value of each variable
    return (np.mean(result[0][-100:]), np.mean(result[1][-100:]),
            np.mean(result[2][-100:]), np.mean(result[3][-100:]))



def main(N, filename, pool=None):
    '''Runs parameter sensitivity on the reduced opioid model'''

    ### Define the parameter space ###

    problem = {
        'num_vars': 9, #number of parameters
        'names': ['alpha', 'beta', 'delta', 'epsilon', 'zeta', 'nu',
                    'mu', 'mu_star', 'sigma'],
        'bounds': [[0,2], [0,2], [0,2], [0,2], [0,2], [0,2],
                    [0,0.1], [0,0.5], [0,2]]
    }

    ### Create an N by num_var matrix of parameter values ###
    param_values = saltelli.sample(problem, N, calc_second_order=True)

    ### Run model ###
    print('Examining the parameter space.')
    poolsize = os.cpu_count()
    chunksize = param_values.shape[0]//poolsize
    output = pool.starmap(run_reduced_model, param_values, chunksize=chunksize)

    ### Parse and save the output ###
    print('Saving and reviewing the results...')
    param_values = pd.DataFrame(param_values, columns=problem['names'])
    # write data to temporary location in case of errors
    with open("raw_result_data.pickle", "wb") as f:
        result = {'output':output, 'param_values':param_values}
        pickle.dump(result, f)
    # Look for errors
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
    # Save results in HDF5 as dataframe
    print('Parsing the results...')
    output = np.array(output)
    # Resave as dataframe in hdf5
    store = pd.HDFStore(filename+'.h5')
    store['param_values'] = param_values
    store['raw_output'] = pd.DataFrame(output, columns=['S', 'P', 'A', 'R'])
    os.remove('raw_result_data.pickle')

    ### Analyze the results and view using Pandas ###
    # Conduct the sobol analysis and pop out the S2 results to a dict
    S2 = {}
    S_sens = sobol.analyze(problem, output[:,0], calc_second_order=True)
    S2['S'] = pd.DataFrame(S_sens.pop('S2'), index=problem['names'],
                           columns=problem['names'])
    S2['S_conf'] = pd.DataFrame(S_sens.pop('S2_conf'), index=problem['names'],
                           columns=problem['names'])
    P_sens = sobol.analyze(problem, output[:,1], calc_second_order=True)
    S2['P'] = pd.DataFrame(P_sens.pop('S2'), index=problem['names'],
                           columns=problem['names'])
    S2['P_conf'] = pd.DataFrame(P_sens.pop('S2_conf'), index=problem['names'],
                           columns=problem['names'])
    A_sens = sobol.analyze(problem, output[:,2], calc_second_order=True)
    S2['A'] = pd.DataFrame(A_sens.pop('S2'), index=problem['names'],
                           columns=problem['names'])
    S2['A_conf'] = pd.DataFrame(A_sens.pop('S2_conf'), index=problem['names'],
                           columns=problem['names'])
    R_sens = sobol.analyze(problem, output[:,3], calc_second_order=True)
    S2['R'] = pd.DataFrame(R_sens.pop('S2'), index=problem['names'],
                           columns=problem['names'])
    S2['R_conf'] = pd.DataFrame(R_sens.pop('S2_conf'), index=problem['names'],
                           columns=problem['names'])
    # Convert the rest to a pandas dataframe
    S_sens = pd.DataFrame(S_sens,index=problem['names'])
    P_sens = pd.DataFrame(P_sens,index=problem['names'])
    A_sens = pd.DataFrame(A_sens,index=problem['names'])
    R_sens = pd.DataFrame(R_sens,index=problem['names'])
    # Plot
    plot_S1_ST(S_sens, P_sens, A_sens, R_sens, True)

    ### Save the analysis ###
    print('Saving...')
    store['S_sens'] = S_sens
    store['P_sens'] = P_sens
    store['A_sens'] = A_sens
    store['R_sens'] = R_sens
    for key in S2.keys():
        store['S2/'+key] = S2[key]
    store.close()



def load_data(filename):
    '''Load analysis data from previous run and return for examination'''

    return pd.HDFStore(filename)



def plot_S1_ST(S_sens, P_sens, A_sens, R_sens, show=True):
    # Gather the S1 and ST results
    S1 = pd.concat({'S':S_sens['S1'], 'P':P_sens['S1'], 'A':A_sens['S1'], 
                   'R':R_sens['S1']}, axis=1)
    ST = pd.concat({'S':S_sens['ST'], 'P':P_sens['ST'], 'A':A_sens['ST'], 
                   'R':R_sens['ST']}, axis=1)
    # Plot
    fig, axes = plt.subplots(ncols=2, figsize=(8, 4.5))
    S1.plot.bar(stacked=True, ax=axes[0], title='First-order indices')
    ST.plot.bar(stacked=True, ax=axes[1], title='Total-order indices')
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

