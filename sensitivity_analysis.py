'''Use SALib to perform Sobol sensitivity analysis on the model parameters'''

import sys, os, time
from collections import OrderedDict
import pickle
import argparse
from multiprocessing import Pool
from SALib.sample import saltelli
from SALib.analyze import sobol
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import heroin_model

default_N = os.cpu_count() - 4
parser = argparse.ArgumentParser()
parser.add_argument("-N", type=int, default=1000,
                    help="obtain N*(2D+2) samples from parameter space")
parser.add_argument("-n", "--ncores", type=int,
                    help="number of cores, defaults to {}".format(default_N))
parser.add_argument("-o", "--filename", type=str, 
                    help="filename to write output to, no extension",
                    default='analysis_{}'.format(time.strftime("%m_%d_%H%M")))
parser.add_argument("--reduced_model", action="store_true",
                    help="run over reduced parameters")
parser.add_argument("--noplot", action="store_true",
                    help="do not plot results at the end")


# Ignore this for now.
def run_reduced_model(alpha,beta,delta,epsilon,zeta,nu,mu,mu_star,sigma):
    '''Defines a model wrapper based on the parameter space in main()'''
    raise NotImplementedError("We haven't found reduced model yet!")

    # Length to run each model
    # tstart = 0
    # tstop =10000
    # # Copy default parameter dict
    # params = dict(heroin_model.params)
    # # Set gamma and xi equal to zero
    # params['gamma'] = 0
    # params['xi'] = 0
    # # Replace other parameter values
    # params['alpha'] = alpha
    # params['beta'] = beta
    # params['delta'] = delta
    # params['epsilon'] = epsilon
    # params['zeta'] = zeta
    # params['nu'] = nu
    # params['mu'] = mu
    # params['mu_star'] = mu_star
    # params['sigma'] = sigma
    # # Get initial conditions
    # S_0 = 0.897
    # P_0 = 0.1
    # A_0 = 0.002
    # R_0 = 0.001
    # # Run model
    # try:
    #     result = heroin_model.solve_odes(S_0,P_0,A_0,R_0,tstart,tstop,params)
    # except:
    #     return (sys.exc_info()[1],None,None,None)
    # # Return just the mean end value of each variable
    # return (np.mean(result[0][-100:]), np.mean(result[1][-100:]),
    #         np.mean(result[2][-100:]), np.mean(result[3][-100:]))



def run_full_model(alpha,beta,delta,epsilon,gamma,xi,zeta,nu,mu,
                   theta_1,theta_2,theta_3,mu_A,mu_H,sigma_A,sigma_H):
    '''Defines a model wrapper based on the parameter space in main()'''
    # Length to run each model
    tstart = 0
    tstop = 10,000
    # Copy default parameter dict
    params = dict(heroin_model.params)
    # Replace other parameter values
    params['alpha'] = alpha
    params['beta'] = beta
    params['delta'] = delta
    params['epsilon'] = epsilon
    params['gamma'] = gamma
    params['xi'] = xi
    params['zeta'] = zeta
    params['nu'] = nu
    params['mu'] = mu
    params['theta_1'] = theta_1
    params['theta_2'] = theta_2
    params['theta_3'] = theta_3
    params['mu_A'] = mu_A
    params['mu_H'] = mu_H
    params['sigma_A'] = sigma_A
    params['sigma_H'] = sigma_H
    
    # Get initial conditions
    S_0 = 0.6221
    P_0 = 0.37
    A_0 = 0.0062
    H_0 = 0.0014
    R_0 = 0.0003
    # Run model
    try:
        result = heroin_model.solve_odes(S_0,P_0,A_0,H_0,R_0,tstart,tstop,params)
    except:
        return (sys.exc_info()[1],None,None,None,None)
    # Return just the mean end value of each variable
    return (np.mean(result[0][-100:]), np.mean(result[1][-100:]),
            np.mean(result[2][-100:]), np.mean(result[3][-100:]),
            np.mean(result[4][-100:]))



def main(N, filename, reduced, pool=None, no_plot=False):
    '''Runs parameter sensitivity on the reduced opioid model'''

    ##### Define the parameter space ######
    if reduced:
        raise NotImplementedError("Reduced model not currently implemented!")
        problem = {
            'num_vars': 9, #number of parameters
            'names': ['alpha', 'beta', 'delta', 'epsilon', 'zeta', 'nu',
                        'mu', 'mu_star', 'sigma'],
            'bounds': [[0,1], [0,1], [0,1], [0,1], [0,1], [0,1],
                        [0,0.1], [0,0.5], [0,1]]
        }
    else:
        problem = {
            'num_vars': 16, #number of parameters
            'names': ['alpha', 'beta', 'delta', 'epsilon', 'gamma', 'xi',
                      'zeta', 'nu', 'mu', 'theta_1', 'theta_2', 'theta_3',
                      'mu_A', 'mu_H', 'sigma_A', 'sigma_H'],
            'bounds': [[0,1], [0,1], [0,1], [0,1], [0,1], [0,1],
                       [0,1], [0,1], [0,0.1], [0,1], [0,1], [0,1],
                       [0,.1], [0,.1], [0,1], [0,1]] #xi was always 0,1
        }

    ### Create an N by num_var matrix of parameter values ###
    param_values = saltelli.sample(problem, N, calc_second_order=True)

    ### Run model ###
    print('Examining the parameter space.')
    if args.ncores is None:
        poolsize = os.cpu_count() - 4
    else:
        poolsize = args.ncores
    chunksize = param_values.shape[0]//poolsize
    if reduced:
        output = pool.starmap(run_reduced_model, param_values, chunksize=chunksize)
    else:
        output = pool.starmap(run_full_model, param_values, chunksize=chunksize)

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
    store['raw_output'] = pd.DataFrame(output, columns=['S', 'P', 'A', 'H', 'R'])
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
    # Convert the rest to a pandas dataframe
    S_sens = pd.DataFrame(S_sens,index=problem['names'])
    P_sens = pd.DataFrame(P_sens,index=problem['names'])
    A_sens = pd.DataFrame(A_sens,index=problem['names'])
    H_sens = pd.DataFrame(H_sens,index=problem['names'])
    R_sens = pd.DataFrame(R_sens,index=problem['names'])

    ### Save the analysis ###
    print('Saving...')
    store['S_sens'] = S_sens
    store['P_sens'] = P_sens
    store['A_sens'] = A_sens
    store['H_sens'] = H_sens
    store['R_sens'] = R_sens
    for key in S2.keys():
        store['S2/'+key] = S2[key]
    store.close()

    # Plot
    if not no_plot:
        plot_S1_ST(S_sens, P_sens, A_sens, H_sens, R_sens, False)



def load_data(filename):
    '''Load analysis data from previous run and return for examination'''

    return pd.HDFStore(filename)



def plot_S1_ST_from_store(store, show=True):
    '''Extract and plot S1 and ST sensitivity data directly from a store object'''

    plot_S1_ST(store['S_sens'], store['P_sens'], store['A_sens'],
               store['H_sens'], store['R_sens'], show)



def print_max_conf(store):
    '''Print off the max confidence interval for each variable in the store,
    for both first-order and total-order indices'''
    for var in ['S_sens', 'P_sens', 'A_sens', 'H_sens', 'R_sens']:
        print('----------- '+var+' -----------')
        print('S1_conf_max: {}'.format(store[var]['S1_conf'].max()))
        print('ST_conf_max: {}'.format(store[var]['ST_conf'].max()))
        print(' ')



def plot_S1_ST(S_sens, P_sens, A_sens, H_sens, R_sens, show=True):
    # Gather the S1 and ST results
    S1 = pd.concat([S_sens['S1'], P_sens['S1'], A_sens['S1'], H_sens['S1'],
                   R_sens['S1']], keys=['S','P','A','H','R'], axis=1)
    ST = pd.concat([S_sens['ST'], P_sens['ST'], A_sens['ST'], H_sens['ST'],
                   R_sens['ST']], keys=['S','P','A','H','R'], axis=1)
    # Plot
    fig, axes = plt.subplots(ncols=2, figsize=(12, 6))
    S1.plot.bar(stacked=True, ax=axes[0])
    ST.plot.bar(stacked=True, ax=axes[1])
    for ax in axes:
        ax.tick_params(labelsize=18)
        ax.legend(fontsize=16)
    axes[0].set_title('First-order indices', fontsize=26)
    axes[1].set_title('Total-order indices', fontsize=26)
    plt.tight_layout()
    if show:
        plt.show()
    else:
        fig.savefig("param_sens_{}.pdf".format(time.strftime("%m_%d_%H%M")))
    return (fig, axes)



### This is left-over from journal article plotting for opioid only ###
def plot_sens_data_double():
    store = load_data("Sensitivity_R0_150000.h5")
    S_sens_02 = store['S_sens']
    A_sens_02 = store['A_sens']
    P_sens_02 = store['P_sens']
    R_sens_02 = store['R_sens']
    store01 = load_data("Sens_R0_small0_01_150000.h5")
    A_sens_01 = store01['A_sens']
    P_sens_01 = store01['P_sens']
    S_sens_01 = store01['S_sens']
    R_sens_01 = store01['R_sens']
    plot_S1_ST_double(S_sens_01, P_sens_01, A_sens_01, R_sens_01,
                      S_sens_02, P_sens_02, A_sens_02, R_sens_02)
    store.close()
    store01.close()
    store = load_data("Sens_full_150000.h5")
    S_sens_02 = store['S_sens']
    A_sens_02 = store['A_sens']
    P_sens_02 = store['P_sens']
    R_sens_02 = store['R_sens']
    store01 = load_data("Sens_full_small0_01_150000.h5")
    A_sens_01 = store01['A_sens']
    P_sens_01 = store01['P_sens']
    S_sens_01 = store01['S_sens']
    R_sens_01 = store01['R_sens']
    plot_S1_ST_double(S_sens_01, P_sens_01, A_sens_01, R_sens_01,
                      S_sens_02, P_sens_02, A_sens_02, R_sens_02)



### This is left-over from journal article plotting for opioid only ###                      
def plot_S1_ST_double(S_sens_01, P_sens_01, A_sens_01, R_sens_01,
                      S_sens_02, P_sens_02, A_sens_02, R_sens_02, H="x"):
    '''Plot unit hypercube with double hypercube'''
    plt.rc('font',family='Arial')
    S1_01 = pd.concat([S_sens_01['S1'], P_sens_01['S1'], A_sens_01['S1'],
                      R_sens_01['S1']], keys=['S','P','A','R'], axis=1)
    ST_01 = pd.concat([S_sens_01['ST'], P_sens_01['ST'], A_sens_01['ST'],
                      R_sens_01['ST']], keys=['S','P','A','R'], axis=1)
    S1_02 = pd.concat([S_sens_02['S1'], P_sens_02['S1'], A_sens_02['S1'],
                      R_sens_02['S1']], keys=['S','P','A','R'], axis=1)
    ST_02 = pd.concat([S_sens_02['ST'], P_sens_02['ST'], A_sens_02['ST'],
                      R_sens_02['ST']], keys=['S','P','A','R'], axis=1)
    print(S1_01)
    # Plot
    fig, axes = plt.subplots(ncols=2, figsize=(12, 6))
    S1_01.plot.bar(stacked=True, ax=axes[0], linewidth=0, legend=False, grid=False)
    ST_01.plot.bar(stacked=True, ax=axes[1], linewidth=0, legend=False, grid=False)
    S1_02.plot.bar(stacked=True, ax=axes[0], linewidth=0, legend=False, grid=False)
    ST_02.plot.bar(stacked=True, ax=axes[1], linewidth=0, legend=False, grid=False)
    for ax in axes:
        h,l = ax.get_legend_handles_labels()
        for i in range(0,8,4):
            for j, pa in enumerate(h[i:i+4]):
                for rect in pa.patches:
                    rect.set_x(rect.get_x()+5/12.*i/4.)
                    if i==0:
                        rect.set_hatch(H)
                    rect.set_width(5/12.)
        ax.set_xticks((np.arange(0, 2*len(S1_01.index), 2) + 5/ 12.)/ 2.)
        ax.tick_params(labelsize=18)
        ax.set_xticklabels(S1_01.index, fontsize=21)
        # add invisible data to add legend
        n = []
        n.append(ax.bar(0, 0, color="gray", hatch=H))
        n.append(ax.bar(0, 0, color="gray", hatch=""))
        l1 = ax.legend(h[4:8], l[4:8], loc=[0.8, 0.64], fontsize=16)
        l2 = ax.legend(n, ["[0,1]","[0,2]"], loc=[0.52, 0.805], fontsize=16)
        ax.add_artist(l1)
    axes[0].set_title('First-order indices', fontsize=26)
    axes[1].set_title('Total-order indices', fontsize=26)
    plt.tight_layout()
    plt.show()



if __name__ == "__main__":
    args = parser.parse_args()
    red = args.reduced_model
    noplot = args.noplot
    if args.ncores is None:
        with Pool() as pool:
            main(args.N, filename=args.filename, reduced=red, pool=pool,
                 no_plot=noplot)
    else:
        with Pool(args.ncores) as pool:
            main(args.N, filename=args.filename, reduced=red, pool=pool,
                 no_plot=noplot)

