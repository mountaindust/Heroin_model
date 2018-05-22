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
from PIL import Image
from matplotlib import gridspec
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
parser.add_argument("--red_model", action="store_true",
                    help="run over only R0 parameters")


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
    S_0 = 0.9435
    P_0 = 0.05 #Study: Boudreau et al +  time increase??
    A_0 = 0.0062 #SAMHSA: https://www.samhsa.gov/data/sites/default/files/NSDUH-FFR2-2015/NSDUH-FFR2-2015.pdf
    R_0 = 0.0003 #HSS Treatment Episode Data Set https://www.samhsa.gov/data/sites/default/files/2014_Treatment_Episode_Data_Set_National_Admissions_9_19_16.pdf
    # Run model
    try:
        result = opioid_model.solve_odes(S_0,P_0,A_0,R_0,tstart,tstop,params)
    except:
        return (sys.exc_info()[1],None,None,None)
    # Return just the mean end value of each variable
    return (result[0][-1], result[1][-1], result[2][-1], result[3][-1])



def run_full_model(alpha,beta,delta,epsilon,gamma,xi,zeta,nu,mu,mu_star,sigma):
    '''Defines a model wrapper based on the parameter space in main()'''
    # Length to run each model
    tstart = 0
    tstop = 10
    #tstop =10000
    # Copy default parameter dict
    params = dict(opioid_model.params)
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
    params['mu_star'] = mu_star
    params['sigma'] = sigma
    # Get initial conditions
    S_0 = 0.9435
    P_0 = 0.05 #Study: Boudreau et al +  time increase??
    A_0 = 0.0062 #SAMHSA: https://www.samhsa.gov/data/sites/default/files/NSDUH-FFR2-2015/NSDUH-FFR2-2015.pdf
    R_0 = 0.0003 #HSS Treatment Episode Data Set https://www.samhsa.gov/data/sites/default/files/2014_Treatment_Episode_Data_Set_National_Admissions_9_19_16.pdf
    # Run model
    try:
        result = opioid_model.solve_odes(S_0,P_0,A_0,R_0,tstart,tstop,params)
    except:
        return (sys.exc_info()[1],None,None,None)
    # Return just the mean end value of each variable
    return (result[0][-1], result[1][-1], result[2][-1], result[3][-1])



def main(N, filename, reduced, pool=None):
    '''Runs parameter sensitivity on the reduced opioid model'''

    ### Define the parameter space ###
    if reduced:
        problem = {
            'num_vars': 9, #number of parameters
            'names': ['alpha', 'beta', 'delta', 'epsilon', 'zeta', 'nu',
                        'mu', 'mu_star', 'sigma'],
            'bounds': [[0.03,0.3], [0.0003,0.03], [0.01,1], [0.8,8], [0.1,2], [0.01,1],
                       [0.001,0.01], [0.005,0.1], [0.01,1]]
        }
    else:
        problem = {
            'num_vars': 11, #number of parameters
            'names': ['alpha', 'beta', 'delta', 'epsilon', 'gamma', 'xi',
                      'zeta', 'nu', 'mu', 'mu_star', 'sigma'],
            'bounds': [[0.02,0.2], [0.00114,0.0114], [0,1], [0.8,8], [0.00235,0.0235], [0,1],
                       [0.2,2], [0,1], [0.002305,0.02305], [0.003652,0.03652], [0,1]] #xi is always 0,1
        }

    ### Create an N by num_var matrix of parameter values ###
    param_values = saltelli.sample(problem, N, calc_second_order=True)

    ### Run model ###
    print('Examining the parameter space.')
    if args.ncores is None:
        poolsize = os.cpu_count()
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

    ### Save the analysis ###
    print('Saving...')
    store['S_sens'] = S_sens
    store['P_sens'] = P_sens
    store['A_sens'] = A_sens
    store['R_sens'] = R_sens
    for key in S2.keys():
        store['S2/'+key] = S2[key]
    store.close()

    # Plot
    plot_S1_ST(S_sens, P_sens, A_sens, R_sens, True)



def load_data(filename):
    '''Load analysis data from previous run and return for examination'''

    return pd.HDFStore(filename)



def plot_S1_ST_from_store(store, show=True):
    '''Extract and plot S1 and ST sensitivity data directly from a store object'''

    plot_S1_ST(store['S_sens'], store['P_sens'], store['A_sens'],
               store['R_sens'], show)



def plot_S1_ST_double_from_stores(store_01, store_02):
    '''Plot [0,1] and [0,2] from two loaded stores'''

    plot_S1_ST_double(store_01['S_sens'], store_01['P_sens'], store_01['A_sens'],
                      store_01['R_sens'], store_02['S_sens'], store_02['P_sens'],
                      store_02['A_sens'], store_02['R_sens'],)



def print_max_conf(store):
    '''Print off the max confidence interval for each variable in the store,
    for both first-order and total-order indices'''
    for var in ['S_sens', 'P_sens', 'A_sens', 'R_sens']:
        print('----------- '+var+' -----------')
        print('S1_conf_max: {}'.format(store[var]['S1_conf'].max()))
        print('ST_conf_max: {}'.format(store[var]['ST_conf'].max()))
        print(' ')



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



def plot_S1_ST(S_sens, P_sens, A_sens, R_sens, show=True):
    # Gather the S1 and ST results
    S1 = pd.concat([S_sens['S1'], P_sens['S1'], A_sens['S1'], 
                   R_sens['S1']], keys=['S','P','A','R'], axis=1) #produces copy
    ST = pd.concat([S_sens['ST'], P_sens['ST'], A_sens['ST'], 
                   R_sens['ST']], keys=['S','P','A','R'], axis=1)
    # Change to greek
    for id in S1.index:
        if id != 'mu_star':
            S1.rename(index={id: r'$\{}$'.format(id)}, inplace=True)
        else:
            S1.rename(index={id: r'$\mu^*$'}, inplace=True)
    for id in ST.index:
        if id != 'mu_star':
            ST.rename(index={id: r'$\{}$'.format(id)}, inplace=True)
        else:
            ST.rename(index={id: r'$\mu^*$'}, inplace=True)
    # Plot
    fig, axes = plt.subplots(ncols=2, figsize=(12, 6))
    S1.plot.bar(stacked=True, ax=axes[0], rot=0, width=0.8)
    ST.plot.bar(stacked=True, ax=axes[1], rot=0, width=0.8)
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



def plot_S1_ST_tbl_from_store(store, show=True):
    S_sens = store['S_sens']
    P_sens = store['P_sens']
    A_sens = store['A_sens']
    R_sens = store['R_sens']
    # Gather the S1 and ST results
    S1 = pd.concat([S_sens['S1'], P_sens['S1'], A_sens['S1'], 
                   R_sens['S1']], keys=['S','P','A','R'], axis=1) #produces copy
    ST = pd.concat([S_sens['ST'], P_sens['ST'], A_sens['ST'], 
                   R_sens['ST']], keys=['S','P','A','R'], axis=1)
    # Change to greek
    for id in S1.index:
        if id != 'mu_star':
            S1.rename(index={id: r'$\{}$'.format(id)}, inplace=True)
        else:
            S1.rename(index={id: r'$\mu^*$'}, inplace=True)
    for id in ST.index:
        if id != 'mu_star':
            ST.rename(index={id: r'$\{}$'.format(id)}, inplace=True)
        else:
            ST.rename(index={id: r'$\mu^*$'}, inplace=True)
    # Plot
    fig = plt.figure(figsize=(12, 6))
    gs = gridspec.GridSpec(1, 3, width_ratios=[2.2,2.2,0.8])
    axes = []
    for ii in range(2):
        axes.append(plt.subplot(gs[ii]))
    S1.plot.bar(stacked=True, ax=axes[0], rot=0, width=0.8)
    ST.plot.bar(stacked=True, ax=axes[1], rot=0, width=0.8)
    for ax in axes:
        ax.tick_params(axis='x', labelsize=18)
        ax.tick_params(axis='y', labelsize=14)
        #ax.get_yaxis().set_visible(False)
        ax.legend(fontsize=16)
    axes[0].set_title('First-order indices', fontsize=26)
    axes[1].set_title('Total-order indices', fontsize=26)
    # Create table
    columns = ('Value Range',)
    rows = list(S1.index)
    # alpha, beta, delta, epsilon, gamma, xi, zeta, nu, mu, mu*, sigma
    cell_text = [['0.03-0.3'], ['0.0003-0.03'], ['0.01-1'], ['0.8-8'], ['0.001-0.1'],
                 ['0-1'], ['0.1-2'], ['0.01-1'], ['0.001-0.01'], ['0.005-0.1'], ['0.01-1']]
    # alpha, beta, delta, epsilon, zeta, nu, mu, mu*, sigma
    #cell_text = [['0.03-0.3'], ['0.0003-0.03'], ['0.01-1'], ['0.8-8'],
    #             ['0.1-2'], ['0.01-1'], ['0.001-0.01'], ['0.005-0.1'], ['0.01-1']]
    tbl_ax = plt.subplot(gs[2])
    the_table = tbl_ax.table(cellText=cell_text, rowLabels=rows, colLabels=columns,
                 loc='center')
    the_table.set_fontsize(18)
    the_table.scale(1,2.3)
    the_table.auto_set_column_width(0)
    tbl_ax.axis('off')
    plt.tight_layout()
    if show:
        plt.show()
    else:
        fig.savefig("param_sens_{}.pdf".format(time.strftime("%m_%d_%H%M")))
    return (fig, axes)



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
    red = args.red_model
    if args.ncores is None:
        with Pool() as pool:
            main(args.N, filename=args.filename, reduced=red, pool=pool)
    else:
        with Pool(args.ncores) as pool:
            main(args.N, filename=args.filename, reduced=red, pool=pool)

