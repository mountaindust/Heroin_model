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
from matplotlib import gridspec
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
def run_reduced_model(alpha,beta_A,delta,epsilon,zeta,nu,mu,mu_star,sigma):
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



def run_full_model(m,beta_A,beta_P,theta_1,epsilon,gamma,sigma,mu,mu_H,
                   theta_2,zeta,theta_3,nu,omega,b,c,d,e,P_0,A_0,H_0,R_0):
    '''Defines a model wrapper based on the parameter space in main()'''
    # Length to run each model
    tstart = 0
    tstop = 6
    # Copy default parameter dict
    params = dict(heroin_model.params)
    # Replace other parameter values
    params['m'] = m                   
    params['beta_A'] = beta_A
    params['beta_P'] = beta_P
    params['theta_1'] = theta_1
    params['epsilon'] = epsilon
    params['gamma'] = gamma
    params['sigma'] = sigma
    params['mu'] = mu
    params['mu_H'] = mu_H
    params['theta_2'] = theta_2
    params['zeta'] = zeta
    params['theta_3'] = theta_3
    params['nu'] = nu
    params['omega'] = omega 
    params['b'] = b
    params['c'] = c
    params['d'] = d
    params['e'] = e

    # Get initial conditions
    S_0 = 1-P_0-A_0-H_0-R_0 
    #P_0 = 0.0710 
    #A_0 = 0.00760 
    #H_0 = 0.00121 
    #R_0 = 0.000443
    
    
    # Run model
    try:
        result = heroin_model.solve_odes(S_0,P_0,A_0,H_0,R_0,tstart,tstop,params)
    except:
        return (sys.exc_info()[1],None,None,None,None)
    # Return just the sensitivity at the end time for each state 
    return (result[0][-1], result[1][-1],
            result[2][-1], result[3][-1],
            result[4][-1])



def main(N, filename, reduced, pool=None, no_plot=False):
    '''Runs parameter sensitivity on the reduced opioid model'''

    ##### Define the parameter space ######
    if reduced:
        raise NotImplementedError("Reduced model not currently implemented!")
        problem = {
            'num_vars': 9, #number of parameters
            'names': ['alpha', 'beta_A', 'delta', 'epsilon', 'zeta', 'nu',
                        'mu', 'mu_star', 'sigma'],
            'bounds': [[0,1], [0,1], [0,1], [0,1], [0,1], [0,1],
                        [0,0.1], [0,0.5], [0,1]]
        }
    else:
        problem = {
            'num_vars': 22, #number of parameters/initial conditions with bounds +/-50%, alpha=m*t+b for first set of bounds and alpha=piecewise linear for second set of bounds
            'names': ['m', 'beta_A', 'beta_P', 'theta_1', 'epsilon', 
                      'gamma', 'sigma', 'mu', 'mu_H', 'theta_2',
                      'zeta','theta_3', 'nu', 'omega', 
                      'b','c', 'd', 'e', 'P_0', 'A_0', 'H_0', 'R_0'],
            'bounds': [[-0.0075,-0.0035], [0.000439,0.001317], [0.0000327,0.0000981], [0.111,0.333], [1.265,3.795], 
                       [0.002525,0.007575], [0.051,0.153], [0.00355,0.01065], [0.0233,0.0699], [0.118,0.354],
                       [0.099,0.297],  [9.85,29.55], [0.0002655,0.0007965], [0.00000000005,0.00000000015], 
                       [0.2,0.3427], [-0.0315,-0.0235], [0.0004885,0.0014655], [0.004415,0.013245], [0.0475,0.1425], [0.00355,0.01065], [0.0002325,0.0006975], [0.002535,0.007605]]  
        }   # for alpha piecewise linear above

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
        plot_S1_ST_tbl(S_sens, P_sens, A_sens, H_sens, R_sens, False,'eps')



def load_data(filename):
    '''Load analysis data from previous run and return for examination'''

    return pd.HDFStore(filename)



def plot_S1_ST_from_store(store, show=True):
    '''Extract and plot S1 and ST sensitivity data directly from a store object'''

    plot_S1_ST(store['S_sens'], store['P_sens'], store['A_sens'],
               store['H_sens'], store['R_sens'], show)


def plot_S1_ST_tbl_from_store(store, show=True, ext='pdf'):
    '''Extract and plot S1 and ST sensitivity data directly from a store object'''

    plot_S1_ST_tbl(store['S_sens'], store['P_sens'], store['A_sens'],
               store['H_sens'], store['R_sens'], show, ext)


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
        ax.set_ylim(bottom=0)
    axes[0].set_title('First-order indices', fontsize=26)
    axes[1].set_title('Total-order indices', fontsize=26)
    plt.tight_layout()
    if show:
        plt.show()
    else:
        fig.savefig("param_sens_{}.pdf".format(time.strftime("%m_%d_%H%M")))
    return (fig, axes)


def plot_S1_ST_tbl(S_sens, P_sens, A_sens, H_sens, R_sens, show=True, ext='eps'):
    '''This is going to be hard coded for the moment because the saved data
    is not as robust as in the locust Run_Sobol. Thus:
    MAKE SURE THE HARD CODED PARAMETER NAMES/INTERVALS MATCH WHAT IS IN MAIN()!'''

    names = ['m', 'beta_A', 'beta_P', 'theta_1', 'epsilon', 
             'gamma', 'sigma', 'mu', 'mu_H', 'theta_2',
             'zeta','theta_3', 'nu', 'omega', 
             'b','c', 'd', 'e', 'P_0', 'A_0', 'H_0', 'R_0']

    bounds = [[-0.0075,-0.0035], [0.000439,0.001317], [0.0000327,0.0000981], [0.111,0.333], [1.265,3.795], 
              [0.002525,0.007575], [0.051,0.153], [0.00355,0.01065], [0.0233,0.0699], [0.118,0.354],
              [0.099,0.297],  [9.85,29.55], [0.0002655,0.0007965], [0.00000000005,0.00000000015], 
              [0.2,0.3427], [-0.0315,-0.0235], [0.0004885,0.0014655], [0.004415,0.013245], [0.0475,0.1425], [0.00355,0.01065], [0.0002325,0.0006975], [0.002535,0.007605]]

    # Gather the S1 and ST results
    S1 = pd.concat([S_sens['S1'], P_sens['S1'], A_sens['S1'], H_sens['S1'],
                   R_sens['S1']], keys=['S','P','A','H','R'], axis=1)
    ST = pd.concat([S_sens['ST'], P_sens['ST'], A_sens['ST'], H_sens['ST'],
                   R_sens['ST']], keys=['S','P','A','H','R'], axis=1)

    ###### Change to greek, LaTeX #####
    for id in S1.index:
        # handle greek ratios
        if id not in ['m', 'b', 'c', 'd', 'e', 'P_0', 'A_0', 'H_0', 'R_0']:
            S1.rename(index={id: r'$\{}$'.format(id)}, inplace=True)
        if id in ['m', 'b', 'c', 'd', 'e']:
            S1.rename(index={id: r'$\widetilde{{{}}}$'.format(id)}, inplace=True)
        if id in ['P_0', 'A_0', 'H_0', 'R_0']:
            S1.rename(index={id: id[0]+r'$_0$'}, inplace=True)

    for id in ST.index:
        if id not in ['m', 'b', 'c', 'd', 'e', 'P_0', 'A_0', 'H_0', 'R_0']:
            ST.rename(index={id: r'$\{}$'.format(id)}, inplace=True)
        if id in ['m', 'b', 'c', 'd', 'e']:
            ST.rename(index={id: r'$\widetilde{{{}}}$'.format(id)}, inplace=True)
        if id in ['P_0', 'A_0', 'H_0', 'R_0']:
            ST.rename(index={id: id[0]+r'$_0$'}, inplace=True)

    # setup for table
    fig = plt.figure(figsize=(13, 6))
    gs = gridspec.GridSpec(1, 3, width_ratios=[3,3,1], wspace=.15, left=0.04,
                            right=0.995, bottom=0.15, top=0.915)
    axes = []
    for ii in range(2):
        axes.append(plt.subplot(gs[ii]))

    s1bars = S1.plot.bar(stacked=True, ax=axes[0], rot=0, width=0.8)
    s2bars = ST.plot.bar(stacked=True, ax=axes[1], rot=0, width=0.8, legend=False)
    for ax in axes:
        ax.tick_params(axis='x', labelsize=14, rotation=0)
        ax.tick_params(axis='y', labelsize=14)
        #ax.get_yaxis().set_visible(False)
        ax.set_ylim(bottom=0)
    axes[0].set_title('First-order indices', fontsize=26)
    axes[1].set_title('Total-order indices', fontsize=26)
    handles, labels = s1bars.get_legend_handles_labels()
    s1bars.legend(reversed(handles), reversed(labels), loc='upper left', fontsize=16)

    # Create table
    columns = ('Value Range',)
    rows = list(S1.index)
    # turn bounds into strings of ranges
    cell_text = []
    for bnd in bounds:
        low = str(bnd[0])
        high = str(bnd[1])
        # concatenate, remove leading zeros
        if low != "0" and low != "0.0":
            low = low.lstrip("0")
        if high != "0" and high != "0.0":
            high = high.lstrip("0")
        # raise any minus signs
        low = low.replace('-','\u00AF')
        high = high.replace('-','\u00AF')
        cell_text.append([low+"-"+high])
    tbl_ax = plt.subplot(gs[2])
    the_table = tbl_ax.table(cellText=cell_text, rowLabels=rows, colLabels=columns,
                loc='left')
    the_table.set_fontsize(10)
    the_table.scale(3.1,1.3)
    the_table.auto_set_column_width(0)
    tbl_ax.axis('off')

    # reposition table
    pos = tbl_ax.get_position()
    newpos = [pos.x0 + 0.11, pos.y0, pos.width, pos.height]
    tbl_ax.set_position(newpos)
    if show:
        plt.show()
    else:
        fig.savefig("param_sens_{}.{}".format(time.strftime("%m_%d_%H%M"), ext))
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
### This is also stale. Could be adapted to compare two different intervals ###                   
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

