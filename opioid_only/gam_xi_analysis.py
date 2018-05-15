'''
Examine the effect of moving gamma and xi away from zero when R0<1.
'''

import sys, os
import pickle
import argparse
from multiprocessing import Pool
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import opioid_model as model

default_N = os.cpu_count()
parser = argparse.ArgumentParser()
parser.add_argument("-n", "--ncores", type=int,
                    help="number of cores, defaults to {}".format(default_N))
parser.add_argument("-o", "--filename", type=str, 
                    help="filename to write output to, no extension",
                    default='gam_xi_analysis')

S_0 = 0.9435
P_0 = 0.05
A_0 = 0.0062
R_0 = 0.0003

tstart = 0
tstop = 10000

model.params['gamma'] = 0
model.params['xi'] = 0

# define the gamma/xi space to explore
gamma_step = 0.0005
gamma_end = 0.05
xi_step = 0.005
xi_end = 0.4



def run_model(gamma, xi):
    '''wrapper around the opioid model solver'''
    # Copy default parameter dict
    params = dict(model.params)
    params['gamma'] = gamma
    params['xi'] = xi
    # Run model
    try:
        result = model.solve_odes(S_0,P_0,A_0,R_0,tstart,tstop,params)
    except:
        return (sys.exc_info()[1],None,None,None)
    # Return just the end value of each variable
    return (result[0][-1],result[1][-1],result[2][-1],result[3][-1])



def main(filename, pool=None):
    # check R0 value of gamma=xi=0 and print to screen
    print("R0 = {}".format(model.compute_R0()))

    # create array of parameter combinations
    param_values = []
    for xi in np.arange(0,xi_end+xi_step,xi_step):
        for gam in np.arange(0,gamma_end+gamma_step,gamma_step):
            param_values.append([gam,xi])
    param_values = np.array(param_values)

    print('Running model...')
    if args.ncores is None:
        poolsize = os.cpu_count()
    else:
        poolsize = args.ncores
    chunksize = param_values.shape[0]//poolsize
    output = pool.starmap(run_model, param_values, chunksize=chunksize)

    ### Parse and save the output ###
    print('Saving and reviewing the results...')
    param_values = pd.DataFrame(param_values, columns=['gamma','xi'])
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
    store['output'] = pd.DataFrame(output, columns=['S', 'P', 'A', 'R'])
    os.remove('raw_result_data.pickle')
    store.close()

    ### Plot the parameter space ###
    plot_param_space(output[:,2], output[:,3], gamma_end, gamma_step, 
                     xi_end, xi_step)



def load_data(filename):
    '''Load analysis data from previous run and return for examination'''

    return pd.HDFStore(filename)



def plot_data():
    store = load_data("gam_xi_alpha02.h5")
    A = store['output']['A'].as_matrix()
    R = store['output']['R'].as_matrix()
    plot_param_space(A, R, gamma_end, gamma_step, xi_end, xi_step)


def plot_param_space(A, R, gam_end, gam_step, xi_end, xi_step, show=True):
    '''Plot the parameter space of R'''
    plt.rc('font',family='Arial')
    # mesh the parameter values
    gam_mesh = np.arange(0, gam_end+gam_step, gam_step)
    gam_len = gam_mesh.size
    xi_mesh = np.arange(0, xi_end+xi_step, xi_step)
    xi_len = xi_mesh.size
    gam_grid, xi_grid = np.meshgrid(gam_mesh, xi_mesh)

    fig, axes = plt.subplots(ncols=2, figsize=(11, 4.5))
    A_plot = axes[0].pcolormesh(gam_grid, xi_grid, A.reshape(xi_len, gam_len))
    A_plot.set_edgecolor('face')
    axes[0].set_title('Addicted @ equilibrium',fontsize=19)
    axes[0].set_xlabel(r'$\gamma$', fontsize=18)
    axes[0].set_ylabel(r'$\xi$', fontsize=18)
    axes[0].tick_params(labelsize=16)
    cbar = plt.colorbar(A_plot, ax=axes[0])
    cbar.ax.tick_params(labelsize=16)
    R_plot = axes[1].pcolormesh(gam_grid, xi_grid, R.reshape(xi_len, gam_len))
    R_plot.set_edgecolor('face')
    axes[1].set_title('Rehabilitating @ equilibrium',fontsize=19)
    axes[1].set_xlabel(r'$\gamma$', fontsize=18)
    axes[1].set_ylabel(r'$\xi$', fontsize=18)
    axes[1].tick_params(labelsize=16)
    cbar = plt.colorbar(R_plot, ax=axes[1])
    cbar.ax.tick_params(labelsize=16)
    plt.tight_layout()
    if show:
        plt.show()
    return (fig, axes)



if __name__ == "__main__":
    args = parser.parse_args()
    if args.ncores is None:
        with Pool() as pool:
            main(filename=args.filename, pool=pool)
    else:
        with Pool(args.ncores) as pool:
            main(filename=args.filename, pool=pool)

