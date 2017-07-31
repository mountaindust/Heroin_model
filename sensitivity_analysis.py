'''Use SALib to perform Sobol sensitivity analysis on the model parameters'''

import sys, os, pickle
from multiprocessing import Pool
from SALib.sample import saltelli
from SALib.analyze import sobol
import numpy as np
import heroin_model

### Define a model wrapper based on the parameter space in main() ###
def run_model(alpha, beta, delta1, delta2, mu1, mu2, irate, hrate,
                epsilon, xi, sigma, gamma_0, dr_add):
    # Length to run each model realization
    endt = 200
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
    params['dr'] = params['de'] + dr_add
    # Get initial conditions
    s0 = heroin_model.s0
    e0 = heroin_model.e0
    i0 = heroin_model.i0
    h0 = heroin_model.h0
    r0 = heroin_model.r0
    # Run model
    try:
        result = heroin_model.heroin_func_ode(s0,e0,i0,h0,r0,endt,params)
    except:
        return (sys.exc_info()[1], 
                None, None, None, None, None, None, None, None, None)
    # Return just the end value and std of each variable
    return (np.mean(result[0][-1]), np.mean(result[1][-1]), 
            np.mean(result[2][-1]), np.mean(result[3][-1]), 
            np.mean(result[4][-1]), np.std(result[0][-50:]),
            np.std(result[1][-50:]), np.std(result[2][-50:]),
            np.std(result[3][-50:]), np.std(result[4][-50:]))



def main(N, pool=None):
    '''Uses the ODE version of the heroin model and conducts parameter sensitivity'''

    ### Define the parameter space ###

    problem = {
        'num_vars': 13, #number of parameters
        'names': ['alpha', 'beta', 'delta1', 'delta2',
                  'mu1', 'mu2', 'irate', 'hrate', 'epsilon',
                  'xi', 'sigma', 'gamma_0', 'dr_add'],
        'bounds': [[0.01, 0.5], [0.01, 0.5], [0, 0.2], [0, 0.8],
                   [0, 0.15], [0, 0.3], [0, 1], [0, 1], [0.01, 0.8],
                   [0, 1], [0, 0.95], [0.1, 10], [0.00001, 0.95]]
    }

    ### Create an N by num_var matrix of parameter values ###
    param_values = saltelli.sample(problem, N, calc_second_order=True)

    ### Run model ###
    print('Examining the parameter space.')
    poolsize = os.cpu_count()
    chunksize = param_values.shape[0]//poolsize
    output = pool.starmap(run_model, param_values, chunksize=chunksize)

    ### Parse and save the output ###
    print('Saving the results...')
    with open("result_data.pickle", "wb") as f:
        pickle.dump(output, f)
    with open("param_values.pickle", "wb") as f:
        pickle.dump(param_values, f)
    print('Reviewing the results...')
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
    std_sum = output[:,5:].sum(1)

    ### Analyze the results ###
    print("S sensitivity")
    print("-------------")
    S_sens = sobol.analyze(problem, output[:,0], calc_second_order=True)
    sobol.print_indices(S_sens, problem, False)
    print("E sensitivity:")
    print("-------------")
    E_sens = sobol.analyze(problem, output[:,1], calc_second_order=True)
    sobol.print_indices(E_sens, problem, False)
    print("I sensitivity:")
    print("-------------")
    I_sens = sobol.analyze(problem, output[:,2], calc_second_order=True)
    sobol.print_indices(I_sens, problem, False)
    print("H sensitivity:")
    print("-------------")
    H_sens = sobol.analyze(problem, output[:,3], calc_second_order=True)
    sobol.print_indices(H_sens, problem, False)
    print("R sensitivity:")
    print("-------------")
    R_sens = sobol.analyze(problem, output[:,4], calc_second_order=True)
    sobol.print_indices(R_sens, problem, False)
    print("Var sensitivity:")
    print("-------------")
    var_sens = sobol.analyze(problem, std_sum, calc_second_order=True)
    sobol.print_indices(var_sens, problem, False)
    
    ### Save the analysis ###
    with open("analysis.pickle", "wb") as f:
        pickle.dump([S_sens, E_sens, I_sens, H_sens, R_sens, var_sens, problem], f)
    print("Analysis (including second order) saved as analysis.pickle.")



if __name__ == "__main__":
    N = 100
    with Pool() as pool:
        main(N, pool)