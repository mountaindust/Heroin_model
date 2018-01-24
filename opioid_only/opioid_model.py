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
params['epsilon'] = 0.74                #P->S rate
params['gamma'] = 1 - params['epsilon'] #P->A
params['xi'] = 0.505                    #fraction of beta due to P
params['zeta'] = 0.25                   #A->R rate of starting treatment
params['nu'] = 0.293*(1-params['delta'])#R->A treatment relapse due to A
params['mu'] = 0.008237                 #nomral death rate
params['mu_star'] = 0.00834             #addiction death rate
params['sigma'] = 0.707*(1-params['delta'])  #R->A natural treatment relapse

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
    Y = np.empty((4)) # S,P,A,R
    S = X[0]
    P = X[1]
    A = X[2]
    R = X[3]
    Y[0] = -params['alpha']*S-params['beta']*(1-params['xi'])*S*A-\
        params['beta']*params['xi']*S*P+params['epsilon']*P+\
        params['delta']*R+params['mu']*(P+R)+params['mu_star']*A
    Y[1] = params['alpha']*S-(params['epsilon']+params['gamma']+params['mu'])*P
    Y[2] = params['gamma']*P+params['sigma']*R+params['beta']*(1-params['xi'])*S*A+\
        params['beta']*params['xi']*S*P+params['nu']*R*A-\
        (params['zeta']+params['mu_star'])*A
    Y[3] = params['zeta']*A-params['nu']*R*A-\
        (params['delta']+params['sigma']+params['mu'])*R
    return Y



def solve_odes(S0=S_0,P0=P_0,A0=A_0,R0=R_0,tstart=tstart,tstop=tstop,p=None):
    '''Solve opioid_odes with given initial conditions, time, and params.'''
    if p is None:
        p = params
    S = [S0]
    P = [P0]
    A = [A0]
    R = [R0]
    # setup solver
    solver = ode(opioid_odes).set_integrator('dopri5')
    solver.set_initial_value([S0,P0,A0,R0], tstart).set_f_params(p)
    # solve
    while solver.successful() and solver.t < tstop:
        solver.integrate(solver.t+1) #integrate up to next integer time
        # record solution at that time
        S.append(solver.y[0])
        P.append(solver.y[1])
        A.append(solver.y[2])
        R.append(solver.y[3])
    return (S,P,A,R)



def compute_R0(p=None):
    '''Check that AFE exists. If so, compute and return R0'''
    if p is None:
        p = params
    if p['gamma'] != 0 or p['xi'] != 0:
        raise ValueError('AFE does not exist with these parameters.')
    else:
        Lambda = 1 - p['sigma']/(p['delta']+p['mu']+p['sigma'])
        S_star = 1 - p['alpha']/(p['alpha']+p['epsilon']+p['mu'])
        return (p['beta']*S_star)/(p['mu']+p['zeta']*Lambda)



def plot_solution(S,P,A,R,tstart=tstart,tstop=tstop,show=True):
    '''Plot a solution set and either show it or return the plot object'''
    t = np.arange(tstart,tstop+1,1)

    fig = plt.figure(figsize=(8, 4.5))
    plt.plot(t, S, label='Susceptible')
    plt.plot(t, P, label="Prescription Users")
    plt.plot(t, A, label="Addicts")
    plt.plot(t, R, label="Recovering Addicts")
    plt.legend()
    plt.xlabel('Time (years)')
    plt.ylabel('Population fraction')
    plt.tight_layout()
    if show:
        plt.show()
    else:
        return fig

if __name__ == "__main__":

    #run model with default values and plot

    sol = solve_odes()
    plot_solution(*sol)