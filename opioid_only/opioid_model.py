import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode

#initial population values
S_0 = 0.9435
P_0 = 0.05 #Study: Boudreau et al +  time increase??
A_0 = 0.0062 #SAMHSA: https://www.samhsa.gov/data/sites/default/files/NSDUH-FFR2-2015/NSDUH-FFR2-2015.pdf
R_0 = 0.0003 #HSS Treatment Episode Data Set https://www.samhsa.gov/data/sites/default/files/2014_Treatment_Episode_Data_Set_National_Admissions_9_19_16.pdf

#temporal info
tstart = 0
tstop = 10

#parameters
params = {}
params['alpha'] = 0.15                  #S->P: prescription rate
params['beta_P'] = 0.00266              #S->A due to P
params['beta_A'] = 0.00094              #S->A due to A
params['delta'] = 0.1                   #R->S: finish recovery
params['epsilon'] = 3.0                 #P->S rate
params['gamma'] = 0.00744               #P->A
params['zeta'] = 0.25                   #A->R rate of starting treatment
params['nu_1'] = 0                      #RP relapse rate
params['nu_2'] = 0                      #RA relapse rate
params['mu'] = 0.007288                 #nomral death rate
params['mu_star'] = 0.01155             #addiction death rate
params['sigma'] = 0.7                   #R->A natural treatment relapse

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
    Y[0] = -params['alpha']*S-params['beta_A']**S*A-\
        params['beta_P']*S*P+params['epsilon']*P+\
        params['delta']*R+params['mu']*(P+R)+params['mu_star']*A
    Y[1] = params['alpha']*S-(params['epsilon']+params['gamma']+params['mu'])*P
    Y[2] = params['gamma']*P+params['sigma']*R+params['beta_A']*S*A+\
        params['beta_P']*S*P-(params['zeta']+params['mu_star'])*A+\
        params['nu_1']*R*P+params['nu_2']*R*A
    Y[3] = params['zeta']*A-(params['delta']+params['sigma']+params['mu'])*R-\
        params['nu_1']*R*P-params['nu_2']*R*A
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
    if p['gamma'] != 0 or p['beta_P'] != 0:
        raise ValueError('AFE does not exist with these parameters.')
    else:
        Lambda = 1 - p['sigma']/(p['delta']+p['mu']+p['sigma'])
        S_star = 1 - p['alpha']/(p['alpha']+p['epsilon']+p['mu'])
        return (p['beta_A']*S_star)/(p['mu']+p['zeta']*Lambda)



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
