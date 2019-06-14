import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode

#initial population values, should add to 1
P_0 = 0.095#0.0835 
A_0 = 0.00647#0.00671
H_0 = 0.000843#0.000874 
R_0 = 0.0584#0.0509
S_0 = 1-P_0-A_0-H_0-R_0

#temporal info, assigning default values
tstart = 0
tstop =  6
#If change tstop and get error, sometimes have to add +1 to part of t linspace like this:
#10*(tstart+tstop+.1)+1

#parameters
params = {}
params['m'] = -0.00483#-0.0156                       #slope of time-dependent alpha: S->P the rate at which people are prescribed opioids 
params['b'] = 0.283#0.303                         #y-intercept of time-dependent alpha: S->P the rate at which people are prescribed opioids
params['beta_A'] = 0.0044#0.00235                  #S->A total probability of becoming addicted to opioids other than by prescription 
params['beta_P'] =  0.000469#0.000141                # S->A proportion of susceptibles that obtain extra prescription opioids OR black market drugs and becomes addicted (Note: MUST BE ZERO FOR AFE)
params['theta_1'] = 0.000502#0.000507                #S->H rate susceptible population becomes addicted to heroin by black market drugs and other addicts 
params['mu'] = 0.00868                      #P,A,H,R->S natural death rate  
params['mu_A'] = 0.00870                    #A->S overdose death rate for opioid addicts 
params['mu_H'] = 0.0507                     #H->S overdose death rate for heroin addicts 
params['gamma'] = 0.00146#0.00115         # P->A rate at which prescribed opioid users become addicted (Note: MUST BE ZERO FOR AFE)  
params['epsilon'] = 2.49#2.54                    #P->S rate at which people come back to the susceptible class after being prescribed opioids (i.e. not addicted)
params['theta_2'] = 0.148#0.0370                  #P->H rate at which opioid prescribed user population becomes addicted to heroin 
params['sigma'] = 0.0283#0.0284                    #R->A rate at which people relapse from treatment into the opioid addicted class 
params['zeta'] = 0.318#0.265                      #A->R rate at which addicted opioid users enter treatment/rehabilitation 
params['theta_3'] = 2.38#3.51                    #A->H rate at which the opioid addicted population becomes addicted to heroin 
params['nu'] = 0.0482#0.00657                      #H->R rate at which heroin users enter treatment/rehabilitation 
params['omega'] = 0.0000000001              #perturbation term for relapse rates
params['c']= -0.0313                        #only for piecewise linear alpha, slope of alpha after Quarter 2 2016 

def alpha(t):
    if t <= 3.25:
       return params['m']*t+params['b']
    else:
       return params['m']*3.25+params['b']-params['c']*3.25+params['c']*t


def update_params(new_params):
    '''Update the params dict with new values contained in new_params'''
    global params
    for key,val in new_params.items():
        if key in params.keys():
            params[key] = val
        else:
            print('Could not find parameter {}.'.format(key))



#think of S, P, A, H, R as coordinates of a vector X
def heroin_odes(t, X, params):
    '''Specifies the heroin model as a system of ODEs.'''
    #Return a new array of given shape and type, without initializing entries.
    Y = np.empty((5)) # S,P,A,H,R
    S = X[0]
    P = X[1]
    A = X[2]
    H = X[3]
    R = X[4]
    
    # For alpha linear
    #Y[0] = -(params['m']*t+params['b'])*S-params['beta_A']*S*A-\
        #params['beta_P']*S*P-params['theta_1']*S*H+\
        #params['epsilon']*P+params['mu']*(P+R)+\
        #(params['mu']+params['mu_A'])*A+(params['mu']+params['mu_H'])*H
   # Y[1] = (params['m']*t+params['b'])*S-params['epsilon']*P-params['gamma']*P-params['theta_2']*P*H-params['mu']*P

   # For alpha piecewise linear
    Y[0] = -alpha(t)*S-params['beta_A']*S*A-\
        params['beta_P']*S*P-params['theta_1']*S*H+\
        params['epsilon']*P+params['mu']*(P+R)+\
        (params['mu']+params['mu_A'])*A+(params['mu']+params['mu_H'])*H
    Y[1] = alpha(t)*S-params['epsilon']*P-params['gamma']*P-params['theta_2']*P*H-params['mu']*P
    Y[2] = params['gamma']*P+params['sigma']*R*A/(A+H+params['omega'])+params['beta_A']*S*A+\
        params['beta_P']*S*P-params['zeta']*A-params['theta_3']*A*H-params['mu']*A-params['mu_A']*A
    Y[3] = params['theta_1']*S*H+params['theta_2']*P*H+params['theta_3']*A*H+\
        params['sigma']*R*H/(A+H+params['omega'])-params['nu']*H-params['mu']*H-params['mu_H']*H
    Y[4] = params['zeta']*A+params['nu']*H-\
        params['sigma']*R*A/(A+H+params['omega'])-params['sigma']*R*H/(A+H+params['omega'])-\
        params['mu']*R
    return Y


#already assigned tstart above, so won't do again here. 
#must do p=None because can never pass a mutable as a default parameter (or else won't be updated when change it)
def solve_odes(S0=S_0,P0=P_0,A0=A_0,H0=H_0,R0=R_0,tstart=tstart,tstop=tstop,p=None):
    '''Solve heroin_odes with given initial conditions, time, and params.'''
    if p is None:
        p = params
    #S, P, A, H, R solution lists
    S = [S0]
    P = [P0]
    A = [A0]
    H = [H0]
    R = [R0]
    # setup solver, dopri15: explicit Runge-Kutta method of order (4)5
    solver = ode(heroin_odes).set_integrator('dopri5')
    #solver = ode(heroin_odes).set_integrator('vode', method='BDF')
    #pass parameters onto function that represents the ODEs
    solver.set_initial_value([S0,P0,A0,H0,R0], tstart).set_f_params(p)
    # solve, solve.t is the current time state of the solver; solver_successful() means no problems while integrating 
    while solver.successful() and solver.t < tstop:
        solver.integrate(solver.t+1) #to integrate up to next integer time, just do solver.t+1, but if we want smoother
        # plots below, so going to solve at more times t so can plot them so do solver.t+.1 for example; solve many times, not just once
        # record solution at that time and append: adds to the list any time the data type is growing and then can convert to data array later 
        #.y are the current solutions states; integrate up to next integer (or next .1) time--only need to record solution at each time
        S.append(solver.y[0])
        P.append(solver.y[1])
        A.append(solver.y[2])
        H.append(solver.y[3])
        R.append(solver.y[4])
    #get list 
    return (S,P,A,H,R)
    #alpha=params['m']*t+params['b']
    #return (alpha) 

# def compute_R0(p=None):
#     '''Check that AFE exists. If so, compute and return R0'''
#     #raise NotImplementedError("Still working on R0 for heroin model!")
#     if p is None:
#          p = params
#     if p['gamma'] != 0 or p['xi'] != 0:
#          raise ValueError('AFE does not exist with these parameters.')
#     else:
#          S_star = (params['epsilon']+params['mu'])/(params['alpha']+params['epsilon']+params['mu'])
#          P_star = params['alpha']/(params['alpha']+params['epsilon']+params['mu'])
#          a= params['zeta']+params['mu']+params['mu_A']
#          b= params['nu']+params['mu']+params['mu_H']
#          c= params['sigma_A']+params['sigma_H']+params['mu']
#          z= params['theta_1']*S_star+params['theta_2']*P_star   
#          r = params['beta']*S_star*(b*c-params['sigma_H']*params['nu'])
#          s = z*(a*c-params['sigma_A']*params['zeta'])
#          detV= a*(b*c-params['sigma_H']*params['nu'])-params['sigma_A']*params['zeta']*b
#          return (((r+s)+((r-s)**(2) + 4*params['beta']*S_star*z*params['sigma_A']*\
#          params['zeta']*params['sigma_H']*params['nu'])**(.5))/(2*detV))
         
    



def plot_solution(S,P,A,H,R,tstart=tstart,tstop=tstop,show=True):
    '''Plot a solution set and either show it or return the plot object'''
    #np.linspace returns evenly spaced values within a given interval; if want more, change solver.integrate(solver.t+.1) above
    # and change t = np.linspace(tstart, tstop+.1, 10*(tstart+tstop+.1) or maybe add +1 at the end here, for example)
    t = np.linspace(tstart, tstop, (tstart+tstop+1))

    fig = plt.figure(figsize=(8, 4.5))
  #  plt.plot(t, S, label='Susceptibles')
  #  plt.plot(t, P, label="Prescription Users")
    plt.plot(t, A, label="Opioid Addicts")
    plt.plot(t, H, label="Heroin and Fentanyl Addicts")
   # plt.plot(t, R, label="Stably Recovered Addicts")
    plt.legend()
    plt.xlabel('Time (years)')
    plt.ylabel('Population proportion')
    #get rid of white space with tight_layout
    plt.tight_layout()
    #in order to get a plot to show 
    if show:
        plt.show()
    else:
        return fig


def plot_addiction_totaled(S,P,A,H,R,tstart=tstart,tstop=tstop,show=True):
    '''Plot a solution set and either show it or return the plot object'''
    #np.linspace returns evenly spaced values within a given interval (0.1 apart in this case)
    t = np.linspace(tstart, tstop, (tstart+tstop+1))
    print(t)
#10*(tstart+tstop+.1)+1
    total = []
    for i in range(len(A)):
            total.append(A[i] + H[i])

    sum_of_values = []
    for i in range(len(S)):
            sum_of_values.append(S[i] + P[i] + A[i] + H[i] + R[i])
    print(sum_of_values)

    fig = plt.figure(figsize=(8, 4.5))
    plt.plot(t, total, label="Total Addicts", color = 'black')
    #plt.plot(t, R, label="Stably Recovered Addicts")
    plt.legend()
    plt.xlabel('Time (years)')
    plt.ylabel('Population proportion')
    #get rid of white space with tight_layout
    plt.tight_layout()
    #in order to get a plot to show 
    if show:
        plt.show()
    else:
        return fig



#comes at end of ODE file (or especially any file you want to run from terminal):
if __name__ == "__main__":

    #run model with default values and plot
    #no parameters because assigned all defaults above; gives tuple 
    sol = solve_odes()
    #need to unpack, use *
    plot_solution(*sol) 
    plot_addiction_totaled(*sol)
    #compute R_0
    #print(compute_R0(p=None))
    print(alpha(0))
    print(alpha(1))
    print(alpha(2))
    print(alpha(3))
    print(alpha(4))
    print(alpha(5))
    print(alpha(6))

    