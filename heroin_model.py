import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode

#initial population values, should add to 1
S_0 = 0.6068 # 1 - proportions in P_0, A_0, H_0, R_0
P_0 = 0.38 # [4] Study: http://annals.org/aim/fullarticle/2646632/prescription-opioid-use-misuse-use-disorders-u-s-adults-2015 (page 293)
A_0 = 0.0062 # [21] https://www.samhsa.gov/data/sites/default/files/NSDUH-FFR1-2015/NSDUH-FFR1-2015/NSDUH-FFR1-2015.pdf (page 25)
H_0 = 0.0026  # [21] SAMSHA: https://www.samhsa.gov/data/sites/default/files/NSDUH-FFR1-2015/NSDUH-FFR1-2015/NSDUH-FFR1-2015.pdf (page 11)
R_0 = 0.0044 # [14] https://d14rmgtrwzf5a.cloudfront.net/sites/default/files/19774-prescription-opioids-and-heroin.pdf (page 17) #((772000+618000)/total population in 2014 = 319,000,000)

#temporal info, assigning default values
tstart = 0
tstop = 10

#parameters
params = {}
params['alpha'] = 0.207                       #S->P the rate at which people are prescribed opioids #[20] https://www.cdc.gov/drugoverdose/pdf/pubs/2017-cdc-drug-surveillance-report.pdf (page 39, table 1A)
params['beta'] = 0.006                     #S->A total probability of becoming addicted to opioids other than by prescription #assume okay for now until find updated/source
params['xi'] =  0.505                      # S->A proportion of susceptibles that obtain extra prescription opioids OR black market drugs and becomes addicted (Note: MUST BE ZERO FOR AFE) #assume okay for now until find updated/source
params['theta_1'] = 0.0003                  #S->H rate susceptible population becomes addicted to heroin by black market drugs and other addicts #ESTIMATED [30] https://www.drugabuse.gov/publications/drugfacts/heroin#ref
params['delta'] = 0.15                     #R->S rate at which people come back to the susceptible class after successfully finishing treatment #[19] https://ac.els-cdn.com/S0740547213000779/1-s2.0-S0740547213000779-main.pdf?_tid=0a47e661-2ac3-45fb-a7da-9fa8c43271af&acdnat=1525189437_134c2f5932fa797a0725faa7c950a0f9 (page 1, 10% successfully treated for opioids + ESTIMATED 5% successfully treated for heroin)
params['mu'] = 0.00844                      #P,A,H,R->S natural death rate  # [24] https://www.cdc.gov/nchs/data/nvsr/nvsr66/nvsr66_05.pdf (page 1)
params['mu_A'] = 0.00004                  #A->S enhanced death rate for opioid addicts (only overdose rate=4/100,000) # [25] https://www.cdc.gov/drugoverdose/data/analysis.html (Opioid Data Analysis, figure)
params['mu_H'] = 0.00004                 #H->S enhanced death rate for heroin addicts (only overdose rate=4/100,000) # [25] https://www.cdc.gov/drugoverdose/data/analysis.html  (Opioid Data Analysis, figure)
params['gamma'] = 0.26          # P->A rate at which prescribed opioid users become addicted (Note: MUST BE ZERO FOR AFE) # [28] https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1360-0443.2010.03052.x
params['epsilon'] = 1-params['gamma']         #P->S rate at which people come back to the susceptible class after being prescribed opioids (i.e. not addicted) #[28] https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1360-0443.2010.03052.x (1-0.26=0.74)
params['theta_2'] = 0.0003                       #P->H rate at which opioid prescribed user population becomes addicted to heroin #[30] https://www.drugabuse.gov/publications/drugfacts/heroin#ref
params['sigma_A'] = 0.9  #R->A rate at which people relapse from treatment into the opioid addicted class #[19] https://ac.els-cdn.com/S0740547213000779/1-s2.0-S0740547213000779-main.pdf?_tid=0a47e661-2ac3-45fb-a7da-9fa8c43271af&acdnat=1525189437_134c2f5932fa797a0725faa7c950a0f9 (page 1, 90% relapse rate by end of one year)
params['zeta'] = 0.25                        #A->R rate at which addicted opioid users enter treatment/rehabilitation #assume okay for now until find updated/source
params['theta_3'] = 0.04                    #A->H rate at which the opioid addicted population becomes addicted to heroin #[14] https://d14rmgtrwzf5a.cloudfront.net/sites/default/files/19774-prescription-opioids-and-heroin.pdf (page 7)
params['sigma_H'] = 0.95    #R->H rate at which people relapse from treatment back into the heroin addicted class #ESTIMATED for now because cannot find source yet
params['nu'] = .15                          #H->R rate at which heroin users enter treatment/rehabilitation #[14] https://d14rmgtrwzf5a.cloudfront.net/sites/default/files/19774-prescription-opioids-and-heroin.pdf (page 17)


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

    Y[0] = -params['alpha']*S-params['beta']*(1-params['xi'])*S*A-\
        params['beta']*params['xi']*S*P-params['theta_1']*S*H+\
        params['epsilon']*P+params['delta']*R+params['mu']*(P+R)+\
        (params['mu']+params['mu_A'])*A+(params['mu']+params['mu_H'])*H
    Y[1] = params['alpha']*S-params['epsilon']*P-params['gamma']*P-params['theta_2']*P*H-params['mu']*P
    Y[2] = params['gamma']*P+params['sigma_A']*R+params['beta']*(1-params['xi'])*S*A+\
        params['beta']*params['xi']*S*P-params['zeta']*A-params['theta_3']*A*H-params['mu']*A-params['mu_A']*A
    Y[3] = params['theta_1']*S*H+params['theta_2']*P*H+params['theta_3']*A*H+params['sigma_H']*R-\
        params['nu']*H-params['mu']*H-params['mu_H']*H
    Y[4] = params['zeta']*A+params['nu']*H-\
        (params['delta']+params['sigma_A']+params['sigma_H']+params['mu'])*R
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
    #pass parameters onto function that represents the ODEs
    solver.set_initial_value([S0,P0,A0,H0,R0], tstart).set_f_params(p)
    # solve, solve.t is the current time state of the solver; solver_successful() means no problems while integrating 
    while solver.successful() and solver.t < tstop:
        solver.integrate(solver.t+1) #integrate up to next integer time; solve many times, not just once
        # record solution at that time and append: adds to the list any time the data type is growing and then can convert to data array later 
        #.y are the current solutions states; integrate up to next integer time--only need to record solution at each time
        S.append(solver.y[0])
        P.append(solver.y[1])
        A.append(solver.y[2])
        H.append(solver.y[3])
        R.append(solver.y[4])
    #get list 
    return (S,P,A,H,R)



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
#          c= params['delta']+params['sigma_A']+params['sigma_H']+params['mu']
#          z= params['theta_1']*S_star+params['theta_2']*P_star   
#          r = params['beta']*S_star*(b*c-params['sigma_H']*params['nu'])
#          s = z*(a*c-params['sigma_A']*params['zeta'])
#          detV= a*(b*c-params['sigma_H']*params['nu'])-params['sigma_A']*params['zeta']*b
#          return (((r+s)+((r-s)**(2) + 4*params['beta']*S_star*z*params['sigma_A']*\
#          params['zeta']*params['sigma_H']*params['nu'])**(.5))/(2*detV))
         

def plot_solution(S,P,A,H,R,tstart=tstart,tstop=tstop,show=True):
    '''Plot a solution set and either show it or return the plot object'''
    #np.arangereturn evenly spaced values within a given interval 
    t = np.arange(tstart,tstop+1,1)

    fig = plt.figure(figsize=(8, 4.5))
    plt.plot(t, S, label='Susceptible')
    plt.plot(t, P, label="Prescription Users")
    plt.plot(t, A, label="Addicts")
    plt.plot(t, H, label="Heroin Addicts")
    plt.plot(t, R, label="Recovering Addicts")
    plt.legend()
    plt.xlabel('Time (years)')
    plt.ylabel('Population fraction')
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

    #no paramaters because assigned all defaults above; gives tuple 
    sol = solve_odes()
    #need to unpack, use *
    plot_solution(*sol) 
    #compute R_0
   # print(compute_R0(p=None))
