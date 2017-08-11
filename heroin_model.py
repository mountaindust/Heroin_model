import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode

CHILD = 48600000 # number of children to remove from US population

#initial population values
##make up h and i populations pretend its the beginning of epiDEmic 
TOTAL01 = 317000000 - CHILD
e01 = 10000000
i01 = 400000
h01 = 100000
r01 = 60000000
s01 = TOTAL01 - (e01 + i01 + h01 + r01)

#inital pop. as fractions
s0 = s01/TOTAL01
e0 = e01/TOTAL01
i0 = i01/TOTAL01
h0 = h01/TOTAL01
r0 = r01/TOTAL01

#parameters:
params = {}
params['alpha'] = .05 #S->E
params['beta'] = .1 #E->S
#delta1&2 exist in [0,1) d2>d1
params['delta1'] = .05 #S->I due to E
params['delta2'] = .25 #S->I due to I
# Heroin adoption rates. mu2>mu1. Also affects I->H through gamma.
params['mu1'] = .00075 #S+E->H linear
params['mu2'] = .1 #Literature? S+E->H due to H

# rehab rates
params['irate'] = .25 #rate from I into R
params['hrate'] = .1 #rate from H into R

#1/3 of people prescribed opioids get hooked/DEpenDEnt on them
params['epsilon'] = .333

params['xi'] = .02 #S+E->R due to H
params['sigma'] = .01 #S+E->R background

#4/5 of new heroin users come from prescription opioids 
#gamma/(mu1+mu2*H_n+gamma) = .8
params['gamma_0'] = 4 #gamma = GAMMA_0*(mu1+mu2*H_n) I->H

#recetivism
params['r'] = 1/1000

#Death rates
params['ds'] = 12305000/1245911000 #.009876 Literature
params['di'] = 6.3/100000 * (TOTAL01/i01) # Literature
params['dh'] = 2.7/100000 * (TOTAL01/h01) # Literature

endt = 200

def heroin_dtds(S,E,I,H,R,tsteps=endt,params=params):
    '''Append arrays S,E,I,H,R with year time-step solutions'''
    mu1 = params['mu1']
    mu2 = params['mu2']
    alpha = params['alpha']
    beta = params['beta']
    sigma = params['sigma']
    delta1 = params['delta1']
    delta2 = params['delta2']
    zeta = params['xi']
    EPSILON = params['epsilon']
    hrate = params['hrate']
    irate = params['irate']
    DS = params['ds']
    DI = params['di']
    DH = params['dh']
    r = params['r']
    for i in range(tsteps):
        S_n = S[-1]
        E_n = E[-1]
        I_n = I[-1]
        H_n = H[-1]
        R_n = R[-1]
        GAMMA = params['gamma_0'] *(mu1+mu2*H_n)
        S.append(S_n + (beta*E_n) - (alpha*S_n) - S_n*(delta1*E_n + delta2*I_n) - S_n*(sigma+zeta*(I_n+H_n)) - S_n*(mu1+mu2*H_n) + DI*I_n + DH*H_n + DS*E_n + (DS+r)*R_n)
        E.append(E_n - (beta*E_n) + (alpha*S_n) - E_n*(sigma+zeta*(I_n+H_n)) - E_n*(mu1+mu2*H_n) - (EPSILON*E_n) - DS*E_n)
        I.append(I_n + (EPSILON*E_n) - (GAMMA*I_n) + S_n*(delta1*E_n + delta2*I_n) - irate*I_n - DI*I_n)
        H.append(H_n + (GAMMA*I_n) + (S_n+E_n)*(mu1+mu2*H_n) - hrate*H_n -DH*H_n)
        R.append(R_n + ((E_n+S_n)*(sigma + zeta*(I_n+H_n))) + irate*I_n + hrate*H_n - (DS+r)*R_n)
        
    assert np.isclose(1, S[-1]+E[-1]+I[-1]+H[-1]+R[-1]), "Final sum not equal to one."
    assert min([min(S), min(E), min(I), min(H), min(R)]) >= 0, "Negative population fraction."



def heroin_solve_odes(S0,E0,I0,H0,R0,version=1,tsteps=endt,params=params):
    '''Solve the heroin model as a system of ODEs'''
    S = [S0]
    E = [E0]
    I = [I0]
    H = [H0]
    R = [R0]
    if version == 1:
        solver = ode(heroin_odes).set_integrator('dopri5')
        solver.set_initial_value([S0,E0,I0,H0,R0], 0).set_f_params(params)
    elif version == 2:
        solver = ode(heroin_odes_v2).set_integrator('dopri5')
        solver.set_initial_value([S0,E0,I0,H0,R0], 0).set_f_params(params)
    elif version == 3:
        solver = ode(heroin_odes_noR).set_integrator('dopri5')
        solver.set_initial_value([S0,E0,I0,H0], 0).set_f_params(params)
        R = None
    while solver.successful() and solver.t < endt:
        solver.integrate(solver.t+1)
        S.append(solver.y[0])
        E.append(solver.y[1])
        I.append(solver.y[2])
        H.append(solver.y[3])
        if version != 3:
            R.append(solver.y[4])
    return (S,E,I,H,R)



def heroin_odes(t, X, params):
    '''Specify the heroin model as a system of ODEs. version=1'''
    Y = []
    params['gamma'] = params['gamma_0'] *(params['mu1']+params['mu2']*X[3])
    Y.append((params['beta']*X[1]) - (params['alpha']*X[0]) - 
        X[0]*(params['delta1']*X[1] + params['delta2']*X[2]) - 
        X[0]*(params['sigma']+params['xi']*(X[2]+X[3])) - 
        X[0]*(params['mu1']+params['mu2']*X[3]) + params['di']*X[2] + 
        params['dh']*X[3] + params['ds']*X[1] + params['ds']*X[4] + params['r']*X[4])
    Y.append(-(params['beta']*X[1]) + (params['alpha']*X[0]) - 
        X[1]*(params['sigma']+params['xi']*(X[2]+X[3])) - 
        X[1]*(params['mu1']+params['mu2']*X[3]) - (params['epsilon']*X[1]) - 
        params['ds']*X[1])
    Y.append((params['epsilon']*X[1]) - (params['gamma']*X[2]) + 
        X[0]*(params['delta1']*X[1] + params['delta2']*X[2]) - 
        params['irate']*X[2] - params['di']*X[2])
    Y.append((params['gamma']*X[2]) + 
        (X[0]+X[1])*(params['mu1']+params['mu2']*X[3]) - 
        params['hrate']*X[3] - params['dh']*X[3])
    Y.append(((X[1]+X[0])*(params['sigma'] + params['xi']*(X[2]+X[3]))) + 
        params['irate']*X[2] + params['hrate']*X[3] - params['ds']*X[4] - params['r']*X[4])
    return Y



def heroin_odes_v2(t, X, params):
    '''Send I and H to S instead of R. version=2'''
    Y = []
    params['gamma'] = params['gamma_0'] *(params['mu1']+params['mu2']*X[3])
    Y.append((params['beta']*X[1]) - (params['alpha']*X[0]) - 
        X[0]*(params['delta1']*X[1] + params['delta2']*X[2]) - 
        X[0]*(params['sigma']+params['xi']*(X[2]+X[3])) - 
        X[0]*(params['mu1']+params['mu2']*X[3]) + 
        params['irate']*X[2] + params['hrate']*X[3] + params['di']*X[2] + 
        params['dh']*X[3] + params['ds']*X[1] + params['ds']*X[4] + params['r']*X[4])
    Y.append(-(params['beta']*X[1]) + (params['alpha']*X[0]) - 
        X[1]*(params['sigma']+params['xi']*(X[2]+X[3])) - 
        X[1]*(params['mu1']+params['mu2']*X[3]) - (params['epsilon']*X[1]) - 
        params['ds']*X[1])
    Y.append((params['epsilon']*X[1]) - (params['gamma']*X[2]) + 
        X[0]*(params['delta1']*X[1] + params['delta2']*X[2]) - 
        params['irate']*X[2] - params['di']*X[2])
    Y.append((params['gamma']*X[2]) + 
        (X[0]+X[1])*(params['mu1']+params['mu2']*X[3]) - 
        params['hrate']*X[3] - params['dh']*X[3])
    Y.append(((X[1]+X[0])*(params['sigma'] + params['xi']*(X[2]+X[3])))
        - params['ds']*X[4] - params['r']*X[4])
    return Y



def heroin_odes_noR(t, X, params):
    '''Same as v2 but with no R. version=3'''
    Y = []
    params['gamma'] = params['gamma_0'] *(params['mu1']+params['mu2']*X[3])
    Y.append((params['beta']*X[1]) - (params['alpha']*X[0]) - 
        X[0]*(params['delta1']*X[1] + params['delta2']*X[2]) -  
        X[0]*(params['mu1']+params['mu2']*X[3]) + 
        params['irate']*X[2] + params['hrate']*X[3] + 
        params['di']*X[2] + params['dh']*X[3] + params['ds']*X[1])
    Y.append(-(params['beta']*X[1]) + (params['alpha']*X[0]) -  
        X[1]*(params['mu1']+params['mu2']*X[3]) - (params['epsilon']*X[1]) - 
        params['ds']*X[1])
    Y.append((params['epsilon']*X[1]) - (params['gamma']*X[2]) + 
        X[0]*(params['delta1']*X[1] + params['delta2']*X[2]) - 
        params['irate']*X[2] - params['di']*X[2])
    Y.append((params['gamma']*X[2]) + 
        (X[0]+X[1])*(params['mu1']+params['mu2']*X[3]) - 
        params['hrate']*X[3] - params['dh']*X[3])
    return Y



def plot_solution(S,E,I,H,R,endt=endt,show=True):
    '''Pass R = None if no R'''

    x  = np.arange(0,endt+1,1)
    
    fig = plt.figure(figsize=(8, 4.5))
    plt.plot(x,S, label = "Suscpetible")
    plt.plot(x,E, label = "Prescription Opioid Users")
    plt.plot(x,I, label = "Prescription Opioid Addicts")
    plt.plot(x,H, label = "Heroin Addicts")
    if R is not None:
        plt.plot(x,R, label = "Resistant")
    plt.legend()
    plt.title('Population of Heroin Epidemic')
    plt.xlabel('Time (year)')
    plt.ylabel('Population')
    
    print("susceptible equals", np.round(S[-1]*TOTAL01))
    print("exposed equals", np.round(E[-1]*TOTAL01))
    print("infected equals", np.round(I[-1]*TOTAL01))
    print("heroin equals", np.round(H[-1]*TOTAL01)) #Goal: between 60,000 and 1 million
    if R is not None:
        print("resistant equals", np.round(R[-1]*TOTAL01))
    plt.tight_layout()
    if show:
        plt.show()
    else:
        return fig
    
    
if __name__ == "__main__":
    
    #establish the list incluDIng initial pops in orDEr to be appeneDEd in the function
    S = [s0]
    E = [e0]
    I = [i0]
    H = [h0]
    R = [r0]
    
    # heroin_dtds(S,E,I,H,R)

    # plot_solution(S,E,I,H,R,endt)

    S,E,I,H,R = heroin_solve_odes(s0,e0,i0,h0,r0,tsteps=endt)
    plot_solution(S,E,I,H,R,endt,False)
    S = [s0]
    E = [e0]
    I = [i0]
    H = [h0]
    R = [r0]
    S,E,I,H,R = heroin_solve_odes(s0,e0,i0,h0,r0,version=2,tsteps=endt)
    plot_solution(S,E,I,H,R,endt,False)
    S = [s0]
    E = [e0]
    I = [i0]
    H = [h0]
    R = [r0]
    S,E,I,H,R = heroin_solve_odes(s0,e0,i0,h0,r0,version=3,tsteps=endt)
    plot_solution(S,E,I,H,R,endt,False)
    plt.show()

    
