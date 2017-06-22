import numpy as np
import matplotlib.pyplot as plt

CHILD = 48600000 # number of children to remove from US population

#initial population values
##make up h and i populations pretend its the beginning of epiDEmic 
TOTAL01 = 317000000 - CHILD
e01 = 10000000
i01 = 400000
h01 = 100000
r01 = 1000000
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
# Heroin adoption rates. mu2>mu1
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
params['gamma_0'] = 4 #gamma = GAMMA_0*(mu1+mu2*H_n)

#Death rates
params['de'] = 12305000/1245911000 #.009876 Literature
params['di'] = 6.3/100000 * (TOTAL01/i01) # Literature
params['dh'] = 2.7/100000 * (TOTAL01/h01) # Literature
params['dr'] = params['de'] + 1/1000 # add recetivism

endt = 100

def heroinfunc(S,E,I,H,R,tsteps=endt,params=params):
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
    DE = params['de']
    DI = params['di']
    DH = params['dh']
    dr = params['dr']
    for i in range(tsteps):
        S_n = S[-1]
        E_n = E[-1]
        I_n = I[-1]
        H_n = H[-1]
        R_n = R[-1]
        GAMMA = params['gamma_0'] *(mu1+mu2*H_n)
        S.append(S_n + (beta*E_n) - (alpha*S_n) - S_n*(delta1*E_n + delta2*I_n) - S_n*(sigma+zeta*(I_n+H_n)) - S_n*(mu1+mu2*H_n) + DI*I_n + DH*H_n + DE*E_n + dr*R_n)
        E.append(E_n - (beta*E_n) + (alpha*S_n) - E_n*(sigma+zeta*(I_n+H_n)) - E_n*(mu1+mu2*H_n) - (EPSILON*E_n) - DE*E_n)
        I.append(I_n + (EPSILON*E_n) - (GAMMA*I_n) + S_n*(delta1*E_n + delta2*I_n) - irate*I_n - DI*I_n)
        H.append(H_n + (GAMMA*I_n) + (S_n+E_n)*(mu1+mu2*H_n) - hrate*H_n -DH*H_n)
        R.append(R_n + ((E_n+S_n)*(sigma + zeta*(I_n+H_n))) + irate*I_n + hrate*H_n - dr*R_n)
        
    assert np.isclose(1, S[-1]+E[-1]+I[-1]+H[-1]+R[-1]), "Final sum not equal to one."


def plot_solution(S,E,I,H,R,endt):

    x  = np.arange(0,endt+1,1)
    
    plt.plot(x,S, label = "Suscpetible")
    plt.plot(x,E, label = "Prescription Opioid Users")
    plt.plot(x,I, label = "Prescription Opioid Addicts")
    plt.plot(x,H, label = "Heroin Addicts")
    plt.plot(x,R, label = "Recovered")
    plt.legend()
    plt.title('Population of Heroin Epidemic JUN15')
    plt.xlabel('Time (year)')
    plt.ylabel('Population')
    
    print("susceptible equals", np.round(S[-1]*TOTAL01))
    print("exposed equals", np.round(E[-1]*TOTAL01))
    print("infected equals", np.round(I[-1]*TOTAL01))
    print("heroin equals", np.round(H[-1]*TOTAL01)) #Goal: between 60,000 and 1 million
    print("recovered equals", np.round(R[-1]*TOTAL01))
    plt.show()
    
    
if __name__ == "__main__":
    
    #establish the list incluDIng initial pops in orDEr to be appeneDEd in the function
    S = [s0]
    E = [e0]
    I = [i0]
    H = [h0]
    R = [r0]
    
    heroinfunc(S,E,I,H,R)

    plot_solution(S,E,I,H,R,endt)
