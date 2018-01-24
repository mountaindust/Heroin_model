import numpy as np
import matplotlib.pyplot as plt

CHILD = 48600000

#initial population values
##make up h and i populations pretend its the beginning of epiDEmic 
TOTAL01 = 317000000 - CHILD
e01 = 10000000
i01 = 400000
h01 = 100000
r01 = 1000000
s01 = TOTAL01 - (e01 + i01 + h01 + r01)

#inital pop. percents
s0 = s01/TOTAL01
e0 = e01/TOTAL01
i0 = i01/TOTAL01
h0 = h01/TOTAL01
r0 = r01/TOTAL01

#Just a check: total inital pop. as a fraction (=1)
total0 = s0 + e0 + i0 + h0 + r0


#parameters:
alpha = .05
beta = .1 
delta1 = .05 #delta1&2 exist in [0,1) d2>d1
delta2 = .25
mu1 = .00075 #mu2>m1
mu2 = .1 #Literature?

irate = .25 #rate from I into R
hrate = .1 #rate from H into R

#1/3 of people prescribed opioids get hooked/DEpenDEnt on them
EPSILON = .333
zeta = .02 #order of magnitude increase
sigma = .01


#4/5 of new heroin users come from prescription opioids 
#gamma/(mu1+mu2*H_n+gamma) = .8
GAMMA_0 = 4 #gamma = GAMMA_0*(mu1+mu2*H_n)


#Death rates
DE = 12305000/1245911000 #.009876
DI = 6.3/100000 * (TOTAL01/i01)
DH = 2.7/100000 * (TOTAL01/h01) 
dr = DE + 1/1000

endt = 100
x  = np.arange(0,endt+1,1)

def heroinfunc(S,E,I,H,R):
    for i in range(endt):
        S_n = S[-1]
        E_n = E[-1]
        I_n = I[-1]
        H_n = H[-1]
        R_n = R[-1]
        GAMMA = GAMMA_0 *(mu1+mu2*H_n)
        S.append(S_n + (beta*E_n) - (alpha*S_n) - S_n*(delta1*E_n + delta2*I_n) - S_n*(sigma+zeta*(I_n+H_n)) - S_n*(mu1+mu2*H_n) + DI*I_n + DH*H_n + DE*E_n + dr*R_n)
        E.append(E_n - (beta*E_n) + (alpha*S_n) - E_n*(sigma+zeta*(I_n+H_n)) - E_n*(mu1+mu2*H_n) - (EPSILON*E_n) - DE*E_n)
        I.append(I_n + (EPSILON*E_n) - (GAMMA*I_n) + S_n*(delta1*E_n + delta2*I_n) - irate*I_n - DI*I_n)
        H.append(H_n + (GAMMA*I_n) + (S_n+E_n)*(mu1+mu2*H_n) - hrate*H_n -DH*H_n)
        R.append(R_n + ((E_n+S_n)*(sigma + zeta*(I_n+H_n))) + irate*I_n + hrate*H_n - dr*R_n)
    plt.plot(x,S, label = "Suscpetible")
    plt.plot(x,E, label = "Prescription Opioid Users")
    plt.plot(x,I, label = "Prescription Opioid Addicts")
    plt.plot(x,H, label = "Heroin Addicts")
    plt.plot(x,R, label = "Recovered")
    plt.legend()
    plt.title('Population of Heroin Epidemic JUN15')
    plt.xlabel('Time (year)')
    plt.ylabel('Population')
    
    print("susceptible equals", S[-1]*TOTAL01)
    print("exposed equals", E[-1]*TOTAL01)
    print("infected equals", I[-1]*TOTAL01)
    print("heroin equals", H[-1]*TOTAL01) #Goal: between 60,000 and 1 million
    print("recovered equals",R[-1]*TOTAL01)
    plt.show()
    return


#establish the list incluDIng initial pops in orDEr to be appeneDEd in the function
S = [s0]
E = [e0]
I = [i0]
H = [h0]
R = [r0]

heroinfunc(S,E,I,H,R)

newpop = S[-1]+E[-1]+I[-1]+H[-1]+R[-1]

#print(newpop)
#print(total0)

assert np.isclose(total0, newpop), "Not equal"


