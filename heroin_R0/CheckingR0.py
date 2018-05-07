import numpy as np 
from numpy import linalg as LA
from numpy.linalg import inv
import matplotlib.pyplot as plt
from scipy.integrate import ode

params = {}
params['alpha'] = 0.207                       #S->P the rate at which people are prescribed opioids #[20] https://www.cdc.gov/drugoverdose/pdf/pubs/2017-cdc-drug-surveillance-report.pdf (page 39, table 1A)
params['beta'] = 0.006                     #S->A total probability of becoming addicted to opioids other than by prescription #assume okay for now until find updated/source
params['xi'] =  0                     # S->A proportion of susceptibles that obtain extra prescription opioids OR black market drugs and becomes addicted (Note: MUST BE ZERO FOR AFE) #assume okay for now until find updated/source
params['theta_1'] = 0.0003                  #S->H rate susceptible population becomes addicted to heroin by black market drugs and other addicts #ESTIMATED [30] https://www.drugabuse.gov/publications/drugfacts/heroin#ref
params['delta'] = 0.15                     #R->S rate at which people come back to the susceptible class after successfully finishing treatment #[19] https://ac.els-cdn.com/S0740547213000779/1-s2.0-S0740547213000779-main.pdf?_tid=0a47e661-2ac3-45fb-a7da-9fa8c43271af&acdnat=1525189437_134c2f5932fa797a0725faa7c950a0f9 (page 1, 10% successfully treated for opioids + ESTIMATED 5% successfully treated for heroin)
params['mu'] = 0.00844                      #P,A,H,R->S natural death rate  # [24] https://www.cdc.gov/nchs/data/nvsr/nvsr66/nvsr66_05.pdf (page 1)
params['mu_A'] = 0.00004                  #A->S enhanced death rate for opioid addicts (only overdose rate=4/100,000) # [25] https://www.cdc.gov/drugoverdose/data/analysis.html (Opioid Data Analysis, figure)
params['mu_H'] = 0.00004                 #H->S enhanced death rate for heroin addicts (only overdose rate=4/100,000) # [25] https://www.cdc.gov/drugoverdose/data/analysis.html  (Opioid Data Analysis, figure)
params['gamma'] = 0          # P->A rate at which prescribed opioid users become addicted (Note: MUST BE ZERO FOR AFE) # [28] https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1360-0443.2010.03052.x
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

         
F = np.matrix( ((params['beta']*(params['epsilon']+params['mu'])/(params['alpha']+params['epsilon']+params['mu']),0, 0),
                 (0, params['theta_1']*(params['epsilon']+params['mu'])/(params['alpha']+params['epsilon']+params['mu'])\
                 +params['theta_2']*params['alpha']/(params['alpha']+params['epsilon']+params['mu']),0), 
                 (0, 0, 0)) ) 

V = np.matrix( ((params['zeta']+params['mu']+params['mu_A'],0,-params['sigma_A']), 
                (0, params['nu']+params['mu']+params['mu_H'], -params['sigma_H']),
                (-params['zeta'],-params['nu'],params['delta']+params['sigma_A']+params['sigma_H']+params['mu'])) )

Vinv = inv(np.matrix(V))

M = np.matmul(F, Vinv)

#print(LA.eigvals(M)) #R0 is the spectral radus of F*Vinv; compare to value computed below
# print(det(V))

def compute_R0(p=None):
      '''Check that AFE exists. If so, compute and return R0'''
      if p is None:
         p = params
      if p['gamma'] != 0 or p['xi'] != 0: #i.e. if gamma not =0 or xi not =0, then AFE DNE. Need at least one to be 0. 
         raise ValueError('AFE does not exist with these parameters.')
      else:
         S_star = (params['epsilon']+params['mu'])/(params['alpha']+params['epsilon']+params['mu'])
         P_star = params['alpha']/(params['alpha']+params['epsilon']+params['mu'])
         a= params['zeta']+params['mu']+params['mu_A']
         b= params['nu']+params['mu']+params['mu_H']
         c= params['delta']+params['sigma_A']+params['sigma_H']+params['mu']
         z= params['theta_1']*S_star+params['theta_2']*P_star   
         detV= a*(b*c-params['sigma_H']*params['nu'])-params['sigma_A']*params['zeta']*b
         print((1/detV)*(1/2)*\
               (params['beta']*S_star*(b*c-params['sigma_H']*params['nu'])+z*(a*c-params['zeta']*params['sigma_A'])+\
               ((params['beta']*S_star*(b*c-params['sigma_H']*params['nu'])+z*(a*c-params['zeta']*params['sigma_A']))**2-\
               4*(params['beta']*S_star*z*(b*c-params['sigma_H']*params['nu'])*(a*c-params['zeta']*params['sigma_A'])-\
               params['beta']*S_star*z*params['sigma_A']*params['nu']\
               *params['sigma_H']*params['zeta']))**.5))

# #do "return" without extra ( ) if don't want to print when run 
#in order to see this value printed in ipython: do "import CheckingR0 as cr" and then "cr.compute_R0()" or if just
# want to see it when run, do the "if_name == "__main__":
#                                       compute_R0
#at the bottom of the program 

if __name__ == "__main__":
    compute_R0()
