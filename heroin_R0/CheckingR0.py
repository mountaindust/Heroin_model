import numpy as np 
from numpy import linalg as LA
from numpy.linalg import inv
import matplotlib.pyplot as plt
from scipy.integrate import ode

params = {}
params['alpha'] = 0.283 #0.303                     #S->P the rate at which people are prescribed opioids #[20] https://www.cdc.gov/drugoverdose/pdf/pubs/2017-cdc-drug-surveillance-report.pdf (page 39, table 1A)
params['beta_A'] = 0.0044 #0.00235                   #S->A total probability of becoming addicted to opioids other than by prescription #assume okay for now until find updated/source
params['beta_P'] =  0                    # S->A proportion of susceptibles that obtain extra prescription opioids OR black market drugs and becomes addicted (Note: MUST BE ZERO FOR AFE) #assume okay for now until find updated/source
params['theta_1'] = 0.000502 #0.000507                  #S->H rate susceptible population becomes addicted to heroin by black market drugs and other addicts #ESTIMATED [30] https://www.drugabuse.gov/publications/drugfacts/heroin#ref
params['mu'] = 0.00868                      #P,A,H,R->S natural death rate  # [24] https://www.cdc.gov/nchs/data/nvsr/nvsr66/nvsr66_05.pdf (page 1)
params['mu_A'] = 0.00870                  #A->S enhanced death rate for opioid addicts (only overdose rate=4/100,000) # [25] https://www.cdc.gov/drugoverdose/data/analysis.html (Opioid Data Analysis, figure)
params['mu_H'] = 0.0507                #H->S enhanced death rate for heroin addicts (only overdose rate=4/100,000) # [25] https://www.cdc.gov/drugoverdose/data/analysis.html  (Opioid Data Analysis, figure)
params['gamma'] = 0          # P->A rate at which prescribed opioid users become addicted (Note: MUST BE ZERO FOR AFE) # [28] https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1360-0443.2010.03052.x
params['epsilon'] = 2.49 #2.54        #P->S rate at which people come back to the susceptible class after being prescribed opioids (i.e. not addicted) #[28] https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1360-0443.2010.03052.x (1-0.26=0.74)
params['theta_2'] = 0.148 #0.0370                       #P->H rate at which opioid prescribed user population becomes addicted to heroin #[30] https://www.drugabuse.gov/publications/drugfacts/heroin#ref
params['sigma'] = 0.0283 #0.0284  #R->A rate at which people relapse from treatment into the opioid addicted class #[19] https://ac.els-cdn.com/S0740547213000779/1-s2.0-S0740547213000779-main.pdf?_tid=0a47e661-2ac3-45fb-a7da-9fa8c43271af&acdnat=1525189437_134c2f5932fa797a0725faa7c950a0f9 (page 1, 90% relapse rate by end of one year)
params['zeta'] = 0.318 #0.265                       #A->R rate at which addicted opioid users enter treatment/rehabilitation #assume okay for now until find updated/source
params['theta_3'] = 2.38 #3.51                   #A->H rate at which the opioid addicted population becomes addicted to heroin #[14] https://d14rmgtrwzf5a.cloudfront.net/sites/default/files/19774-prescription-opioids-and-heroin.pdf (page 7)
params['nu'] = 0.0482 #0.00657                         #H->R rate at which heroin users enter treatment/rehabilitation #[14] https://d14rmgtrwzf5a.cloudfront.net/sites/default/files/19774-prescription-opioids-and-heroin.pdf (page 17)
params['omega']=0.0000000001


def update_params(new_params):
    '''Update the params dict with new values contained in new_params'''
    global params
    for key,val in new_params.items():
        if key in params.keys():
            params[key] = val
        else:
            print('Could not find parameter {}.'.format(key))

         
F = np.matrix( ((params['beta_A']*(params['epsilon']+params['mu'])/(params['alpha']+params['epsilon']+params['mu']),0),
                 (0, params['theta_1']*(params['epsilon']+params['mu'])/(params['alpha']+params['epsilon']+params['mu'])\
                 +params['theta_2']*params['alpha']/(params['alpha']+params['epsilon']+params['mu']))) ) 

V = np.matrix( ((params['zeta']+params['mu']+params['mu_A'],0), 
                (0, params['nu']+params['mu']+params['mu_H'])) )

Vinv = inv(np.matrix(V))

M = np.matmul(F, Vinv)

print(LA.eigvals(M)) #R0 is the spectral radus of F*Vinv; compare to value computed below
#print(det(V))

def compute_R0(p=None):
      '''Check that AFE exists. If so, compute and return R0'''
      if p is None:
         p = params
      if p['gamma'] != 0 or p['beta_P'] != 0: #i.e. if gamma not =0 or beta_P not =0, then AFE DNE. Need at least one to be 0. 
         raise ValueError('AFE does not exist with these parameters.')
      else:
         print(params['beta_A']*(params['epsilon']+params['mu'])/((params['alpha']+params['epsilon']+params['mu'])*(params['zeta']+params['mu']+params['mu_A'])),
         (params['theta_1']*(params['epsilon']+params['mu'])+params['theta_2']*params['alpha'])/((params['alpha']+params['epsilon']+params['mu'])*(params['nu']+params['mu']+params['mu_H'])))

# #do "return" without extra ( ) if don't want to print when run 
#in order to see this value printed in ipython: do "import CheckingR0 as cr" and then "cr.compute_R0()" or if just
# want to see it when run, do the "if_name == "__main__":
#                                       compute_R0
#at the bottom of the program 

if __name__ == "__main__":
    compute_R0()
