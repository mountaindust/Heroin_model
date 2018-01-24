#! /usr/bin/env python3

'''Loop over a set of parameters in heroin model to observe changes'''

import numpy as np
import matplotlib.pyplot as plt
import heroin_model

def main():
    '''Examine steady states across a single parameter'''
    
    # Set up variables to store solutions
    avg_final = []
    std_final = []
    
    # Set up parameters
    params = heroin_model.params
    ### transition rates ###
    #alpha = .05 #S->E
    #beta = .1 #E->S
    #epsilon = .333 # E->I Literature
    #gamma_0 = 4 #gamma = GAMMA_0*(mu1+mu2*H_n) I->H
    #sigma = .01 #S+E->R background
    ### nonlinear transitions rates ###
    #xi = .02 #S+E->R due to H
    #delta1 = .05 #S->I due to E
    #delta2 = .25 #S->I due to I d2>d1. in [0,1)
    ### Heroin adoption rates. mu2>mu1. Affects I->H through gamma ###
    #mu1 = .00075 #S+E->H linear
    #mu2 = .1 #Literature? S+E->H due to H
    ### rehab rates ###
    #irate = .25 #rate from I into R
    #hrate = .1 #rate from H into R
    ### Death rates ###
    #de = 12305000/1245911000 #.009876 Literature
    #di = 6.3/100000 * (TOTAL01/i01) # Literature
    #dh = 2.7/100000 * (TOTAL01/h01) # Literature
    #dr = de + 1/1000 # add recetivism
    params['mu1'] = 0.005
    params['dr'] = params['de'] + 0.01
    params['irate'] = .05
    params['alpha'] = .1
    params['beta'] = .09
    
    # Parameter to loop over
    param_name = 'xi'
    param_values = np.arange(0,15+1,1)
    
    # Time to run each simulation
    simt = 300
    
    # Amount of end-time to avg over
    endt = 50
    
    # Loop
    for val in param_values:
        # Replace parameter
        params[param_name] = val
        
        # Use default initial conditions
        S = [heroin_model.s0]
        E = [heroin_model.e0]
        I = [heroin_model.i0]
        H = [heroin_model.h0]
        R = [heroin_model.r0]
        
        heroin_model.heroinfunc(S,E,I,H,R,simt,params)
        
        # Record result
        avg_S = np.mean(S[-endt:]); std_S = np.std(S[-endt:])
        avg_E = np.mean(E[-endt:]); std_E = np.std(E[-endt:])
        avg_I = np.mean(I[-endt:]); std_I = np.std(I[-endt:])
        avg_H = np.mean(H[-endt:]); std_H = np.std(H[-endt:])
        avg_R = np.mean(R[-endt:]); std_R = np.std(R[-endt:])
        
        avg_final.append((avg_S, avg_E, avg_I, avg_H, avg_R))
        std_final.append((std_S, std_E, std_I, std_H, std_R))
        
    # Plot result
    
    plt.figure()
    plt.title('Avg final value')
    plt.plot(param_values, [avg[0] for avg in avg_final], label='S')
    plt.plot(param_values, [avg[1] for avg in avg_final], label='E')
    plt.plot(param_values, [avg[2] for avg in avg_final], label='I')
    plt.plot(param_values, [avg[3] for avg in avg_final], label='H')
    plt.plot(param_values, [avg[4] for avg in avg_final], label='R')
    plt.legend()
    plt.xlabel(param_name)
    plt.xlim([param_values[0], param_values[-1]])
    plt.ylabel('Final avg value')
    
    print('Possible oscillations:')
    THRESHOLD = 1000/heroin_model.TOTAL01
    std_final = np.array(std_final)
    cycles = std_final > THRESHOLD
    for n,var in enumerate(['S','E','I','H','R']):
        if not any(cycles[:,n]):
            print(var + ': None.')
        else:
            print(var + ': ' + cycles[:,n])
    plt.show()
    
    
    
def main2D():
    '''Examine steady states across 2 parameters simultaneously'''
    
    plt.rcParams['image.cmap'] = 'cool'
    
    # Set up variables to store solutions
    avg_final = []
    std_final = []
    
    # Set up parameters
    params = heroin_model.params
    
    # Parameters to loop over
    param_name_1 = 'xi'
    param_name_2 = 'dr'
    p_val_1_start = 0
    p_val_1_end = 15
    p_val_1_step = 1
    p_val_2_start = 0.001
    p_val_2_end = 0.05
    p_val_2_step = 0.001
    param_values_1 = np.arange(p_val_1_start,p_val_1_end+p_val_1_step,
                               p_val_1_step)
    param_values_2 = params['de'] + np.arange(p_val_2_start,p_val_2_end+p_val_2_step,
                                              p_val_2_step)
    
    # Time to run each simulation
    simt = 300
    
    # Amount of end-time to avg over
    endt = 50
    
    # Loop
    for val1 in param_values_1:
        avg_final.append([])
        std_final.append([])
        # replace parameter
        params[param_name_1] = val1
        for val2 in param_values_2:
            # replace parameter
            params[param_name_2] = val2
            
            # Use default initial conditions
            S = [heroin_model.s0]
            E = [heroin_model.e0]
            I = [heroin_model.i0]
            H = [heroin_model.h0]
            R = [heroin_model.r0]
            
            heroin_model.heroinfunc(S,E,I,H,R,simt,params)
        
            # Record result
            avg_S = np.mean(S[-endt:]); std_S = np.std(S[-endt:])
            avg_E = np.mean(E[-endt:]); std_E = np.std(E[-endt:])
            avg_I = np.mean(I[-endt:]); std_I = np.std(I[-endt:])
            avg_H = np.mean(H[-endt:]); std_H = np.std(H[-endt:])
            avg_R = np.mean(R[-endt:]); std_R = np.std(R[-endt:])
            
            avg_final[-1].append((avg_S, avg_E, avg_I, avg_H, avg_R))
            std_final[-1].append((std_S, std_E, std_I, std_H, std_R))
            
    # Collect results
    avg_final = np.array(avg_final)
    std_final = np.array(std_final)
    
    # Plot
    param_values_1 = np.arange(p_val_1_start,p_val_1_end+2*p_val_1_step,
                               p_val_1_step)
    param_values_2 = params['de'] + np.arange(p_val_2_start,p_val_2_end+p_val_2_step,
                                              p_val_2_step)
    
    plt.figure()
    ax = plt.subplot(231)
    plt.title('S')
    plt.pcolormesh(param_values_2, param_values_1, avg_final[:,:,0])
    plt.xlabel(param_name_2)
    plt.ylabel(param_name_1)
    plt.colorbar()
    ticks = ax.get_xticks()
    ax.set_xticks(ticks[:-1:2])
    ax = plt.subplot(232)
    plt.title('E')
    plt.pcolormesh(param_values_2, param_values_1, avg_final[:,:,1], cmap='hot')
    plt.xlabel(param_name_2)
    plt.ylabel(param_name_1)
    plt.colorbar()
    ticks = ax.get_xticks()
    ax.set_xticks(ticks[:-1:2])
    ax = plt.subplot(233)
    plt.title('H')
    plt.pcolormesh(param_values_2, param_values_1, avg_final[:,:,3], cmap='hot')
    plt.xlabel(param_name_2)
    plt.ylabel(param_name_1)
    plt.colorbar()
    ticks = ax.get_xticks()
    ax.set_xticks(ticks[:-1:2])
    ax = plt.subplot(234)
    plt.title('R')
    plt.pcolormesh(param_values_2, param_values_1, avg_final[:,:,4])
    plt.xlabel(param_name_2)
    plt.ylabel(param_name_1)
    plt.colorbar()
    ticks = ax.get_xticks()
    ax.set_xticks(ticks[:-1:2])
    ax = plt.subplot(235)
    plt.title('I')
    plt.pcolormesh(param_values_2, param_values_1, avg_final[:,:,2], cmap='hot')
    plt.xlabel(param_name_2)
    plt.ylabel(param_name_1)
    plt.colorbar()
    ticks = ax.get_xticks()
    ax.set_xticks(ticks[:-1:2])
    plt.subplot(236)
    plt.title('Cycles')
    THRESHOLD = 1000/heroin_model.TOTAL01
    cycles = std_final > THRESHOLD
    cycles_plot = np.empty(cycles.shape[:2])
    for ii in range(cycles.shape[0]):
        for jj in range(cycles.shape[1]):
            cycles_plot[ii,jj] = int(any(cycles[ii,jj,:]))
            
    plt.pcolormesh(param_values_2, param_values_1, cycles_plot, cmap='Greys')
    plt.xlabel(param_name_2)
    plt.ylabel(param_name_1)
    
    plt.tight_layout()
    plt.show()
            
    
if __name__ == "__main__":
    main2D()
