import csv
import numpy as np
import pandas as pd
from collections import Counter

# SPAHR_ABM class:

# Arguments to constructor:
#   1. parameters = dict
#       - Contains parameters for model to run
#   2. constants = array of length 2
#       - specifies if alpha, mu_a is constant, respectively

class SPAHR_ABM:

    class_names = ["[run number]", '[step]', 
                "count turtles with [class = \"S\"] / turtle-count",
                "count turtles with [class = \"P\"] / turtle-count",
                "count turtles with [class = \"A\"] / turtle-count",
                "count turtles with [class = \"H\"] / turtle-count",
                "count turtles with [class = \"R\"] / turtle-count"]

    # Constructor
    def __init__(self, parameters, constants = [1, 1]):

        self.n_turtles  = parameters['turtle-count']
        self.n_ticks    = parameters['num-ticks']
        self.delta_t    = parameters['delta_t']

        self.m_tilde    = parameters['m_tilde']         
        self.b_tilde    = parameters['b_tilde']  
        self.c_tilde    = parameters['c_tilde']       
        self.beta_A     = parameters['beta_A']   
        self.beta_P     = parameters['beta_P']    
        self.theta_1    = parameters['theta_1']  
        self.mu         = parameters['mu']        

        self.gamma      = parameters['gamma']                
        self.epsilon    = parameters['epsilon'] 
        self.theta_2    = parameters['theta_2']
    
        self.zeta       = parameters['zeta']      
        self.theta_3    = parameters['theta_3']
        self.d_tilde    = parameters['d_tilde']         
        self.e_tilde    = parameters['e_tilde']     

        self.nu         = parameters['nu']
        self.mu_h       = parameters['mu_H']

        self.sigma      = parameters['sigma']      

        self.S0         = parameters['S0']
        self.P0         = parameters['P0']
        self.A0         = parameters['A0']
        self.H0         = parameters['H0']
        self.R0         = parameters['R0']

        self.name       = str(self.n_turtles) + "_turtles_ABM"

        if (constants[0] == 1):
            self.alpha_const = True
            self.alpha = parameters["alpha"]
        else:
            self.alpha_const = False
        
        if (constants[1] == 1):
            self.mu_a_const = True
            self.mu_a = parameters["mu_a"]
        else:
            self.mu_a_const = False

    # Piecewise linear version of alpha_t
    # For internal use only
    def __alpha_t(self, ticks):
        if (not self.alpha_const):
            t = (ticks * self.delta_t)

            if (t <= 3.25):
                self.alpha = ( self.m_tilde * t + self.b_tilde )
            else:
                self.alpha = ( (self.m_tilde * 3.25) + self.b_tilde + self.c_tilde * (t - 3.25) )

    # Linear version of mu_a_t
    # For internal use only
    def __mu_a_t(self, ticks):
        if (not self.mu_a_const):
            t = (ticks * self.delta_t)
            self.mu_a = (self.d_tilde * t + self.e_tilde)

    def robust_par_toDF(self):

        rcount = self.model_data.shape[0]

        self.model_data["turtle-count"]   = [self.n_turtles] * rcount
        self.model_data["mu"]             = [self.mu] * rcount
        self.model_data["delta_t"]        = [self.delta_t] * rcount
        self.model_data["initial_S"]      = [self.S0] * rcount
        self.model_data["beta_A"]         = [self.beta_A] * rcount
        self.model_data["beta_P"]         = [self.beta_P] * rcount
        self.model_data["theta_1"]        = [self.theta_1] * rcount
        self.model_data["constant_alpha"] = [int(self.alpha_const)] * rcount
        self.model_data["alpha_c"]        = [self.alpha] * rcount
        self.model_data["m_tilde"]        = [self.m_tilde] * rcount
        self.model_data["b_tilde"]        = [self.b_tilde] * rcount
        self.model_data["c_tilde"]        = [self.c_tilde] * rcount
        self.model_data["initial_P"]      = [self.P0] * rcount
        self.model_data["gamma"]          = [self.gamma] * rcount
        self.model_data["epsilon"]        = [self.epsilon] * rcount
        self.model_data["theta_2"]        = [self.theta_2] * rcount
        self.model_data["initial_A"]      = [self.A0] * rcount
        self.model_data["zeta"]           = [self.zeta] * rcount
        self.model_data["theta_3"]        = [self.theta_3] * rcount
        self.model_data["constant_mu_a"]  = [self.mu_a_const] * rcount
        self.model_data["mu_a_c"]         = [self.mu_a] * rcount
        self.model_data["d_tilde"]        = [self.d_tilde] * rcount
        self.model_data["e_tilde"]        = [self.e_tilde] * rcount
        self.model_data["initial_H"]      = [self.H0] * rcount
        self.model_data["v_nu"]           = [self.nu] * rcount
        self.model_data["mu_h"]           = [self.mu_h] * rcount
        self.model_data["initial_R"]      = [self.R0] * rcount
        self.model_data["sigma"]          = [self.sigma] * rcount

    # Function for internal use only
    # Runs the full model once
    # Main implementation of ABM
    def __run_model(self, alpha_constant = True, mu_a_constant = True, relapse_mem = False):

        # Initialize the number of turtles

        S0  = round(self.S0 * self.n_turtles)
        P0  = round(self.P0 * self.n_turtles)
        A0  = round(self.A0 * self.n_turtles)
        H0  = round(self.H0 * self.n_turtles)
        R0  = round(self.R0 * self.n_turtles)

        turtles = ['S' for ii in range(S0)]
        turtles += ['P' for ii in range(P0)]
        turtles += ['A' for ii in range(A0)]
        turtles += ['H' for ii in range(H0)]
        turtles += ['R' for ii in range(R0)]

        # Implementing addiction history in model runs
        addict_hist = np.zeros(self.n_turtles)

        self.turtle_hist = []

        # running loop

        for tick in range(self.n_ticks):
            self.turtle_hist.append(list(turtles))      # Add new list of turtles to back
            cnt_turtles = Counter(self.turtle_hist[-1]) # Get count of list in back

            self.__alpha_t(tick)
            self.__mu_a_t(tick)

            for n, turtle in enumerate(self.turtle_hist[-1]):

                # S case
                if turtle == 'S':
                    lam = (self.alpha + 
                           self.beta_A * (cnt_turtles['A'] / self.n_turtles) + 
                           self.beta_P * (cnt_turtles['P'] / self.n_turtles) +
                           self.theta_1 * (cnt_turtles['H'] / self.n_turtles) + self.mu)
                    chng = 1 - np.exp(-lam * self.delta_t)
                    rf1 = np.random.rand()
                    if rf1 < chng:                  # Change block
                        rf2 = np.random.rand()
                        s_to_p = self.alpha / lam
                        s_to_a = (self.beta_A * (cnt_turtles['A'] / self.n_turtles) + 
                                  self.beta_P * (cnt_turtles['P'] / self.n_turtles)) / lam
                        s_to_h = self.theta_1 * (cnt_turtles['H'] / self.n_turtles) / lam

                        if rf2 < s_to_p:
                            turtles[n] = 'P'
                        elif rf2 < (s_to_a + s_to_p):
                            turtles[n] = 'A'
                        elif rf2 < (s_to_h + s_to_a + s_to_p):
                            turtles[n] = 'H'

                # PS Case
                elif turtle == 'P':

                    lam = (self.epsilon + self.gamma + 
                          (self.theta_2 * (cnt_turtles['H'] / self.n_turtles)) + self.mu)

                    chng = 1 - np.exp(-lam * self.delta_t)
                    rf1 = np.random.rand()

                    if rf1 < chng:                  # Change block
                        rf2 = np.random.rand()

                        p_to_s = self.epsilon / lam
                        p_to_a = self.gamma / lam
                        p_to_h = self.theta_2 * (cnt_turtles['H'] / self.n_turtles) / lam

                        if rf2 < p_to_s: 
                            turtles[n] = 'S'
                        elif rf2 < (p_to_a + p_to_s):
                            turtles[n] = 'A'
                        elif rf2 < (p_to_h + p_to_a + p_to_s):
                            turtles[n] = 'H'
                        else:                # DEATH CASE
                            turtles[n] = 'S'

                # A Case
                elif turtle == 'A':
                    lam = ( self.zeta + self.theta_3 * (cnt_turtles['H'] / self.n_turtles) + 
                            self.mu + self.mu_a )

                    chng = 1 - np.exp(-lam * self.delta_t)
                    rf1 = np.random.rand()

                    if rf1 < chng:                  # Change Block
                        rf2 = np.random.rand()
                        a_to_r = self.zeta / lam
                        a_to_h = self.theta_3 * (cnt_turtles['H'] / self.n_turtles) / lam
                        if rf2 < a_to_r:
                            turtles[n] = 'R'

                            if(relapse_mem):          # Addiction history
                                addict_hist[n] = 'A'

                        elif rf2 < (a_to_h + a_to_r):
                            turtles[n] = 'H'
                        else:
                            turtles[n] = 'S'

                # H Case
                elif turtle == 'H':
                    lam = (self.nu + self.mu + self.mu_h)
                    chng = 1 - np.exp(-lam * self.delta_t)
                    rf1 = np.random.rand()
                    if rf1 < chng:                  # Change Block
                        rf2 = np.random.rand()
                        h_to_r = self.nu / lam
                        if rf2 < h_to_r:
                            turtles[n] = 'R'

                            if(relapse_mem):
                                addict_hist[n] = 'H'

                        else:
                            turtles[n] = 'S'

                # R Case
                else:
                    
                    if(relapse_mem):
                        lam = self.sigma + self.mu

                    else:
                        lam = (self.sigma * cnt_turtles['A'] / (cnt_turtles['A']+cnt_turtles['H']+0.0000001) +
                               self.sigma * cnt_turtles['H'] / (cnt_turtles['A']+cnt_turtles['H']+0.0000001) + self.mu)

                    chng = 1 - np.exp(-lam * self.delta_t)
                    rf1 = np.random.rand()

                    if rf1 < chng:                # Change block
                        rf2 = np.random.rand()

                        if(relapse_mem):  # Relapse memory code
                            r_to_aorh = (self.sigma) / lam 

                            if (rf2 < r_to_aorh):
                                turtles[n] = addict_hist[n]
                            else: # Death
                                turtles[n] = 'S'

                        else: # Non relapse memory code
                            r_to_a = (self.sigma * cnt_turtles['A'] / (cnt_turtles['A']+cnt_turtles['H']+0.0000001)) / lam
                            r_to_h = (self.sigma * cnt_turtles['H'] / (cnt_turtles['A']+cnt_turtles['H']+0.0000001)) / lam

                            if rf2 < r_to_a:
                                turtles[n] = 'A'
                            elif rf2 < (r_to_h + r_to_a):
                                turtles[n] = 'H'
                            else: # Death
                                turtles[n] = 'S'

        # include last result and return
        self.turtle_hist.append(list(turtles)) # Append final list
        # This value is now store in the class as turtle_hist


    # Creates dataframe from running the model for a given number of runs
    #   Arguments:
    #       n_runs: this is the number of runs the model should run
    #       relapse_mem: decides whether to run relapse memory code (=True) or not (=False)
    #       progress_bar: print a progress bar (=True) or not (=False)

    # See Jupyter Notebook for more info on running
    def run_model(self, n_runs, relapse_mem = False, progress_bar = True):

        np.random.seed()

        self.model_data = pd.DataFrame(columns = self.class_names)

        print("Running SPAHR ABM: " + str(n_runs) + " runs")

        row_counter = 0

        for run in range(n_runs + 1):
            self.__run_model()
            data = []
            for j, i in enumerate(self.turtle_hist):
                cnt_turtles = Counter(i)
                data.append(run)
                data.append(j)
                        
                data.append( cnt_turtles['S'] / self.n_turtles )
                data.append( cnt_turtles['P'] / self.n_turtles )
                data.append( cnt_turtles['A'] / self.n_turtles )
                data.append( cnt_turtles['H'] / self.n_turtles )
                data.append( cnt_turtles['R'] / self.n_turtles )
                self.model_data.loc[row_counter] = data
                row_counter += 1
                data = []

            # Progress bar established here
            if ((( run / n_runs * 100) % 10 == 0) and progress_bar):
                bar = ""
                pct = (run / n_runs * 100)
                for i in range(0, int(run / n_runs * 10)):
                    bar += '-'
                print(str(pct) + "% " + bar)

        if(progress_bar):
            print("Finished")

    # Writes the run_model data to a csv file
    # Arguments:
    #       file_name: specifies file name to output results to
    #       robust: includes parameters in output (=True) or not (=False)
    def data_to_csv(self, file_name = "NONE", robust = True):

        try:
            self.model_data
        except NameError:
            print("Model has not yet been ran")
            print("Run run_model function")
            return -1
        else:

            file_name = str(file_name)

            if (file_name == "NONE"): # file_name will have custom name for class if none specified
                file_name = str(self.name) + ".csv"

            if ((file_name[-4:-1] + file_name[-1]) != ".csv"): # If declared file_name doesn't have ".csv" suffix
                file_name = str(file_name) + ".csv"

            if (robust):  # Write parameters if robust
                self.robust_par_toDF()

            print("Writing to " + str(file_name))
            self.model_data.to_csv(file_name, index = False)

    # Writes the run_model data to an excel file
    # Arguments:
    #       file_name: specifies file name to output results to
    #       robust: includes parameters in output (=True) or not (=False)
    def data_to_excel(self, file_name = "NONE", robust = True):

        try:
            self.model_data
        except NameError:
            print("Model has not yet been run")
            print("Run run_model function")
            return -1
        else:

            file_name = str(file_name)

            if (file_name == "NONE"): # file_name will have custom name for class if none specified
                file_name = str(self.name) + ".csv"

            if ((file_name[-4:-1] + file_name[-1]) != ".csv"): # If declared file_name doesn't have ".csv" suffix
                file_name = str(file_name) + ".csv"

            if (robust): # Write parameters to df if robust
                self.robust_par_toDF()

            print("Writing to " + str(file_name))
            self.model_data.to_excel(file_name, index = False)
