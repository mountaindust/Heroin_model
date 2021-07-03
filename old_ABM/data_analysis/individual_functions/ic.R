# Generates initial conditions dataframes for ABM

#   Arguments:
#     1. abm_data: character name of abm
#     2. option:   which ic dataframe to get
#       - "overall", "S", "P", "A", "H", or "R"
#     3. N: tells whether your data was run on SNPAHR (TRUE) or SPAHR (FALSE)


generate_ic = function(abm_data, option = "overall", N = FALSE){
  
  ab = read.csv(abm_data) # Read in our data
  nm = colnames(ab)
  
  # Makes vector of names of every parameter to record
  if (N){
    nm = subset(nm, !(nm %in% c("X.run.number.", "X.step.", "count.turtles.with..class....S.....turtle.count", "count.turtles.with..class....N.....turtle.count", "count.turtles.with..class....P_S.....turtle.count", "count.turtles.with..class....P_N.....turtle.count","count.turtles.with..class....A.....turtle.count", "count.turtles.with..class....R.....turtle.count", "count.turtles.with..class....H.....turtle.count")))
  }
  else {
    nm = subset(nm, !(nm %in% c("X.run.number.", "X.step.", "count.turtles.with..class....S.....turtle.count", "count.turtles.with..class....P.....turtle.count", "count.turtles.with..class....A.....turtle.count", "count.turtles.with..class....R.....turtle.count", "count.turtles.with..class....H.....turtle.count")))
  }
  
  vals = ab[1,nm]       # Values for Initial Conditions
  nv = as.numeric(vals) # Convert to numeric
  
  icdf = data.frame(Initial_Value = nv) # Stands for initial condition dataframe
  row.names(icdf) = nm                  # Assign row names
  
  max_ticks = max(ab$X.step.)       # Get max number of ticks
  runs      = max(ab$X.run.number.) # Get max number of runs
  
  const_a    = icdf["constant_alpha", 1] # 1 if alpha is constant
  const_mu_a = icdf["constant_mu_a", 1]  # 1 if mu_a is constant
  
  # Overall statstics/parameters:
  
  if (option == "overall") {  
    
    if (N) {
      row_o        = c("turtle.count", "mu", "delta_t", "m", "gamma", "r", "constant_alpha", "constant_mu_a")
    } else {
      row_o        = c("turtle.count", "mu", "delta_t", "constant_alpha", "constant_mu_a")
    }
      
    icdf_overall = icdf[row_o, ]
    icdf_overall = data.frame(Initial_Value = icdf_overall)
    icdf_overall = rbind(icdf_overall, c(max_ticks))
    icdf_overall = rbind(icdf_overall, c(runs))
      
    row.names(icdf_overall) = c(row_o, "Max Ticks", "Runs")
    df_result = icdf_overall
  }
  
  if (option == "S"){
    # S parameters:
    if (N){
      row_s  = c("initial_S", "beta_A", "beta_P", "theta_1", "rho_S")
    } else {
      row_s  = c("initial_S", "beta_A", "beta_P", "theta_1")
    }
    
    icdf_s = icdf[row_s,]
    icdf_s = data.frame(Initial_Value = icdf_s)
    
    if (const_a == 1){ # If alpha is constant, add on alpha_c
      
      alpha_c = icdf["alpha_c",1]
      icdf_s = rbind(icdf_s, c(alpha_c))
      row_s  = c(row_s, "alpha")
      
    } else { # If alpha piecewise linear, get tilde values
      
      m_tilde = icdf["m_tilde", 1]        # Gets value
      icdf_s  = rbind(icdf_s, c(m_tilde)) # Adds on value to dataframe
      
      b_tilde = icdf["b_tilde", 1]
      icdf_s  = rbind(icdf_s, c(m_tilde))
      
      c_tilde = icdf["c_tilde", 1]
      icdf_s  = rbind(icdf_s, c(m_tilde))
      
      row_s   = c(row_s, "m_tilde", "b_tilde", "c_tilde")
    }
    
    row.names(icdf_s) = row_s # Set new row names
    df_result = icdf_s
  }
  
  if (option == "N"){
    row_s = c("initial_N", "rho_N")
    
    icdf_s = icdf[row_s,]
    icdf_s = data.frame(Initial_Value = icdf_s)     
    
    if (const_a == 1){ # If alpha is constant, add on alpha_c
      
      alpha_c = icdf["alpha_c",1]
      icdf_s = rbind(icdf_s, c(alpha_c))
      row_s  = c(row_s, "alpha")
      
    } else { # If alpha piecewise linear, get tilde values
      
      m_tilde = icdf["m_tilde", 1]        # Gets value
      icdf_s  = rbind(icdf_s, c(m_tilde)) # Adds on value to dataframe
      
      b_tilde = icdf["b_tilde", 1]
      icdf_s  = rbind(icdf_s, c(m_tilde))
      
      c_tilde = icdf["c_tilde", 1]
      icdf_s  = rbind(icdf_s, c(m_tilde))
      
      row_s   = c(row_s, "m_tilde", "b_tilde", "c_tilde")
    }
    
    row.names(icdf_s) = row_s # Set new row names
    df_result = icdf_s
  }
  
  if (option == "P_S"){
    row_p  = c("initial_P_S", "gamma_S", "theta_2", "epsilon_S")
    icdf_p = icdf[row_p,]
    icdf_p = data.frame(Initial_Value = icdf_p)
    
    row.names(icdf_p) = row_p
    df_result = icdf_p
  }
  if (option == "P_N"){
    row_p  = c("initial_P_N", "gamma_N", "epsilon_N")
    icdf_p = icdf[row_p,]
    icdf_p = data.frame(Initial_Value = icdf_p)
    
    row.names(icdf_p) = row_p
    df_result = icdf_p
  }
  if (option == "P"){
    # P parameters:
    row_p  = c("initial_P", "gamma", "theta_2", "epsilon")
    icdf_p = icdf[row_p,]
    icdf_p = data.frame(Initial_Value = icdf_p)
    row.names(icdf_p) = row_p
    df_result = icdf_p
  }
  
  if (option == "A"){
    # A parameters:
    row_a  = c("initial_A", "zeta", "theta_3")
    icdf_a = icdf[row_a, ]
    icdf_a = data.frame(Initial_Value = icdf_a)
    
    if (const_mu_a == 1){ # If mu_a constant, add on mu_a_c
      
      mu_a_c = icdf["mu_a_c", 1]
      icdf_a = rbind(icdf_a, c(mu_a_c))
      row_a  = c(row_a, "mu_a")
      
    } else { # If mu_a linear, add on the d_tilde, e_tilde terms
      
      d_tilde = icdf["d_tilde", 1]
      icdf_a  = rbind(icdf_a, c(d_tilde))
      
      e_tilde = icdf["e_tilde", 1]
      icdf_a  = rbind(icdf_a, c(e_tilde))
      
      row_a   = c(row_a, "d_tilde", "e_tilde")
    }
    
    row.names(icdf_a) = row_a # Assign row names
    df_result = icdf_a
  }
  
  if (option == "H"){
    # H parameters:
    row_h  = c("initial_H", "mu_h", "v_nu")
    icdf_h = icdf[row_h,]
    icdf_h = data.frame(Initial_Value = icdf_h)
    row.names(icdf_h) = row_h
    df_result = icdf_h
    
  }
  
  if (option == "R"){
    # R parameters:
    row_r  = c("initial_R", "sigma")
    icdf_r = icdf[row_r,]
    icdf_r = data.frame(Initial_Value = icdf_r)
    row.names(icdf_r) = row_r 
    df_result = icdf_r
  }
  
  return (df_result)
}