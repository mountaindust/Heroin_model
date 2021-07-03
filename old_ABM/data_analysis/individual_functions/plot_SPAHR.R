
# Function to contain all of the plotting for SPAHR
#   - Plots every ABM run, ABM mean, and ODE on one plot

#  Arguments:
#     1. abm: name of abm data
#       - Two options:
#         1) abm is the name of a file
#           - If so, this file is converted and retrieved
#         2) abm is a pre-made dataframe that assumes ABM standard naming
#     2. ode: name of ode data
#       - Two options:
#         1) ode is the name of a file 
#           - If so, this file is converted and retrieved
#         2) ode is a pre-made dataframe
#     3. plot_class: name of class to plot (default "S")
#       - current options: "S","P","A","H","R"
#     4. N: specifies whether you are plotting the SPAHR or SNPAHR model
#       - if SPAHR, N = FALSE
#       - if SNPAHR, N = TRUE


plotSPAHR = function(abm, ode, plot_class = "S"){
  
  # If "adm" is a name of a file, we need to retrieve that file
  abm_data = ABM_df(abm)
  
  # If "ode" is a name of a file, retrieve the file
  if (is.character(ode)){
    ode_data = read.csv(as.character(ode))
  }
  
  steps = max(abm_data$step)
  runs  = max(abm_data$run)
  
  if(plot_class == "S"){
    plt_c   = abm_data$S
    plt_ode = ode_data$S
    
    # Plot parameters
    plot_head = "S Proportions"
    clr       = "steelblue"
  }
  if(plot_class == "N"){
    plt_c   = abm_data$N
    plt_ode = ode_data$N
    
    plot_head = "N Proportions"
    clr       = "turquoise"
  }
  if(plot_class == "P"){
    plt_c   = abm_data$P
    plt_ode = ode_data$P
    
    plot_head = "P Proportions"
    clr       = "yellow3"
  }
  if(plot_class == "P_S"){
    plt_c    = abm_data$P_S
    plt_ode  = ode_data$PS
    
    plot_head = "P_S Proportions"
    clr       = "yellow3"
  }
  if(plot_class == "P_N"){
    plt_c    = abm_data$P_N
    plt_ode  = ode_data$PN
    
    plot_head = "P_N Proportoins"
    clr       = "orange"
  }
  if(plot_class == "A"){
    plt_c   = abm_data$A
    plt_ode = ode_data$A 
    
    plot_head = "A Proportions"
    clr       = "tomato"
  }
  if(plot_class == "H"){
    plt_c    = abm_data$H
    plt_ode = ode_data$H 
    
    plot_head = "H Proportions"
    clr       = "violet"
  }
  if(plot_class == "R"){
    plt_c   = abm_data$R
    plt_ode = ode_data$R
    
    plot_head = "R Proportions"
    clr       = "palegreen3"
  }

  
  # Find range we need for plotting:
  
  mx = max(plt_c)
  mn = min(plt_c)
  
  r_vec = find_range(upper = mx, lower = mn, range_pct = 0.35) # 0.35 works best w/ legend
  
  # Plotting the Counts individually
  
  num_vec = seq(0, steps, by = 1)
  
  for (i in seq(1, as.integer(runs), by = 1)){
    
    counts = plt_c[abm_data$run == i]
    
    if (i == 1){
      plot(num_vec, counts, type = "l", col = clr, ylim = c(r_vec[2], r_vec[1]), main = plot_head, xlab = "Ticks", ylab = "Proportion")
    }
    else {
      lines(num_vec, counts, type = "l", col = clr)
      
    }
  }
  
  # Plot the mean line on the same graph
  mval = plt_c[abm_data$step == 0 & abm_data$run == 1]
  for (i in seq(1, as.integer(steps), by = 1)){  
    mval = c(mval, mean(plt_c[abm_data$step == i]))
  }
  
  lines(num_vec, mval,type = "l", col = "black")     # Plot the mean line
  lines(num_vec, plt_ode, type = "l", col = "blue") # Plot the ODE data
  
  #  legend(0, 0.07, title = "Legend",legend = c("ABM", "ABM Mean", "ODE"), col = c("orange", "black", "blue"), lty = c(1,1,1),ncol = 1)
  legend("topright", cex = 0.8, title = "Legend",legend = c("ABM", "ABM Mean", "ODE"), col = c(clr, "black", "blue"), lty = c(1,1,1),ncol = 1)
  
}