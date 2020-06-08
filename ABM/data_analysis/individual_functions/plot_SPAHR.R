
# Function to contain all of the plotting for SPAHR
#   - Plots every ABM run, ABM mean, and ODE on one plot

#  Arguments:
#     1. abm: name of abm data
#     2. ode: name of ode data
#     3. plot_class: name of class to plot (default "S")
#       - current options: "S","P","A","H","R"

plotSPAHR = function(abm, ode, plot_class = "S", N = FALSE){
  
  abm_data = read.csv(as.character(abm))
  ode_data = read.csv(as.character(ode))
  
  abm_data = abm_data %>% rename(run = "X.run.number.", step = "X.step.")

  steps = max(abm_data$step)
  runs  = max(abm_data$run)
  
  if(plot_class == "S"){
    abm_data = abm_data %>% rename(S = "count.turtles.with..class....S.....turtle.count")
    plt_c   = abm_data$S
    plt_ode = ode_data$S
    
    # Plot parameters
    plot_head = "S Proportions"
    clr       = "steelblue"
  }
  if(plot_class == "N"){
    abm_data = abm_data %>% rename(N = "count.turtles.with..class....N.....turtle.count")
    plt_c   = abm_data$N
    plt_ode = ode_data$N
    
    plot_head = "N Proportions"
    clr       = "turquoise"
  }
  if(plot_class == "P"){
    abm_data = abm_data %>% rename(P = "count.turtles.with..class....P.....turtle.count")
    plt_c   = abm_data$P
    plt_ode = ode_data$P
    
    plot_head = "P Proportions"
    clr       = "yellow2"
  }
  if(plot_class == "P_S"){
    abm_data = abm_data %>% rename(P_S = "count.turtles.with..class....P_S.....turtle.count")
    plt_c    = abm_data$P_S
    #plt_ode = <what>
    
    plot_head = "P_S Proportions"
    clr       = "yellow2"
    
  }
  if(plot_class == "P_N"){
    abm_data = abm_data %>% rename(P_N = "count.turtles.with..class....P_N.....turtle.count")
    plt_c    = abm_data$P_N
    #plt_ode = <what>
    
    plot_head = "P_N Proportoins"
    clr       = "orange"
  }
  if(plot_class == "A"){
    abm_data = abm_data %>% rename(A = "count.turtles.with..class....A.....turtle.count")
    plt_c   = abm_data$A
    plt_ode = ode_data$A 
    
    plot_head = "A Proportions"
    clr       = "tomato"
  }
  if(plot_class == "H"){
    abm_data = abm_data %>% rename(H = "count.turtles.with..class....H.....turtle.count")
    plt_c    = abm_data$H
    plt_ode = ode_data$H 
    
    plot_head = "H Proportions"
    clr       = "violet"
  }
  if(plot_class == "R"){
    abm_data = abm_data %>% rename(R = "count.turtles.with..class....R.....turtle.count")
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
  
  mval = plt_c[abm_data$step == 0 & abm_data$run == 1]
  for (i in seq(1, as.integer(steps), by = 1)){  
    mval = c(mval, mean(plt_c[abm_data$step == i]))
  }
  
  lines(num_vec, mval,type = "l", col = "black")     # Plot the mean line
  lines(num_vec, plt_ode, type = "l", col = "blue") # Plot the ODE data
  
  #  legend(0, 0.07, title = "Legend",legend = c("ABM", "ABM Mean", "ODE"), col = c("orange", "black", "blue"), lty = c(1,1,1),ncol = 1)
  legend("topright", cex = 0.8, title = "Legend",legend = c("ABM", "ABM Mean", "ODE"), col = c(clr, "black", "blue"), lty = c(1,1,1),ncol = 1)
  
  
}

# find_range: **********************
# Gets the range of values to include for the plots
#     Arguments:
#       1. upper     : maximum value of the data being considered
#       2. lower     : minimum value of the data
#       3. range_pct : percent of the range to go up/down for range
find_range = function(upper, lower, range_pct){
  
  ran = upper - lower # Range to consider
  
  ran_add = ran * range_pct # Add to max/min
  
  u = upper + ran_add # Max and min vals
  d = lower - ran_add
  
  if (d < 0){ d = 0 } # Reset lower to 0 if needed
  return (c(u, d))
}