# Plot H counts in the BehaviorSpace csv file

plotH <- function(abm, ode){
  
  abm_data = read.csv(as.character(abm))
  ode_data = read.csv(as.character(ode))
  
#  abm_data = abm_data %>% rename(run = "X.run.number.", step = "X.step.", S = "count.turtles.with...class....S...", P = "count.turtles.with...class....P...", A = "count.turtles.with...class....A...", R = "count.turtles.with...class....R...", H = "count.turtles.with...class....H...")
#  abm_data = abm_data %>% rename(run = "[run number]", step = "[step]", S = "count turtles with [class = \"S\"] / turtle-count", P = "count turtles with [class = \"P\"] / turtle-count", A = "count turtles with [class = \"A\"] / turtle-count", R = "count turtles with [class = \"R\"] / turtle-count", H = "count turtles with [class = \"H\"] / turtle-count")
  abm_data = abm_data %>% rename(run = "X.run.number.", step = "X.step.", S = "count.turtles.with..class....S.....turtle.count", P = "count.turtles.with..class....P.....turtle.count", A = "count.turtles.with..class....A.....turtle.count", R = "count.turtles.with..class....R.....turtle.count", H = "count.turtles.with..class....H.....turtle.count")
  
  steps = max(abm_data$step)
  runs  = max(abm_data$run)
  
  # Find range we need for plotting:
  
  mx = max(abm_data$H)
  mn = min(abm_data$H)
  
  r_vec = find_range(upper = mx, lower = mn, range_pct = 0.35)
  
  # Plotting the H counts individually **************************************** 
  
  num_vec = seq(0, steps, by = 1)
  
  for (i in seq(1, as.integer(runs), by = 1)){
    
    H_counts = abm_data$H[abm_data$run == i] # Get counts for each run
    
    if (i == 1){ # Initial plot
      plot(num_vec, H_counts, type = "l", col = "orange", ylim = c(r_vec[2], r_vec[1]), main = "H Proportions", xlab = "Ticks", ylab = "Proportion")
    }
    else { # Adding lines thereafter
      lines(num_vec, H_counts, type = "l", col = "orange")
      
    }
  }
  
  mean_H = abm_data$H[abm_data$step == 0 & abm_data$run == 1] # Intialize mean value
  
  for (i in seq(1, as.integer(steps), by = 1)){ # Creating the mean line for ABM
    mean_H = c(mean_H, mean(abm_data$H[abm_data$step == i]))
  }
  
  lines(num_vec, mean_H,type = "l", col = "black")     # Plot the mean line
  lines(num_vec, ode_data$H, type = "l", col = "blue") # Plot the ODE data
  
  # Gives legend to plot
  legend("topright", cex = 0.8, title = "Legend",legend = c("ABM", "ABM Mean", "ODE"), col = c("orange", "black", "blue"), lty = c(1,1,1),ncol = 1)
  
  
  end_mean = tail(mean_H, n = 1)     # Gets the ending mean value of the ABM
  end_ODE  = tail(ode_data$H, n = 1) # Gets the ending value of the ODE
  ret = (c(end_ODE, end_mean))
}