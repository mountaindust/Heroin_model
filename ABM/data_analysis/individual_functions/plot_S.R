# Plot S counts in the BehaviorSpace csv file


plotS <- function(abm, ode){
  abm_data = read.csv(as.character(abm))
  ode_data = read.csv(as.character(ode))

  # Rename our columns:
 # abm_data = abm_data %>% rename(run = "X.run.number.", step = "X.step.", S = "count.turtles.with...class....S...", P = "count.turtles.with...class....P...", A = "count.turtles.with...class....A...", R = "count.turtles.with...class....R...", H = "count.turtles.with...class....H...")
#  abm_data = abm_data %>% rename(run = "[run number]", step = "[step]", S = "count turtles with [class = \"S\"] / turtle-count", P = "count turtles with [class = \"P\"] / turtle-count", A = "count turtles with [class = \"A\"] / turtle-count", R = "count turtles with [class = \"R\"] / turtle-count", H = "count turtles with [class = \"H\"] / turtle-count")
  abm_data = abm_data %>% rename(run = "X.run.number.", step = "X.step.", S = "count.turtles.with..class....S.....turtle.count", P = "count.turtles.with..class....P.....turtle.count", A = "count.turtles.with..class....A.....turtle.count", R = "count.turtles.with..class....R.....turtle.count", H = "count.turtles.with..class....H.....turtle.count")
  
  steps = max(abm_data$step)
  runs  = max(abm_data$run)
  
  # Find range we need for plotting:
  
  mx = max(abm_data$S)
  mn = min(abm_data$S)
  
  r_vec = find_range(upper = mx, lower = mn, range_pct = 0.35)
  
  # Plotting the S counts individually ****************************************
  mean_S = abm_data$S[abm_data$step == 0 & abm_data$run == 1]
  
  num_vec = seq(0, steps, by = 1)

  for (i in seq(1, as.integer(runs), by = 1)){
    S_counts = abm_data$S[abm_data$run == i]
  
    if (i == 1){
      plot(num_vec, S_counts, type = "l", col = "steelblue", ylim = c(r_vec[2], r_vec[1]), main = "S Proportions", xlab = "Ticks", ylab = "Proportion")
    }
    else {
      lines(num_vec, S_counts, type = "l", col = "steelblue")
    }
  
  }
  
  
  mean_S = abm_data$S[abm_data$step == 0 & abm_data$run == 1]
  for (i in seq(1, as.integer(steps), by = 1)){  
    mean_S = c(mean_S, mean(abm_data$S[abm_data$step == i]))
  }

  lines(num_vec, mean_S,type = "l", col = "black")
  lines(num_vec, ode_data$S, type = "l", col = "blue")

#  legend(-1, 0.757, title = "Legend",legend = c("ABM", "ABM Mean", "ODE"), col = c("steelblue", "black", "blue"), lty = c(1,1,1),ncol = 1)
  legend("topright", cex = 0.8, title = "Legend",legend = c("ABM", "ABM Mean", "ODE"), col = c("steelblue", "black", "blue"), lty = c(1,1,1),ncol = 1)
  
  end_mean = tail(mean_S, n = 1)
  end_ODE  = tail(ode_data$S, n = 1)
  ret = (c(end_ODE, end_mean))
}