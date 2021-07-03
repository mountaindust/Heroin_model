# File containing helper functions for the plotting functions

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


# meanLineABM:
#   Returns a vector corresponding to the mean line for ABM data

# Arguments:
#   1. abm = two options:
#       1) character value corresponding to filename
#       2) dataframe

#   2. plot_class = character value corresponding to class of interest
#       - Options: 'S', 'N', 'P', 'P_S', 'P_N', 'A', 'H', 'R'

#   3. need_name = boolean corresponding to whether dataframe needs to be renamed
#       - If dataframe needs to be renamed: TRUE
#       - If dataframe does not need renaming: FALSE

meanLineABM = function(abm, plot_class, need_name = TRUE){
  
  # Get runs data
  if (is.character(abm)){
    abm = read.csv(as.character(abm))
  }
  
  if (need_name){
    abm = abm %>% rename(run = "X.run.number.", step = "X.step.")
  }
  
  # Get steps and runs for each
  steps = max(abm$step)
  runs  = max(abm$run)
  
  # Option for class
  
  if(plot_class == "S"){
    if (need_name){
      abm = abm %>% rename(S = "count.turtles.with..class....S.....turtle.count")
    }
    
    interest = abm$S
  }
  if(plot_class == "N"){
    if (need_name){
      abm = abm %>% rename(N = "count.turtles.with..class....N.....turtle.count")
    }
    
    interest = abm$N
  }
  if(plot_class == "P"){
    if(need_name){
      abm = abm %>% rename(P = "count.turtles.with..class....P.....turtle.count")
    }
    
    interest = abm$P
  }
  if(plot_class == "P_S"){
    if (need_name){
      abm = abm %>% rename(P_S = "count.turtles.with..class....P_S.....turtle.count")
    }
    
    interest = abm$P_S
  }
  if(plot_class == "P_N"){
    if (need_name){
      abm = abm %>% rename(P_N = "count.turtles.with..class....P_N.....turtle.count")
    }
    
    interest = abm$P_N
  }
  if(plot_class == "A"){
    if (need_name){
      abm = abm %>% rename(A = "count.turtles.with..class....A.....turtle.count")
    }
    
    interest = abm$A
  }
  if(plot_class == "H"){
    if (need_name){
      abm = abm %>% rename(H = "count.turtles.with..class....H.....turtle.count")
    }
    
    interest = abm$H
  }
  if(plot_class == "R"){
    if (need_name){
      abm = abm %>% rename(R = "count.turtles.with..class....R.....turtle.count")
    }
    
    interest = abm$R
  }
  
  mval = c()
  for (i in seq(1, as.integer(steps), by = 1)){  
    mval = c(mval, mean(interest[abm$step == i]))
  }
  
  return(mval)
}


# ABM_df: generates a csv file name for ABM data based off of the name of the file

# Arguments:
#   1. abm: name of file to be retrieved, converted
#   2. N: specifies whether abm is SPAHR or SNPAHR
#     - If abm is SPAHR data, N = FALSE
#     - If abm is SNPAHR data, N = TRUE

ABM_df = function(abm){
  
  if (!(is.character(abm))){ # If abm is already a dataframe, don't need conversion
    return(abm)
  }
  
  abm_data = read.csv(as.character(abm))
  
  abm_data = abm_data %>% rename(run = "X.run.number.", step = "X.step.")
  
  abm_data = abm_data %>% rename(S = "count.turtles.with..class....S.....turtle.count")
  abm_data = abm_data %>% rename(A = "count.turtles.with..class....A.....turtle.count")
  abm_data = abm_data %>% rename(H = "count.turtles.with..class....H.....turtle.count")
  abm_data = abm_data %>% rename(R = "count.turtles.with..class....R.....turtle.count")
  
  N = FALSE
  if("count.turtles.with..class....N.....turtle.count" %in% colnames(abm)){
    N = TRUE
  }
  
  if (N){
    abm_data = abm_data %>% rename(N = "count.turtles.with..class....N.....turtle.count")
    abm_data = abm_data %>% rename(P_S = "count.turtles.with..class....P_S.....turtle.count")
    abm_data = abm_data %>% rename(P_N = "count.turtles.with..class....P_N.....turtle.count")
  } else {
    abm_data = abm_data %>% rename(P = "count.turtles.with..class....P.....turtle.count")
  }
  
  return(abm_data)
}


# abm_data MUST BE A DATAFRAME already read in

# downsample_ABM: downsamples ABM steps by specified factor

# Arguments:
#   1. adm_data: dataframe on which to downsample
#   2. every: specifies the factor which you want to downsample by
#     - i.e. your final df will have (nrow(abm_data) / every) rows
#     - Gets steps that are multiples of this number
#   3. verbose: option to print steps to screen during runtime
#

downsample_ABM = function(abm_data, every = 10, verbose = TRUE){
  
  runs = max(abm_data$run)
  steps = max(abm_data$step)
  
  ret_abm = data.frame()
  
  for (j in seq(0, runs)){
    
    # Access all of the necessary steps within each run
    for (i in seq(from = 0, to = steps, by = every)){
    
      if ((i == 0) && (j == 0)){ # Initialize new dataframe
        
        ret_abm = data.frame(subset(abm_data, step == i & run == j))
        colnames(ret_abm) = colnames(abm_data)                       # Make sure we rename the df
        
      } else {
        ret_abm = rbind(ret_abm, subset(abm_data, step == i & run == j))
      }
      
      if (verbose){
        print(paste( paste("run: ", as.character(j)), paste("--- step: ", as.character(i)))) # print progress
      }
    }
  }

  # Convert values in step col
  ret_abm$step = ret_abm$step / every
  
  return (ret_abm)
}