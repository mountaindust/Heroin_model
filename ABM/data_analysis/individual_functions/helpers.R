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