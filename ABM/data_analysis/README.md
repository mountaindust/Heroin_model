# Data Analysis Tools for ABM
This directory contains all the tools needed to analyze the data outputted by BehaviorSpace from NetLogo. 

### Contents
#### 1. [data](https://github.com/mountaindust/Heroin_model/tree/master/ABM/data_analysis/data)
- This folder contains all of the active data being analyzed
- This is meant to organize output from BehaviorSpace, so the user should specify this as the output path when running BehaviorSpace
#### 2. [individual_functions](https://github.com/mountaindust/Heroin_model/tree/master/ABM/data_analysis/individual_functions)
- This folder contains the plotting functions for each class of the model (S, P, A, H, R) as well as a helpers.R file
- These files are all sourced in SPARH_dashboard.Rmd, so they must remain in this folder unaltered
#### 3. [SPARH_dashboard.Rmd](https://github.com/mountaindust/Heroin_model/blob/master/ABM/data_analysis/SPARH_dashboard.Rmd)
- This is the Rmarkdown file used to output a dashboard that allows us to visualize our data
- A detailed guide on how to use this file is located below

## Using the Analysis Tools
As stated above, the main tool for analysis in this directory is located in [SPARH_dashboard.Rmd](https://github.com/mountaindust/Heroin_model/blob/master/ABM/data_analysis/SPARH_dashboard.Rmd). You can view our current working dashboard at SPARH_dashboard.html. 
