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
#### 4. [SPARH_dashboard.html](https://github.com/mountaindust/Heroin_model/blob/master/ABM/data_analysis/SPARH_dashboard.html)
- This is the output from knitting the SPARH_dashboard.Rmd file. 

## Using the Analysis Tools
As stated above, the main tool for analysis in this directory is located in [SPARH_dashboard.Rmd](https://github.com/mountaindust/Heroin_model/blob/master/ABM/data_analysis/SPARH_dashboard.Rmd). You can view our current working dashboard at SPARH_dashboard.html. </br> </br>

Before trying to run the dashboard code, I would highly recommend downloading [Rstudio](https://rstudio.com), a free IDE that makes editing and knitting RMD files very simple. This tutorial will show the user how to use Rstudio as an IDE for creating the dashboard.

Running this code is very simple, and as a user, you should not have to do much to get an output. I will walk through the steps to changing the RMD file so that you can run it.

### 1. Install needed packages
- The RMD code uses these packages:
1. tidyverse
2. knitr

To install these packages, simply run the following command in the console (bottom of the page in Rstudio).
```
install.packages("tidyverse")
install.packages("knitr")
```

### 2. Set working directory
- In order to ensure that your program is reading data and function files from the correct location, you must set your working directory in Rstudio as so:
```
setwd(<directory containing SPARH_dashboard.Rmd>)
```
