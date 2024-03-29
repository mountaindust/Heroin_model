# Data Analysis Tools for ABM
This directory contains all the tools needed to analyze the data outputted by BehaviorSpace from NetLogo. 

### Contents
#### 1. [individual_functions](https://github.com/mountaindust/Heroin_model/tree/master/ABM/data_analysis/individual_functions)
- This folder contains the plotting function for the model as well as some formatting functions associated with the dashboard
- These files are all sourced in SPAHR_dashboard.Rmd, so they must remain unaltered in the individual_functions folder
#### 2. [SPAHR_dashboard.Rmd](https://github.com/mountaindust/Heroin_model/blob/master/ABM/data_analysis/SPAHR_dashboard.Rmd)
- This is the Rmarkdown file used to output a dashboard that allows us to visualize our data
- A detailed guide on how to use this file is located below

## Setup
Before trying to run the dashboard code, I would highly recommend downloading [Rstudio](https://rstudio.com), a free IDE that makes editing and knitting RMD files very simple. This tutorial will show the user how to use Rstudio as an IDE for creating the dashboard. </br>

There are a few modifications that will have to be made to the [SPAHR_dashboard.Rmd](https://github.com/mountaindust/Heroin_model/blob/master/ABM/data_analysis/SPAHR_dashboard.Rmd) file before you can begin analyzing data. 

### 1. Install needed packages
The RMD code uses these packages:
1. tidyverse
2. knitr
3. kableExtra

To install these packages, simply run the following command in the console (bottom of the window in Rstudio).
```
install.packages("tidyverse")
install.packages("knitr")
```

### 2. Set working directory
In order to ensure that your program is reading data and function files from the correct location, you must set your working directory in Rstudio as so:
```
setwd(<directory containing SPAHR_dashboard.Rmd>)
```
I have included a line (line 12) in the RMD file that does this step for you. Your only job is to ensure that you replace the contents of the setwd function with your own working directory. This line of code is shown below. </br>
![Setwd line](https://github.com/mountaindust/Heroin_model/blob/master/ABM/supporting_docs/setwd_pic.png) </br>
Note: It may be beneficial to change this line to set your working directory automatically when you "Knit" the RMD code.

## Using the Analysis Tools
As stated above, the main tool for analysis in this directory is located in [SPAHR_dashboard.Rmd](https://github.com/mountaindust/Heroin_model/blob/master/ABM/data_analysis/SPAHR_dashboard.Rmd).</br> </br>

Running this code should require very little work on the users part after setup is complete. I will walk through the steps to analyzing BehaviorSpace output with this RMD code:

### 1. Obtain data

First, you must obtain data from BehaviorSpace in NetLogo. You can read [this tutorial](https://ccl.northwestern.edu/netlogo/docs/behaviorspace.html) about BehaviorSpace and how to run an experiment. The RMD code only reads the "Table" output from BehaviorSpace, so you only need to output this format. 

### 2. Modify BehaviorSpace CSV Output

This next step is absolutely crucial as the RMD code will not work if this step is skipped. The BehaviorSpace output will include a small heading at the top of the CSV file. This must be removed in order for the RMD code to read the CSV file. In your BehaviorSpace output file, you must delete all the rows highlighted below: </br>

![Rows to delete](https://github.com/mountaindust/Heroin_model/blob/master/ABM/supporting_docs/rows_to_delete.png) </br>

Then your file will be presented as a normal CSV file. If you do not have access to a spreadsheet editor such as Microsoft Excel or Apple's Numbers, then you can delete this top content from a standard text editor. This step can be seen below: </br>

![Deleting rows from a text file](https://github.com/mountaindust/Heroin_model/blob/master/ABM/supporting_docs/text_delete_rows.png)

### 3. Change name of file variables
You must also change the names of the input files inside of the RMD code. The names of the ABM data file and the ODE data file are named "abm_data" and "ode_data", respectively. They are shown below: </br>

![Names of files](https://github.com/mountaindust/Heroin_model/blob/master/ABM/supporting_docs/name_of_files.png) </br>
Note: This screenshot is from an older version of this code, so the line numbers are displayed differently in this image than they will be in your downloaded RMD file.</br>

Remember that the names of these files must be directory-dependent. For example, if SPAHR_dashboard.Rmd was contained in the "Desktop" directory (your working directory), and your BehaviorSpace output data  (named "data_abm.csv") was located in "ABM" directory, where "ABM" is inside of your "Desktop" directory, then you would enter the following on line 21:
```
abm_data = "ABM/data_abm.csv"
```

### 4. Knit the RMD file
If you have completed all of these steps, you should be ready to knit the RMD file without any issues. The code will output a .html file that is the intercative dashboard built by the code. </br>

In the toolbar at the top of the Rstudio editor, you will see an icon labeled knit, as seen below in the orange oval:
![Knit](https://github.com/mountaindust/Heroin_model/blob/master/ABM/supporting_docs/knit_pic.png) </br>

After selecting this button, you should see a dashboard immediately output to your computer's display. This code creates a dashboard called "SPAHR_dashboard.html" (same as file name, just with an .html suffix) in the same directory location as SPAHR_dashboard.Rmd on your local machine. You can view this output via a web browser. </br>

Note: You may experience some errors when Knitting this code. These will be displayed below the editor under the title "R Markdown". These errors may be associated to issues involving the steps above.


