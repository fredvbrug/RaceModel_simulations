# RaceModel_simulations
Here we provide the source code for the simulations that are reported in "Capturing the ability to inhibit actions and impulsive behaviors: A consensus guide to the stop-signal task". Detailed information about the simulations can be found in this manuscript. There are two main folders in this repository.

### The 'tau_x_omission' folder
In a first set of simulations ('tau_x_omission' folder), we compared different SSRT estimation methods. Performance in the stop-signal task was simulated based on assumptions of the independent race model, and we explored how SSRT estimation bias and reliability were influenced by (1) the number of trials in a task and (2) changes in go performance (in particular, skew (or tau) of the go RT distribution, and the proportion of go omissions). In the main 'tau_x_omission' folder, three R scripts can be found:

- *simulation.R*: the R code used for the simulations.
- *estimator.R*: the R code used to estimate SSRTs for the simulated 'participants'.
- *analyser.R*: the R code used to process the SSRT estimates and to create the figures in the manuscript.

The folder also contains three subfolders:

- */simulated_data*: this folder contains the simulated data reported in the manuscript. When new simulations are run, the data will appear in this folder (the old files will be overwritten).
- */processed_data*: this folder contains the relevant (averaged) go and stop estimates for each simulated 'participant'
- */summary_data*: this folder contains the plots with the data summaries (as reported in the manuscript).

### The 'power_tests' folder
The second set of simulations ('power_tests' folder) explored how different parameters affected the power to detect SSRT differences. This code can be adjusted to determine the required sample size under varying conditions, or acceptable levels of go omissions and RT distribution skew.

In the main 'power_tests' folder, three R scripts can be found:

- *simulation.R*: the R code used to simulate the data. To reduce the number of files, SSRT is immediately estimated.
- *effectsize.R*: the R code used to determine true and observed effect size for t-tests.
- *power.R*: the R code used to determine observed power for each condition and to create the figures for the manuscript.

The folder also contains three subfolders:

- */simulated_data*: this folder contains the simulated data reported in the manuscript. When new simulations are run, the data will appear in this folder (and old files are overwritten).
- */processed_data*: this folder contains the achieved power estimates, as determined by the first part of the power.R script. These data are saved in case the figures (second part of the power.R script) have to be updated (without having to run the power calculations again).
- */summary_data*: this folder contains the plots with the data summaries (as reported in the manuscript).
