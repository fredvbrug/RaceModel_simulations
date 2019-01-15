# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# script to process the SSRT estimates
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

# clear everything
rm(list = ls()) 

# load libraries
library (ggplot2) # to plot the data
library (doBy) # for data summaries
library (reshape) # to melt and cast data
library (reshape2) # to melt and cast data
library (plyr)
library (Hmisc)

# -----------------------------------------------------------------------------
# STEP: load all the relevant output files
# -----------------------------------------------------------------------------
# create a new data frame to combine all the data
data <- data.frame() 

# determine how many simulated data files we have 
files <- dir(path='./processed_data', pattern='*.Rdata') 

# load the files and combine them
for (i in files) {
  name <- paste("./processed_data/", i, sep = "") # path to the file
  load(name) # load the file

  # add the content to the data frame
  data <- rbind (data, signal.cast)
  rm(signal.cast) 
}

# reorder levels pmissReq and tau
data$pmissReq = factor(data$pmissReq,levels(data$pmissReq)[c(1,5,2,3,4)])
data$tau = factor(data$tau,levels(data$tau)[c(1,5,2,3,4)])

# determine if signal-respond RT is sometimes longer than no-signal RT
check1.violations <- subset(data, raceCheck < 0)
nrow(check1.violations)
table(check1.violations$tau, check1.violations$pmissReq, check1.violations$Ntrials)

# does this influence the SSRT estimates? 
check2 <- subset(data, Ntrials == 100 & tau %in% c('1', '50') & raceCheck > 0)
check2.violations <- subset(data, Ntrials == 100 & tau  %in% c('1', '50') & raceCheck < 0)
check2.all <- subset(data, Ntrials == 100 & tau  %in% c('1', '50'))
  
mean(check2$SSRTall_diff)
mean(check2.violations$SSRTall_diff)
mean(check2.all$SSRTall_diff)