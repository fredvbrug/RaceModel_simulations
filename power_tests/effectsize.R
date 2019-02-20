# -----------------------------------------------------------------------------
# script to determine effect size for t-tests
# -----------------------------------------------------------------------------
rm(list = ls()) # clear everything
library (plyr) # library for ddply
library(pwr) # library for some initial power tests

# --- process data of all simulations ---

# get the file names
files <- dir(path='./simulated_data_test/', pattern='simulation_')

# process all files
signal.cast <- data.frame() # create empty data frame

for (i in files) {
  name <- paste("./simulated_data_test/", i, sep = "")
  load (name)
  
  signal.cast <- simulation
  
  print(name)
  print(ddply(signal.cast, .(group), summarise, mean.true =mean(SSRTtrue), sd.true = sd(SSRTtrue), mean.est =mean(SSRTall), sd.est = sd(SSRTall), sd.go = sd(goRT)))
  test <- t.test(SSRTall ~ group, data = signal.cast, paired = F)
  
  d <- abs(test$statistic[[1]]*sqrt(1/nrow(subset(signal.cast,group==1))+1/nrow(subset(signal.cast,group==2)))) # Calculate d
  print(pwr.t.test(n = c(1:6)*32, d = d, sig.level = 0.05, type = "two.sample", alternative = "two.sided")$power)
  print(d)
  
  test <- t.test(SSRTtrue ~ group, data = signal.cast, paired = F)
  dtrue <- abs(test$statistic[[1]]*sqrt(1/nrow(subset(signal.cast,group==1))+1/nrow(subset(signal.cast,group==2)))) # Calculate d
  print(dtrue)
}