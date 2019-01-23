# -----------------------------------------------------------------------------
# script to estimate the SSRTs etc.
# -----------------------------------------------------------------------------
rm(list = ls()) # clear everything
library (plyr) # library for ddply
library(pwr) # library for some initial power tests

# --- auxiliary functions ---

# function to analyse stop-signal data at once...
  funcSignal <- function(data){

    # signal data: prespond, ssd, true SSRT, and signal-respond RT
    signal <- subset(data, signals == 'signal')
    presp <-  mean(signal$resp)
    ssd <- mean(signal$used.ssd)
    SSRTtrue <- mean(signal$SSRT.true)

    signal.resp <- subset(signal, resp == '1')
    signal.resp.rt <- mean(signal.resp$RT.true)

    # no signal data: with and without missed responses
    # for the missed responses, set RT to max RT of the subject/condition
    nosignal <- subset(data, signals == 'nosignal')
    nosignal_resp <- subset(nosignal, resp == '1')
    nosignal$RT.true <- ifelse(nosignal$RT.true == 9999, max(nosignal_resp$RT.true), nosignal$RT.true)
    pmiss <- 1 - mean(nosignal$resp) # determine the probability of a missed go response

    # --- estimate 1 --- all no-signal trials are INcluded when the nth RT is determined
    ## determine nth RT
    nthRT_all <- quantile(nosignal$RT.true, probs = presp, type = 6)

    ## SSRT(integration) = nthRT - ssd
    SSRTall <- nthRT_all - ssd
    SSRTall_diff <- SSRTall - SSRTtrue # calculate the difference with the true SSRT

    #  --- estimate 4 ---  estimate SSRT with the mean method
    mRT <- mean(nosignal_resp$RT.true)
    SSRTmean <- mRT - ssd
    SSRTmean_diff <- SSRTmean - SSRTtrue # calculate the difference with the true SSRT

    # Also calculate no-signal RT for go trials with a response,
    # and the difference with signal-respond RT
    nosignal.resp.rt <- mean(nosignal_resp$RT.true)
    race.check <- nosignal.resp.rt - signal.resp.rt

    # Return all data
    return(data.frame(presp = presp, pmiss = pmiss, ssd = ssd, SSRTtrue = SSRTtrue,
                      nthRT_all = nthRT_all, SSRTall = SSRTall, SSRTmean = SSRTmean,
                      SSRTall_diff = SSRTall_diff, SSRTmean_diff = SSRTmean_diff,
                      sRT = signal.resp.rt, goRT = nosignal.resp.rt, raceCheck = race.check))
  }


# --- process data of all simulations ---

# get the file names
files <- dir(path='./simulated_data/', pattern='simulation_')

# process all files
for (i in files) {
  name <- paste("./simulated_data/", i, sep = "")
  load (name)

  #change factors/variables
  simulation$signals <- factor (simulation$signals, levels = 0:1, labels = c("nosignal", "signal"))
  simulation$subject <- factor(simulation$subject)

  # determine if there was a 'response' or not
  simulation$resp <- ifelse(simulation$outcome == "s-inhibit", 0, 1) # obviously, there is no response on a signal-inhibit
  simulation$resp <- ifelse(simulation$RT.true == 9999, 0, simulation$resp) # when RT = 9999, trial = go failure

  # analyse the stop-signal data
  signal.cast <- ddply(simulation, .(subject, group), funcSignal)

  # save the processed data for the simulation
  file.name <- paste('./processed_data/', unlist(strsplit(i, "\\."))[1],'_processed.Rdata', sep="")
  save (signal.cast, file = file.name)

  # print progress (and immediately check if the simulations make sense)
  print(unlist(strsplit(i, "\\."))[1])
}

# --- perform some checks to check quality of the simulations etc. ---

# get the file names
files <- dir(path='./processed_data/', pattern='simulation_')

# load the files, and calculate for each simulation mean, sd, and power.
for (i in files) {
  name <- paste("./processed_data/", i, sep = "")
  load (name)

  print('-------------------')
  print(name)
  print(ddply(signal.cast, .(group), summarise, mean.true =mean(SSRTtrue), sd.true = sd(SSRTtrue), mean.est =mean(SSRTall), sd.est = sd(SSRTall), sd.go = sd(goRT)))
  test <- t.test(SSRTall ~ group, data = signal.cast, paired = F)
  d <- abs(test$statistic[[1]]*sqrt(1/nrow(subset(signal.cast,group==1))+1/nrow(subset(signal.cast,group==2)))) # Calculate d
  print(pwr.t.test(n = c(1:6)*32, d = d, sig.level = 0.05, type = "two.sample", alternative = "two.sided")$power)
}
