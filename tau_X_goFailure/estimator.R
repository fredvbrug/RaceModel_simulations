# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# script to estimate the SSRTs
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

library (plyr)

# --- STEP:  analyse stop-signal data ---

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
    
    # -----------------------------------------------------------------------------
    # estimate 1: all no-signal trials are INcluded when the nth RT is determined
    # -----------------------------------------------------------------------------
    # nth RT with missed go's included 
    nthRT_all <- quantile(nosignal$RT.true, probs = presp, type = 6) 
    
    # SSRT estimated with the integration method = nthRT - ssd
    SSRTall <- nthRT_all - ssd
    SSRTall_diff <- SSRTall - SSRTtrue # calculate the difference with the true SSRT
    
    # ------------------------------------------------------------------------------
    # estimate 2: missed no-signal trials are EXcluded when the nth RT is determined
    # ------------------------------------------------------------------------------
    # nth RT with missed responses excluded
    nthRT_resp <- quantile(nosignal_resp$RT.true, probs = presp, type = 6) 
    
    # SSRT estimated with the integration method = nthRT - ssd
    SSRTresp <- nthRT_resp - ssd
    SSRTresp_diff <- SSRTresp - SSRTtrue # calculate the difference with the true SSRT
    
    # ------------------------------------------------------------------------------
    # estimate 3: p(respond|signal) is adjused using the Tannock et al. method
    #             missed no-signal trials are INcluded in the go distribution
    # ------------------------------------------------------------------------------
    # adjust p(respond) using the Tannock et al. method 
    pinhibit <- 1 - presp # the adjustment method works with p(inhibit|signal)
    pmiss <- 1 - mean(nosignal$resp) # determine the probability of a missed go response
    pinhibit_adj <- (pinhibit - pmiss)/(1 - pmiss)
    presp_adj <- 1 - pinhibit_adj
    
    # SSRT estimated with the integration method = nthRT - ssd
    nthRT_adj <- quantile(nosignal$RT.true, probs = presp_adj, type = 6) 
    SSRTadj <- nthRT_adj - ssd
    SSRTadj_diff <- SSRTadj - SSRTtrue # calculate the difference with the true SSRT
    
    # ------------------------------------------------------------------------------
    # estimate 4: for the sake of completeness, estimate SSRT with the mean method
    # ------------------------------------------------------------------------------
    mRT <- mean(nosignal_resp$RT.true)
    SSRTmean <- mRT - ssd
    SSRTmean_diff <- SSRTmean - SSRTtrue # calculate the difference with the true SSRT
    
    # ------------------------------------------------------------------------------
    # Also calculate no-signal RT for go trials with a response, 
    # and the difference with signal-respond RT
    # ------------------------------------------------------------------------------
    nosignal.resp.rt <- mean(nosignal_resp$RT.true)
    race.check <- nosignal.resp.rt - signal.resp.rt
    
    # ------------------------------------------------------------------------------
    # Return all data
    # ------------------------------------------------------------------------------
    return(data.frame(presp = presp, pmiss = pmiss, presp_adj = presp_adj, ssd = ssd, SSRTtrue = SSRTtrue, 
                      nthRT_all = nthRT_all, nthRT_resp = nthRT_resp, nthRT_adj = nthRT_adj,
                      SSRTall = SSRTall, SSRTresp = SSRTresp, SSRTadj = SSRTadj, SSRTmean = SSRTmean,
                      SSRTall_diff = SSRTall_diff, SSRTresp_diff = SSRTresp_diff, SSRTadj_diff = SSRTadj_diff, SSRTmean_diff = SSRTmean_diff,
                      sRT = signal.resp.rt, goRT = nosignal.resp.rt, raceCheck = race.check))
  }

  
# determine how many simulated data files we have 
files <- dir(path='./simulated_data', pattern='*.Rdata')

# process the data of all simulations
for (i in files) {
  name <- paste("./simulated_data/", i, sep = "")
  load (name)
  
  #extract some information about the simulation
  labels <-  unlist(strsplit(files[match(i, files)], "_")) 
  N <- as.numeric(labels[2]) # the number of trials
  Tau <- as.numeric(labels[4]) # the number of trials
  Pm <- as.numeric(unlist(strsplit(labels[6], "\\."))[1]) # the percentage of missed responses

  #change factors/variables
  simulation$signals <- factor (simulation$signals, levels = 0:1, labels = c("nosignal", "signal"))
  simulation$subject <- factor(simulation$subject)
  simulation$Ntrials <- factor(N)
  simulation$pmissReq <- factor(Pm)
  simulation$tau <- factor(Tau)
  
  # determine if there was a 'response' or not
  simulation$resp <- ifelse(simulation$outcome == "s-inhibit", 0, 1) # obviously, there is no response on a signal-inhibit
  simulation$resp <- ifelse(simulation$RT.true == 9999, 0, simulation$resp) # when RT = 9999, trial = go failure
  
  # analyse the stop-signal data  
  signal.cast <- ddply(simulation, .(subject, Ntrials, pmissReq, tau), funcSignal)
  
  # save the outcome (Rdata for R-users; csv copy for others)
  filename <- sprintf ("./processed_data/output_Ntrials_%d_RTtau_%d_RTmiss_%d.Rdata", N, Tau, Pm)
  save (signal.cast, file = filename)

  filename <- sprintf ("./processed_data/output_Ntrials_%d_RTtau_%d_RTmiss_%d.csv", N, Tau, Pm)
  write.csv (signal.cast, file = filename, row.names=FALSE)
}


# create a new data frame to combine all the data
data <- data.frame() 

# determine how many processed data files we have 
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

# create numeric version of the subject factor; will be used in some analyser functions
data$subjectNum <- as.numeric(data$subject)

# save the combined data file
save (data, file = './processed_data/combined_individual_data.Rdata')

