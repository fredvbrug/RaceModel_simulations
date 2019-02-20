# -----------------------------------------------------------------------------
# script to generate the simulated data for the power analyses
# -----------------------------------------------------------------------------

rm(list = ls()) # clear everything
library (gamlss) # library to create random ex-gaussian distributions
library(truncnorm) # library to create truncated normal distribution

# ---- auxiliary functions ----
# function to run the actual simulation
simulator <- function (pSIM, pRT, pSSRT, set) {
  simulation <- data.frame() # create empty data frame
  
  for (s in 1:pSIM$NSIM){
    for (g in 1:nrow(pRT)){
      #reset the start SSD value at the beginning of each simulation
      SSD = pSIM$SSD.start +  pRT$tau[g]
      
      # get an ex-gaussian RT distribution
      ## step 1: determine mean of the subject; this mean is allowed to differ from population mean
      ## but the lowest possible subject mean = 300 ms
      sub.mean <- rtruncnorm(1, a=300, b=Inf, mean = pRT$mu[g], sd = pRT$sdpop)
      
      ## step 2: get for this subject an ex-gaussian go-RT distribution
      RT.true <- rexGAUS(pSIM$NTRIALS, mu=sub.mean, sigma=pRT$sigma[g], nu= pRT$tau[g]) # ex-gaussian RT distribution
      
      ## step 3a: check for missed trials...
      RT.true <- ifelse(RT.true > pSIM$cutoff, 9999, RT.true) # label these with 9999
      
      ## step 3b: ...and insert extra ones if needed
      RT.true <- sort(RT.true)  # start with putting the slow trials ('original' go failures) at the end
      missed <- sum(RT.true>pSIM$cutoff) # determine how many 'go failures' occur
      extra_missed <- (pSIM$NTRIAL * pRT$miss[g]/100) - missed # determine how many extra 'go failures' we need
      
      if (extra_missed > 0){
        trial_missed <- sample(1:(pSIM$NTRIAL-missed), extra_missed) # get a random selection of trials that we will turn into a 'go failure'
        RT.true[trial_missed] <- 9999 # replace the values for these trials
      }
      
      # step 4: randomize the RTs again
      RT.true <- sample(RT.true)
      
      # get an ex-gaussian SSRT distribution
      ## step 1: determine mean of the subject; this mean is allowed to differ from population mean
      ## but the lowest possible subject mean = 100 ms
      sub.mean <- rtruncnorm(1, a=100, b=Inf, mean = pSSRT$mu[g], sd = pSSRT$sdpop)
      
      # step 2: get a ex-gaussian distribution
      SSRT.true <- rexGAUS(pSIM$NTRIALS, mu=sub.mean, sigma=pSSRT$sigma[g], nu= pSSRT$tau[g]) # ex-gaussian SSRT distribution
      
      # make a stop-signal list & adjust SSRTs
      signals <- sample(1:(pSIM$NTRIAL)) # make a stop-signal list
      signals <- ifelse(signals > (pSIM$NTRIAL * pSIM$PSIGNALS), 0, 1)
      SSRT.true = SSRT.true * signals
      
      # do a very simple race simulation for NTRIALS
      # if RT > SSRT + SSD, then signal-inhibit; else signal-respond
      # SSD will be adjusted accordingly
      used.ssd <- rep(0, pSIM$NTRIALS) #define/reset variable 'used.ssd'
      outcome <-  rep(0, pSIM$NTRIALS) #define/reset variable 'outcome'
      race <- rep(0, pSIM$NTRIALS)  #define/reset variable 'race'' (difference relative finishing times)
      
      # start a 'race' loop with tracking
      for (i in 1:pSIM$NTRIALS){
        used.ssd[i] <- signals[i] * SSD
        if (signals[i] == 0)
          outcome[i] <- "no-signal"
        else {
          race[i] <- RT.true[i] - (SSRT.true[i] + SSD)
          if (race[i] > 0){
            outcome[i] <- "s-inhibit"
            SSD <- SSD + 50
          } else{
            outcome[i] <- "s-respond"
            SSD <- SSD - 50
          }
        }
      }
      
      # put everything in a single data frame, and then immediately process these data
      tmp <- data.frame(signals, RT.true, SSRT.true, used.ssd, outcome)
      tmp$signals <- factor (tmp$signals, levels = 0:1, labels = c("nosignal", "signal"))
      tmp$resp <- ifelse(tmp$outcome == "s-inhibit", 0, 1) # obviously, there is no response on a signal-inhibit
      tmp$resp <- ifelse(tmp$RT.true == 9999, 0, tmp$resp) # when RT = 9999, trial = go failure
      
      tmp.processed <- funcSignal(tmp)
      tmp.processed$group <- g
      tmp.processed$subject <- s
      
      # combine the simulation with the other simulations
      simulation <- rbind(simulation, tmp.processed)
    }
    
    # print the progress
    text <- sprintf ("simulation %d, pair nmbr. %d", set, s)
    print(text)
  }
  return(simulation)
}

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



# ---- define simulation parameters ----

# define some variables that are fixed across groups
# NSIM = number of simulated 'subjects'
# NTRIALS = total number of trials (no-signal + signal) per simulation
# SSD.start = start value of SSD (tau of GO RT will be added though)
# PSIGNALS = # proportion of stop signals
# cut-off = MAX GO RT
pSIM <- data.frame(NSIM = 5000, NTRIALS = 100, SSD.start = 300, PSIGNALS = .25, cutoff = 1500)

# determine the shape of RT distribution for the two groups:
# mu: population mean of normal part of RT distribution
# sdpop: mu (subject) is a random value from a normal population distribution with mean = mu(pop) and SD = sdpop
# sigma: standard deviation of normal part of RT distribution
# tau: mean of exponential part of the RT distribution
# miss: percentage missed go responses
pRT <- data.frame (mu = c(500,500),
                   sdpop = c(100,100),
                   sigma = c(50,50),
                   tau = c(50,50),
                   miss = c(0,0))


# determine shape of SSRT distrbution
# mu: population mean of normal part of SSRT distribution
# sdpop: mu (subject) is a random value from a normal population distribution with mean = mu(pop) and SD = sdpop
# sigma: standard deviation of normal part of SSRT distribution
# tau: mean of exponential part of the SSRT distribution
pSSRT <- data.frame (mu = c(200,200),
                     sdpop = c(40,40),
                     sigma = c(20,20),
                     tau = c(10,10))

# do the simulations for two groups using a couple of
# parameter combinations and write data file
# Adjust some values. Note that mean(subject) = mu(subject) + tau(subject)
# Thus, if both mu and tau are increased with 5 ms, mean(subject) increases with +10 ms
ntrials.adj <- c(0, 100, 300)
go.adjust <-  data.frame(mu=c(0,25,75),sdpop=0, sigma=c(0,2.5,7.5),tau=c(0,25,75),miss=c(0,5,10))
stop.adjust <- data.frame(mu=c(10,15),sdpop=0, sigma=c(1,1.5),tau=c(5,10))

counter <- 1
sim.info <- data.frame()

for (n in 1:length(ntrials.adj)){
  for (g in 1:nrow(go.adjust)){
    for (s in 1:nrow(stop.adjust)){
      # adjust ntrials
      pSIM.adj <- pSIM
      pSIM.adj$NTRIALS <- pSIM$NTRIALS + ntrials.adj[n]
      
      # adjust parameters for the go distribution
      # mean RT and p(miss) can differ between groups
      pRT.adj <- pRT
      pRT.adj[2,] <- pRT.adj[2,]+go.adjust[g,]
      
      # adjust parameters for the stop distribution
      # mean SSRT can differ between groups
      pSSRT.adj <- pSSRT
      pSSRT.adj[2,] <- pSSRT.adj[2,]+stop.adjust[s,]
      
      # run the simulation and save the data
      simulation <- simulator(pSIM.adj, pRT.adj, pSSRT.adj, counter)
      file.name <- sprintf('./simulated_data/simulation_set_%d.Rdata', counter)
      save(simulation, file = file.name)
      
      # save the information about the simulation
      sim.current <- data.frame(simulation = counter,
                                ntrials = pSIM$NTRIALS + ntrials.adj[n],
                                G1_Gmu = pRT.adj$mu[1],
                                G1_Gsdpop = pRT.adj$sdpop[1],
                                G1_Gsigma = pRT.adj$sigma[1],
                                G1_Gtau = pRT.adj$tau[1],
                                G1_Gmiss = pRT.adj$miss[1],
                                G1_Smu = pSSRT.adj$mu[1],
                                G1_Ssdpop = pSSRT.adj$sdpop[1],
                                G1_Ssigma = pSSRT.adj$sigma[1],
                                G1_Stau = pSSRT.adj$tau[1],
                                G2_Gmu = pRT.adj$mu[2],
                                G2_Gsdpop = pRT.adj$sdpop[2],
                                G2_Gsigma = pRT.adj$sigma[2],
                                G2_Gtau = pRT.adj$tau[2],
                                G2_Gmiss = pRT.adj$miss[2],
                                G2_Smu = pSSRT.adj$mu[2],
                                G2_Ssdpop = pSSRT.adj$sdpop[2],
                                G2_Ssigma = pSSRT.adj$sigma[2],
                                G2_Stau = pSSRT.adj$tau[2]
      )
      sim.info <- rbind(sim.info, sim.current)
      counter <- counter+1
    }
  }
}

# save sim.info file
save(sim.info, file = './simulated_data/information_simulations.Rdata')
